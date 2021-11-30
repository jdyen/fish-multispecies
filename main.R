# script to fit Poisson-log Matrix Normal model to fish assemblage abundance 
#   data using variational and MCMC sampling with cmdstan

# Author: Jian Yen
# Date created: 2021/12/01

# packages for data access and manipulation
library(aae.hydro)
library(qs)
library(RPostgres)
library(DBI)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

# packages for analyses
library(cmdstanr)

# load data

# check if fish data exist
if ("vefmap-compiled.qs" %in% dir("data")) {
  
  # load from disk if possible
  vefmap <- qread("data/vefmap-compiled.qs")
  
} else {
  
  # connect to db
  con <- dbConnect(
    RPostgres::Postgres(),
    dbname = "arispatialdb",
    host = "ari-spatial-poc-db.cluster-custom-cepp1cnsvaah.ap-southeast-2.rds.amazonaws.com", 
    port = "5432", 
    user = "aae_user",
    password = rstudioapi::askForPassword("Database password")
  )
  
  # view flat vefmap database
  vefmap <- tbl(con, in_schema(sql("aquatic_data"), sql("v_vefmap_flat_data")))
  
  # remove species that won't be included in analyses, filter to target waterbodies,
  #   and grab EF data only
  vefmap <- vefmap %>% filter(
    gear_type %in% c("EFB", "EF_BM", "EF_BP", "EF_LB", "EF_MB"),
    !scientific_name %in% c("Anura spp.", "Maccullochella", "Hyriidae spp.", "unidentified turtle"),
    waterbody %in% c("GOULBURN", "WIMMERA", "LODDON", "CAMPASPE", "BROKEN", "PYRAMID CK", "LITTLE MURRAY", "MacKenzie", "Burnt Creek", "Mt William Creek")
  )
  
  # and collect
  vefmap <- vefmap %>% collect()
  
  # find nearest flow gauges for each site
  vefmap_coords <- vefmap %>% distinct(site_name, waterbody, lon, lat)
  gauges <- vector("list", length = nrow(vefmap_coords))
  for (i in seq_len(nrow(vefmap_coords))) {
    gauges[[i]] <- tbl(con, sql(
      paste0(
        "SELECT * FROM stream_network.get_closest_gauges(",
        vefmap_coords$lon[i], " , ", vefmap_coords$lat[i],
        ")"
      )
    )) %>%
      filter(
        field == "Y"
      ) %>%
      collect()
  }
  
  # disconnect from db
  dbDisconnect(con)
  
  # clean up species and waterbody names
  vefmap <- vefmap %>% mutate(
    scientific_name = clean_scinames(scientific_name),
    waterbody = clean_waterbody(waterbody),
    reach = add_reach(notes),
    reach = check_reach(reach, waterbody, site_name),
    region = add_region(waterbody)
  )
  
  # add manually identified gauges
  vefmap <- add_gauge(vefmap)
  
  # save to disk
  qsave(vefmap, file = "data/vefmap-compiled.qs")
  
}

# check if flow data exist
if ("flow-compiled.qs" %in% dir("data")) {
  
  # load from disk if possible
  flow <- qread("data/flow-compiled.qs")
  
} else {
  
  # find unique stations so we can download data just for those
  all_gauges <- vefmap %>%
    group_by(best_station, reach, waterbody) %>%
    summarise(
      min_date = min(sdate),
      max_date = max(sdate)
    ) %>%
    mutate(
      min_date = floor_date(min_date, unit = "years") - years(2),
      max_date = ymd(paste0(year(max_date), "-12-31"))
    )
  
  # go through each and grab data one gauge at a time
  flow <- vector("list", length = nrow(all_gauges))
  names(flow) <- all_gauges$best_station
  for (i in seq_along(flow)) {
    flow[[i]] <- fetch_hydro(
      sites = all_gauges$best_station[i],
      start = all_gauges$min_date[i],
      end = all_gauges$max_date[i],
      include_missing = TRUE,
      options = list(
        varfrom = c("100.00", "450.00"),
        varto = c("141.00", "450.00")
      )
    )
  }
  
  # save to disk
  qsave(flow, file = "data/flow-compiled.qs")
  
}

# impute missing flow values
for (i in seq_along(flow)) {
  flow[[i]] <- flow[[i]] %>% filter(variable_code == "141.00")
  flow[[i]]$value <- impute_year(flow[[i]]$value, flow[[i]]$date_formatted)
  flow[[i]]$value <- impute_rolling(flow[[i]]$value, recursive = TRUE, max_iter = 50)
}

# define size class bins
size_breaks <- c(0, 30, 50, 100, 200, 400, 1600)
vefmap <- vefmap %>% mutate(
  size_class = cut(length_mm, breaks = size_breaks, labels = FALSE)
)

# calculate electro effort
effort <- vefmap %>% 
  distinct(sdate, site_name, waterbody, reach, seconds, id_survey, id_surveyevent) %>%
  group_by(sdate, site_name, waterbody, reach, id_survey) %>%
  summarise(
    effort = sum(seconds)
  )

# compile abundances by species and size class for records with
#   length data only
size_abundance <- vefmap %>%
  filter(!is.na(length_mm)) %>%
  group_by(waterbody, reach, site_name, id_survey, sdate, scientific_name, size_class) %>%
  summarise(
    count = n()
  ) %>%
  left_join(effort, by = c("sdate", "waterbody", "reach", "site_name", "id_survey")) %>%
  ungroup %>%
  arrange(waterbody, reach, sdate)

size_abundance <- size_abundance %>%
  pivot_wider(
    id_cols = c("waterbody", "reach", "site_name", "id_survey", "sdate", "effort"),
    names_from = c(scientific_name, size_class),
    values_from = count,
    values_fill = 0
  )
size_abundance_counts <- size_abundance %>% select(7:ncol(size_abundance)) %>% as.matrix()
size_abundance_offset <- size_abundance %>% select(effort) %>% as.matrix()

# add some covariates
size_abundance_covariates <- get_metrics(flow, size_abundance)
size_abundance_covariates <- size_abundance_covariates %>%
  mutate(
    waterbody_id = rebase_factor(waterbody),
    reach_id = rebase_factor(paste(waterbody, reach, sep = "_")),
    site_id = rebase_factor(site_name),
    year_id = rebase_factor(year(sdate))
  ) %>%
  select(-effort)
size_abundance_covariates_clean <- size_abundance_covariates %>% 
  select(-waterbody, -reach, -site_name, -id_survey, -sdate) %>% 
  as.matrix

# standardise rownames to match everything up
rownames(size_abundance_counts) <- rownames(size_abundance_offset) <- rownames(size_abundance_covariates_clean) <- 
  paste(size_abundance$site_name, size_abundance$id_survey, sep = "_")
size_abundance_collated <- prepare_data(
  counts = size_abundance_counts,
  covariates = size_abundance_covariates_clean,
  offset = size_abundance_offset
)

# settings for variational Bayes and MCMC
seed <- 5623256

# settings for MCMC
iter <- 400
warmup <- 400
chains <- 2
refresh <- floor(iter / 5)

# initial matrix normal model (AR1)
y_time_flat <- aperm(array(c(aperm(y_time, c(3, 4, 1, 2))), dim = c(nsp * nclass, nsite, nyear)), c(2, 3, 1))
dat <- list(
  N = N_time, K = K,
  nsp = nsp, nclass = nclass,
  ntime = nyear,
  ar_order = 1L,
  nriver = nriver,
  nreach = nreach,
  river = river[seq_len(N_time)],
  reach = reach[seq_len(N_time)],
  X = X_sim,
  y = y_time_flat
)
file <- file.path("src/matrix_normal_ar.stan")
mod <- cmdstan_model(file)
# fit <- mod$sample(
#   data = dat, 
#   iter_warmup = warmup,
#   iter_sampling = iter,
#   chains = chains, 
#   parallel_chains = chains,
#   refresh = refresh
# )
fit <- mod$variational(
  data = dat, 
  seed = seed,
  output_samples = 1000
)
fit$save_object(file = "outputs/fitted/fit_struc.rds")

# unstructured multivariate normal model (AR1)
dat <- list(
  N = N_time, K = K,
  nn = nsp, nm = nclass, nq = nsp * nclass,
  ntime = nyear,
  ar_order = 1L,
  nriver = nriver,
  nreach = nreach,
  river = river[seq_len(N_time)],
  reach = reach[seq_len(N_time)],
  X = X_sim,
  y = y_time_flat,
  eta = 1
)
init_set <- lapply(
  seq_len(chains), 
  function(x) list(
    L = diag(nsp * nclass),
    eps = array(1, dim = c(N_time, K, nsp * nclass)),
    sigma = rep(1, nsp * nclass)
  )
)
file <- file.path("src/multi_normal_unstructured_ar.stan")
mod <- cmdstan_model(file)q
# fit <- mod$sample(
#   data = dat, 
#   init = init_set,
#   iter_warmup = iter,
#   iter_sampling = warmup,
#   chains = chains, 
#   parallel_chains = chains,
#   refresh = refresh
# )
fit <- mod$variational(
  data = dat, 
  seed = seed,
  output_samples = 1000
)
fit$save_object(file = "outputs/fitted/fit_unstruc.rds")

# MVN species, IID size classes
dat <- list(
  N = N_time, K = K,
  nsp = nsp, nclass = nclass,
  ntime = nyear,
  ar_order = 1L,
  nriver = nriver,
  nreach = nreach,
  river = river[seq_len(N_time)],
  reach = reach[seq_len(N_time)],
  X = X_sim,
  y = y_time_flat
)
file <- file.path("src/matrix_normal_ar_single_dim.stan")
mod <- cmdstan_model(file)
# fit <- mod$sample(
#   data = dat, 
#   iter_warmup = iter,
#   iter_sampling = warmup,
#   chains = chains, 
#   parallel_chains = chains,
#   refresh = refresh
# )
fit <- mod$variational(
  data = dat, 
  seed = seed,
  output_samples = 1000
)
fit$save_object(file = "outputs/fitted/fit_mvn_spp.rds")

# MVN size classes, IID species
y_time_flat <- aperm(array(c(aperm(aperm(y_time, c(1:2, 4, 3)), c(3, 4, 1, 2))), dim = c(15 * 8, 120, 5)), c(2, 3, 1))
dat <- list(
  N = N_time, K = K,
  nsp = nclass, nclass = nsp,  ## NOTE SWITCH HERE (_sp params are actually _class params)
  ntime = nyear,
  ar_order = 1L,
  nriver = nriver,
  nreach = nreach,
  river = river[seq_len(N_time)],
  reach = reach[seq_len(N_time)],
  X = X_sim,
  y = y_time_flat
)
# fit <- mod$sample(
#   data = dat, 
#   iter_warmup = iter,
#   iter_sampling = warmup,
#   chains = chains, 
#   parallel_chains = chains,
#   refresh = refresh
# )
fit <- mod$variational(
  data = dat, 
  seed = seed,
  output_samples = 1000
)

fit$save_object(file = "outputs/fitted/fit_mvn_size.rds")
