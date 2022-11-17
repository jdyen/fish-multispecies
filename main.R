# script to fit Poisson-log Matrix Normal model to fish assemblage abundance 
#   data using MCMC sampling with rstan

# Author: Jian Yen
# Date created: 2021/12/01
# Updated: 2022/06/29

# packages for data access and manipulation
library(aae.hydro)
library(qs)
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)

# packages for analyses
library(rstan)

# some helper functions
source("R/load-data.R")

# load fish data
vefmap <- qread("data/vefmap-compiled.qs")
vefmap <- vefmap %>% mutate(
  reach = check_reach(reach, waterbody, site_name)
)

# load flow data
flow <- qread("data/flow-compiled.qs")

# impute missing flow values
for (i in seq_along(flow)) {
  flow[[i]] <- flow[[i]] %>% filter(variable_code == "141.00")
  flow[[i]]$value <- impute_year(flow[[i]]$value, flow[[i]]$date_formatted)
  flow[[i]]$value <- impute_rolling(flow[[i]]$value, recursive = TRUE, max_iter = 50)
}

# define size class bins
## DO LIFE STAGES FOR EACH SPECIES? c(early, juv, adult)?
# size_breaks <- c(0, 30, 50, 100, 200, 400, 1600)
# vefmap <- vefmap %>% mutate(
#   size_class = cut(length_mm, breaks = size_breaks, labels = FALSE)
# )
# VERSION USING QUANTILES TO DEFINE 3 CATEGORIES 
#   CAN EASILY REPLACE WITH SPECIES-SPECIFIC LENGTH THRESHOLDS
vefmap <- vefmap %>% 
  group_by(scientific_name) %>%
  summarise(
    low_cutoff = quantile(length_mm, probs = 0.25, na.rm = TRUE),
    high_cutoff = quantile(length_mm, probs = 0.75, na.rm = TRUE)
  ) %>%
  right_join(vefmap, by = "scientific_name") %>%
  mutate(size_class = ifelse(length_mm < low_cutoff, 1, ifelse(length_mm > high_cutoff, 3, 2)))

# calculate electro effort
effort <- vefmap %>% 
  distinct(sdate, site_name, waterbody, reach, seconds, id_survey, id_surveyevent) %>%
  group_by(sdate, site_name, waterbody, reach, id_survey) %>%
  summarise(
    effort = sum(seconds)
  )

# define priority species
priority_spp <- c(
  "cyprinus_carpio",
  "maccullochella_peelii", 
  "macquaria_ambigua", 
  "perca_fluviatilis", 
  "gambusia_holbrooki", 
  "gadopsis_marmoratus", 
  "melanotaenia_fluviatilis", 
  "retropinna_semoni", 
  "bidyanus_bidyanus", 
  "maccullochella_macquariensis"
)

# compile abundances by species and size class for records with
#   length data only
size_abundance <- vefmap %>%
  filter(!is.na(length_mm), scientific_name %in% priority_spp) %>%
  group_by(waterbody, reach, site_name, id_survey, sdate, scientific_name, size_class, gear_type) %>%
  summarise(
    count = n()
  ) %>%
  left_join(effort, by = c("sdate", "waterbody", "reach", "site_name", "id_survey")) %>%
  ungroup %>%
  arrange(waterbody, reach, sdate)

# grab a list of unique sites/years and fill an array with counts for each
#   site, year, species, and size class
info <- size_abundance %>% 
  select(waterbody, reach, site_name, id_survey, sdate, effort, gear_type) %>%
  mutate(syear = year(sdate)) %>%
  distinct()
info <- add_gauge(info)
nsite <- length(unique(info$site_name))
nyear <- length(unique(info$syear))
nsp <- length(unique(size_abundance$scientific_name))
nclass <- length(unique(size_abundance$size_class))
y <- array(0, dim = c(nsite, nyear, nsp, nclass))
dimnames(y) <- list(
  sort(unique(info$site_name)),
  sort(unique(info$syear)),
  sort(unique(size_abundance$scientific_name)),
  sort(unique(size_abundance$size_class))
)
for (i in seq_len(nrow(size_abundance))) {
  tmp <- size_abundance[i, ]
  y[tmp$site_name, as.character(year(tmp$sdate)), tmp$scientific_name, as.character(tmp$size_class)] <- tmp$count
}

# which sites were visited?
visited <- apply(y, c(1, 2), sum) > 0

# calculate covariates
threshold_list <- c(
  "goulburn_river_r4" = 600,
  "goulburn_river_r5" = 600,
  "broken_river_r2" = 15,
  "broken_river_r3" = 15,
  "broken_creek_r4" = 200,
  "broken_creek_r5" = 200,
  "campaspe_river_r2" = 10,
  "campaspe_river_r3" = 10,
  "campaspe_river_r4" = 10,
  "loddon_river_r2" = 25,
  "loddon_river_r3" = 25,
  "loddon_river_r4" = 25,
  "loddon_river_r5" = 125,
  "pyramid_creek_rNA" = 90,
  "wimmera_river_r2" = 15,
  "wimmera_river_r3" = 15,
  "burnt_creek_r1" = 1,
  "burnt_creek_r2" = 1,
  "little_murray_river_rNA" = 1,
  "mackenzie_river_r2" = 10,
  "mackenzie_river_r3" = 10,
  "mt_william_creek_rNA" = 5
)
size_abundance_covariates <- get_metrics(flow, info, threshold_list)
covars <- c("ave_winter", "ave_spring", "ave_summer", "ave_antecedent", "low_flow", "cv_flow", "survey_flow")
ncovar <- length(covars)
X <- array(NA, dim = c(nsite, nyear, ncovar))
dimnames(X) <- list(
  sort(unique(info$site_name)),
  sort(unique(info$syear)),
  covars
)
for (i in seq_len(nrow(size_abundance_covariates))) {
  tmp <- size_abundance_covariates[i, ]
  X[tmp$site_name, as.character(year(tmp$sdate)), covars] <- unlist(tmp[covars])
}
X <- lapply(seq_len(nyear), function(i, .x) .x[, i, ], .x = X)
X <- do.call(rbind, X)
X <- apply(X, 2, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x))
X_std <- sweep(X, 2, colMeans(X), "-")
X_std <- sweep(X_std, 2, apply(X_std, 2, sd), "/")

# calculate effort data
eff <- array(0, dim = c(nsite, nyear))
dimnames(eff) <- list(
  sort(unique(info$site_name)),
  sort(unique(info$syear))
)
gear <- array(0, dim = c(nsite, nyear))
dimnames(gear) <- dimnames(eff)
for (i in seq_len(nrow(info))) {
  tmp <- info[i, ]
  eff[tmp$site_name, as.character(year(tmp$sdate))] <- tmp$effort
  gear[tmp$site_name, as.character(year(tmp$sdate))] <- tmp$gear_type
}
eff <- ifelse(eff == 0, median(eff[eff > 0]), eff)

# define some extra vars
idx <- match(dimnames(y)[[1]], info$site_name)
river <- info$waterbody[idx]
reach <- paste(river, info$reach[idx], sep = "_")
river <- rebase_factor(river)
reach <- rebase_factor(reach)

# expand over years to match response matrix dims
river <- rep(river, times = nyear)
reach <- rep(reach, times = nyear)
site <- rep(seq_len(nsite), times = nyear)
year <- rep(seq_len(nyear), each = nsite)

# initial matrix normal model (AR1)
y_flat <- aperm(array(c(aperm(y, c(3, 4, 1, 2))), dim = c(nsp * nclass, nsite, nyear)), c(2, 3, 1))
y_diminished <- do.call(rbind, lapply(seq_len(nyear), function(i, .x) .x[, i, ], .x = y_flat))
y_diminished <- y_diminished[c(visited), ]
X_diminished <- X_std[c(visited), ]
eff_diminished <- c(eff)[c(visited)]
river <- rebase_factor(river[c(visited)])
reach <- rebase_factor(reach[c(visited)])
site <- rebase_factor(site[c(visited)])
year <- rebase_factor(year[c(visited)])
gear <- rebase_factor(gear[c(visited)])
dat <- list(
  N = nrow(y_diminished), K = ncovar,
  nsp = nsp, nclass = nclass, nq = nsp * nclass,
  X = X_diminished,
  y = y_diminished,
  log_effort = log(eff_diminished),
  nriver = max(river),
  nreach = max(reach),
  nsite = max(site),
  nyear = max(year),
  ngear = max(gear),
  river = river,
  reach = reach,
  site = site,
  year = year,
  gear = gear,
  main_scale = 3.,
  sigma_scale = 3.,
  sigma_scale2 = 3.,
  ar_scale = 1.
)

# add a flattened version of y
## TRANSPOSED BECAUSE MODEL DEFINES mu[Q, N] and then flattens this
dat$yflat <- c(t(dat$y))

## TRANSPOSE X FOR CURRENT VERSION TOO
dat$X <- t(dat$X)

# and some extra terms for the zero-inflation part of the model
dat$nflat <- length(dat$yflat)
dat$zero_idx <- which(dat$yflat == 0)
dat$nzero <- length(dat$zero_idx)
dat$nonzero_idx <- which(dat$yflat > 0)
dat$notzero <- length(dat$nonzero_idx)

# add spline settings
# dat$hs_df <- rep(2, times = dat$K)
# dat$hs_df_global <- rep(2, times = dat$K)
# dat$hs_df_slab <- rep(2, times = dat$K)
# dat$hs_scale_global <- rep(3., times = dat$K)
# dat$hs_scale_slab <- rep(3., times = dat$K)

# add extra info on missing surveys
y_previous_idx <- matrix(0, nrow = nsite, ncol = nyear)
y_previous_idx[visited] <- seq_len(sum(visited))
y_previous_idx <- cbind(rep(0, nsite), y_previous_idx[, -nyear])
visited_previous <- y_previous_idx > 0
y_previous_idx <- c(y_previous_idx)[visited]  # filter on visited (not visited_previous) to match y
visited_previous <- c(visited_previous)[visited]  # filter on visited (not visited_previous) to match y

# add these to data set
dat$prev_idx <- y_previous_idx
dat$visited_prev <- visited_previous

# create a variable with y + 1 to save converting arrays to matrices in Stan
dat$yp1 <- dat$y + 1

# add index for missing observations and update prev_idx to match these correctly
missing_idx <- !dat$visited_prev
dat$nmissing <- sum(missing_idx)
dat$prev_idx[missing_idx] <- seq_len(dat$nmissing)

# settings for MCMC
seed <- 1125343481

# settings for MCMC
iter <- 2000
warmup <- 1000
chains <- 4
cores <- 4
thin <- 2

# compile and sample from matrix normal model
stan_file <- file.path("src/matrix_normal_nocovar.stan")
mod <- stan_model(file = stan_file)

# define some initial conditions
empirical_corr_sp <- cor(apply(y, c(1, 3), sum))
empirical_corr_class <- cor(apply(y, c(1, 4), sum))
init <- lapply(
  seq_len(chains),
  function(x) list(
    L_sp = t(chol(empirical_corr_sp)),
    L_class = t(chol(empirical_corr_class))
  )
)
draws_mat_normal <- sampling(
  object = mod,
  data = dat,
  init = init,
  pars = c(
    "Sigma_sp",
    "Sigma_class",
    "alpha",
    "sigma_river",
    "sigma_reach",
    "sigma_site",
    "sigma_year",
    "sigma_gear",
    "mu"
  ),
  chains = chains,
  iter = iter,
  warmup = warmup,
  thin = thin,
  cores = cores,
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed = seed
)

# save fitted
qsave(draws_mat_normal, file = "outputs/fitted/mat-normal-draws.qs")
  
# remove fitted model to free up space for other models
rm(draws_mat_normal)

# compile and sample from species-grouped model
stan_file <- file.path("src/matrix_normal_single_nocovar.stan")
mod <- stan_model(file = stan_file)

# need a new data file (species are target here, repeat for class below)
dat_spp <- dat
dat_spp$ntarget <- dat_spp$nsp
dat_spp$nother <- dat_spp$nclass

# define some initial conditions
init <- lapply(
  seq_len(chains),
  function(x) list(
    L_target = t(chol(empirical_corr_sp))
  )
)

# sample from posterior
draws_mat_normal_species <- sampling(
  object = mod,
  data = dat_spp,
  init = init,
  pars = c(
    "Sigma_target",
    "sigma_other",
    "alpha",
    "sigma_river",
    "sigma_reach",
    "sigma_site",
    "sigma_year",
    "sigma_gear",
    "mu"
  ),
  chains = chains,
  iter = iter,
  warmup = warmup,
  thin = thin,
  cores = cores,
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed = seed
)

# save fitted
qsave(draws_mat_normal_species, file = "outputs/fitted/mat-normal-species-draws.qs")

# remove fitted model to free up space for other models
rm(draws_mat_normal_species)

# repeat for class-grouped model
dat_class <- dat
dat_class$ntarget <- dat_spp$nclass
dat_class$nother <- dat_spp$nsp

# define some initial conditions
init <- lapply(
  seq_len(chains),
  function(x) list(
    L_target = t(chol(empirical_corr_class))
  )
)

# sample from posterior
draws_mat_normal_class <- sampling(
  object = mod,
  data = dat_class,
  init = init,
  pars = c(
    "Sigma_target",
    "sigma_other",
    "alpha",
    "sigma_river",
    "sigma_reach",
    "sigma_site",
    "sigma_year",
    "sigma_gear",
    "mu"
  ),
  chains = chains,
  iter = iter,
  warmup = warmup,
  thin = thin,
  cores = cores,
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed = seed
)

# save fitted
qsave(draws_mat_normal_class, file = "outputs/fitted/mat-normal-class-draws.qs")

# remove fitted model to free up space for other models
rm(draws_mat_normal_class)

# compile and sample from species-grouped model
stan_file <- file.path("src/multi_normal_nocovar.stan")
mod <- stan_model(file = stan_file)

# define some initial conditions
empirical_corr <- kronecker(empirical_corr_sp, empirical_corr_class)
init <- lapply(
  seq_len(chains),
  function(x) list(
    L = t(chol(empirical_corr))
  )
)

# sample from posterior
draws_multi_normal <- sampling(
  object = mod,
  data = dat,
  init = init,
  pars = c(
    "Sigma",
    "alpha",
    "sigma_river",
    "sigma_reach",
    "sigma_site",
    "sigma_year",
    "sigma_gear",
    "mu"
  ),
  chains = chains,
  iter = iter,
  warmup = warmup,
  thin = thin,
  cores = cores,
  control = list(adapt_delta = 0.9, max_treedepth = 15),
  seed = seed
)

# save fitted
qsave(draws_multi_normal, file = "outputs/fitted/multi-normal-draws.qs")

# remove fitted model to free up space for other models
rm(draws_multi_normal)

# summarise fitted models
file_names <- dir("outputs/fitted")
draws <- vector("list", length = length(file_names))
for (i in seq_along(draws))
  draws[[i]] <- qread(paste0("outputs/fitted/", file_names[i]))

# summarise each and extract diagnostics
for (i in seq_along(draws)) {
  summary_tmp <- summary(draws[[i]])$summary
  diagnostics[[i]] <- data.frame(
    model = file_names[i],
    par = rownames(summary_tmp),
    rhat = summary_tmp[, "Rhat"],
    n_eff = summary_tmp[, "n_eff"]
  )
  rm(summary_tmp)
}
diagnostics <- do.call(rbind, diagnostics)

# summarise the covariance matrices
# draws_mat <- as.matrix(draws_mat_normal)
# Sigma_sp <- draws_mat[, grepl("Sigma_sp", colnames(draws_mat))]
# Sigma_class <- draws_mat[, grepl("Sigma_class", colnames(draws_mat))]
# cov_sp <- matrix(apply(Sigma_sp, 2, median), ncol = dat$nsp)
# cov_cl <- matrix(apply(Sigma_class, 2, median), ncol = dat$nclass)

# and convert these to correlations for plotting
# cor_sp <- cov2cor(cov_sp)
# diag(cor_sp) <- 0
# cor_class <- cov2cor(cov_cl)
# diag(cor_class) <- 0
