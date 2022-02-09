# script to fit Poisson-log Matrix Normal model to fish assemblage abundance 
#   data using variational and MCMC sampling with cmdstan

# Author: Jian Yen
# Date created: 2021/12/01

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

# define priority species
priority_spp <- c(
  'cyprinus_carpio',
  'maccullochella_peelii', 
  'macquaria_ambigua', 
  'perca_fluviatilis', 
  'gambusia_holbrooki', 
  'gadopsis_marmoratus', 
  'melanotaenia_fluviatilis', 
  'retropinna_semoni', 
  'bidyanus_bidyanus', 
  'maccullochella_macquariensis'
)

# compile abundances by species and size class for records with
#   length data only
size_abundance <- vefmap %>%
  filter(!is.na(length_mm), scientific_name %in% priority_spp) %>%
  group_by(waterbody, reach, site_name, id_survey, sdate, scientific_name, size_class) %>%
  summarise(
    count = n()
  ) %>%
  left_join(effort, by = c("sdate", "waterbody", "reach", "site_name", "id_survey")) %>%
  ungroup %>%
  arrange(waterbody, reach, sdate)

# grab a list of unique sites/years and fill an array with counts for each
#   site, year, species, and size class
info <- size_abundance %>% 
  select(waterbody, reach, site_name, id_survey, sdate, effort) %>%
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

## NOTE: some surveys at a single site occur twice in a year; these are currently
##   omitted from analyses due to data format (site x year)

# which sites were visited?
visited <- apply(y, c(1, 2), sum) > 0

# calculate covariates
size_abundance_covariates <- get_metrics(flow, info)
covars <- c("ave_spring", "ave_summer", "ave_antecedent", "low_flow", "cv_flow")
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
for (i in seq_len(nrow(info))) {
  tmp <- info[i, ]
  eff[tmp$site_name, as.character(year(tmp$sdate))] <- tmp$effort
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
  river = river,
  reach = reach,
  site = site,
  year = year
)

# think through this, should be identifying each survey based on
#  its position in diminished data set, then linking each year to
#  the previous year. Records 0 if site not visited in previous
#  year, which can be used in the Stan model to check for form of
#  ar_term.
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

# settings for variational Bayes and MCMC
seed <- 353532

# settings for MCMC
iter <- 4000
warmup <- 2000
chains <- 4
cores <- 4

# fit model with alternative parameterisation of matrix normal
file_alt <- file.path("src/matrix_normal_alt_param.stan")
mod_alt <- stan_model(file = file_alt)

# define some initial conditions
empirical_corr_sp <- cor(apply(y, c(1, 3), sum))
empirical_corr_class <- cor(apply(y, c(1, 4), sum))
init_alt <- lapply(
  seq_len(chains),
  function(x) list(
    L_sp = chol(empirical_corr_sp),
    L_class = chol(empirical_corr_class)
  )
)
draws_alt <- sampling(
  object = mod_alt,
  data = dat,
  init = init_alt,
  pars = c("Sigma_sp", "Sigma_class", "alpha", "beta", "rho", "tau", "sigma_river", "sigma_reach", "sigma_site", "sigma_year", "log_y_missing", "log_y_shift_missing"),
  chains = chains,
  iter = iter,
  warmup = warmup,
  cores = cores,
  control = list(adapt_delta = 0.75, max_treedepth = 15),
  seed = seed
)
vb_alt <- vb(
  object = mod_alt,
  data = dat,
  init = init_alt[[1]],
  pars = c("Sigma_sp", "Sigma_class", "alpha", "beta", "rho", "tau", "sigma_river", "sigma_reach", "sigma_site", "sigma_year", "log_y_missing", "log_y_shift_missing"),
  algorithm = "fullrank",
  adapt_iter = 1000,
  iter = 1000,
  seed = seed
)

# save fitted
qsave(draws_alt, file = "outputs/fitted/mat-normal-draws.qs")

# load and summarise fitted
draws_alt <- qread("outputs/fitted/mat-normal-draws.qs")
sum_alt <- summary(draws_alt)

cov_sp <- matrix(
  sum_alt$summary[grepl("Sigma_sp", rownames(sum_alt$summary)), "mean"],
  nrow = nsp,
  byrow = TRUE
)
rownames(cov_sp) <- colnames(cov_sp) <- sort(unique(size_abundance$scientific_name))
cov_class <- matrix(
  sum_alt$summary[grepl("Sigma_class", rownames(sum_alt$summary)), "mean"],
  nrow = nclass,
  byrow = TRUE
)
cor_sp <- cov2cor(cov_sp)
cor_class <- cov2cor(cov_class)

# think about variances and covariances and correlations
#   Size classes seem to have bigger absolute covars and corrs

# extra other effects
sum_alt$summary[grepl("sigma", rownames(sum_alt$summary)), ]
sum_alt$summary[grepl("alpha", rownames(sum_alt$summary)), ]
sum_alt$summary[grepl("beta", rownames(sum_alt$summary)), ]
sum_alt$summary[grepl("tau", rownames(sum_alt$summary)), ]
sum_alt$summary[grepl("rho", rownames(sum_alt$summary)), ]

# plot all against priors

# design main output plots

## TODO: decide if singledim or unstruc models are feasible or worthwhile
##   should be easy by setting one of the mats to I(nsp/nclass) or just writing
##   an unstruc version
