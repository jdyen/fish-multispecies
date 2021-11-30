# functions to tidy fish and flow data

# function to clean species names
clean_scinames <- function(x) {
  
  # clean up common inconsistencies
  x <- tolower(x)
  x <- gsub("\\.", "", x)
  x <- gsub(" ", "_", x)
  x <- gsub("spp", "sp", x)
  
  # and rename some species to simplify
  x[grepl("galaxias", x)] <- "galaxias_sp"
  x[grepl("philypnodon", x)] <- "philypnodon_sp"
  x[grepl("craterocephalus", x)] <- "craterocephalus_sp"
  x[grepl("hypseleotris", x)] <- "hypseleotris_sp"
  
  # return
  x
  
}

# function to clean waterbody names
clean_waterbody <- function(x) {
  
  # tweak formatting
  x <- tolower(x)
  x <- gsub(" ", "_", x)
  x <- gsub("_ck$", "_creek", x)
  
  # return
  x
  
}

# function to extract reach info from notes column
add_reach <- function(x) {
  
  # extract reach id from notes (a string)
  x <- substr(x, 8, 8)
  x <- gsub(";", NA, x)
  
  # return numeric version of this
  as.numeric(x)
  
}

# function to check reaches and add missing reach info for Wimmera sites
check_reach <- function(x, system, site) {
  
  # check mackenzie sites (1-5 are in reach 2, reach 3 otherwise)
  idx <- system == "mackenzie"
  x[idx] <- ifelse(site[idx] %in% paste0("MR_", 1:5), 2, 3)
  
  # and burnt ck sites (1-5 and 11 are in reach 1, 2 otherwise)
  idx <- system == "burnt_creek"
  x[idx] <- ifelse(site[idx] %in% paste0("BC_", 6:10), 2, 1)
  
  # wimmera sites
  idx <- system == "wimmera"
  x[idx] <- ifelse(is.na(x[idx]), 2, x[idx])
  
  # return
  x
  
}

# define a lookup table and function to add region (catchment)
#   info based on waterbody
region_lookup <- c(
  "goulburn" = 405,
  "wimmera" = 415,
  "yarra" = 229,
  "loddon" = 407,
  "campaspe" = 406,
  "broken" = 404,
  "glenelg" = 238,
  "cardinia_creek" = 228,
  "barwon" = 233,
  "werribee" = 231,
  "tarwin" = 227, 
  "thomson" = 225,
  "pyramid_creek" = 407,
  "little_murray" = 407,
  "taylors_creek" = 230,
  "mackenzie" = 415,
  "burnt_creek" = 415,
  "mt_william_creek" = 415,
  "moorabool" = 232,
  "macalister" = 225 
)
add_region <- function(x) {
  region_lookup[x]
}

# function to find nearest station based on output nearest gauges from spatial database
find_nearest <- function(x, y) {
  
  best <- y %>% filter(grepl(!!x$waterbody, stname, ignore.case = TRUE))
  if (nrow(best) == 0)
    best <- y
  
  best <- best %>% arrange(min_distance_m)
  
  best <- best[1, ]
  
  as.data.frame(c(x, best))
  
}

# gauges
gauge_lookup <- c(
  "broken_r2" = 404216,            # Casey's weir head gauge, backup: 404241 Casey's Weir
  "broken_r3" = 404224,            # Gowangardie
  "broken_r4" = 404224,            # Gowangardie
  "broken_r5" = 404214,            # Katamite
  "burnt_creek_r1" = 415223,       # Wonwondah
  "burnt_creek_r2" = 415223,       # Wonwondah
  "campaspe_r2" = 406201,          # Barnadown
  "campaspe_r3" = 406202,          # Rochester
  "campaspe_r4" = 406202,          # Rochester, closest missing data (406276 Fehrings)
  "goulburn_r4" = 405200,          # Murchison, backup: 405276 Loch Garry
  "goulburn_r5" = 405232,          # McCoys
  "little_murray_rNA" = 409399,    # Little Murray Weir
  "loddon_r2" = 407229,            # Serpentine Weir, backup: 407240 Laanecoorie head gauge
  "loddon_r3" = 407229,            # Serpentine Weir, backup: 407240 laanecoorie head gauge  
  "loddon_r4" = 407205,            # Appin South, backup: 407323 Yando Rd
  "loddon_r5" = 407202,            # Kerang
  "mackenzie_r2" = 415202,         # Wartook Reservoir, backup: 415251 Mck Ck
  "mackenzie_r3" = 415251,         # McK river at McK Ck
  "mt_william_creek_rNA" = 415203, # Lake Lonsdale tail gauge
  "pyramid_creek_rNA" = 407295,    # Box Ck @ Mansfields Bridge, backup: 407294 Flannery's Bridge 
  "wimmera_r2" = 415200,           # Horsham, backup 415261 Quantong
  "wimmera_r3" = 415256            # U/S Dimboola, backup: 415261 Quantong
)

# function to check and replace flow gauges with appropriate ones
add_gauge <- function(x) {
  
  x %>% mutate(
    sysreach = paste(waterbody, reach, sep = "_r"),
    best_station = gauge_lookup[sysreach]
  ) %>% select(-sysreach)
  
}

# function to rebase characters or factors as integers
rebase_factor <- function(x) {
  as.integer(as.factor(x))
}

# function to calculate flow metrics
get_metrics <- function(x, y) {
  
  # work out relevant gauge for each reach
  reaches <- paste(y$waterbody, y$reach, sep = "_r")
  y <- y %>% mutate(gauge = gauge_lookup[reaches])
  
  # go through each gauge and survey date and grab relevant metrics
  gauge_date <- y %>% distinct(gauge, sdate)
  out <- list(length = nrow(gauge_date))
  for (i in seq_len(nrow(gauge_date))) {
    xsub <- x[[as.character(gauge_date$gauge[i])]]
    out[[i]] <- get_metrics_internal(xsub, gauge_date$sdate[i])
  }
  out <- do.call(rbind, out)
  
  # standardise all
  out_std <- sweep(out, 2, colMeans(out), "-")
  out_std <- sweep(out_std, 2, apply(out, 2, sd), "/")
  colnames(out_std) <- paste0(colnames(out), "_std")
  out <- cbind(out, out_std)
  
  # link with gauge and date info
  out <- cbind(gauge_date, out)
  
  # add metrics back into fish data, drop gauge column
  #  and return
  y %>% 
    select(waterbody, reach, site_name, id_survey, sdate, effort, gauge) %>%
    left_join(out, by = c("gauge", "sdate")) %>%
    select(-gauge)
  
}
  
# function to calculate specific flow metrics
get_metrics_internal <- function(x, date) {
  
  # calculate 10th percentile for low flow calcs
  q10 <- quantile(x$value, probs = 0.1, na.rm = TRUE)
  
  # work out start and end dates for each season
  start_year <- year(date)
  end_year <- year(date)
  idx <- month(date) %in% c(1:6)
  start_year[idx] <- start_year - 1L
  end_year[!idx] <- end_year + 1L
  
  out <- data.frame(
    ave_spring = calculate(
      x$value,
      x$date_formatted,
      resolution = survey(
        season = 9:11,
        start = dmy(paste("01-09-", start_year)),
        end = dmy(paste("30-11-", start_year))
      ),
      standardise = by_mean(range(year(x$date_formatted)))
    )$metric,
    ave_summer = calculate(
      x$value,
      x$date_formatted,
      resolution = survey(
        season = 12:15,
        start = dmy(paste("01-12-", start_year)),
        end = dmy(paste("31-03-", end_year))
      ),
      standardise = by_mean(range(year(x$date_formatted)))
    )$metric,
    ave_antecedent = calculate(
      x$value,
      x$date_formatted,
      resolution = survey(
        season = 7:18,
        lag = 1,
        start = dmy(paste("01-12-", start_year)),
        end = dmy(paste("31-03-", end_year))
      ),
      standardise = by_mean(range(year(x$date_formatted)))
    )$metric,
    low_flow = calculate(
      x$value,
      x$date_formatted,
      resolution = survey(
        season = 7:18,
        start = dmy(paste("01-07-", start_year)),
        end = dmy(paste("30-06-", end_year))
      ),
      fun = days_below,
      threshold = q10
    )$metric,
    cv_flow = calculate(
      x$value,
      x$date_formatted,
      resolution = survey(
        season = 7:18,
        start = dmy(paste("01-07-", start_year)),
        end = dmy(paste("30-06-", end_year))
      ),
      fun = cv_fun
    )$metric
  )
  
  # return
  out

}
  
# function to calculate coef of variation accounting for 0 SD variables
cv_fun <- function(x, ...) {
  out <- sd(x, ...) / mean(x, ...)
  if (anyNA(out))
    out[is.na(out)] <- sd(x, ...)
  out
}