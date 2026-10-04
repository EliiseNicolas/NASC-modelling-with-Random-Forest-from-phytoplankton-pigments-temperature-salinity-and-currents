# ==============================================================================
# 01_data_prep.R
#
# Random-forest pipeline: loading and cleaning of the dataset
#
# Learning dataset (01_build_learning_dataset.R; learning_dataset_file() in
# config.R), one row per ESU: time_nasc, lat_nasc, lon_nasc, nasc, fod
# (character, the string "NA" when missing), ftle, Chla_total and
# <pigment>_totpig. Only the ESU with day == `day_code` (config.R) are kept,
# as in the data-validation scripts.
#
# Defines : load_and_clean()
# ==============================================================================

# Load the dataset of one frequency and keep the daytime, complete rows.
# The NASC is converted to log10 here (it is not in the source file).
#   freq : frequency (kHz)
# Returns a list:
#   df         : tibble with time, lat, lon, year, NASC (log10), fod (factor),
#                the numeric covariates, and x_km, y_km (km from the mean
#                position, used by the spatial blocking)
#   fod_levels : levels of fod
load_and_clean <- function(freq) {
  in_file <- learning_dataset_file(freq)
  if (!file.exists(in_file)) {
    stop("File not found: ", in_file,
         " -- run 01_build_learning_dataset.R first")
  }
  ds <- readRDS(in_file)

  n_total <- nrow(ds)
  ds      <- ds[ds$day == day_code, ]
  n_day   <- nrow(ds)

  df <- tibble(
    time = ds$time_nasc,
    lat  = ds$lat_nasc,
    lon  = ds$lon_nasc,
    year = lubridate::year(ds$time_nasc),
    NASC = log10(pmax(ds[[response_raw]], nasc_log_floor)),
    fod  = factor(ifelse(ds$fod == "NA", NA_character_, ds$fod))
  )
  df[covariates_num] <- ds[covariates_num]

  df <- df[stats::complete.cases(df[, covariates_all]) & is.finite(df$NASC), ]

  # Coordinates in km, centred on the mean position (local flat approximation)
  km_per_deg_lat <- 111.32
  km_per_deg_lon <- 111.32 * cos(mean(df$lat, na.rm = TRUE) * pi / 180)
  df <- df %>% mutate(
    x_km = (lon - mean(lon)) * km_per_deg_lon,
    y_km = (lat - mean(lat)) * km_per_deg_lat
  )

  cat(sprintf(
    "[%d kHz] %d ESU -> %d with day == %d -> %d kept (complete rows)\n",
    freq, n_total, n_day, day_code, nrow(df)
  ))
  cat(sprintf("  Years: %s\n", paste(sort(unique(df$year)), collapse = ", ")))

  list(df = df, fod_levels = levels(df$fod))
}
