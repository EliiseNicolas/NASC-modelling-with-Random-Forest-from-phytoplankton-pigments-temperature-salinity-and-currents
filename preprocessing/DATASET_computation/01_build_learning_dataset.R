# ==============================================================================
# 01_build_learning_dataset.R
#
# Learning dataset: NASC and colocated covariates per ESU
#
# Input  : per frequency, four tables with the same ESU in the same order
#          NASC per ESU (03_compute_nasc.R)
#          <nasc_dir>/NASC_per_esu_<years_tag>_<freq>kHz.rds
#          FOD classes colocated with the ESU (05_colocate_fod_nasc.R)
#          <fod_nasc_dir>/fod_colocated_NASC_per_esu_<years_tag>_<freq>kHz.rds
#          pigments colocated with the ESU (02_colocate_pigments_nasc.R)
#          <pigments_nasc_dir>/
#            pigments_colocated_NASC_per_esu_<n>x<n>_<years_tag>_<freq>kHz.rds
#          FTLE colocated with the ESU (02_colocate_ftle_nasc.R)
#          <ftle_nasc_dir>/ftle_colocated_NASC_per_esu_<years_tag>_<freq>kHz.rds
#
# Steps  : for each frequency
#          1) check that the four tables hold the same ESU in the same order
#             (the script stops otherwise)
#          2) gather the NASC, the FOD class, the pigments and the FTLE in one
#             data frame
#          3) checks: years, distance between each ESU and its matched pixels
#
# Output : <learning_dataset_dir>/learning_dataset_<years_tag>_<freq>kHz.rds,
#          a data frame with one row per ESU:
#            time_nasc, lat_nasc, lon_nasc : date and position of the ESU
#            day                           : 3 = day, 1 = night
#            nasc                          : NASC
#            lat_fod, lon_fod, fod         : matched FOD grid point and FOD
#                                            class (text, see below)
#            lat_pig, lon_pig, ...         : matched pigment pixel, then all
#                                            the variables of the pigment table
#                                            (Chla, ..., Chla_total,
#                                            <pigment>_totpig)
#            lat_ftle, lon_ftle, ftle      : matched FTLE pixel and FTLE
#
# The FOD class is stored as text of fixed width (" 1", ..., "13", and "NA"
# when missing): as a factor, its levels then sort in the numeric order of
# the classes.
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # freqs, nasc_file(), fod_nasc_file(), ftle_nasc_file(),
                     # pigments_nasc_file(), learning_dataset_file()

# Columns identifying the ESU in the colocated tables
esu_cols <- c("time", "lat_sv", "lon_sv")


# ---- Functions ---------------------------------------------------------------

# Stop if the rows of a colocated table are not the ESU of the NASC table, in
# the same order.
#   nasc      : NASC table (time, lat, lon)
#   colocated : colocated table (time, lat_sv, lon_sv)
#   name      : name of the colocated table, for the error message
check_same_esu <- function(nasc, colocated, name) {
  same <- nrow(colocated) == nrow(nasc) &&
    all(colocated$time == nasc$time) &&
    all(colocated$lat_sv == nasc$lat) &&
    all(colocated$lon_sv == nasc$lon)
  if (!isTRUE(same)) {
    stop("The ", name, " table does not hold the same ESU as the NASC table.")
  }
}


# ---- Learning dataset, per frequency -----------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  nasc <- readRDS(nasc_file(freq))
  fod  <- readRDS(fod_nasc_file(freq))
  pig  <- readRDS(pigments_nasc_file(freq))
  ftle <- readRDS(ftle_nasc_file(freq))

  # ---- 1) Same ESU in the four tables ----

  check_same_esu(nasc, fod, "FOD")
  check_same_esu(nasc, pig, "pigment")
  check_same_esu(nasc, ftle, "FTLE")

  # ---- 2) Assembly ----

  ds <- data.frame(
    time_nasc = nasc$time,
    lat_nasc  = nasc$lat,
    lon_nasc  = nasc$lon,
    day       = nasc$day,
    nasc      = nasc$NASC
  )

  # FOD: matched grid point and class
  ds$lat_fod <- fod$lat_fod
  ds$lon_fod <- fod$lon_fod
  ds$fod     <- format(fod$fod_cluster)

  # Pigments: matched pixel and all the pigment variables
  pig_vars <- setdiff(names(pig), esu_cols)
  ds[pig_vars] <- pig[pig_vars]

  # FTLE: matched pixel and value
  ds$lat_ftle <- ftle$lat_ftle
  ds$lon_ftle <- ftle$lon_ftle
  ds$ftle     <- ftle$ftle

  # ---- 3) Checks ----

  str(ds)
  cat("Years:", unique(format(ds$time_nasc, "%Y")), "\n")
  cat("Days :", length(unique(as.Date(ds$time_nasc))), "\n")
  cat("ESU with a FOD grid point:", sum(!is.na(ds$lat_fod)), "/", nrow(ds),
      "\n")

  # Distance between each ESU and its matched pixels (degrees)
  for (covariate in c("fod", "pig", "ftle")) {
    cat("Matched", covariate, "pixel - ESU, latitude then longitude:\n")
    print(summary(ds[[paste0("lat_", covariate)]] - ds$lat_nasc))
    print(summary(ds[[paste0("lon_", covariate)]] - ds$lon_nasc))
  }

  # ---- Save ----

  out_file <- learning_dataset_file(freq)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ds, out_file)
  cat("File saved:", out_file, "\n")
}
