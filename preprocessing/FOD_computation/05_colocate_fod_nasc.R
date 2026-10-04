# ==============================================================================
# 05_colocate_fod_nasc.R
#
# Colocation of the FOD classes with the acoustic ESU (NASC)
#
# Input  : FOD classes on the grid (04_fod_clustering.R) and FOD grid
#          (02_fod_bspline.R), in <fod_dir>
#            cluster_transition_map_renamed.rds   (lon x lat x time)
#            lon.rds, lat.rds, time.rds
#          NASC per ESU (03_compute_nasc.R)
#          <nasc_dir>/NASC_per_esu_<years_tag>_<freq>kHz.rds
#
# Steps  : for each frequency and each ESU, take the FOD class (clusters and
#          transitions, final numbering) of the nearest grid point, on the
#          same day
#
# Output : <fod_nasc_dir>/
#            fod_colocated_NASC_per_esu_<years_tag>_<freq>kHz.rds
#          a data frame with one row per ESU:
#            time, lat_sv, lon_sv : date and position of the ESU
#            lat_fod, lon_fod     : matched FOD grid point
#            fod_cluster          : FOD class (NA if the day is missing from
#                                   the FOD or the grid point is masked)
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # freqs, directories, nasc_file(), fod_nasc_file()

out_dir <- fod_nasc_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---- FOD classes -------------------------------------------------------------

fod_map  <- readRDS(file.path(fod_dir, "cluster_transition_map_renamed.rds"))
lon_fod  <- as.vector(readRDS(file.path(fod_dir, "lon.rds")))
lat_fod  <- as.vector(readRDS(file.path(fod_dir, "lat.rds")))
time_fod <- readRDS(file.path(fod_dir, "time.rds"))

fod_dates <- format(as.Date(time_fod), "%Y-%m-%d")

stopifnot(
  all(dim(fod_map) == c(length(lon_fod), length(lat_fod), length(time_fod)))
)

# Classes present in the FOD maps (to compare with the colocated classes)
cat("Classes in the FOD maps:", sort(unique(as.vector(fod_map))), "\n")


# ---- Colocation, per frequency -----------------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  nasc      <- readRDS(nasc_file(freq))
  esu_dates <- format(as.Date(nasc$time), "%Y-%m-%d")

  colocated <- data.frame(
    time        = nasc$time,
    lat_sv      = nasc$lat,
    lon_sv      = nasc$lon,
    lat_fod     = NA_real_,
    lon_fod     = NA_real_,
    fod_cluster = NA_integer_
  )

  for (date_i in unique(esu_dates)) {
    ind      <- which(esu_dates == date_i)   # ESU of this day
    idx_date <- match(date_i, fod_dates)     # this day in the FOD maps

    if (is.na(idx_date)) {
      cat("Day not found in the FOD maps:", date_i, "\n")
      next
    }

    # Nearest FOD grid point of each ESU
    idx_lon <- sapply(nasc$lon[ind], function(x) which.min(abs(lon_fod - x)))
    idx_lat <- sapply(nasc$lat[ind], function(x) which.min(abs(lat_fod - x)))

    colocated$lat_fod[ind]     <- lat_fod[idx_lat]
    colocated$lon_fod[ind]     <- lon_fod[idx_lon]
    colocated$fod_cluster[ind] <- fod_map[cbind(idx_lon, idx_lat, idx_date)]
  }

  cat("Classes in the colocated data:",
      sort(unique(colocated$fod_cluster)), "\n")

  # Check: distance between each ESU and its matched grid point (degrees)
  print(summary(colocated$lat_fod - colocated$lat_sv))
  print(summary(colocated$lon_fod - colocated$lon_sv))
  str(colocated)

  out_file <- fod_nasc_file(freq)
  saveRDS(colocated, out_file)
  cat("File saved:", out_file, "\n")
}
