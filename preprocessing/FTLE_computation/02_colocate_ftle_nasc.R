# ==============================================================================
# 02_colocate_ftle_nasc.R
#
# Colocation of the FTLE with the acoustic ESU (NASC)
#
# Input  : daily FTLE NetCDF files in <raw_ftle_dir>, named
#          map_<YYYY-MM-DD>*.nc, variable FTLE (lon, lat, time) with a single
#          time step
#          NASC per ESU (03_compute_nasc.R)
#          <nasc_dir>/NASC_per_esu_<years_tag>_<freq>kHz.rds
#
# Steps  : for each frequency and each ESU, take the FTLE of the nearest pixel
#          in the FTLE map of the same day
#
# Output : <ftle_nasc_dir>/
#            ftle_colocated_NASC_per_esu_<years_tag>_<freq>kHz.rds
#          a data frame with one row per ESU:
#            time, lat_sv, lon_sv : date and position of the ESU
#            lat_ftle, lon_ftle   : matched FTLE pixel
#            ftle                 : FTLE value (NA if there is no map on that
#                                   day or no data in the pixel)
# ==============================================================================

library(ncdf4)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # years_tag, freqs, directories, nasc_file()

in_dir  <- raw_ftle_dir
out_dir <- ftle_nasc_dir

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---- FTLE files --------------------------------------------------------------

files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)
if (length(files) == 0) stop("No FTLE file in ", in_dir)


# ---- Colocation, per frequency -----------------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  nasc      <- readRDS(nasc_file(freq))
  esu_dates <- format(as.Date(nasc$time), "%Y-%m-%d")

  colocated <- data.frame(
    time     = nasc$time,
    lat_sv   = nasc$lat,
    lon_sv   = nasc$lon,
    lat_ftle = NA_real_,
    lon_ftle = NA_real_,
    ftle     = NA_real_
  )

  for (date_i in unique(esu_dates)) {
    ind <- which(esu_dates == date_i)   # ESU of this day

    # FTLE map of this day
    file_i <- files[grepl(paste0("map_", date_i, ".*\\.nc$"), basename(files))]

    if (length(file_i) == 0) {
      cat("Day not found in the FTLE files:", date_i, "\n")
      next
    }

    ds       <- nc_open(file_i[1])
    lon_ftle <- ds$dim$lon$vals
    lat_ftle <- ds$dim$lat$vals
    ftle_map <- ncvar_get(ds, "FTLE")   # (lon x lat)
    nc_close(ds)

    # Nearest FTLE pixel of each ESU
    idx_lon <- sapply(nasc$lon[ind], function(x) which.min(abs(lon_ftle - x)))
    idx_lat <- sapply(nasc$lat[ind], function(x) which.min(abs(lat_ftle - x)))

    colocated$lat_ftle[ind] <- lat_ftle[idx_lat]
    colocated$lon_ftle[ind] <- lon_ftle[idx_lon]
    colocated$ftle[ind]     <- ftle_map[cbind(idx_lon, idx_lat)]
  }

  # Check: distance between each ESU and its matched pixel (degrees)
  print(summary(colocated$lat_ftle - colocated$lat_sv))
  print(summary(colocated$lon_ftle - colocated$lon_sv))
  str(colocated)

  out_file <- file.path(
    out_dir,
    paste0("ftle_colocated_NASC_per_esu_", years_tag, "_", freq, "kHz.rds")
  )
  saveRDS(colocated, out_file)
  cat("File saved:", out_file, "\n")
}
