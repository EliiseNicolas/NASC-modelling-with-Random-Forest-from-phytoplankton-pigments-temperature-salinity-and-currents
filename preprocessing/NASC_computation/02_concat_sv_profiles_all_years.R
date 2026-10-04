# ==============================================================================
# 02_concat_sv_profiles_all_years.R
#
# NASC computation, step 2/3: concatenation of the years
#
# Input  : cropped Sv profiles per year and frequency
#          (01_filter_crop_sv_profiles.R)
#          <sv_per_year_dir>/<freq>kHz/Sv_<year>_<freq>kHz.rds
#
# Steps  : for each frequency, check that all years share the same depth grid,
#          then stack the profiles of all years.
#
# Output : <sv_all_years_dir>/Sv_<years_tag>_<freq>kHz.rds, a list with
#            profiles : matrix (n_profiles x n_depth), Sv in dB
#            lat, lon, time, day : one value per profile
#            depth    : common depth grid (m)
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # years, freqs, sv_year_file(), sv_all_years_file()


# ---- Concatenation -----------------------------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  # Read the yearly files
  sv_years <- lapply(years, function(year) {
    readRDS(sv_year_file(year, freq))
  })
  names(sv_years) <- years

  for (year in years) {
    cat(year, ":", nrow(sv_years[[year]]$profiles), "profiles\n")
  }

  # The depth grid must be identical for all years
  depth <- sv_years[[1]]$depth
  for (year in years) {
    if (!identical(sv_years[[year]]$depth, depth)) {
      stop("Depth grid of ", year, " differs from ", years[1],
           " at ", freq, " kHz")
    }
  }

  # Stack profiles and metadata
  profiles <- do.call(rbind, lapply(unname(sv_years), `[[`, "profiles"))
  metadata <- do.call(rbind, lapply(unname(sv_years), function(sv) {
    data.frame(lat = sv$lat, lon = sv$lon, time = sv$time, day = sv$day)
  }))
  stopifnot(nrow(profiles) == nrow(metadata))

  sv_concat <- list(
    profiles = profiles,
    lat      = metadata$lat,
    lon      = metadata$lon,
    time     = metadata$time,
    day      = metadata$day,
    depth    = depth
  )
  str(sv_concat)

  out_file <- sv_all_years_file(freq)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sv_concat, out_file)
}
