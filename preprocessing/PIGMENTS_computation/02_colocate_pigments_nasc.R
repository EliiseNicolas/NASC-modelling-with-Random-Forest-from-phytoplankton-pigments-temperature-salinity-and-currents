# ==============================================================================
# 02_colocate_pigments_nasc.R
#
# Colocation of the pigments with the acoustic ESU (NASC)
#
# Input  : pigment dataset of all years, already cropped to the FOD area and
#          filtered (01_create_ds_pigments_filtered_all_years.R)
#          <pigments_file>
#          NASC per ESU (03_compute_nasc.R)
#          <nasc_dir>/NASC_per_esu_<years_tag>_<freq>kHz.rds
#
# Steps  : for each frequency and each ESU, average every variable of the
#          pigment dataset (concentrations and ratios) over a window of
#          `pigments_window_size` x `pigments_window_size` pixels (config.R)
#          centred on the nearest pixel, on the same day
#
# Output : <pigments_nasc_dir>/
#            pigments_colocated_NASC_per_esu_<n>x<n>_<years_tag>_<freq>kHz.rds
#          (n = pigments_window_size), a data frame with one row per ESU:
#            time, lat_sv, lon_sv : date and position of the ESU
#            lat_pig, lon_pig     : matched pigment pixel
#            Chla, Per, ...       : mean concentration of each pigment
#            Chla_total           : Chla alone
#            <pigment>_totpig     : mean ratio of the pigment to the sum of the
#                                   pigments other than Chla
#          Pigment columns are NA if the day is missing from the pigment
#          dataset, NaN if the window has no valid value.
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # freqs, directories, pigments_file, nasc_file(),
                     # pigments_nasc_file(), pigments_window_size

in_file <- pigments_file
out_dir <- pigments_nasc_dir

# Width of the averaging window, in pixels (odd)
window_size <- pigments_window_size

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Indices of a window of `size` pixels centred on the pixel `index`, limited
# to the grid 1..n.
get_window <- function(index, n, size) {
  half <- floor(size / 2)
  max(1, index - half):min(n, index + half)
}


# ---- Pigment dataset ---------------------------------------------------------

pigments <- readRDS(in_file)

lon_pig   <- pigments$lon
lat_pig   <- pigments$lat
pig_dates <- format(pigments$date, "%Y-%m-%d")
n_lon     <- length(lon_pig)
n_lat     <- length(lat_pig)

# All the variables of the dataset (pigments and ratios), (date x lon x lat)
var_names <- setdiff(names(pigments), c("lon", "lat", "date"))

# Column names in the output: "c_cond_Chla" -> "Chla", ratios unchanged
col_names <- sub("^c_cond_", "", var_names)


# ---- Colocation, per frequency -----------------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  nasc      <- readRDS(nasc_file(freq))
  esu_dates <- format(as.Date(nasc$time), "%Y-%m-%d")

  colocated <- data.frame(
    time    = nasc$time,
    lat_sv  = nasc$lat,
    lon_sv  = nasc$lon,
    lat_pig = NA_real_,
    lon_pig = NA_real_
  )
  colocated[col_names] <- NA_real_

  for (date_i in unique(esu_dates)) {
    ind      <- which(esu_dates == date_i)   # ESU of this day
    idx_date <- match(date_i, pig_dates)     # this day in the pigment dataset

    if (is.na(idx_date)) {
      cat("Day not found in the pigment dataset:", date_i, "\n")
      next
    }

    # Nearest pigment pixel of each ESU
    idx_lon <- sapply(nasc$lon[ind], function(x) which.min(abs(lon_pig - x)))
    idx_lat <- sapply(nasc$lat[ind], function(x) which.min(abs(lat_pig - x)))

    colocated$lat_pig[ind] <- lat_pig[idx_lat]
    colocated$lon_pig[ind] <- lon_pig[idx_lon]

    for (v in seq_along(var_names)) {
      # Map (lon x lat) of the variable on this day
      map_day <- pigments[[var_names[v]]][idx_date, , ]

      # Mean over the window around each ESU
      colocated[ind, col_names[v]] <- sapply(seq_along(ind), function(j) {
        lon_win <- get_window(idx_lon[j], n_lon, window_size)
        lat_win <- get_window(idx_lat[j], n_lat, window_size)
        mean(map_day[lon_win, lat_win], na.rm = TRUE)
      })
    }
  }

  # Check: distance between each ESU and its matched pixel (degrees)
  print(summary(colocated$lat_pig - colocated$lat_sv))
  print(summary(colocated$lon_pig - colocated$lon_sv))
  str(colocated)

  out_file <- pigments_nasc_file(freq)
  saveRDS(colocated, out_file)
  cat("File saved:", out_file, "\n")
}
