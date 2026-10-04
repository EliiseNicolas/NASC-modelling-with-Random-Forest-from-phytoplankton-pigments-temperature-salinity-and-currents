# ==============================================================================
# 01_filter_crop_sv_profiles.R
#
# NASC computation, step 1/3: NA diagnostic and crop of the Sv profiles
#
# Input  : OBSAUSTRAL echo-integration NetCDF transects (2018, 2021, 2022, 2023)
#          Sv(depth, time, frequency) at 18, 38, 70, 120 and 200 kHz.
#
# Steps  : I)  NA diagnostic per year and frequency, used to choose the depth
#              below which Sv is mostly missing (-> `depth_max`).
#          II) Crop each year / frequency to the study area, the selected day
#              codes and the depth range ]depth_min, depth_max[.
#
# Output : <sv_per_year_dir>/<freq>kHz/Sv_<year>_<freq>kHz.rds, a list with
#            profiles : matrix (n_profiles x n_depth), Sv in dB
#            lat, lon : position of each profile
#            depth    : depth grid (m)
#            time     : POSIXct (UTC)
#            day      : day code of each profile
# ==============================================================================

library(ncdf4)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # years, freqs, study area, directories, sv_year_file()

in_dir   <- raw_acoustic_dir
in_files <- c(
  "2018" = "LOCEAN_SOOP-BA_A_20180105T121559Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20180201T103636Z_C-20260522T153853Z.nc",
  "2021" = "LOCEAN_SOOP-BA_A_20210122T143044Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20210307T041919Z_C-20260609T160140Z.nc",
  "2022" = "LOCEAN_SOOP-BA_A_20220201T075726Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20220305T041935Z_C-20260827T064134Z.nc",
  "2023" = "LOCEAN_SOOP-BA_A_20230123T103153Z_MARIONDUFRESNE_FV02_EchointegrationAcoustic-18-38-70-120-200_END-20230227T021804Z_C-20260728T105027Z.nc"
)
fig_dir  <- file.path(fig_root, "NASC_computation", "01_diag_NA")

stopifnot(all(years %in% names(in_files)))

# What to run
run_diagnostic <- TRUE   # part I (figures + console summary)
run_filtering  <- TRUE   # part II (writes the .rds files)

# Day codes kept (variable `day` of the NetCDF files)
day_codes <- c(1, 3)

# Depth range (m). Maximum depth per frequency, chosen from part I.
depth_min <- 25
depth_max <- c("18" = 700, "38" = 700, "70" = 500, "120" = 300, "200" = 150)

# % of NA left in ]depth_min, depth_max[ (part I, whole transect):
#
#   freq      depth_max   2018     2021      2023
#   18 kHz    700 m       4.03 %   17.79 %   4.14 %
#   38 kHz    700 m       2.94 %   17.66 %   2.07 %
#   70 kHz    500 m       3.28 %   14.08 %   1.51 %
#   120 kHz   300 m       3.96 %    9.47 %   1.20 %
#   200 kHz   150 m       3.69 %    3.83 %   0.69 %

# `time` is stored as days since this origin
time_origin <- as.POSIXct("1950-01-01", tz = "UTC")


# ---- Functions ---------------------------------------------------------------

# Read one frequency of a transect file.
# Returns a list: Sv (depth x time), depth, time (POSIXct), lat, lon, day.
read_sv <- function(path, freq) {
  ds <- nc_open(path)
  on.exit(nc_close(ds))

  idx_freq <- which(ncvar_get(ds, "instrument_frequency") == freq)
  stopifnot(length(idx_freq) == 1)

  list(
    # Sv is (depth, time, frequency): read only the requested frequency
    Sv    = ncvar_get(ds, "Sv", start = c(1, 1, idx_freq), count = c(-1, -1, 1)),
    depth = ncvar_get(ds, "depth"),
    time  = time_origin + ncvar_get(ds, "time") * 86400,
    lat   = ncvar_get(ds, "latitude"),
    lon   = ncvar_get(ds, "longitude"),
    day   = ncvar_get(ds, "day")
  )
}


# Part I: NA diagnostic for one year and one frequency, on the whole transect
# (no filtering on position or day code). Prints a summary and saves a figure.
diagnose_na <- function(path, year, freq, depth_min, depth_max, fig_dir) {
  d     <- read_sv(path, freq)
  is_na <- is.na(d$Sv)

  cat(sprintf("NA, all depths       : %d (%.2f %%)\n",
              sum(is_na), 100 * mean(is_na)))

  in_range   <- d$depth > depth_min & d$depth < depth_max
  is_na_crop <- is_na[in_range, , drop = FALSE]

  cat(sprintf("NA, %g-%g m       : %d (%.2f %%)\n",
              depth_min, depth_max, sum(is_na_crop), 100 * mean(is_na_crop)))
  cat("Depths with no data at all (m) :",
      d$depth[rowSums(is_na) == ncol(is_na)], "\n")
  cat("Profiles with no data at all   :",
      sum(colSums(is_na) == nrow(is_na)), "\n")

  fig_subdir <- file.path(fig_dir, paste0(freq, "kHz"))
  dir.create(fig_subdir, recursive = TRUE, showWarnings = FALSE)
  filename <- file.path(
    fig_subdir,
    paste0("diagnostic_NA_", year, "_", freq, "kHz_transect_all_dataset.png")
  )

  png(filename, width = 1600, height = 900, res = 150)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))

  plot(
    x = d$depth, y = 100 * rowMeans(is_na), type = "l",
    xlab = "Depth (m)", ylab = "% NA",
    main = "NA per depth", cex.main = 0.8
  )
  plot(
    x = d$time, y = colSums(is_na_crop), type = "l",
    xlab = "Time", ylab = "Number of NA",
    main = paste0("NA per profile\n(cropped ", depth_min, "-", depth_max, " m)"),
    cex.main = 0.8
  )
  mtext(
    paste0("NA diagnostic ", year, " transect dataset ", freq, " kHz"),
    outer = TRUE, cex = 1.5
  )
  dev.off()

  cat("Figure saved:", filename, "\n")
  invisible(filename)
}


# Part II: crop one year and one frequency in depth, position and day code.
# Bounds are excluded, in depth and in position.
crop_sv <- function(path, year, freq, depth_min, depth_max) {
  d <- read_sv(path, freq)

  depth_idx <- which(d$depth > depth_min & d$depth < depth_max)
  time_idx  <- which(
    d$lat > lat_min & d$lat < lat_max &
      d$lon > lon_min & d$lon < lon_max &
      d$day %in% day_codes
  )

  Sv    <- d$Sv[depth_idx, time_idx, drop = FALSE]
  depth <- d$depth[depth_idx]

  # 2022 only: drop the first depth level and round the depths, so that the
  # grid starts at 27 m and matches the other years (checked in script 02).
  if (year == "2022") {
    Sv    <- Sv[-1, , drop = FALSE]
    depth <- round(depth[-1])
  }

  list(
    profiles = t(Sv),   # (n_profiles x n_depth)
    lat      = d$lat[time_idx],
    lon      = d$lon[time_idx],
    depth    = depth,
    time     = d$time[time_idx],
    day      = d$day[time_idx]
  )
}


# ---- Part I: NA diagnostic ---------------------------------------------------

if (run_diagnostic) {
  for (freq in freqs) {
    for (year in years) {
      cat("\n--- NA diagnostic:", year, "-", freq, "kHz ---\n")
      diagnose_na(
        path      = file.path(in_dir, in_files[[year]]),
        year      = year,
        freq      = freq,
        depth_min = depth_min,
        depth_max = depth_max[[as.character(freq)]],
        fig_dir   = fig_dir
      )
    }
  }
}


# ---- Part II: filtering ------------------------------------------------------

if (run_filtering) {
  for (freq in freqs) {
    for (year in years) {
      sv_cropped <- crop_sv(
        path      = file.path(in_dir, in_files[[year]]),
        year      = year,
        freq      = freq,
        depth_min = depth_min,
        depth_max = depth_max[[as.character(freq)]]
      )

      cat(sprintf(
        "%s - %g kHz: %d profiles x %d depths (%g-%g m)\n",
        year, freq, nrow(sv_cropped$profiles), ncol(sv_cropped$profiles),
        min(sv_cropped$depth), max(sv_cropped$depth)
      ))

      out_file <- sv_year_file(year, freq)
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      saveRDS(sv_cropped, out_file)
    }
  }
}
