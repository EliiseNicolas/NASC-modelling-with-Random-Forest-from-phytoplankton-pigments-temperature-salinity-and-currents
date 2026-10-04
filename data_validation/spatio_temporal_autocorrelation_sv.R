# ==============================================================================
# spatio_temporal_autocorrelation_sv.R
#
# Data validation: spatio-temporal dependence of Sv
#
# Input  : Sv profiles of all years, per frequency
#          (02_concat_sv_profiles_all_years.R)
#          <sv_all_years_dir>/Sv_<years_tag>_<freq>kHz.rds, a list with
#          profiles (pings x depths, dB), depth, lon, lat, time and day
#          (3 = day, 1 = night). One ping = one Sv profile.
#
# Steps  : for each frequency, and for all the pings, the daytime pings then
#          the night-time pings
#          1) mean Sv of each ping over depth
#          2) vertical coverage: number of valid depth bins per ping
#          3) omnidirectional variogram
#          4) omnidirectional spatial correlogram
#          5) temporal correlogram (by pairs: the time step is irregular)
#
# Output : in <fig_dir>
#            diag_autocorrelation_spatio_temp_Sv_<freq>kHz_<all|day|night>.png
# ==============================================================================

library(gstat)
library(ggplot2)
library(dplyr)
library(patchwork)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # freqs, directories, sv_all_years_file()

fig_dir <- file.path(fig_root_validation, "spatio_temporal_autocorrelation")

# Subsets of pings analysed: values of `day` kept (NULL = all the pings)
day_subsets <- list(all = NULL, day = 3, night = 1)

# Layer analysed: c(min, max) in m, e.g. c(50, 200); NULL = whole water column
depth_range <- NULL

min_depth_obs <- 5         # minimum number of valid depth bins to keep a ping
n_sample      <- 3000      # maximum number of pings kept (random sample)
time_unit     <- "hours"   # unit of the time lags
seed          <- 123

# Directions of the directional variogram and correlogram, with one colour
# each (same labels and colours in every figure)
direction_levels <- c("Meridional (N-S)", "Diagonal NE-SW", "Zonal (E-W)",
                      "Diagonal NW-SE")
direction_colors <- c(
  "Meridional (N-S)" = "#1b9e77",
  "Diagonal NE-SW"   = "#7570b3",
  "Zonal (E-W)"      = "#d95f02",
  "Diagonal NW-SE"   = "#e7298a"
)

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions: helpers ------------------------------------------------------

# Add coordinates in km (local flat approximation): x_km (longitude scaled at
# the mean latitude) and y_km. Valid over the extent of the transects; a
# proper projection (package sf) would be better over a wider latitude range.
lonlat_to_km <- function(dat) {
  lat0 <- mean(dat$lat, na.rm = TRUE)
  dat %>%
    mutate(x_km = lon * 111.32 * cos(lat0 * pi / 180), y_km = lat * 111.32)
}

# Default maximum distance of the variograms: half the diagonal of the area.
default_cutoff <- function(dat) {
  dx <- diff(range(dat$x_km))
  dy <- diff(range(dat$y_km))
  0.5 * sqrt(dx^2 + dy^2)
}

# Placeholder panel, when a diagnostic cannot be computed.
empty_plot <- function(message) {
  patchwork::wrap_elements(grid::textGrob(message))
}

# Draw `n_pairs` random pairs among `n` observations (with replacement, fixed
# seed). Returns the indices i and j of the two members of each pair.
sample_pair_indices <- function(n, n_pairs) {
  set.seed(seed)
  i <- sample(seq_len(n), n_pairs, replace = TRUE)
  j <- sample(seq_len(n), n_pairs, replace = TRUE)
  list(i = i, j = j)
}

# Random pairs of observations of `dat` (x_km, y_km, value).
# Returns a data frame with one row per valid pair (non-zero distance, two
# finite values): distance (km), angle (degrees in [0, 180), 0 = East-West,
# 90 = North-South) and values z1, z2 of the two members.
sample_pairs <- function(dat, n_pairs) {
  idx <- sample_pair_indices(nrow(dat), n_pairs)

  dx <- dat$x_km[idx$i] - dat$x_km[idx$j]
  dy <- dat$y_km[idx$i] - dat$y_km[idx$j]
  distance <- sqrt(dx^2 + dy^2)
  angle    <- (atan2(dy, dx) * 180 / pi) %% 180

  z1 <- dat$value[idx$i]
  z2 <- dat$value[idx$j]

  keep <- distance > 0 & is.finite(z1) & is.finite(z2)
  data.frame(distance = distance[keep], angle = angle[keep],
             z1 = z1[keep], z2 = z2[keep])
}

# Direction class of an angle (degrees in [0, 180), see sample_pairs()).
angle_to_direction <- function(angle) {
  direction <- cut(
    angle, breaks = c(-Inf, 22.5, 67.5, 112.5, 157.5, Inf),
    labels = c("Zonal (E-W)", "Diagonal NE-SW", "Meridional (N-S)",
               "Diagonal NW-SE", "Zonal (E-W)")
  )
  factor(direction, levels = direction_levels)
}


# ---- Functions: ping data ----------------------------------------------------

# Mean Sv of each ping over depth.
# Sv is logarithmic: the mean is taken in the linear domain (10^(Sv/10)) and
# converted back to dB, otherwise it is biased. Pings are processed by chunks
# to limit the memory peak (the matrix has ~230 million cells).
#   profiles_mat : Sv matrix (pings x depths, dB)
#   chunk_size   : number of pings per chunk
# Returns a list: value (mean Sv per ping, dB; NaN without any valid bin) and
# n_valid (number of valid depth bins per ping).
mean_sv_by_ping <- function(profiles_mat, chunk_size = 50000) {
  n <- nrow(profiles_mat)
  out_mean   <- numeric(n)
  out_nvalid <- integer(n)

  for (start in seq(1, n, by = chunk_size)) {
    end   <- min(start + chunk_size - 1, n)
    block <- profiles_mat[start:end, , drop = FALSE]
    out_mean[start:end]   <- 10 * log10(rowMeans(10^(block / 10), na.rm = TRUE))
    out_nvalid[start:end] <- rowSums(!is.na(block))
  }
  list(value = out_mean, n_valid = out_nvalid)
}

# Mean Sv per ping, with the vertical coverage.
#   profiles      : Sv dataset (list: profiles, depth, lon, lat, time, day)
#   depth_range   : c(min, max) to analyse one layer; NULL = whole water column
#   min_depth_obs : minimum number of valid depth bins to keep a ping (avoids
#                   means computed on 1-2 values)
#   n_sample      : maximum number of pings kept (random sample, fixed seed)
#   day_value     : values of profiles$day kept; NULL = all the pings
# Returns a data frame with one row per ping kept: lon, lat, time, value (mean
# Sv, dB), n_valid, frac_valid. The attribute "coverage" summarises the effect
# of the NA (share of NA, pings in total / kept / sampled).
make_ping_data <- function(profiles, depth_range = NULL, min_depth_obs = 5,
                           n_sample = 3000, day_value = NULL) {
  # Day / night filter first, to keep the matrix and the metadata aligned
  # (without filter, the matrix is not copied)
  if (!is.null(day_value)) {
    keep_rows <- profiles$day %in% day_value
    mat  <- profiles$profiles[keep_rows, , drop = FALSE]
    lon  <- profiles$lon[keep_rows]
    lat  <- profiles$lat[keep_rows]
    time <- profiles$time[keep_rows]
  } else {
    mat  <- profiles$profiles
    lon  <- profiles$lon
    lat  <- profiles$lat
    time <- profiles$time
  }

  if (!is.null(depth_range)) {
    keep_cols <- profiles$depth >= depth_range[1] &
      profiles$depth <= depth_range[2]
    mat <- mat[, keep_cols, drop = FALSE]
  }
  n_depth <- ncol(mat)

  sv <- mean_sv_by_ping(mat)

  dat <- data.frame(
    lon = lon, lat = lat, time = time,
    value = sv$value, n_valid = sv$n_valid, frac_valid = sv$n_valid / n_depth
  )

  pct_na_global <- 1 - sum(dat$n_valid) / (nrow(dat) * n_depth)
  n_pings_total <- nrow(dat)

  dat <- dat %>% filter(is.finite(value), n_valid >= min_depth_obs)
  n_pings_kept <- nrow(dat)

  if (nrow(dat) > n_sample) {
    set.seed(seed)
    dat <- dat %>% slice_sample(n = n_sample)
  }

  attr(dat, "coverage") <- data.frame(
    n_depth_bins = n_depth, pct_na_global = pct_na_global,
    n_pings_total = n_pings_total, n_pings_kept = n_pings_kept,
    n_pings_sampled = nrow(dat)
  )
  dat
}


# ---- Functions: figures ------------------------------------------------------

# Map of the number of valid depth bins per ping.
#   dat   : output of make_ping_data()
#   label : name of the variable, for the title
make_coverage_plot <- function(dat, label) {
  cov <- attr(dat, "coverage")
  subtitle <- sprintf(
    "%.0f%% NA (depth) | %d/%d pings kept (min_depth_obs) | %d sampled",
    100 * cov$pct_na_global, cov$n_pings_kept, cov$n_pings_total,
    cov$n_pings_sampled
  )
  ggplot(dat, aes(x = lon, y = lat, color = n_valid)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_viridis_c(name = "Valid\ndepth bins") +
    labs(title = paste("Vertical coverage -", label), subtitle = subtitle,
         x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# Variogram of the mean Sv per ping, omnidirectional or directional.
# Caution: a transect mostly samples the direction of the ship track. The
# directions the ship rarely follows have few pairs and must be read with
# care, unlike on a regular grid covering the whole area.
#   cutoff      : maximum distance (km); NULL = default_cutoff()
#   alpha       : directions (degrees clockwise from North, gstat convention:
#                 0 = North-South, 90 = East-West) among 0, 45, 90, 135;
#                 ignored if directional = FALSE
#   directional : TRUE = one curve per direction; FALSE = a single curve, all
#                 directions pooled
# Returns a ggplot (a placeholder panel with fewer than 30 pings).
make_variogram <- function(dat, label, cutoff = NULL,
                           alpha = c(0, 45, 90, 135), directional = TRUE) {
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(empty_plot(paste("Not enough data -", label)))

  if (is.null(cutoff)) cutoff <- default_cutoff(dat)

  if (!directional) {
    v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat,
                   cutoff = cutoff, width = cutoff / 15)
    return(
      ggplot(v, aes(x = dist, y = gamma)) +
        geom_point(size = 1.8) +
        geom_line() +
        labs(title = paste("Omnidirectional variogram -", label),
             x = "Distance (km)", y = "Semivariance") +
        theme_minimal()
    )
  }

  v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat,
                 cutoff = cutoff, width = cutoff / 15, alpha = alpha)

  alpha_labels <- c("0" = "Meridional (N-S)", "45" = "Diagonal NE-SW",
                    "90" = "Zonal (E-W)", "135" = "Diagonal NW-SE")
  v$direction <- factor(alpha_labels[as.character(v$dir.hor)],
                        levels = direction_levels)

  ggplot(v, aes(x = dist, y = gamma, color = direction)) +
    geom_point(size = 1.6) +
    geom_line() +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    labs(title = paste("Directional variogram -", label),
         x = "Distance (km)", y = "Semivariance", color = "Direction") +
    theme_minimal()
}

# Spatial correlogram of the mean Sv per ping, omnidirectional or directional:
# correlation between the two members of random pairs of pings, per distance
# class.
#   n_pairs           : number of pairs drawn
#   n_bins            : number of distance classes
#   min_pairs_per_bin : minimum number of pairs to keep a class
#   directional       : TRUE = one curve per direction; FALSE = a single
#                       curve, all directions pooled
# Returns a ggplot (a placeholder panel with fewer than 30 pings).
make_correlogram_spatial <- function(dat, label, n_pairs = 50000, n_bins = 20,
                                     min_pairs_per_bin = 30,
                                     directional = TRUE) {
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(empty_plot(paste("Not enough data -", label)))

  pairs <- sample_pairs(dat, n_pairs)
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)

  if (!directional) {
    cor_df <- pairs %>%
      group_by(bin) %>%
      summarise(distance = mean(distance), correlation = cor(z1, z2),
                n = n(), .groups = "drop") %>%
      filter(n >= min_pairs_per_bin)

    return(
      ggplot(cor_df, aes(x = distance, y = correlation)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        geom_point(size = 1.8) +
        geom_line() +
        labs(title = paste("Omnidirectional spatial correlogram -", label),
             x = "Distance (km)", y = "Correlation") +
        ylim(-1, 1) +
        theme_minimal()
    )
  }

  pairs$direction <- angle_to_direction(pairs$angle)

  cor_df <- pairs %>%
    group_by(direction, bin) %>%
    summarise(distance = mean(distance), correlation = cor(z1, z2),
              n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)

  ggplot(cor_df, aes(x = distance, y = correlation, color = direction)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.6) +
    geom_line() +
    scale_color_manual(values = direction_colors, drop = FALSE) +
    labs(title = paste("Directional spatial correlogram -", label),
         x = "Distance (km)", y = "Correlation", color = "Direction") +
    ylim(-1, 1) +
    theme_minimal()
}

# Temporal correlogram of the mean Sv per ping.
# A correlogram by pairs, as in space, replaces the usual ACF, which needs a
# regular time step: pings are a few seconds apart, with gaps. The time-lag
# classes are quantiles (not of fixed width), because the lags between random
# pairs are very skewed (many long lags, few short ones) over transects that
# last several days.
# CAUTION: the ship moves continuously, so a short time lag almost always
# means a short distance. This correlogram mixes temporal and spatial
# decorrelation; it is not a pure effect of time.
#   n_pairs           : number of pairs drawn
#   n_bins            : number of time-lag classes
#   min_pairs_per_bin : minimum number of pairs to keep a class
#   time_unit         : unit of the time lags (`units` argument of difftime)
# Returns a ggplot (a placeholder panel with fewer than 30 pings).
make_correlogram_temporal <- function(dat, label, n_pairs = 50000, n_bins = 20,
                                      min_pairs_per_bin = 30,
                                      time_unit = "hours") {
  n <- nrow(dat)
  if (n < 30) return(empty_plot(paste("Not enough data -", label)))

  idx <- sample_pair_indices(n, n_pairs)

  dt <- abs(as.numeric(
    difftime(dat$time[idx$i], dat$time[idx$j], units = time_unit)
  ))
  z1 <- dat$value[idx$i]
  z2 <- dat$value[idx$j]

  keep  <- dt > 0 & is.finite(z1) & is.finite(z2)
  pairs <- data.frame(dt = dt[keep], z1 = z1[keep], z2 = z2[keep])

  breaks <- unique(
    quantile(pairs$dt, probs = seq(0, 1, length.out = n_bins + 1))
  )
  pairs$bin <- cut(pairs$dt, breaks = breaks, include.lowest = TRUE)

  cor_df <- pairs %>%
    group_by(bin) %>%
    summarise(dt = mean(dt), correlation = cor(z1, z2), n = n(),
              .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)

  ggplot(cor_df, aes(x = dt, y = correlation)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.8) +
    geom_line() +
    labs(title = paste("Temporal correlogram -", label),
         x = paste0("Time lag (", time_unit, ")"), y = "Correlation") +
    ylim(-1, 1) +
    theme_minimal()
}

# Four-panel figure (coverage, variogram, spatial and temporal correlograms)
# of one frequency and one subset of pings, saved in <fig_dir>.
#   freq        : frequency (kHz), for the title and the file name
#   subset_name : name of the subset of pings ("all", "day" or "night")
#   day_value   : values of profiles$day kept; NULL = all the pings
#   label       : name of the variable, for the titles
# Returns the patchwork figure.
make_profile_plots <- function(profiles, freq, subset_name, day_value,
                               label = "Sv") {
  cat("Processing:", label, freq, "kHz -", subset_name, "\n")

  dat <- make_ping_data(profiles, depth_range = depth_range,
                        min_depth_obs = min_depth_obs, n_sample = n_sample,
                        day_value = day_value)
  cov <- attr(dat, "coverage")
  cat(sprintf(
    "  -> %.0f%% NA (depth) | %d/%d pings kept | %d sampled\n",
    100 * cov$pct_na_global, cov$n_pings_kept, cov$n_pings_total,
    cov$n_pings_sampled
  ))

  p_coverage  <- make_coverage_plot(dat, label)
  p_variogram <- make_variogram(dat, label, directional = FALSE)
  p_spatial   <- make_correlogram_spatial(dat, label, directional = FALSE)
  p_temporal  <- make_correlogram_temporal(dat, label, time_unit = time_unit)

  p <- (p_coverage + p_variogram) / (p_spatial + p_temporal) +
    patchwork::plot_annotation(
      title = paste0("Spatio-temporal dependence - ", label, " ", freq,
                     " kHz (", subset_name, " pings)")
    )

  ggsave(
    filename = file.path(
      fig_dir,
      paste0("diag_autocorrelation_spatio_temp_Sv_", freq, "kHz_",
             subset_name, ".png")
    ),
    plot = p, width = 14, height = 10, dpi = 300, units = "in"
  )
  p
}


# ---- Diagnostics, per frequency and subset of pings --------------------------

for (freq in freqs) {
  profiles <- readRDS(sv_all_years_file(freq))

  for (subset_name in names(day_subsets)) {
    p <- make_profile_plots(profiles, freq, subset_name,
                            day_value = day_subsets[[subset_name]])
    print(p)
  }
}

cat("Figures saved in:", fig_dir, "\n")
