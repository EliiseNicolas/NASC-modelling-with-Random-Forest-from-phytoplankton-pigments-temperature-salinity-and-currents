# ==============================================================================
# spatio_temporal_autocorrelation_pigments.R
#
# Data validation: spatio-temporal dependence of the pigments
#
# Input  : pigment dataset of all years, cropped to the FOD area and filtered
#          (01_create_ds_pigments_filtered_all_years.R)
#          <pigments_file>
#
# Steps  : for each pigment
#          1) time mean of each pixel, with the coverage: number of valid
#             observations per pixel
#          2) directional variogram of the time mean (anisotropy)
#          3) directional spatial correlogram
#          4) temporal autocorrelation (ACF) of the spatial mean
#          then, for all the pigments together (standardised values)
#          5) variogram and correlogram along the latitude axis and along the
#             longitude axis
#
# Output : in <fig_dir>
#            diag_autocorrelation_spatio_temp_<pigment>.png
#            diag_autocorrelation_lat_lon_all_pigments.png
# ==============================================================================

library(gstat)
library(ggplot2)
library(dplyr)
library(patchwork)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, pigments_file

in_file <- pigments_file
fig_dir <- file.path(fig_root_validation, "spatio_temporal_autocorrelation")

pig_names <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb",
               "DvChla")

# Figure of each pigment
n_sample_pigment <- 5000   # maximum number of pixels kept (random sample)
min_obs_pigment  <- 10     # minimum number of valid dates to keep a pixel

# Figure of all the pigments, per axis
n_sample_axis <- 3000
min_obs_axis  <- 3
tol_hor       <- 22.5      # angular tolerance around the axis (degrees)

seed <- 123

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

# One colour per pigment, the same in the four panels of the last figure
# (otherwise a pigment missing from one panel would shift the colours)
pigment_colors <- setNames(scales::hue_pal()(length(pig_names)), pig_names)

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

# Put the array of a pigment in the order (time, lon, lat).
# With an array in (time, lat, lon), expand.grid(lon, lat) + as.vector() would
# silently misalign the values and the coordinates: the order is checked and
# fixed. Error if the dimensions match neither (lon, lat) nor (lat, lon).
#   x    : 3-D array of a pigment
#   pigs : pigment dataset (for the length of lon and lat)
fix_lonlat_order <- function(x, pigs) {
  d <- dim(x)
  n_lon <- length(pigs$lon)
  n_lat <- length(pigs$lat)
  if (d[2] == n_lon && d[3] == n_lat) return(x)
  if (d[2] == n_lat && d[3] == n_lon) return(aperm(x, c(1, 3, 2)))
  stop("Array dimensions do not match pigs$lon / pigs$lat.")
}


# ---- Functions: pixel data ---------------------------------------------------

# Time mean of a pigment per pixel, with the coverage.
# With ~80-90 % of NA (clouds), the number of valid observations used in the
# mean of each pixel (n_valid) is counted, and a minimum (min_obs) is required
# to keep the pixel: this avoids unstable means computed on 1 or 2 values.
#   pigment    : name of the pigment (e.g. "Chla")
#   pigs       : pigment dataset
#   date_index : indices of the dates kept (NULL = all)
#   n_sample   : maximum number of pixels kept (random sample, fixed seed)
#   min_obs    : minimum number of valid observations to keep a pixel
# Returns a data frame with one row per pixel kept: lon, lat, value (time
# mean), n_valid, frac_valid. The attribute "coverage" summarises the effect
# of the NA (share of NA, pixels in total / kept / sampled).
make_spatial_data <- function(pigment, pigs, date_index = NULL,
                              n_sample = 3000, min_obs = 3) {
  x <- pigs[[paste0("c_cond_", pigment)]]
  x <- fix_lonlat_order(x, pigs)
  if (!is.null(date_index)) x <- x[date_index, , , drop = FALSE]

  n_time  <- dim(x)[1]
  n_valid <- apply(x, c(2, 3), function(v) sum(!is.na(v)))
  z <- apply(x, c(2, 3), mean, na.rm = TRUE)
  z[is.nan(z)] <- NA

  dat <- expand.grid(lon = pigs$lon, lat = pigs$lat)
  dat$value      <- as.vector(z)
  dat$n_valid    <- as.vector(n_valid)
  dat$frac_valid <- dat$n_valid / n_time

  pct_na_global  <- 1 - sum(dat$n_valid) / (n_time * nrow(dat))
  n_pixels_total <- nrow(dat)

  dat <- dat %>% filter(is.finite(value), n_valid >= min_obs)
  n_pixels_kept <- nrow(dat)

  if (nrow(dat) > n_sample) {
    set.seed(seed)
    dat <- dat %>% slice_sample(n = n_sample)
  }

  attr(dat, "coverage") <- data.frame(
    pigment = pigment, n_time = n_time, pct_na_global = pct_na_global,
    n_pixels_total = n_pixels_total, n_pixels_kept = n_pixels_kept,
    n_pixels_sampled = nrow(dat)
  )
  dat
}


# ---- Functions: figure of one pigment ----------------------------------------

# Map of the number of valid observations per pixel.
#   dat   : output of make_spatial_data()
#   label : name of the pigment, for the title
make_coverage_plot <- function(dat, label) {
  cov <- attr(dat, "coverage")
  subtitle <- sprintf(
    "%.0f%% NA (clouds) | %d/%d pixels kept (min_obs) | %d sampled",
    100 * cov$pct_na_global, cov$n_pixels_kept, cov$n_pixels_total,
    cov$n_pixels_sampled
  )
  ggplot(dat, aes(x = lon, y = lat, color = n_valid)) +
    geom_point(size = 1.2) +
    scale_color_viridis_c(name = "Valid\nobservations") +
    labs(title = paste("Temporal coverage -", label), subtitle = subtitle,
         x = "Longitude", y = "Latitude") +
    theme_minimal()
}

# Directional variogram of a pigment. A strong lon / lat anisotropy gives very
# different curves from one direction to another.
#   cutoff : maximum distance (km); NULL = default_cutoff()
#   alpha  : directions (degrees clockwise from North, gstat convention:
#            0 = North-South, 90 = East-West) among 0, 45, 90, 135
# Returns a ggplot (a placeholder panel with fewer than 30 pixels).
make_variogram <- function(dat, label, cutoff = NULL,
                           alpha = c(0, 45, 90, 135)) {
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(empty_plot(paste("Not enough data -", label)))

  if (is.null(cutoff)) cutoff <- default_cutoff(dat)

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

# Directional spatial correlogram of a pigment: correlation between the two
# members of random pairs of pixels, per direction (the four classes of the
# variogram) and per distance class.
#   n_pairs           : number of pairs drawn
#   n_bins            : number of distance classes
#   min_pairs_per_bin : minimum number of pairs to keep a class
# Returns a ggplot (a placeholder panel with fewer than 30 pixels).
make_correlogram_spatial <- function(dat, label, n_pairs = 50000, n_bins = 20,
                                     min_pairs_per_bin = 30) {
  dat <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(empty_plot(paste("Not enough data -", label)))

  pairs <- sample_pairs(dat, n_pairs)
  pairs$direction <- angle_to_direction(pairs$angle)
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)

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

# Temporal autocorrelation (ACF) of the spatial mean of a pigment.
# The lag is a number of successive dates of the dataset: it is a number of
# days only within one year without missing date.
#   lag_max : maximum lag
# Returns a ggplot (a placeholder panel if the series is too short).
make_acf <- function(pigment, pigs, lag_max = 60) {
  x <- pigs[[paste0("c_cond_", pigment)]]
  x <- fix_lonlat_order(x, pigs)

  spatial_mean <- apply(x, 1, mean, na.rm = TRUE)
  spatial_mean[is.nan(spatial_mean)] <- NA

  if (sum(!is.na(spatial_mean)) < lag_max + 5) {
    return(empty_plot(paste("Series too short -", pigment)))
  }

  acf_obj <- acf(spatial_mean, lag.max = lag_max, na.action = na.pass,
                 plot = FALSE)
  acf_df <- data.frame(lag = acf_obj$lag[, 1, 1], acf = acf_obj$acf[, 1, 1])

  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(title = paste("Temporal autocorrelation -", pigment),
         x = "Lag (successive dates)", y = "ACF") +
    theme_minimal()
}

# Four-panel figure of one pigment (coverage, variogram, spatial correlogram
# and ACF), saved in <fig_dir>.
#   n_sample, min_obs : see make_spatial_data()
# Returns a list: plot (patchwork figure) and coverage (effect of the NA).
make_pigment_plots <- function(pigment, pigs, n_sample = 3000, min_obs = 3) {
  cat("Processing:", pigment, "\n")

  dat <- make_spatial_data(pigment, pigs, n_sample = n_sample,
                           min_obs = min_obs)
  cov <- attr(dat, "coverage")
  cat(sprintf(
    "  -> %.0f%% NA | %d/%d pixels kept (min_obs = %d) | %d sampled\n",
    100 * cov$pct_na_global, cov$n_pixels_kept, cov$n_pixels_total, min_obs,
    cov$n_pixels_sampled
  ))

  p_coverage  <- make_coverage_plot(dat, pigment)
  p_variogram <- make_variogram(dat, pigment)
  p_spatial   <- make_correlogram_spatial(dat, pigment)
  p_acf       <- make_acf(pigment, pigs)

  p <- (p_coverage + p_variogram) / (p_spatial + p_acf) +
    patchwork::plot_annotation(
      title = paste("Spatio-temporal dependence -", pigment)
    )

  ggsave(
    filename = file.path(
      fig_dir, paste0("diag_autocorrelation_spatio_temp_", pigment, ".png")
    ),
    plot = p, width = 14, height = 10, dpi = 300, units = "in"
  )

  list(plot = p, coverage = cov)
}


# ---- Functions: figure of all the pigments, per axis -------------------------

# Variogram of a pigment along one axis.
#   axis    : "lat" (North-South, gstat alpha = 0) or "lon" (East-West,
#             alpha = 90)
#   cutoff  : maximum distance (km); NULL = default_cutoff()
#   tol_hor : angular tolerance around the axis (degrees)
# Returns a data frame (pigment, dist, gamma), or NULL with fewer than 30
# pixels.
variogram_axis <- function(dat, label, axis = c("lat", "lon"), cutoff = NULL,
                           tol_hor = 22.5) {
  axis <- match.arg(axis)
  dat  <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(NULL)
  if (is.null(cutoff)) cutoff <- default_cutoff(dat)

  v <- variogram(value ~ 1, locations = ~x_km + y_km, data = dat,
                 cutoff = cutoff, width = cutoff / 15,
                 alpha = if (axis == "lat") 0 else 90, tol.hor = tol_hor)
  data.frame(pigment = label, dist = v$dist, gamma = v$gamma)
}

# Correlogram of a pigment along one axis.
#   axis              : "lat" (North-South) or "lon" (East-West)
#   n_pairs           : number of pairs drawn
#   n_bins            : number of distance classes
#   min_pairs_per_bin : minimum number of pairs to keep a class
#   tol_hor           : angular tolerance around the axis (degrees)
# Returns a data frame (pigment, distance, correlation), or NULL with too few
# pixels or pairs.
correlogram_axis <- function(dat, label, axis = c("lat", "lon"),
                             n_pairs = 50000, n_bins = 20,
                             min_pairs_per_bin = 30, tol_hor = 22.5) {
  axis <- match.arg(axis)
  dat  <- lonlat_to_km(dat)
  if (nrow(dat) < 30) return(NULL)

  pairs <- sample_pairs(dat, n_pairs)

  in_axis <- if (axis == "lat") {
    # Meridional (N-S)
    pairs$angle >= (90 - tol_hor) & pairs$angle <= (90 + tol_hor)
  } else {
    # Zonal (E-W)
    pairs$angle <= tol_hor | pairs$angle >= (180 - tol_hor)
  }

  pairs <- pairs[in_axis, ]
  if (nrow(pairs) < min_pairs_per_bin) return(NULL)
  pairs$bin <- cut(pairs$distance, breaks = n_bins, include.lowest = TRUE)

  cor_df <- pairs %>%
    group_by(bin) %>%
    summarise(distance = mean(distance), correlation = cor(z1, z2),
              n = n(), .groups = "drop") %>%
    filter(n >= min_pairs_per_bin)

  data.frame(pigment = label, distance = cor_df$distance,
             correlation = cor_df$correlation)
}

# Variograms of all the pigments along one axis, with the mean trend (loess)
# in black.
#   df : data frame (pigment, dist, gamma)
plot_variogram_axis <- function(df, title) {
  ggplot(df, aes(x = dist, y = gamma, color = pigment)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2,
                method = "loess") +
    scale_color_manual(values = pigment_colors, limits = pig_names) +
    labs(title = title, x = "Distance (km)", y = "Semivariance (z-score)",
         color = "Pigment") +
    theme_minimal()
}

# Correlograms of all the pigments along one axis, with the mean trend (loess)
# in black.
#   df : data frame (pigment, distance, correlation)
plot_correlogram_axis <- function(df, title) {
  ggplot(df, aes(x = distance, y = correlation, color = pigment)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(alpha = 0.6) +
    geom_point(size = 1, alpha = 0.6) +
    geom_smooth(aes(group = 1), color = "black", se = FALSE, linewidth = 1.2,
                method = "loess") +
    scale_color_manual(values = pigment_colors, limits = pig_names) +
    labs(title = title, x = "Distance (km)", y = "Correlation",
         color = "Pigment") +
    ylim(-1, 1) +
    theme_minimal()
}


# ---- Read the dataset --------------------------------------------------------

pigs <- readRDS(in_file)
cat("Latitude range :", range(pigs$lat), "\n")
cat("Longitude range:", range(pigs$lon), "\n")


# ---- 1) to 4) Figure of each pigment -----------------------------------------

plots_pigments   <- list()
coverage_summary <- list()

# An error on one pigment does not stop the others
for (pig in pig_names) {
  res <- tryCatch(
    make_pigment_plots(pig, pigs, n_sample = n_sample_pigment,
                       min_obs = min_obs_pigment),
    error = function(e) {
      cat("  ERROR for", pig, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (!is.null(res)) {
    plots_pigments[[pig]]   <- res$plot
    coverage_summary[[pig]] <- res$coverage
  }
}

# Effect of the NA, all pigments
print(bind_rows(coverage_summary))

for (pig in names(plots_pigments)) print(plots_pigments[[pig]])


# ---- 5) Variogram and correlogram per axis, all pigments ---------------------

variogram_lat   <- list()
variogram_lon   <- list()
correlogram_lat <- list()
correlogram_lon <- list()

for (pig in pig_names) {
  dat <- tryCatch(
    make_spatial_data(pig, pigs, n_sample = n_sample_axis,
                      min_obs = min_obs_axis),
    error = function(e) NULL
  )
  if (is.null(dat) || nrow(dat) < 30) next

  # z-score per pigment, to make the curves comparable
  dat$value <- as.numeric(scale(dat$value))

  variogram_lat[[pig]] <- variogram_axis(dat, pig, "lat", tol_hor = tol_hor)
  variogram_lon[[pig]] <- variogram_axis(dat, pig, "lon", tol_hor = tol_hor)
  correlogram_lat[[pig]] <- correlogram_axis(dat, pig, "lat",
                                             tol_hor = tol_hor)
  correlogram_lon[[pig]] <- correlogram_axis(dat, pig, "lon",
                                             tol_hor = tol_hor)
}

p_variogram_lat <- plot_variogram_axis(
  bind_rows(variogram_lat), "Variogram - latitude axis (meridional)"
)
p_variogram_lon <- plot_variogram_axis(
  bind_rows(variogram_lon), "Variogram - longitude axis (zonal)"
)
p_correlogram_lat <- plot_correlogram_axis(
  bind_rows(correlogram_lat), "Correlogram - latitude axis (meridional)"
)
p_correlogram_lon <- plot_correlogram_axis(
  bind_rows(correlogram_lon), "Correlogram - longitude axis (zonal)"
)

# guides = "collect" merges the four identical legends into one
p_axes <- (p_variogram_lat + p_variogram_lon) /
  (p_correlogram_lat + p_correlogram_lon) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Spatial dependence per axis - all pigments"
  ) &
  theme(legend.position = "right")

print(p_axes)

ggsave(
  filename = file.path(
    fig_dir, "diag_autocorrelation_lat_lon_all_pigments.png"
  ),
  plot = p_axes, width = 16, height = 10, dpi = 300, units = "in"
)

cat("Figures saved in:", fig_dir, "\n")
