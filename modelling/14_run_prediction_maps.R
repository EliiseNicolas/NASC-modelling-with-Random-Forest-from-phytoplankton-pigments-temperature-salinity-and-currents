# ==============================================================================
# 14_run_prediction_maps.R
#
# Random-forest pipeline: prediction maps (complete pixels only, no
# imputation of the missing covariates)
#
# Input  : prediction dataset (02_build_prediction_dataset.R)
#          <prediction_dataset_file>, a list with date, lon, lat, ftle, pig
#          (a list of arrays, named like the covariates) and fod; all the
#          arrays are (date x lon x lat)
#          models of the target frequency and scheme
#          (10_run_nested_cv_training.R: models.rds)
#          learning dataset of the target frequency
#          (01_build_learning_dataset.R), for the FOD levels of the models
#
# Steps  : 1) single pass over the dates: predict each day once (mean of the
#             models of the folds), draw the daily map and accumulate the sum
#             and the number of valid values of the month (constant memory)
#          2) monthly composites: map of the mean prediction, with a second
#             map (and colour bar) of the number of daily predictions
#             averaged in each pixel
#
# Output : in <model_out_root>/<freq>kHz/rf/<scheme>/predictions_monthly/
#            monthly_composite_<YYYY-MM>.csv   (lon, lat, pred, n_valid,
#                                               purity)
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/predictions_daily/
#            pred_<YYYYMMDD>.png
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/predictions_monthly/
#            monthly_mean_<YYYY-MM>.png
#            monthly_mean_n_days_<YYYY-MM>.png
#
# A pixel is predicted on a given day only if all its covariates are
# available: the number of daily predictions behind a monthly mean varies
# from one pixel to another, from 1 to the number of dates of the month in
# the prediction dataset.
#
# Usage  : optional target, to set in the console before source()
#            target_freq <- 38 ; target_scheme <- "naive_RS_80_20"
#            target_months <- "2023-02"   # strongly advised for a first run
#            save_daily_maps <- FALSE
#          target_months = NULL processes all the months of the grid (several
#          hours).
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning and prediction dataset files
for (f in c("00_model_config.R", "01_data_prep.R", "03_models.R")) {
  source(file.path(modelling_dir, f))
}

if (!exists("target_freq"))     target_freq     <- default_target_freq
if (!exists("target_scheme"))   target_scheme   <- default_target_scheme
if (!exists("target_months"))   target_months   <- NULL
if (!exists("save_daily_maps")) save_daily_maps <- TRUE

in_file <- prediction_dataset_file
out_dir <- scheme_out_dir(target_freq, target_scheme)
fig_dir <- scheme_fig_dir(target_freq, target_scheme)

monthly_out_dir <- file.path(out_dir, "predictions_monthly")
daily_fig_dir   <- file.path(fig_dir, "predictions_daily")
monthly_fig_dir <- file.path(fig_dir, "predictions_monthly")

models_file <- file.path(out_dir, "models.rds")
if (!file.exists(models_file)) {
  stop("models.rds not found in ", out_dir,
       " -- this frequency / scheme has not been trained.")
}
if (!file.exists(in_file)) {
  stop("Prediction dataset not found: ", in_file,
       " -- run 02_build_prediction_dataset.R first")
}

dir.create(monthly_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(daily_fig_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(monthly_fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Read the models, the observations and the grid --------------------------

models <- readRDS(models_file)

# Levels of fod in the data the models were trained on
fod_levels <- load_and_clean(target_freq)$fod_levels

grid_ds <- readRDS(in_file)

base_grid <- expand.grid(lon = grid_ds$lon, lat = grid_ds$lat,
                         KEEP.OUT.ATTRS = FALSE)
n_pix     <- nrow(base_grid)
lon_range <- range(base_grid$lon)
lat_range <- range(base_grid$lat)

subtitle <- sprintf("%d kHz - %s", target_freq, target_scheme)


# ---- Functions ---------------------------------------------------------------

# Covariates of the grid at one date, one row per pixel.
#   date_idx : index of the date in the grid
extract_grid_for_date <- function(date_idx) {
  grid <- base_grid
  grid$ftle <- as.vector(grid_ds$ftle[date_idx, , ])

  # Pigment covariates of the models
  for (v in setdiff(covariates_num, "ftle")) {
    grid[[v]] <- as.vector(grid_ds$pig[[v]][date_idx, , ])
  }

  # fod: matched on the levels of the model, whatever the spaces (" 1" vs
  # "1"); the string "NA" and the empty string are missing values
  fod_chr <- trimws(as.character(as.vector(grid_ds$fod[date_idx, , ])))
  fod_chr[fod_chr %in% c("NA", "")] <- NA_character_
  grid$fod <- factor(fod_levels[match(fod_chr, trimws(fod_levels))],
                     levels = fod_levels)
  grid
}

# Mean prediction of the models of the folds, on the pixels where all the
# covariates are available (NA elsewhere).
predict_complete_pixels <- function(grid_df) {
  complete <- stats::complete.cases(grid_df[, covariates_all])
  pred <- rep(NA_real_, nrow(grid_df))
  if (any(complete)) {
    sub <- grid_df[complete, , drop = FALSE]
    pred[complete] <- rowMeans(sapply(models, function(m) predict_rf(m, sub)))
  }
  pred
}


# ---- Dates to predict --------------------------------------------------------

dates  <- as.Date(grid_ds$date)
months <- format(dates, "%Y-%m")

months_to_process <- if (is.null(target_months)) {
  sort(unique(months))
} else {
  target_months
}
dates_idx <- which(months %in% months_to_process)
if (length(dates_idx) == 0) {
  stop("No date of the grid matches target_months = ",
       paste(target_months, collapse = ", "))
}
cat(sprintf("%d dates to predict (%s)\n", length(dates_idx),
            paste(months_to_process, collapse = ", ")))


# ---- 1) Daily predictions, daily maps and monthly accumulation ---------------

month_sum  <- list()         # sum of the predictions, per pixel
month_n    <- list()         # number of valid predictions, per pixel
month_days <- integer(0)     # number of days

for (k in seq_along(dates_idx)) {
  i    <- dates_idx[k]
  ym   <- months[i]
  grid <- extract_grid_for_date(i)
  pred <- predict_complete_pixels(grid)

  if (is.null(month_sum[[ym]])) {
    month_sum[[ym]] <- numeric(n_pix)
    month_n[[ym]]   <- integer(n_pix)
    month_days[ym]  <- 0L
  }
  valid <- !is.na(pred)
  month_sum[[ym]][valid] <- month_sum[[ym]][valid] + pred[valid]
  month_n[[ym]]          <- month_n[[ym]] + valid
  month_days[ym]         <- month_days[ym] + 1L

  if (save_daily_maps) {
    grid$pred <- pred
    p <- ggplot(grid[valid, ], aes(x = lon, y = lat, fill = pred)) +
      geom_tile() +
      scale_fill_viridis_c() +
      coord_quickmap(xlim = lon_range, ylim = lat_range) +
      labs(title = "Predicted NASC (RF, complete pixels only)",
           subtitle = sprintf("%s - %s - %.1f%% of pixels predicted",
                              subtitle, format(dates[i]), 100 * mean(valid)),
           x = "Longitude", y = "Latitude", fill = "log10(NASC)") +
      theme_pipeline
    ggsave(
      file.path(daily_fig_dir,
                sprintf("pred_%s.png", format(dates[i], "%Y%m%d"))),
      p, width = 8, height = 6, dpi = 150
    )
  }

  if (k %% 5 == 0 || k == length(dates_idx)) {
    cat(sprintf("  ... %d / %d dates\n", k, length(dates_idx)))
  }
}


# ---- 2) Monthly composites --------------------------------------------------

cat("\nMonthly composites...\n")

for (ym in names(month_sum)) {
  n_days    <- month_days[ym]
  mean_pred <- ifelse(month_n[[ym]] > 0, month_sum[[ym]] / month_n[[ym]],
                      NA_real_)
  composite <- tibble(
    lon     = base_grid$lon,
    lat     = base_grid$lat,
    pred    = mean_pred,
    n_valid = month_n[[ym]],            # daily predictions averaged
    purity  = month_n[[ym]] / n_days    # share of the days with a prediction
  )
  valid_pixels <- composite %>% filter(!is.na(pred))
  write.csv(valid_pixels,
            file.path(monthly_out_dir,
                      sprintf("monthly_composite_%s.csv", ym)),
            row.names = FALSE)

  # Mean prediction
  p_mean <- ggplot(valid_pixels, aes(x = lon, y = lat, fill = pred)) +
    geom_tile() +
    scale_fill_viridis_c(limits = range(valid_pixels$pred)) +
    coord_quickmap(xlim = lon_range, ylim = lat_range) +
    labs(title = "Predicted NASC -- monthly composite (RF)",
         subtitle = sprintf("%s - %s (%d days)", subtitle, ym, n_days),
         x = "Longitude", y = "Latitude", fill = "Mean\nlog10(NASC)") +
    theme_pipeline
  ggsave(file.path(monthly_fig_dir, sprintf("monthly_mean_%s.png", ym)),
         p_mean, width = 8, height = 6, dpi = 150)

  # Number of daily predictions averaged in each pixel, with its own colour
  # bar (same scale for all the months with the same number of days)
  p_n_days <- ggplot(valid_pixels, aes(x = lon, y = lat, fill = n_valid)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = c(1, max(n_days, 2))) +
    coord_quickmap(xlim = lon_range, ylim = lat_range) +
    labs(title = "Number of daily predictions in the monthly mean",
         subtitle = sprintf("%s - %s (%d days)", subtitle, ym, n_days),
         x = "Longitude", y = "Latitude", fill = "Number\nof days") +
    theme_pipeline

  p_side <- (p_mean + p_n_days) +
    plot_annotation(
      title = sprintf("Monthly composite and number of days -- %s", ym)
    )
  ggsave(
    file.path(monthly_fig_dir, sprintf("monthly_mean_n_days_%s.png", ym)),
    p_side, width = 14, height = 6, dpi = 150
  )

  cat(sprintf(
    "  %s: %d days, %.1f%% of pixels predicted, median of %g days per pixel\n",
    ym, n_days, 100 * mean(composite$n_valid > 0),
    median(valid_pixels$n_valid)
  ))
}

cat("\nResults saved in:", out_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
