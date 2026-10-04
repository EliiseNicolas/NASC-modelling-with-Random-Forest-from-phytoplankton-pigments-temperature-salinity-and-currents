# ==============================================================================
# 00_model_config.R
#
# Random-forest pipeline: settings shared by the modelling scripts
#
# Sourced by the 1x_run_*.R scripts, after config.R.
#
# Defines : packages
#           frequencies, response and covariates
#           cross-validation schemes (naive, spatially blocked) and inner CV
#           tuning grid of the random forest
#           noise levels of the robustness test
#           default target of the post-processing scripts (11 to 14)
#           output folders, metrics and plot theme
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ranger)
  library(lubridate)
})


# ---- Scope -------------------------------------------------------------------

model_freqs <- c(38, 120)   # kHz
model_name  <- "rf"         # only model of the pipeline (used in the paths)


# ---- Variables ---------------------------------------------------------------

response_raw   <- "nasc"   # NASC column of the dataset, not in log10 yet
nasc_log_floor <- 1e-6     # avoids log10(0) = -Inf when nasc = 0

covariates_num <- c(
  "ftle", "Chla_total",
  "Per_totpig", "But_totpig", "Fuco_totpig",
  "Hex_totpig", "Allo_totpig", "Zea_totpig", "Chlb_totpig", "DvChla_totpig"
)
covariates_all <- c(covariates_num, "fod")


# ---- Cross-validation schemes ------------------------------------------------

# Naive scheme: repeated random split
naive_n_folds    <- 10
naive_train_frac <- 0.8

# Spatially blocked schemes: one per resolution. There is no buffer: the only
# guard against excessive extrapolation is that a test block is limited to ONE
# year (see 02_folds.R), the training set being all the rest (including the
# same area in the other years).
spatial_resolutions <- list(
  list(cellsize_km = c(1000, 1000), label = "1000x1000km"),
  list(cellsize_km = c(20, 20),     label = "20x20km")
)
block_min_n      <- 50    # minimum observations for a (block, year) to be used
block_folds_frac <- 0.3   # target number of folds = this share of the blocks,
block_min_folds  <- 5     # ... at least this number
block_max_folds  <- 30    # ... and at most this number


# ---- Inner CV (nested, for a tuning without leakage) -------------------------

# Simplification: repeated random subsampling, not a nested spatial blocking
# (which would be much more expensive).
inner_cv_repeats    <- 3
inner_cv_train_frac <- 0.8


# ---- Tuning grid of the random forest ----------------------------------------

rf_tuning_grid <- expand.grid(
  mtry          = c(2, 3, 4),
  min.node.size = c(10, 20, 30),
  num.trees     = c(300, 500),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)


# ---- Robustness to Gaussian noise --------------------------------------------

# Standard deviation of the noise, as a share of the standard deviation of
# each covariate
noise_levels <- c(0, 0.0025, 0.005, 0.01, 0.02, 0.05)


# ---- Default target of the post-processing scripts (11 to 14) ----------------

# Frequency and scheme explained / mapped when `target_freq` and
# `target_scheme` are not set before sourcing the script
default_target_freq   <- 38
default_target_scheme <- "blocked_spatial_1000x1000km"


# ---- Output folders ----------------------------------------------------------

# Results (models, tables) of one frequency and one scheme
scheme_out_dir <- function(freq, scheme) {
  file.path(model_out_root, paste0(freq, "kHz"), model_name, scheme)
}

# Figures of one frequency and one scheme
scheme_fig_dir <- function(freq, scheme) {
  file.path(fig_root_modelling, paste0(freq, "kHz"), model_name, scheme)
}


# ---- Metrics -----------------------------------------------------------------

rmse_fn <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))

# R2 against a constant prediction `baseline_mean`
r2_fn <- function(obs, pred, baseline_mean) {
  1 - sum((obs - pred)^2, na.rm = TRUE) /
    sum((obs - baseline_mean)^2, na.rm = TRUE)
}

# RMSE divided by the standard deviation of the observations
nrmse_fn <- function(obs, pred) {
  rmse_fn(obs, pred) / stats::sd(obs, na.rm = TRUE)
}


# ---- Plots -------------------------------------------------------------------

theme_pipeline <- theme_bw()

# `a` unless it is NULL, then `b`
`%||%` <- function(a, b) if (is.null(a)) b else a
