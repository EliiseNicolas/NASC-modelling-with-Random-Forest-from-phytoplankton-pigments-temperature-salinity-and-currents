# ==============================================================================
# 11_run_shap_treeshap.R
#
# Random-forest pipeline: exact SHAP values (TreeSHAP)
#
# Input  : learning dataset of the target frequency
#          (01_build_learning_dataset.R)
#          models and fold indices of the target frequency and scheme
#          (10_run_nested_cv_training.R: models.rds, fold_indices.rds).
#          Nothing is trained again.
#
# Steps  : for each fold
#          1) unify the ranger model on a sample of the training set
#          2) TreeSHAP on a sample of the test set
#          3) mean |SHAP| per covariate
#          then
#          4) importance: mean |SHAP|, per fold and summarised
#          5) direction of the effect: signed SHAP vs value of the covariate
#
# Output : in <model_out_root>/<freq>kHz/rf/<scheme>/
#            shap_treeshap_by_fold.csv, shap_treeshap_summary.csv
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/
#            08a_shap_treeshap_importance.png
#            08b_shap_treeshap_direction.png
#
# Requires the package `treeshap`:
#   remotes::install_github("ModelOriented/treeshap")
#
# fod is a factor: ranger (respect.unordered.factors = "ignore", the default
# in regression) treats it as ORDERED, in the order of its levels. It is
# therefore given to treeshap as integer codes, which reproduces exactly the
# splits of the model.
#
# Usage  : optional target, to set in the console before source()
#            target_freq <- 38 ; target_scheme <- "naive_RS_80_20"
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), modelling_dir
for (f in c("00_model_config.R", "01_data_prep.R")) {
  source(file.path(modelling_dir, f))
}

if (!requireNamespace("treeshap", quietly = TRUE)) {
  stop("Package `treeshap` not installed: ",
       "remotes::install_github('ModelOriented/treeshap')")
}

if (!exists("target_freq"))   target_freq   <- default_target_freq
if (!exists("target_scheme")) target_scheme <- default_target_scheme

max_test_rows       <- 2000   # test rows explained per fold
max_background_rows <- 5000   # training rows used to unify the model
n_sample_per_fold   <- 500    # rows kept per fold in the direction plot
seed                <- 1

out_dir <- scheme_out_dir(target_freq, target_scheme)
fig_dir <- scheme_fig_dir(target_freq, target_scheme)

in_files <- file.path(out_dir, c("models.rds", "fold_indices.rds"))
if (!all(file.exists(in_files))) {
  stop("Missing files in ", out_dir, ": ",
       paste(basename(in_files[!file.exists(in_files)]), collapse = ", "),
       "\n  -> this frequency / scheme has not been trained ",
       "(see target_freq / target_scheme).")
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Covariates as a data frame, with fod as integer codes.
encode_for_treeshap <- function(x) {
  x <- as.data.frame(x)
  x$fod <- as.integer(x$fod)
  x
}


# ---- Read the data, the models and the folds ---------------------------------

df     <- load_and_clean(target_freq)$df
models <- readRDS(file.path(out_dir, "models.rds"))
folds  <- readRDS(file.path(out_dir, "fold_indices.rds"))


# ---- 1) to 3) TreeSHAP, per fold ---------------------------------------------

set.seed(seed)
shap_per_fold <- list()
shap_samples  <- list()

for (fid in names(folds)) {
  f <- folds[[fid]]
  test_idx <- if (length(f$test) > max_test_rows) {
    sample(f$test, max_test_rows)
  } else {
    f$test
  }
  bg_idx <- if (length(f$train) > max_background_rows) {
    sample(f$train, max_background_rows)
  } else {
    f$train
  }

  x_bg   <- encode_for_treeshap(df[bg_idx, covariates_all])
  x_test <- encode_for_treeshap(df[test_idx, covariates_all])

  cat(sprintf("fold %s: treeshap on %d rows (background: %d rows)...\n",
              fid, nrow(x_test), nrow(x_bg)))

  unified <- tryCatch(
    treeshap::ranger.unify(models[[fid]], x_bg),
    error = function(e) {
      cat("  [!] ranger.unify:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(unified)) next

  shap_res <- tryCatch(
    treeshap::treeshap(unified, x_test, verbose = FALSE),
    error = function(e) {
      cat("  [!] treeshap:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(shap_res)) next

  shap_values <- as.data.frame(shap_res$shaps)
  mean_abs    <- colMeans(abs(shap_values))
  shap_per_fold[[fid]] <- tibble(variable = names(mean_abs),
                                 importance = unname(mean_abs), fold_id = fid)

  keep <- if (nrow(shap_values) > n_sample_per_fold) {
    sample(nrow(shap_values), n_sample_per_fold)
  } else {
    seq_len(nrow(shap_values))
  }
  shap_samples[[fid]] <- bind_rows(lapply(names(shap_values), function(v) {
    tibble(fold_id = fid, variable = v,
           shap = shap_values[[v]][keep], value = x_test[[v]][keep])
  }))
}

if (length(shap_per_fold) == 0) {
  stop("No fold processed -- see the [!] messages above.")
}


# ---- 4) Importance: mean |SHAP| ----------------------------------------------

shap_all <- bind_rows(shap_per_fold)
shap_summary <- shap_all %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance),
            sd_importance = sd(importance), .groups = "drop") %>%
  arrange(desc(mean_importance))

write.csv(shap_all, file.path(out_dir, "shap_treeshap_by_fold.csv"),
          row.names = FALSE)
write.csv(shap_summary, file.path(out_dir, "shap_treeshap_summary.csv"),
          row.names = FALSE)

subtitle <- sprintf("RF - %d kHz - %s", target_freq, target_scheme)

p_importance <- ggplot(
  shap_summary,
  aes(x = reorder(variable, mean_importance), y = mean_importance)
) +
  geom_col(fill = "darkgreen") +
  geom_errorbar(aes(ymin = mean_importance - sd_importance,
                    ymax = mean_importance + sd_importance), width = 0.2) +
  coord_flip() +
  labs(title = "SHAP importance (TreeSHAP), mean |SHAP|", subtitle = subtitle,
       x = NULL, y = "Mean |SHAP|") +
  theme_pipeline
ggsave(file.path(fig_dir, "08a_shap_treeshap_importance.png"), p_importance,
       width = 8, height = 6, dpi = 150)


# ---- 5) Direction of the effect ----------------------------------------------

# Signed SHAP, coloured by the value of the covariate scaled to [0, 1]
# (fod is categorical: no colour)
direction_df <- bind_rows(shap_samples) %>%
  group_by(variable) %>%
  mutate(value_scaled = if (variable[1] == "fod") {
    NA_real_
  } else {
    (value - min(value)) / (max(value) - min(value) + 1e-12)
  }) %>%
  ungroup()

p_direction <- ggplot(
  direction_df,
  aes(x = shap, y = reorder(variable, abs(shap), FUN = mean),
      color = value_scaled)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_jitter(height = 0.22, size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(na.value = "grey60",
                        name = "Covariate\nvalue\n(scaled)") +
  labs(title = "Direction of the effect (signed SHAP per covariate)",
       subtitle = paste0(subtitle, " -- grey: fod (categorical)"),
       x = "SHAP value (effect on log10 NASC)", y = NULL) +
  theme_pipeline
ggsave(file.path(fig_dir, "08b_shap_treeshap_direction.png"), p_direction,
       width = 9, height = 6, dpi = 150)

cat("Results saved in:", out_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
