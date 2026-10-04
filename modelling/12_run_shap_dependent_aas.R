# ==============================================================================
# 12_run_shap_dependent_aas.R
#
# Random-forest pipeline: SHAP values for DEPENDENT covariates
# (Aas, Jullum & Loland 2019), with the package `shapr`
#
# Input  : learning dataset of the target frequency
#          (01_build_learning_dataset.R)
#          models and fold indices of the target frequency and scheme
#          (10_run_nested_cv_training.R: models.rds, fold_indices.rds).
#          Nothing is trained again.
#
# Steps  : for the first `max_folds` folds
#          1) Shapley values of a sample of the test set, with a sample of the
#             training set as background
#          2) mean |SHAP| per covariate
#          then
#          3) importance: mean |SHAP|, per fold and summarised
#
# Output : in <model_out_root>/<freq>kHz/rf/<scheme>/
#            shap_aas_dependent_by_fold.csv, shap_aas_dependent_summary.csv
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/
#            09_shap_aas_dependent.png
#
# Requires the package `shapr` (old API, shapr < 1.0, and new API, >= 1.0).
#
# Why, in addition to TreeSHAP: TreeSHAP assumes independent covariates, and
# the pigments are strongly correlated. shapr estimates the CONDITIONAL
# distribution of the covariates that are left out.
#
# Method: the approaches of Aas et al. 2019 ("empirical", "gaussian",
# "copula") only handle NUMERIC covariates, and fod is categorical. If the
# requested approach fails for that reason, the script FALLS BACK on "ctree"
# (Redelmeier, Jullum & Aas 2020, same team, handles categorical covariates)
# and says so in its outputs. The values are then NOT strictly those of Aas
# 2019 -- to be stated in the report.
#
# Very expensive: there are 2^11 = 2048 coalitions, so the number of folds,
# of rows explained and of background rows is kept small.
#
# Usage  : optional target, to set in the console before source()
#            target_freq <- 38 ; target_scheme <- "naive_RS_80_20"
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), modelling_dir
for (f in c("00_model_config.R", "01_data_prep.R")) {
  source(file.path(modelling_dir, f))
}

if (!requireNamespace("shapr", quietly = TRUE)) {
  stop("Package `shapr` not installed: install.packages('shapr')")
}

if (!exists("target_freq"))   target_freq   <- default_target_freq
if (!exists("target_scheme")) target_scheme <- default_target_scheme

# "empirical", "gaussian", "copula" (Aas 2019) or "ctree"
approach <- "empirical"

n_background  <- 100   # training rows used as background
max_test_rows <- 30    # predictions explained per fold
max_folds     <- 3     # number of folds processed (the computation is heavy)
seed          <- 1

shapr_version <- utils::packageVersion("shapr")
new_api       <- shapr_version >= "1.0.0"
cat(sprintf("shapr %s -> %s API\n", as.character(shapr_version),
            if (new_api) "new" else "old"))

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

# Shapley values of `x_test` with shapr, with either API.
#   model    : ranger model
#   x_train  : background rows
#   x_test   : rows to explain
#   p0       : prediction without any covariate (mean response)
#   approach : shapr approach
# Returns a data frame with one column per covariate.
run_shapr <- function(model, x_train, x_test, p0, approach) {
  if (new_api) {
    explanation <- shapr::explain(model = model, x_explain = x_test,
                                  x_train = x_train, approach = approach,
                                  phi0 = p0, seed = seed)
    shap_values <- explanation$shapley_values_est
    if (is.null(shap_values)) shap_values <- explanation$shapley_values
  } else {
    explainer   <- shapr::shapr(x_train, model)
    explanation <- shapr::explain(x_test, approach = approach,
                                  explainer = explainer, prediction_zero = p0)
    shap_values <- explanation$dt
  }
  as.data.frame(shap_values)
}


# ---- Read the data, the models and the folds ---------------------------------

df     <- load_and_clean(target_freq)$df
models <- readRDS(file.path(out_dir, "models.rds"))
folds  <- readRDS(file.path(out_dir, "fold_indices.rds"))


# ---- 1) and 2) Shapley values, per fold --------------------------------------

set.seed(seed)
fold_ids      <- names(folds)[seq_len(min(max_folds, length(folds)))]
shap_per_fold <- list()
approach_used <- character(0)

for (fid in fold_ids) {
  f <- folds[[fid]]
  bg_idx   <- sample(f$train, min(n_background, length(f$train)))
  test_idx <- sample(f$test, min(max_test_rows, length(f$test)))
  x_train  <- as.data.frame(df[bg_idx, covariates_all])
  x_test   <- as.data.frame(df[test_idx, covariates_all])
  p0       <- mean(df$NASC[bg_idx])

  cat(sprintf("fold %s: shapr (%s) on %d predictions...\n",
              fid, approach, nrow(x_test)))

  used <- approach
  shap_values <- tryCatch(
    run_shapr(models[[fid]], x_train, x_test, p0, approach),
    error = function(e) {
      cat("  [!] approach", approach, "failed:", conditionMessage(e), "\n")
      if (approach == "ctree") return(NULL)

      cat("  -> falling back on 'ctree' (handles categorical covariates;",
          "NOT identical to Aas 2019)\n")
      used <<- "ctree"
      tryCatch(
        run_shapr(models[[fid]], x_train, x_test, p0, "ctree"),
        error = function(e2) {
          cat("  [!] ctree also failed:", conditionMessage(e2), "\n")
          NULL
        }
      )
    }
  )
  if (is.null(shap_values)) next

  # Keep the covariates (drops the column of the baseline)
  shap_values <- shap_values[, intersect(names(shap_values), covariates_all),
                             drop = FALSE]
  mean_abs <- colMeans(abs(shap_values))
  shap_per_fold[[fid]] <- tibble(variable = names(mean_abs),
                                 importance = unname(mean_abs),
                                 fold_id = fid, approach = used)
  approach_used <- c(approach_used, used)
}

if (length(shap_per_fold) == 0) {
  stop("No fold explained -- see the [!] messages above (shapr version?).")
}


# ---- 3) Importance: mean |SHAP| ----------------------------------------------

shap_all <- bind_rows(shap_per_fold)
shap_summary <- shap_all %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance),
            sd_importance = sd(importance), .groups = "drop") %>%
  arrange(desc(mean_importance))

write.csv(shap_all, file.path(out_dir, "shap_aas_dependent_by_fold.csv"),
          row.names = FALSE)
write.csv(shap_summary, file.path(out_dir, "shap_aas_dependent_summary.csv"),
          row.names = FALSE)

approach_label <- paste(unique(approach_used), collapse = "/")

p <- ggplot(
  shap_summary,
  aes(x = reorder(variable, mean_importance), y = mean_importance)
) +
  geom_col(fill = "purple") +
  geom_errorbar(aes(ymin = mean_importance - sd_importance,
                    ymax = mean_importance + sd_importance), width = 0.2) +
  coord_flip() +
  labs(title = "SHAP importance (dependent covariates, shapr)",
       subtitle = sprintf("RF - %d kHz - %s - approach: %s - %d folds",
                          target_freq, target_scheme, approach_label,
                          length(shap_per_fold)),
       x = NULL, y = "Mean |SHAP|") +
  theme_pipeline
ggsave(file.path(fig_dir, "09_shap_aas_dependent.png"), p,
       width = 8, height = 6, dpi = 150)

cat("Approach(es) used:", approach_label, "\n")
cat("Results saved in:", out_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
