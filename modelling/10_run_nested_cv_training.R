# ==============================================================================
# 10_run_nested_cv_training.R
#
# Random-forest pipeline, main script: nested CV, training and diagnostics
#
# Input  : learning dataset (NASC, FOD, pigments and FTLE per ESU), one file
#          per frequency (01_build_learning_dataset.R)
#          <learning_dataset_dir>/learning_dataset_<years_tag>_<freq>kHz.rds
#          restricted to the ESU with day == `day_code` (config.R)
#
# Steps  : for each frequency of `model_freqs` and each CV scheme
#          1) nested CV: inner tuning without leakage, refit and evaluation
#             on the test set of each fold (04_diagnostics.R)
#          2) learning curve
#          3) observed vs predicted NASC
#          4) distribution of the NASC, training vs test set
#          5) map of the residuals
#          6) robustness to Gaussian noise (relative and absolute)
#          7) variable importance
#          8) location of the training and test sets
#          9) overall RMSE, R2 and NRMSE
#          Every figure exists aggregated over the folds and per fold.
#
# Output : in <model_out_root>/<freq>kHz/rf/<scheme>/
#            fold_indices.rds, fold_assignment.csv
#            models.rds, metrics.csv, obs_pred.csv, importance.csv, tuning.csv
#            learning_curve_summary.csv, learning_curve_detail.csv
#            noise_robustness.csv
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/
#            01a_learning_curve_agg.png, 01b_learning_curve_by_fold.png
#            02a_obs_pred_agg.png, 02b_obs_pred_by_fold.png
#            03a_distrib_train_test_agg.png, 03b_distrib_train_test_by_fold.png
#            04a_residual_map_agg.png, 04b_residual_map_by_fold.png
#            05a_noise_relative_agg.png, 05b_noise_relative_by_fold.png
#            05c_noise_absolute_agg.png, 05d_noise_absolute_by_fold.png
#            06_importance.png, 07_spatial_train_test.png
#          <model_out_root>/nrmse_summary_all.csv
#
# Model validation: NRMSE is the criterion used to compare the schemes and
# the frequencies (summary table at the end of the script); RMSE and R2 are
# still written on each plot.
#
# Resuming (`skip_existing`): the nested CV is the expensive part. If
# models.rds, metrics.csv, obs_pred.csv, importance.csv, tuning.csv and
# fold_indices.rds already exist for a frequency and a scheme, they are
# reloaded (fold indices included, so that they are never drawn again) and
# only the post-processing is done again.
#
# Errors (safe_step()): each post-processing step is protected. An error in a
# plot prints the name of the step and the R message, then the script goes on
# with the next steps.
#
# Usage  : optional filters, to set in the console before source()
#            only_freqs <- 38 ; only_schemes <- "naive_RS_80_20"
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), modelling_dir
for (f in c("00_model_config.R", "01_data_prep.R", "02_folds.R", "03_models.R",
            "04_diagnostics.R", "05_plots.R")) {
  source(file.path(modelling_dir, f))
}

skip_existing <- TRUE   # reload an existing nested CV instead of running it

# Optional filters (NULL = everything)
if (!exists("only_freqs"))   only_freqs   <- NULL
if (!exists("only_schemes")) only_schemes <- NULL

dir.create(model_out_root, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Run a post-processing step; on error, print the name of the step and the
# message, and go on.
safe_step <- function(name, expr) {
  tryCatch(expr, error = function(e) {
    cat(sprintf("  [ERROR in step '%s'] %s\n", name, conditionMessage(e)))
    invisible(NULL)
  })
}


# ---- Nested CV and diagnostics, per frequency and scheme ---------------------

summary_list <- list()

# An error that is not protected by safe_step() prints its message, its call
# and the call stack (to locate the faulty instruction), then stops the script.
withCallingHandlers({

  for (freq in model_freqs) {
    if (!is.null(only_freqs) && !(freq %in% only_freqs)) next

    cat("\n============================================================\n")
    cat("NESTED CV --", freq, "kHz\n")
    cat("============================================================\n")

    df      <- load_and_clean(freq)$df
    schemes <- build_all_schemes(df)

    for (scheme_name in names(schemes)) {
      if (!is.null(only_schemes) && !(scheme_name %in% only_schemes)) next
      scheme <- schemes[[scheme_name]]
      if (length(scheme$folds) == 0) next

      cat(sprintf("\n-- %s (%d folds) --\n", scheme_name,
                  length(scheme$folds)))

      out_dir <- scheme_out_dir(freq, scheme_name)
      fig_dir <- scheme_fig_dir(freq, scheme_name)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

      # ---- 1) Nested CV (or reload) ----

      folds_file <- file.path(out_dir, "fold_indices.rds")
      cv_files <- c(
        models     = file.path(out_dir, "models.rds"),
        metrics    = file.path(out_dir, "metrics.csv"),
        obs_pred   = file.path(out_dir, "obs_pred.csv"),
        importance = file.path(out_dir, "importance.csv"),
        tuning     = file.path(out_dir, "tuning.csv")
      )
      resume <- skip_existing && all(file.exists(c(cv_files, folds_file)))

      if (resume) {
        cat("  [done] nested CV reloaded from disk (fold indices included)\n")
        scheme$folds <- readRDS(folds_file)
        cv_res <- list(
          models     = readRDS(cv_files[["models"]]),
          metrics    = read.csv(cv_files[["metrics"]]),
          obs_pred   = read.csv(cv_files[["obs_pred"]]),
          importance = read.csv(cv_files[["importance"]]),
          tuning     = read.csv(cv_files[["tuning"]])
        )
      } else {
        # The fold indices themselves are saved (not only the results):
        # scripts 11 to 13 reload them instead of drawing them again
        saveRDS(scheme$folds, folds_file)
        fold_assignment <- imap_dfr(scheme$folds, function(f, fid) {
          bind_rows(
            tibble(row_id = f$train, fold_id = fid, set = "train"),
            tibble(row_id = f$test,  fold_id = fid, set = "test")
          )
        })
        write.csv(fold_assignment, file.path(out_dir, "fold_assignment.csv"),
                  row.names = FALSE)

        cv_res <- run_nested_cv_scheme(scheme)

        saveRDS(cv_res$models, cv_files[["models"]])
        for (tab in c("metrics", "obs_pred", "importance", "tuning")) {
          write.csv(cv_res[[tab]], cv_files[[tab]], row.names = FALSE)
        }
      }

      label <- sprintf("RF - %d kHz - %s", freq, scheme_name)
      folds <- scheme$folds

      # ---- 2) Learning curve ----

      safe_step("learning curve", {
        m <- cv_res$metrics
        best_params_by_fold <- setNames(
          lapply(seq_len(nrow(m)), function(i) {
            list(mtry          = m$best_mtry[i],
                 min.node.size = m$best_min_node_size[i],
                 num.trees     = m$best_num_trees[i])
          }),
          m$fold_id
        )
        lc <- compute_learning_curve(df, folds, best_params_by_fold)
        write.csv(lc$summary,
                  file.path(out_dir, "learning_curve_summary.csv"),
                  row.names = FALSE)
        write.csv(lc$detail,
                  file.path(out_dir, "learning_curve_detail.csv"),
                  row.names = FALSE)

        rmse_range <- range(c(0, lc$detail$rmse_train, lc$detail$rmse_test),
                            na.rm = TRUE)
        save_fig(plot_learning_curve_agg(lc$summary, label, rmse_range),
                 fig_dir, "01a_learning_curve_agg.png")
        save_fig(plot_learning_curve_by_fold(lc$detail, label, rmse_range),
                 fig_dir, "01b_learning_curve_by_fold.png",
                 width = 10, height = 8)
      })

      # ---- 3) Observed vs predicted ----

      safe_step("observed vs predicted", {
        axis_range <- range(c(cv_res$obs_pred$obs, cv_res$obs_pred$pred),
                            na.rm = TRUE)
        save_fig(plot_obs_pred_agg(cv_res$obs_pred, label, axis_range),
                 fig_dir, "02a_obs_pred_agg.png")
        save_fig(plot_obs_pred_by_fold(cv_res$obs_pred, cv_res$metrics, label,
                                       axis_range),
                 fig_dir, "02b_obs_pred_by_fold.png", width = 10, height = 8)
      })

      # ---- 4) Distribution, training vs test set ----

      safe_step("train / test distribution", {
        save_fig(plot_train_test_distribution(df, folds, label,
                                              by_fold = FALSE),
                 fig_dir, "03a_distrib_train_test_agg.png")
        save_fig(plot_train_test_distribution(df, folds, label,
                                              by_fold = TRUE),
                 fig_dir, "03b_distrib_train_test_by_fold.png",
                 width = 10, height = 8)
      })

      # ---- 5) Map of the residuals ----

      safe_step("map of the residuals", {
        max_abs  <- max(abs(cv_res$obs_pred$residual))
        res_lims <- c(-max_abs, max_abs)
        save_fig(plot_residual_map(cv_res$obs_pred, label, res_lims,
                                   by_fold = FALSE),
                 fig_dir, "04a_residual_map_agg.png")
        save_fig(plot_residual_map(cv_res$obs_pred, label, res_lims,
                                   by_fold = TRUE),
                 fig_dir, "04b_residual_map_by_fold.png",
                 width = 10, height = 8)
      })

      # ---- 6) Robustness to noise ----

      safe_step("robustness to noise", {
        noise_df <- compute_noise_robustness(df, folds, cv_res$models)
        write.csv(noise_df, file.path(out_dir, "noise_robustness.csv"),
                  row.names = FALSE)
        save_fig(plot_noise_robustness(noise_df, relative = TRUE, label,
                                       by_fold = FALSE),
                 fig_dir, "05a_noise_relative_agg.png")
        save_fig(plot_noise_robustness(noise_df, relative = TRUE, label,
                                       by_fold = TRUE),
                 fig_dir, "05b_noise_relative_by_fold.png", width = 10)
        save_fig(plot_noise_robustness(noise_df, relative = FALSE, label,
                                       by_fold = FALSE),
                 fig_dir, "05c_noise_absolute_agg.png")
        save_fig(plot_noise_robustness(noise_df, relative = FALSE, label,
                                       by_fold = TRUE),
                 fig_dir, "05d_noise_absolute_by_fold.png", width = 10)
      })

      # ---- 7) Variable importance ----

      safe_step("importance", {
        save_fig(plot_importance(cv_res$importance, label),
                 fig_dir, "06_importance.png")
      })

      # ---- 8) Location of the training and test sets ----

      safe_step("location of the training and test sets", {
        save_fig(plot_spatial_train_test(scheme, label),
                 fig_dir, "07_spatial_train_test.png", width = 10, height = 8)
      })

      # ---- 9) Overall RMSE, R2 and NRMSE (all the test sets pooled) ----

      obs  <- cv_res$obs_pred$obs
      pred <- cv_res$obs_pred$pred
      scheme_summary <- tibble(
        freq    = freq,
        scheme  = scheme_name,
        n_folds = length(folds),
        rmse    = rmse_fn(obs, pred),
        r2      = r2_fn(obs, pred, mean(obs)),
        nrmse   = nrmse_fn(obs, pred)
      )
      summary_list[[length(summary_list) + 1]] <- scheme_summary
      cat(sprintf("  -> RMSE = %.4f | R2 = %.4f | NRMSE = %.4f\n",
                  scheme_summary$rmse, scheme_summary$r2,
                  scheme_summary$nrmse))
    }
  }

}, error = function(e) {
  cat("\n[UNPROTECTED ERROR]", conditionMessage(e), "\n")
  cat("  Call:", paste(deparse(conditionCall(e)), collapse = " "), "\n")
  cat("  Call stack (deepest last):\n")
  for (cl in tail(sys.calls(), 14)) {
    cat("    ", substr(paste(deparse(cl), collapse = " "), 1, 160), "\n")
  }
})


# ---- Summary -----------------------------------------------------------------

nrmse_summary <- bind_rows(summary_list)
write.csv(nrmse_summary, file.path(model_out_root, "nrmse_summary_all.csv"),
          row.names = FALSE)
cat("\n=== NRMSE summary (validation criterion) ===\n")
print(nrmse_summary)

cat("\nResults saved in:", model_out_root, "\n")
cat("Figures saved in:", fig_root_modelling, "\n")
