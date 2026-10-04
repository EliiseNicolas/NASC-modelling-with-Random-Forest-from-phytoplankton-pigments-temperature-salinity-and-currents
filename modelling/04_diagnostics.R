# ==============================================================================
# 04_diagnostics.R
#
# Random-forest pipeline: nested CV (inner tuning without leakage), metrics,
# learning curve and robustness to noise
#
# Defines : tune_inner(), run_nested_cv_scheme(), compute_learning_curve(),
#           compute_noise_robustness()
# ==============================================================================

# Inner tuning of ONE outer fold: grid search using only subsamples of the
# outer training set. The outer test set is never seen here.
# Criterion: mean RMSE over the inner splits.
#   df        : the observations
#   train_idx : row indices of the outer training set
#   grid      : tuning grid (one row per combination of hyperparameters)
# Returns a list: best_params and tuning_results (idx, mean_rmse, sd_rmse).
tune_inner <- function(df, train_idx, grid = rf_tuning_grid) {
  inner_folds <- build_inner_folds(train_idx)

  results <- map_dfr(seq_len(nrow(grid)), function(i) {
    params <- as.list(grid[i, ])
    rmses <- map_dbl(inner_folds, function(f) {
      model <- fit_rf(df[f$train, ], params)
      rmse_fn(df$NASC[f$test], predict_rf(model, df[f$test, ]))
    })
    tibble(idx = i, mean_rmse = mean(rmses), sd_rmse = stats::sd(rmses))
  })

  best_idx <- results$idx[which.min(results$mean_rmse)]
  list(best_params = as.list(grid[best_idx, ]), tuning_results = results)
}

# Nested CV of ONE scheme. For each outer fold: inner tuning, refit on the
# whole outer training set with the best hyperparameters, evaluation on the
# outer test set.
#   scheme : a scheme (see 02_folds.R)
# Returns NULL for an empty scheme, otherwise a list:
#   models     : final model of each fold
#   metrics    : one row per fold (sizes, RMSE, R2, NRMSE, best hyperparameters)
#   obs_pred   : observed and predicted NASC of the test sets
#   importance : permutation importance, per fold
#   tuning     : results of the inner tuning, per fold
run_nested_cv_scheme <- function(scheme) {
  df    <- scheme$data
  folds <- scheme$folds
  if (length(folds) == 0) return(NULL)

  models          <- list()
  metrics_list    <- list()
  obs_pred_list   <- list()
  importance_list <- list()
  tuning_list     <- list()

  for (fid in names(folds)) {
    f <- folds[[fid]]
    cat(sprintf("    fold %s (n_train = %d, n_test = %d): inner tuning...\n",
                fid, length(f$train), length(f$test)))

    tuned       <- tune_inner(df, f$train)
    final_model <- fit_rf(df[f$train, ], tuned$best_params)

    pred_test  <- predict_rf(final_model, df[f$test, ])
    pred_train <- predict_rf(final_model, df[f$train, ])
    obs_test   <- df$NASC[f$test]
    obs_train  <- df$NASC[f$train]

    metrics_list[[fid]] <- tibble(
      fold_id            = fid,
      n_train            = length(f$train),
      n_test             = length(f$test),
      rmse_train         = rmse_fn(obs_train, pred_train),
      rmse_test          = rmse_fn(obs_test, pred_test),
      r2_test            = r2_fn(obs_test, pred_test, mean(obs_train)),
      nrmse_test         = nrmse_fn(obs_test, pred_test),
      best_mtry          = tuned$best_params$mtry,
      best_min_node_size = tuned$best_params$min.node.size,
      best_num_trees     = tuned$best_params$num.trees
    )

    obs_pred_list[[fid]] <- tibble(
      fold_id = fid, obs = obs_test, pred = pred_test,
      residual = obs_test - pred_test,
      lon = df$lon[f$test], lat = df$lat[f$test], set = "test"
    )

    importance_list[[fid]] <- get_importance_rf(final_model) %>%
      mutate(fold_id = fid)
    tuning_list[[fid]] <- tuned$tuning_results %>% mutate(fold_id = fid)
    models[[fid]] <- final_model
  }

  list(
    models     = models,
    metrics    = bind_rows(metrics_list),
    obs_pred   = bind_rows(obs_pred_list),
    importance = bind_rows(importance_list),
    tuning     = bind_rows(tuning_list)
  )
}

# Learning curve, per fold: the model is refitted on growing shares of the
# training set, with the hyperparameters already selected for the fold (no
# new tuning).
#   best_params_by_fold : list of hyperparameters, named by fold
#   fractions           : shares of the training set
# Returns a list: detail (one row per fold and share) and summary (mean and
# standard deviation over the folds, per share).
compute_learning_curve <- function(df, folds, best_params_by_fold,
                                   fractions = c(0.1, seq(0.2, 1, by = 0.2)),
                                   seed = 1) {
  set.seed(seed)
  detail <- imap_dfr(folds, function(f, fid) {
    params <- best_params_by_fold[[fid]]
    map_dfr(fractions, function(frac) {
      n_sub   <- max(20, round(frac * length(f$train)))
      sub_idx <- sample(f$train, min(n_sub, length(f$train)))
      model   <- fit_rf(df[sub_idx, ], params)
      pred_train <- predict_rf(model, df[sub_idx, ])
      pred_test  <- predict_rf(model, df[f$test, ])
      tibble(
        fold_id    = fid,
        fraction   = frac,
        rmse_train = rmse_fn(df$NASC[sub_idx], pred_train),
        rmse_test  = rmse_fn(df$NASC[f$test], pred_test)
      )
    })
  })

  summary <- detail %>%
    group_by(fraction) %>%
    summarise(mean_rmse_train = mean(rmse_train),
              sd_rmse_train   = sd(rmse_train),
              mean_rmse_test  = mean(rmse_test),
              sd_rmse_test    = sd(rmse_test),
              .groups = "drop")

  list(detail = detail, summary = summary)
}

# Robustness to Gaussian noise, per fold: the model already fitted for the
# fold predicts its test set, with noise added to the numeric covariates.
#   models : final model of each fold, named by fold
#   levels : standard deviation of the noise, as a share of the standard
#            deviation of each covariate
# Returns one row per fold and level: rmse and rmse_relative (% of the RMSE
# without noise).
compute_noise_robustness <- function(df, folds, models, levels = noise_levels,
                                     seed = 1) {
  set.seed(seed)
  sds <- sapply(covariates_num, function(v) stats::sd(df[[v]], na.rm = TRUE))

  imap_dfr(folds, function(f, fid) {
    model   <- models[[fid]]
    test_df <- df[f$test, ]
    obs     <- test_df$NASC
    map_dfr(levels, function(level) {
      noisy <- test_df
      for (v in covariates_num) {
        noisy[[v]] <- noisy[[v]] + rnorm(nrow(noisy), 0, level * sds[[v]])
      }
      tibble(fold_id = fid, noise_level = level,
             rmse = rmse_fn(obs, predict_rf(model, noisy)))
    })
  }) %>%
    group_by(fold_id) %>%
    mutate(rmse_relative = rmse / rmse[noise_level == 0] * 100) %>%
    ungroup()
}
