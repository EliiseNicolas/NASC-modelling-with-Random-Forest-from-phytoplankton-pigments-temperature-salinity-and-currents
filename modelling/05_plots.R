# ==============================================================================
# 05_plots.R
#
# Random-forest pipeline: figures, aggregated over the folds and per fold
#
# RMSE and R2 are always written on the observed vs predicted plots, even
# though NRMSE is the criterion used to validate the model (see the summary
# table of 10_run_nested_cv_training.R).
#
# Defines : save_fig(), annotate_metrics(), plot_learning_curve_agg(),
#           plot_learning_curve_by_fold(), plot_obs_pred_agg(),
#           plot_obs_pred_by_fold(), plot_train_test_distribution(),
#           plot_residual_map(), plot_noise_robustness(),
#           plot_spatial_train_test(), plot_importance()
# ==============================================================================

set_colors <- c(Train = "steelblue", Test = "firebrick")

# Save a figure as PNG in `fig_dir`.
#   width, height : size (inches)
save_fig <- function(plot, fig_dir, filename, width = 8, height = 6) {
  ggsave(file.path(fig_dir, filename), plot,
         width = width, height = height, dpi = 150)
}

# Text of the metrics written on a plot.
annotate_metrics <- function(rmse, r2, nrmse = NULL) {
  label <- sprintf("RMSE = %.3f\nR\u00b2 = %.3f", rmse, r2)
  if (!is.null(nrmse)) label <- paste0(label, sprintf("\nNRMSE = %.3f", nrmse))
  label
}


# ---- 1) Learning curve -------------------------------------------------------

# Aggregated: mean +/- standard deviation over the folds.
#   lc_summary : `summary` of compute_learning_curve()
plot_learning_curve_agg <- function(lc_summary, subtitle = "", ylim = NULL) {
  p <- ggplot(lc_summary, aes(x = fraction)) +
    geom_ribbon(aes(ymin = mean_rmse_train - sd_rmse_train,
                    ymax = mean_rmse_train + sd_rmse_train),
                fill = set_colors[["Train"]], alpha = 0.2) +
    geom_ribbon(aes(ymin = mean_rmse_test - sd_rmse_test,
                    ymax = mean_rmse_test + sd_rmse_test),
                fill = set_colors[["Test"]], alpha = 0.2) +
    geom_line(aes(y = mean_rmse_train, color = "Train")) +
    geom_point(aes(y = mean_rmse_train, color = "Train")) +
    geom_line(aes(y = mean_rmse_test, color = "Test")) +
    geom_point(aes(y = mean_rmse_test, color = "Test")) +
    scale_color_manual(values = set_colors) +
    labs(title = "Learning curve (aggregated)", subtitle = subtitle,
         x = "Share of the training set", y = "RMSE", color = NULL) +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# Per fold.
#   lc_detail : `detail` of compute_learning_curve()
plot_learning_curve_by_fold <- function(lc_detail, subtitle = "", ylim = NULL) {
  long <- lc_detail %>%
    select(fold_id, fraction, rmse_train, rmse_test) %>%
    pivot_longer(c(rmse_train, rmse_test),
                 names_to = "type", values_to = "rmse") %>%
    mutate(type = ifelse(type == "rmse_train", "Train", "Test"))

  p <- ggplot(long, aes(x = fraction, y = rmse, color = type)) +
    geom_line() +
    geom_point(size = 0.8) +
    facet_wrap(~fold_id) +
    scale_color_manual(values = set_colors) +
    labs(title = "Learning curve, per fold", subtitle = subtitle,
         x = "Share of the training set", y = "RMSE", color = NULL) +
    theme_pipeline
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}


# ---- 2) Observed vs predicted NASC -------------------------------------------

# Aggregated: all the folds pooled, with RMSE, R2 and NRMSE.
#   obs_pred : `obs_pred` of run_nested_cv_scheme()
plot_obs_pred_agg <- function(obs_pred, subtitle = "", axis_limits = NULL) {
  rmse  <- rmse_fn(obs_pred$obs, obs_pred$pred)
  r2    <- r2_fn(obs_pred$obs, obs_pred$pred, mean(obs_pred$obs))
  nrmse <- nrmse_fn(obs_pred$obs, obs_pred$pred)

  p <- ggplot(obs_pred, aes(x = obs, y = pred, color = fold_id)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(linetype = "dashed") +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = annotate_metrics(rmse, r2, nrmse), size = 3) +
    labs(title = "Observed vs predicted NASC (out-of-fold, aggregated)",
         subtitle = subtitle,
         x = "Observed NASC (log10)", y = "Predicted NASC (log10)",
         color = "Fold") +
    theme_pipeline
  if (!is.null(axis_limits)) {
    p <- p + coord_cartesian(xlim = axis_limits, ylim = axis_limits)
  }
  p
}

# Per fold.
#   metrics : `metrics` of run_nested_cv_scheme()
plot_obs_pred_by_fold <- function(obs_pred, metrics, subtitle = "",
                                  axis_limits = NULL) {
  labels_df <- metrics %>%
    mutate(label = annotate_metrics(rmse_test, r2_test, nrmse_test))

  p <- ggplot(obs_pred, aes(x = obs, y = pred)) +
    geom_point(alpha = 0.4, size = 0.6) +
    geom_abline(linetype = "dashed", color = "red") +
    geom_text(data = labels_df, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.05, vjust = 1.1, size = 2.5, inherit.aes = FALSE) +
    facet_wrap(~fold_id) +
    labs(title = "Observed vs predicted NASC, per fold", subtitle = subtitle,
         x = "Observed NASC (log10)", y = "Predicted NASC (log10)") +
    theme_pipeline
  if (!is.null(axis_limits)) {
    p <- p + coord_cartesian(xlim = axis_limits, ylim = axis_limits)
  }
  p
}


# ---- 3) Distribution of the NASC, training vs test set -----------------------

#   by_fold : FALSE = all the folds pooled; TRUE = one panel per fold
plot_train_test_distribution <- function(df, folds, subtitle = "",
                                         by_fold = FALSE) {
  if (!by_fold) {
    all_train <- unique(unlist(lapply(folds, `[[`, "train")))
    all_test  <- unique(unlist(lapply(folds, `[[`, "test")))
    long <- bind_rows(
      tibble(NASC = df$NASC[all_train], set = "Train"),
      tibble(NASC = df$NASC[all_test],  set = "Test")
    )
    return(
      ggplot(long, aes(x = NASC, fill = set)) +
        geom_histogram(position = "identity", alpha = 0.5, bins = 40) +
        labs(title = "NASC distribution, train vs test (aggregated)",
             subtitle = subtitle, x = "NASC (log10)", y = "Count",
             fill = NULL) +
        theme_pipeline
    )
  }

  long <- imap_dfr(folds, function(f, fid) {
    bind_rows(
      tibble(fold_id = fid, NASC = df$NASC[f$train], set = "Train"),
      tibble(fold_id = fid, NASC = df$NASC[f$test],  set = "Test")
    )
  })
  ggplot(long, aes(x = NASC, fill = set)) +
    geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
    facet_wrap(~fold_id, scales = "free_y") +
    labs(title = "NASC distribution, train vs test, per fold",
         subtitle = subtitle, x = "NASC (log10)", y = "Count", fill = NULL) +
    theme_pipeline
}


# ---- 4) Map of the residuals -------------------------------------------------

#   color_limits : limits of the colour scale; NULL = symmetric around 0
#   by_fold      : FALSE = all the folds pooled; TRUE = one panel per fold
plot_residual_map <- function(obs_pred, subtitle = "", color_limits = NULL,
                              by_fold = FALSE) {
  max_abs <- max(abs(obs_pred$residual))
  limits  <- color_limits %||% c(-max_abs, max_abs)

  p <- ggplot(obs_pred, aes(x = lon, y = lat, color = residual)) +
    geom_point(size = 1) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, limits = limits) +
    coord_quickmap() +
    labs(title = "Map of the residuals (observed - predicted)",
         subtitle = subtitle, x = "Longitude", y = "Latitude",
         color = "Residual") +
    theme_pipeline
  if (by_fold) p <- p + facet_wrap(~fold_id)
  p
}


# ---- 5) Robustness to noise --------------------------------------------------

#   noise_df : output of compute_noise_robustness()
#   relative : TRUE = RMSE in % of the RMSE without noise; FALSE = RMSE
#   by_fold  : FALSE = mean +/- standard deviation over the folds;
#              TRUE = one curve per fold
plot_noise_robustness <- function(noise_df, relative = TRUE, subtitle = "",
                                  by_fold = FALSE) {
  y_col <- if (relative) "rmse_relative" else "rmse"
  y_lab <- if (relative) {
    "Relative RMSE (% of the RMSE without noise)"
  } else {
    "RMSE"
  }
  kind  <- if (relative) "relative" else "absolute"
  x_lab <- "Noise level (share of the standard deviation)"

  if (!by_fold) {
    agg <- noise_df %>%
      group_by(noise_level) %>%
      summarise(mean_y = mean(.data[[y_col]]), sd_y = sd(.data[[y_col]]),
                .groups = "drop")
    return(
      ggplot(agg, aes(x = noise_level, y = mean_y)) +
        geom_ribbon(aes(ymin = mean_y - sd_y, ymax = mean_y + sd_y),
                    alpha = 0.2) +
        geom_line() +
        geom_point() +
        labs(title = sprintf("Robustness to Gaussian noise (%s)", kind),
             subtitle = subtitle, x = x_lab, y = y_lab) +
        theme_pipeline
    )
  }

  ggplot(noise_df,
         aes(x = noise_level, y = .data[[y_col]], color = fold_id)) +
    geom_line() +
    geom_point(size = 0.8) +
    labs(title = sprintf("Robustness to Gaussian noise (%s), per fold", kind),
         subtitle = subtitle, x = x_lab, y = y_lab, color = "Fold") +
    theme_pipeline
}


# ---- 6) Location of the training and test sets, per fold ---------------------

#   scheme : a scheme (see 02_folds.R)
plot_spatial_train_test <- function(scheme, subtitle = "") {
  status_levels <- c("Not used", "Train", "Test")

  status_by_fold <- imap_dfr(scheme$folds, function(f, fid) {
    status <- scheme$data %>% select(lon, lat)
    status$status <- "Not used"
    status$status[f$train] <- "Train"
    status$status[f$test]  <- "Test"
    status$fold_id <- fid
    status
  }) %>%
    mutate(status = factor(status, levels = status_levels)) %>%
    arrange(fold_id, status)   # test points drawn last, on top

  ggplot(status_by_fold, aes(x = lon, y = lat, color = status)) +
    geom_point(size = 0.1, alpha = 0.6) +
    scale_color_manual(
      values = c("Test" = "red", "Train" = "steelblue", "Not used" = "grey85")
    ) +
    coord_quickmap() +
    facet_wrap(~fold_id) +
    labs(title = "Location of the training and test sets, per fold",
         subtitle = subtitle, x = "Longitude", y = "Latitude", color = NULL) +
    theme_pipeline
}


# ---- 7) Variable importance (permutation) ------------------------------------

# Mean +/- standard deviation over the folds.
#   importance_all : `importance` of run_nested_cv_scheme()
plot_importance <- function(importance_all, subtitle = "") {
  summary_df <- importance_all %>%
    group_by(variable) %>%
    summarise(mean_imp = mean(importance), sd_imp = sd(importance),
              .groups = "drop") %>%
    arrange(desc(mean_imp))

  ggplot(summary_df, aes(x = reorder(variable, mean_imp), y = mean_imp)) +
    geom_col(fill = "orange") +
    geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp),
                  width = 0.2) +
    coord_flip() +
    labs(title = "Mean variable importance (+/- sd over the folds)",
         subtitle = subtitle, x = NULL, y = "Importance (permutation)") +
    theme_pipeline
}
