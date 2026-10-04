# ==============================================================================
# 13_run_tsne.R
#
# Random-forest pipeline: t-SNE of the covariate space
#
# Input  : learning dataset of the target frequency
#          (01_build_learning_dataset.R)
#          fold indices and out-of-fold predictions of the target frequency
#          and scheme (10_run_nested_cv_training.R: fold_indices.rds,
#          obs_pred.csv)
#
# Steps  : 1) t-SNE of a sample of the observations, on the standardised
#             numeric covariates (fod is left out: it is categorical and
#             t-SNE uses Euclidean distances)
#          2) map coloured by NASC
#          3) map coloured by FOD
#          4) map of the test points coloured by residual
#
# Output : in <model_out_root>/<freq>kHz/rf/<scheme>/
#            tsne_coords.csv
#          in <fig_root_modelling>/<freq>kHz/rf/<scheme>/
#            10a_tsne_nasc.png, 10b_tsne_fod.png, 10c_tsne_residual.png
#
# Requires the package `Rtsne`.
#
# t-SNE reduces the covariate space to two dimensions while preserving the
# LOCAL similarities between observations. It shows structure and groups; it
# is NOT a measure of importance or of direction of effect (it complements
# SHAP, it does not replace it). It is used here to check whether the
# prediction errors gather in a particular region of the covariate space.
#
# Usage  : optional target, to set in the console before source()
#            target_freq <- 38 ; target_scheme <- "naive_RS_80_20"
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), modelling_dir
for (f in c("00_model_config.R", "01_data_prep.R")) {
  source(file.path(modelling_dir, f))
}

if (!requireNamespace("Rtsne", quietly = TRUE)) {
  stop("Package `Rtsne` not installed: install.packages('Rtsne')")
}

if (!exists("target_freq"))   target_freq   <- default_target_freq
if (!exists("target_scheme")) target_scheme <- default_target_scheme

max_rows   <- 5000   # t-SNE is expensive: the observations are subsampled
perplexity <- 30
seed       <- 1

out_dir <- scheme_out_dir(target_freq, target_scheme)
fig_dir <- scheme_fig_dir(target_freq, target_scheme)

in_files <- file.path(out_dir, c("fold_indices.rds", "obs_pred.csv"))
if (!all(file.exists(in_files))) {
  stop("Missing files in ", out_dir, ": ",
       paste(basename(in_files[!file.exists(in_files)]), collapse = ", "),
       "\n  -> run 10_run_nested_cv_training.R first.")
}

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Read the data and the predictions ---------------------------------------

df       <- load_and_clean(target_freq)$df
obs_pred <- read.csv(file.path(out_dir, "obs_pred.csv"))

subtitle <- sprintf("%d kHz - %s", target_freq, target_scheme)


# ---- 1) t-SNE ----------------------------------------------------------------

set.seed(seed)
sample_idx <- if (nrow(df) > max_rows) {
  sample(seq_len(nrow(df)), max_rows)
} else {
  seq_len(nrow(df))
}

x_scaled <- scale(df[sample_idx, covariates_num, drop = FALSE])

cat(sprintf("t-SNE on %d rows, %d numeric covariates...\n",
            nrow(x_scaled), ncol(x_scaled)))
set.seed(seed)
tsne_res <- Rtsne::Rtsne(x_scaled, perplexity = perplexity,
                         check_duplicates = FALSE)

tsne_df <- tibble(
  tsne1 = tsne_res$Y[, 1],
  tsne2 = tsne_res$Y[, 2],
  NASC  = df$NASC[sample_idx],
  fod   = df$fod[sample_idx],
  lon   = df$lon[sample_idx],
  lat   = df$lat[sample_idx]
)
write.csv(tsne_df, file.path(out_dir, "tsne_coords.csv"), row.names = FALSE)


# ---- 2) Coloured by NASC -----------------------------------------------------

p_nasc <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = NASC)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_viridis_c() +
  labs(title = "t-SNE of the covariate space, coloured by NASC",
       subtitle = subtitle) +
  theme_pipeline
ggsave(file.path(fig_dir, "10a_tsne_nasc.png"), p_nasc,
       width = 8, height = 6, dpi = 150)


# ---- 3) Coloured by FOD ------------------------------------------------------

p_fod <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = fod)) +
  geom_point(alpha = 0.5, size = 0.8) +
  labs(title = "t-SNE of the covariate space, coloured by FOD",
       subtitle = subtitle) +
  theme_pipeline
ggsave(file.path(fig_dir, "10b_tsne_fod.png"), p_fod,
       width = 8, height = 6, dpi = 150)


# ---- 4) Test points coloured by residual -------------------------------------

# obs_pred.csv only holds the test rows (out-of-fold) and has no row
# identifier: the residuals are joined by lon / lat. This join is approximate
# (an observation tested in several folds keeps one residual only).
if (nrow(obs_pred) > 0) {
  residuals_df <- obs_pred %>%
    select(lon, lat, residual) %>%
    distinct(lon, lat, .keep_all = TRUE)
  tsne_test <- inner_join(tsne_df, residuals_df, by = c("lon", "lat"))

  if (nrow(tsne_test) > 10) {
    res_lim <- max(abs(tsne_test$residual))
    p_residual <- ggplot(tsne_test,
                         aes(x = tsne1, y = tsne2, color = residual)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red",
                            midpoint = 0, limits = c(-res_lim, res_lim)) +
      labs(title = "t-SNE, test points coloured by residual",
           subtitle = subtitle) +
      theme_pipeline
    ggsave(file.path(fig_dir, "10c_tsne_residual.png"), p_residual,
           width = 8, height = 6, dpi = 150)
  } else {
    cat("  [!] too few lon / lat matches for the residual plot -- skipped\n")
  }
}

cat("Results saved in:", out_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
