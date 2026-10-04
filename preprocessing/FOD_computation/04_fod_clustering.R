# ==============================================================================
# 04_fod_clustering.R
#
# Functional Oceanographic Domains (FOD), step 3/3: clustering
#
# Input  : PCA scores (03_fod_pca.R), grid and smoothed profiles
#          (02_fod_bspline.R), all in <fod_dir>
#
# Steps  : 1) Gaussian mixture model (Mclust) on the PCA scores
#          2) renumber the clusters from South to North
#          3) soft classification: profiles whose highest probability is below
#             `prob_threshold` are transition profiles
#          4) transition classes, from the two most probable clusters
#          5) put the classification back on the (lon, lat, time) grid
#          6) mean temperature / salinity profile of each class
#          7) one map per day
#
# Classes: 1-6  = clusters C1..C6 (South -> North)
#          7-13 = transitions between two clusters (see `fod_classes`)
#          0    = transition between two clusters that are not in the list
#
# Output : in <fod_dir> (clusters numbered South -> North in every file)
#            mclust_model.rds                 fitted model (Mclust numbering)
#            cluster_rename.rds               Mclust number -> final number
#            cluster.rds, max_probability.rds, cluster_probabilities.rds,
#            cluster_soft.rds                 per profile
#            cluster_transition_renamed.rds   final class per profile
#            cluster_transition_flat_renamed.rds, cluster_transition_map_renamed.rds
#                                             final class on the grid
#            cluster_probability_maps.rds     probabilities on the grid
#            transition_information.rds, transition_summary_renamed.rds,
#            transition_codes.rds, fod_classes.rds, nclust.rds,
#            softclass_threshold.rds
#            {mean,q1,q3}_{temperature,salinity}_clusters_transitions.rds
#          in <fig_dir>
#            mclust_bic.png, temperature_clusters_transitions.png,
#            salinity_clusters_transitions.png
#          in <map_dir>
#            cluster_transition_map_renamed_<YYYYMMDD>.png
# ==============================================================================

library(mclust)
library(ggplot2)
library(sf)
library(rnaturalearth)   # also requires the `rnaturalearthdata` package


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories

fig_dir <- file.path(fig_root, "FOD_computation", "04_fod_clustering")
map_dir <- file.path(fig_dir, "fod_maps")

# Gaussian mixture model
n_clusters_tested <- 4:6     # the number of clusters is chosen by BIC
gmm_model         <- "VVV"
prob_threshold    <- 0.75    # below: transition profile
seed              <- 1

# FALSE: reuse <fod_dir>/mclust_model.rds if it exists (same data, same model).
# TRUE : fit the model again. The Mclust numbering may then change:
#        `cluster_rename` below must be checked.
refit_gmm <- FALSE

# Renumbering of the clusters from South to North.
# Position = Mclust number, value = final number.
# !! Specific to the fitted model: check it after any new fit.
cluster_rename <- c(1L, 2L, 4L, 3L, 6L, 5L)

# Transition classes (final numbering): pair of clusters -> class code
transition_codes <- c(
  "1-2" = 7L, "1-3" = 8L, "2-3" = 9L, "3-4" = 10L,
  "4-6" = 11L, "4-5" = 12L, "5-6" = 13L
)

# Classes shown in the figures, ordered from South to North
fod_classes <- data.frame(
  code  = c(1L, 7L, 2L, 8L, 9L, 3L, 10L, 4L, 11L, 12L, 5L, 13L, 6L),
  label = c("C1", "T1-2", "C2", "T1-3", "T2-3", "C3", "T3-4",
            "C4", "T4-6", "T4-5", "C5", "T5-6", "C6"),
  col   = c("#2166AC", "#3B73B9", "#67A9CF", "#3FA7B5", "#45B97C", "#1A9850",
            "#C46A00", "#A6D96A", "#E85D04", "#C9184A", "#FDAE61", "#8F1D3F",
            "#D73027"),
  stringsAsFactors = FALSE
)

dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Two most probable clusters of each row of a probability matrix.
# Returns a two-column matrix (first, second).
top2_clusters <- function(prob) {
  n      <- nrow(prob)
  first  <- rep(1L, n)
  second <- rep(NA_integer_, n)
  p1     <- prob[, 1]
  p2     <- rep(-Inf, n)

  for (k in seq_len(ncol(prob))[-1]) {
    x <- prob[, k]

    new_first  <- x > p1
    new_second <- !new_first & x > p2

    second[new_first] <- first[new_first]
    p2[new_first]     <- p1[new_first]
    first[new_first]  <- k
    p1[new_first]     <- x[new_first]

    second[new_second] <- k
    p2[new_second]     <- x[new_second]
  }
  cbind(first = first, second = second)
}

# Put a per-profile vector back on the (lon, lat, time) grid.
# `flat = TRUE` returns the flattened grid instead of the array.
to_grid <- function(x, mask, dims, flat = FALSE) {
  x_flat <- rep(NA, length(mask))
  x_flat[mask] <- x
  if (flat) x_flat else array(x_flat, dim = dims)
}

# Mean, first and third quartile profile of each class.
# Returns three lists (mean, q1, q3) named by class code; classes without any
# profile are skipped.
class_profile_stats <- function(profiles, classes, codes) {
  stats <- list(mean = list(), q1 = list(), q3 = list())

  for (code in codes) {
    ind <- classes == code
    if (!any(ind)) next

    x   <- profiles[ind, , drop = FALSE]
    q   <- apply(x, 2, quantile, probs = c(0.25, 0.75), na.rm = TRUE)
    key <- as.character(code)

    stats$mean[[key]] <- colMeans(x, na.rm = TRUE)
    stats$q1[[key]]   <- q[1, ]
    stats$q3[[key]]   <- q[2, ]
  }
  stats
}

# Mean profile (line) and interquartile range (band) of each class
plot_class_profiles <- function(stats, depth, classes, xlab, out_file) {
  classes <- classes[as.character(classes$code) %in% names(stats$mean), ]
  keys    <- as.character(classes$code)

  png(out_file, width = 2000, height = 2800, res = 300)

  plot(
    NULL,
    xlim = range(unlist(stats$q1[keys]), unlist(stats$q3[keys])),
    ylim = rev(range(depth)),
    xlab = xlab, ylab = "Depth (m)"
  )
  for (i in seq_along(keys)) {
    polygon(
      x = c(stats$q1[[keys[i]]], rev(stats$q3[[keys[i]]])),
      y = c(depth, rev(depth)),
      col = adjustcolor(classes$col[i], alpha.f = 0.25), border = NA
    )
    lines(stats$mean[[keys[i]]], depth, col = classes$col[i], lwd = 3)
  }
  legend("bottomright", legend = classes$label, col = classes$col,
         lwd = 3, bty = "n", cex = 0.85, ncol = 2)

  dev.off()
}


# ---- Load the inputs ---------------------------------------------------------

scores      <- readRDS(file.path(fod_dir, "pca_scores.rds"))   # (n x nharm)
mask_common <- readRDS(file.path(fod_dir, "mask_common.rds"))
lon         <- as.vector(readRDS(file.path(fod_dir, "lon.rds")))
lat         <- as.vector(readRDS(file.path(fod_dir, "lat.rds")))
depth       <- as.vector(readRDS(file.path(fod_dir, "depth.rds")))
time        <- readRDS(file.path(fod_dir, "time.rds"))

grid_dims <- c(length(lon), length(lat), length(time))

stopifnot(
  length(mask_common) == prod(grid_dims),
  sum(mask_common) == nrow(scores)
)


# ---- 1) Gaussian mixture model -----------------------------------------------

model_file <- file.path(fod_dir, "mclust_model.rds")

if (refit_gmm || !file.exists(model_file)) {
  # Mclust initialises the EM on a random subsample when the data are large
  set.seed(seed)
  gmm <- Mclust(scores, G = n_clusters_tested, modelNames = gmm_model)
  saveRDS(gmm, model_file)
} else {
  cat("Reusing the model saved in", model_file, "\n")
  gmm <- readRDS(model_file)
}

stopifnot(gmm$n == nrow(scores), gmm$d == ncol(scores))
rm(scores)

nclust <- gmm$G
cat("Number of clusters chosen by BIC:", nclust, "\n")

if (nclust != length(cluster_rename) ||
    !setequal(cluster_rename, seq_len(nclust))) {
  stop("`cluster_rename` does not match the ", nclust, " clusters of the model: ",
       "update `cluster_rename`, `transition_codes` and `fod_classes`.")
}

png(file.path(fig_dir, "mclust_bic.png"), width = 2000, height = 1600, res = 300)
plot(gmm, what = "BIC")
dev.off()


# ---- 2) Renumber the clusters from South to North ----------------------------

# Column j of `cluster_prob` = probability of the final cluster j
cluster_prob <- gmm$z[, order(cluster_rename), drop = FALSE]
dimnames(cluster_prob) <- NULL
cluster      <- cluster_rename[gmm$classification]


# ---- 3) Soft classification --------------------------------------------------

# Highest probability of each profile
max_prob <- cluster_prob[, 1]
for (k in seq_len(nclust)[-1]) max_prob <- pmax(max_prob, cluster_prob[, k])

cluster_soft <- cluster
cluster_soft[max_prob < prob_threshold] <- 0L   # 0 = transition profile


# ---- 4) Transition classes ---------------------------------------------------

idx_transition  <- which(cluster_soft == 0)
prob_transition <- cluster_prob[idx_transition, , drop = FALSE]

top2 <- top2_clusters(prob_transition)
pair <- paste(pmin(top2[, 1], top2[, 2]), pmax(top2[, 1], top2[, 2]), sep = "-")

# Pairs that are not in `transition_codes` keep the code 0
transition_class <- unname(transition_codes[pair])
transition_class[is.na(transition_class)] <- 0L

fod_class <- cluster_soft
fod_class[idx_transition] <- transition_class

i_row <- seq_len(nrow(prob_transition))
transition_info <- data.frame(
  profile            = idx_transition,
  cluster_1          = top2[, 1],
  cluster_2          = top2[, 2],
  prob_1             = prob_transition[cbind(i_row, top2[, 1])],
  prob_2             = prob_transition[cbind(i_row, top2[, 2])],
  pair               = pair,
  transition_cluster = transition_class,
  stringsAsFactors   = FALSE
)

# Number of profiles per class
class_counts <- data.frame(
  code  = c(0L, fod_classes$code),
  label = c("Unidentified transition", fod_classes$label),
  stringsAsFactors = FALSE
)
class_counts$n <- as.integer(table(factor(fod_class, levels = class_counts$code)))
print(class_counts)

transition_summary <- data.frame(
  code       = unname(transition_codes),
  transition = names(transition_codes),
  n          = as.integer(table(factor(fod_class, levels = transition_codes))),
  stringsAsFactors = FALSE
)

saveRDS(cluster_rename, file.path(fod_dir, "cluster_rename.rds"))
saveRDS(cluster, file.path(fod_dir, "cluster.rds"))
saveRDS(max_prob, file.path(fod_dir, "max_probability.rds"))
saveRDS(cluster_prob, file.path(fod_dir, "cluster_probabilities.rds"))
saveRDS(cluster_soft, file.path(fod_dir, "cluster_soft.rds"))
saveRDS(fod_class, file.path(fod_dir, "cluster_transition_renamed.rds"))
saveRDS(transition_info, file.path(fod_dir, "transition_information.rds"))
saveRDS(transition_summary,
        file.path(fod_dir, "transition_summary_renamed.rds"))
saveRDS(transition_codes, file.path(fod_dir, "transition_codes.rds"))
saveRDS(fod_classes, file.path(fod_dir, "fod_classes.rds"))
saveRDS(nclust, file.path(fod_dir, "nclust.rds"))
saveRDS(prob_threshold, file.path(fod_dir, "softclass_threshold.rds"))

rm(gmm, prob_transition, transition_info, top2, pair)
gc()


# ---- 5) Back on the (lon, lat, time) grid ------------------------------------

fod_map <- to_grid(fod_class, mask_common, grid_dims)

saveRDS(to_grid(fod_class, mask_common, grid_dims, flat = TRUE),
        file.path(fod_dir, "cluster_transition_flat_renamed.rds"))
saveRDS(fod_map, file.path(fod_dir, "cluster_transition_map_renamed.rds"))

prob_maps <- lapply(seq_len(nclust), function(cl) {
  to_grid(cluster_prob[, cl], mask_common, grid_dims)
})
names(prob_maps) <- paste0("cluster_", seq_len(nclust))

saveRDS(prob_maps, file.path(fod_dir, "cluster_probability_maps.rds"))
rm(prob_maps, cluster_prob)
gc()


# ---- 6) Mean temperature / salinity profile of each class --------------------

temp_rec <- readRDS(file.path(fod_dir, "temperature_profiles.rds"))
stopifnot(nrow(temp_rec) == length(fod_class))
stats_temp <- class_profile_stats(temp_rec, fod_class, fod_classes$code)
rm(temp_rec)
gc()

sal_rec <- readRDS(file.path(fod_dir, "salinity_profiles.rds"))
stopifnot(nrow(sal_rec) == length(fod_class))
stats_sal <- class_profile_stats(sal_rec, fod_class, fod_classes$code)
rm(sal_rec)
gc()

for (s in c("mean", "q1", "q3")) {
  saveRDS(stats_temp[[s]], file.path(
    fod_dir, paste0(s, "_temperature_clusters_transitions.rds")
  ))
  saveRDS(stats_sal[[s]], file.path(
    fod_dir, paste0(s, "_salinity_clusters_transitions.rds")
  ))
}

plot_class_profiles(
  stats_temp, depth, fod_classes, xlab = "Temperature (°C)",
  out_file = file.path(fig_dir, "temperature_clusters_transitions.png")
)
plot_class_profiles(
  stats_sal, depth, fod_classes, xlab = "Salinity (PSU)",
  out_file = file.path(fig_dir, "salinity_clusters_transitions.png")
)


# ---- 7) One map per day ------------------------------------------------------

world   <- ne_countries(scale = "medium", returnclass = "sf")
map_df  <- expand.grid(lon = lon, lat = lat)
fod_pal <- setNames(fod_classes$col, fod_classes$label)

for (i in seq_along(time)) {
  # Class 0 (unidentified transition) is not in `fod_classes`: it becomes NA
  # and is drawn in white, like the cells without data.
  map_df$fod <- factor(as.vector(fod_map[, , i]),
                       levels = fod_classes$code, labels = fod_classes$label)

  p <- ggplot() +
    geom_raster(data = map_df, aes(x = lon, y = lat, fill = fod)) +
    geom_sf(data = world, fill = "grey85", color = "grey50", linewidth = 0.2) +
    coord_sf(xlim = range(lon), ylim = range(lat)) +
    scale_fill_manual(
      values = fod_pal, name = "FOD", na.value = "white", drop = FALSE,
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(
      x = "Longitude", y = "Latitude",
      title = paste("Functional Oceanographic Domains (clusters) -",
                    format(time[i], "%Y-%m-%d"))
    ) +
    theme_minimal()

  ggsave(
    filename = file.path(
      map_dir,
      paste0("cluster_transition_map_renamed_", format(time[i], "%Y%m%d"), ".png")
    ),
    plot = p, width = 8, height = 8, units = "in", dpi = 300
  )
}

cat("Results saved in:", fod_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
cat("Maps saved in:", map_dir, "\n")
