# ==============================================================================
# 02_folds.R
#
# Random-forest pipeline: cross-validation schemes
#
# A scheme is a list: scheme (name), data (the observations) and folds (a
# named list of folds, each with the row indices `train` and `test`).
#
# Spatially blocked schemes:
#   - no buffer: the training set is ALL the observations except the test
#     block itself (including the same area in the other years)
#   - a test block is restricted to ONE year (not the whole cell, all years
#     together). This is what limits the extrapolation: in training, the model
#     still sees observations of the SAME area in other years.
#
# Defines : build_naive_scheme(), assign_spatial_block(),
#           stratify_by_latitude(), build_blocked_scheme(),
#           build_all_schemes(), build_inner_folds()
# ==============================================================================


# ---- Naive Monte-Carlo CV ----------------------------------------------------

# Repeated random split of the observations.
#   n_folds    : number of splits
#   train_frac : share of the observations in the training set
build_naive_scheme <- function(df, n_folds = naive_n_folds,
                               train_frac = naive_train_frac, seed = 1) {
  set.seed(seed)
  n <- nrow(df)
  folds <- map(seq_len(n_folds), function(k) {
    train_idx <- sample(seq_len(n), size = floor(train_frac * n))
    list(train = train_idx, test = setdiff(seq_len(n), train_idx))
  })
  names(folds) <- paste0("naive_", seq_len(n_folds))
  list(scheme = "naive_RS_80_20", data = df, folds = folds)
}


# ---- Spatial blocking, without buffer, one year per test block ---------------

# Add the spatial block of each observation (columns block_x, block_y and
# spatial_block).
#   cellsize_km : size of the blocks along x and y (km)
assign_spatial_block <- function(df, cellsize_km) {
  df %>% mutate(
    block_x = floor(x_km / cellsize_km[1]),
    block_y = floor(y_km / cellsize_km[2]),
    spatial_block = paste(block_x, block_y, sep = "_")
  )
}

# Draw `n_target` candidates, balanced between the three latitude tertiles.
#   candidates_df : one row per candidate, with the column `lat_col`
stratify_by_latitude <- function(candidates_df, n_target,
                                 lat_col = "mean_lat", seed = 1) {
  set.seed(seed)
  if (nrow(candidates_df) <= n_target) return(candidates_df)

  breaks <- quantile(candidates_df[[lat_col]],
                     probs = seq(0, 1, length.out = 4), na.rm = TRUE)
  candidates_df$tertile <- cut(candidates_df[[lat_col]], breaks = breaks,
                               include.lowest = TRUE, labels = FALSE)

  n_per_tertile <- ceiling(n_target / 3)
  picked <- candidates_df %>%
    group_by(tertile) %>%
    group_modify(~ slice_sample(.x, n = min(n_per_tertile, nrow(.x)))) %>%
    ungroup()
  if (nrow(picked) > n_target) picked <- slice_sample(picked, n = n_target)
  picked
}

# Spatially blocked scheme: each fold tests one (block, year) and trains on
# all the other observations.
#   cellsize_km : size of the blocks along x and y (km)
#   label       : name of the resolution, used in the name of the scheme
build_blocked_scheme <- function(df, cellsize_km, label, seed = 1) {
  scheme_name <- paste0("blocked_spatial_", label)
  df_blocked  <- assign_spatial_block(df, cellsize_km)

  # Candidates: (block, year) with at least `block_min_n` observations
  candidates <- df_blocked %>%
    group_by(spatial_block, year) %>%
    summarise(n = n(), mean_lat = mean(lat), .groups = "drop") %>%
    filter(n >= block_min_n)

  n_blocks <- nrow(candidates)
  if (n_blocks == 0) {
    cat(sprintf("  [!] %s: no (block, year) with >= %d obs -- empty scheme\n",
                label, block_min_n))
    return(list(scheme = scheme_name, data = df_blocked, folds = list()))
  }

  # Number of folds: a share of the candidates, within fixed bounds
  n_target <- min(block_max_folds,
                  max(block_min_folds, round(block_folds_frac * n_blocks)))
  n_target <- min(n_target, n_blocks)
  selected <- stratify_by_latitude(candidates, n_target, seed = seed)

  cat(sprintf("  %s: %d (block, year) candidates -> %d folds (target: %d)\n",
              label, n_blocks, nrow(selected), n_target))

  folds <- map(seq_len(nrow(selected)), function(i) {
    test_idx <- which(df_blocked$spatial_block == selected$spatial_block[i] &
                        df_blocked$year == selected$year[i])
    # No buffer: everything else is in the training set
    train_idx <- setdiff(seq_len(nrow(df_blocked)), test_idx)
    list(train = train_idx, test = test_idx)
  })
  names(folds) <- paste0("s_", selected$spatial_block, "_", selected$year)

  list(scheme = scheme_name, data = df_blocked, folds = folds)
}

# All the schemes: the naive one, then one blocked scheme per resolution.
# Returns a list of schemes, named by scheme.
build_all_schemes <- function(df) {
  schemes <- list()
  schemes[["naive_RS_80_20"]] <- build_naive_scheme(df)
  for (res in spatial_resolutions) {
    schemes[[paste0("blocked_spatial_", res$label)]] <-
      build_blocked_scheme(df, res$cellsize_km, res$label)
  }
  schemes
}


# ---- Inner CV (nested): repeated random subsampling --------------------------

# Inner folds, drawn inside the training set of an outer fold. They are only
# used for the tuning and never see the outer test set: no leakage between the
# tuning and the final evaluation.
#   train_idx  : row indices of the outer training set
#   repeats    : number of inner splits
#   train_frac : share of `train_idx` in each inner training set
build_inner_folds <- function(train_idx, repeats = inner_cv_repeats,
                              train_frac = inner_cv_train_frac, seed = 1) {
  set.seed(seed)
  n <- length(train_idx)
  map(seq_len(repeats), function(k) {
    inner_train <- sample(train_idx, size = floor(train_frac * n))
    list(train = inner_train, test = setdiff(train_idx, inner_train))
  })
}
