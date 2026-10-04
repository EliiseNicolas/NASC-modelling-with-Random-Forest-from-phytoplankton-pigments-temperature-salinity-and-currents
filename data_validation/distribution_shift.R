# ==============================================================================
# distribution_shift.R
#
# Data validation: distribution shift between years and between FOD zones
#
# Input  : merged dataset (NASC, pigments, FTLE and FOD per ESU) at `freq` kHz
#          <dataset_dir>/
#            NASC_per_esu_pig_ftle_fod_<years_tag>_transect_<freq>kHz.rds
#          restricted to the ESU with day == `day_code`
#
# Steps  : how much the covariates (FTLE, pigments) and the NASC differ from
#          one year to another and from one FOD zone to another
#          1) preparation: partitions (year, zone) and daily blocks
#          2) univariate shift: overlaid densities, then, for each variable
#             and each level (one year or one zone against all the others),
#             Kolmogorov-Smirnov statistic, overlap of the densities and share
#             of the observations outside the range of the others
#          3) multivariate shift (adversarial validation): a random forest
#             tries to recognise the level from the covariates. An AUC close
#             to 0.5 means similar distributions, close to 1 a strong shift;
#             the permutation importance shows which variables carry it.
#
# Output : in <fig_dir> (<part> = year or zone)
#            figures  densities_<part>.png, ks_<part>.png, overlap_<part>.png,
#                     out_of_range_<part>.png, adversarial_auc.png,
#                     adversarial_importance_<part>.png
#            tables   univariate_<part>.csv, adversarial_auc.csv,
#                     adversarial_importance.csv
#            distribution_shift_results.rds   all the results, to reload them
#                                             without computing them again
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ranger)

options(scipen = 999)
set.seed(42)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, dataset_file()

freq     <- 38   # kHz
day_code <- 3    # ESU kept: 3 = day (1 = night)

in_file <- dataset_file(freq)
fig_dir <- file.path(fig_root_validation, "distribution_shift")

# Covariates. `fod` is not one of them: it is the partition.
pig_ratios   <- c("Per_totpig", "But_totpig", "Fuco_totpig", "Hex_totpig",
                  "Allo_totpig", "Zea_totpig", "Chlb_totpig", "DvChla_totpig")
covars       <- c("ftle", "Chla", "total_pig", pig_ratios)
covars_phys  <- c("ftle")   # without the pigments (fewer missing values)

# Univariate shift: minimum number of values in each of the two groups
min_group_size <- 20

# Adversarial validation
n_max     <- 20000   # maximum number of observations drawn per class
n_folds   <- 5       # folds of the cross-validation
num_trees <- 300     # trees of each forest

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Print a figure and save it as PNG in <fig_dir>.
#   name          : file name, without extension
#   width, height : size (cm)
save_plot <- function(p, name, width = 25, height = 18) {
  print(p)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p,
         width = width, height = height, units = "cm", dpi = 300)
}

# Save a table as CSV in <fig_dir>.
#   name : file name, without extension
save_table <- function(tab, name) {
  write.csv(tab, file.path(fig_dir, paste0(name, ".csv")), row.names = FALSE)
}

# Overlaid densities of each variable, one curve per level.
#   long : long table (year, zone, variable, value)
#   part : partition, "year" or "zone"
# Returns a ggplot with one panel per variable.
density_plot <- function(long, part) {
  ggplot(filter(long, !is.na(.data[[part]])),
         aes(x = value, colour = .data[[part]])) +
    geom_density(linewidth = 0.5) +
    facet_wrap(~ variable, scales = "free") +
    labs(x = NULL, y = "Density", colour = part) +
    theme_bw()
}

# Overlap coefficient of two densities: their common area (1 = identical,
# 0 = disjoint).
#   a, b : vectors of values
#   n    : number of points where the densities are evaluated
overlap <- function(a, b, n = 512) {
  rng <- range(c(a, b))
  if (diff(rng) == 0) return(1)
  da <- density(a, from = rng[1], to = rng[2], n = n)
  db <- density(b, from = rng[1], to = rng[2], n = n)
  sum(pmin(da$y, db$y)) * diff(da$x[1:2])
}

# Effect sizes per variable and per level (level against the rest).
# Comparisons where one of the two groups has fewer than `min_size` values are
# skipped.
#   long : long table (year, zone, variable, value)
#   part : partition, "year" or "zone"
# Returns a tibble with one row per variable and level: n (size of the level),
# D_ks (Kolmogorov-Smirnov statistic), overlap (see overlap()) and
# pct_out_of_range (% of values outside the min-max of the rest).
univariate_shift <- function(long, part, min_size = min_group_size) {
  long %>%
    group_by(variable) %>%
    group_modify(function(df, key) {
      df <- df[!is.na(df[[part]]), ]
      g  <- droplevels(df[[part]])
      bind_rows(lapply(levels(g), function(lev) {
        a <- df$value[g == lev]
        b <- df$value[g != lev]
        if (length(a) < min_size || length(b) < min_size) return(NULL)
        tibble(
          level            = lev,
          n                = length(a),
          D_ks             = unname(suppressWarnings(ks.test(a, b)$statistic)),
          overlap          = overlap(a, b),
          pct_out_of_range = 100 * mean(a < min(b) | a > max(b))
        )
      }))
    }) %>%
    ungroup()
}

# Heat map of a metric, variables in rows and levels in columns.
#   tab    : table with the columns level, variable and the metric
#   metric : name of the column to show
heatmap_plot <- function(tab, metric, title) {
  ggplot(tab, aes(x = level, y = variable, fill = .data[[metric]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(.data[[metric]], 2)), size = 3) +
    scale_fill_viridis_c(name = metric) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal()
}

# AUC from the rank formula (Mann-Whitney).
#   y : observed classes (0 or 1)
#   p : predicted scores of class 1
auc <- function(y, p) {
  r  <- rank(p)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Adversarial validation: can each level be recognised from the covariates?
# For each level, a random forest is trained to tell it from the rest (at most
# `n_max` observations per class), with a cross-validation by daily blocks.
# Only the rows where all the covariates are finite are kept.
#   dat       : data frame with the partition, `block` and the covariates
#   part      : partition, "year" or "zone"
#   vars      : names of the covariates
#   n_max     : maximum number of observations drawn per class
#   k         : number of folds
#   num_trees : number of trees of each forest
# Returns a list of two tibbles: auc (level, n, out-of-sample auc) and imp
# (level, variable, permutation importance averaged over the folds).
adversarial_validation <- function(dat, part, vars, n_max = 20000, k = 5,
                                   num_trees = 300) {
  d0 <- dat %>%
    select(all_of(c(part, "block", vars))) %>%
    filter(!is.na(.data[[part]]), if_all(all_of(vars), is.finite))
  part_levels <- levels(droplevels(d0[[part]]))

  res <- lapply(part_levels, function(lev) {
    # Class 1 = the level, class 0 = the rest; balanced subsampling
    d <- d0 %>%
      mutate(y = factor(as.integer(.data[[part]] == lev))) %>%
      group_by(y) %>%
      slice_sample(n = n_max) %>%
      ungroup()
    if (n_distinct(d$y) < 2) return(NULL)

    # Cross-validation by daily blocks
    blocks <- unique(d$block)
    d$fold <- sample(rep_len(1:k, length(blocks)))[match(d$block, blocks)]

    p   <- numeric(nrow(d))
    imp <- 0
    for (f in 1:k) {
      train <- d$fold != f
      rf <- ranger(y ~ ., data = d[train, c("y", vars)], probability = TRUE,
                   num.trees = num_trees, importance = "permutation")
      p[!train] <- predict(rf, d[!train, vars])$predictions[, "1"]
      imp <- imp + rf$variable.importance / k
    }

    list(
      auc = tibble(level = lev, n = nrow(d),
                   auc = auc(as.integer(as.character(d$y)), p)),
      imp = tibble(level = lev, variable = names(imp), importance = imp)
    )
  })

  res <- Filter(Negate(is.null), res)
  list(auc = bind_rows(lapply(res, `[[`, "auc")),
       imp = bind_rows(lapply(res, `[[`, "imp")))
}


# ---- Read the dataset --------------------------------------------------------

ds <- readRDS(in_file)
str(ds)
ds <- ds[ds$day == day_code, ]


# ---- 1) Preparation ----------------------------------------------------------

# In the dataset, a missing FOD class is the string "NA"
dat <- ds %>%
  filter(!is.na(nasc), !is.na(lat_nasc), !is.na(lon_nasc),
         !is.na(time_nasc)) %>%
  mutate(zone     = factor(na_if(fod, "NA")),
         year     = factor(format(time_nasc, "%Y")),
         block    = as.character(as.Date(time_nasc)),   # daily blocks (CV)
         log_nasc = log10(nasc + 1))

# Checks: number of ESU per zone, and per year and zone
print(table(dat$zone, useNA = "ifany"))
print(table(dat$year, dat$zone, useNA = "ifany"))


# ---- 2) Univariate shift (covariates and NASC) -------------------------------

long <- dat %>%
  select(year, zone, all_of(covars), log_nasc) %>%
  pivot_longer(c(all_of(covars), log_nasc),
               names_to = "variable", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(variable = factor(variable, levels = c("log_nasc", covars)))

save_plot(density_plot(long, "year"), "densities_year", width = 30, height = 22)
save_plot(density_plot(long, "zone"), "densities_zone", width = 30, height = 22)

shift_year <- univariate_shift(long, "year")
shift_zone <- univariate_shift(long, "zone")
save_table(shift_year, "univariate_year")
save_table(shift_zone, "univariate_zone")

save_plot(
  heatmap_plot(shift_year, "D_ks",
               "Kolmogorov-Smirnov D (year vs other years)"),
  "ks_year"
)
save_plot(
  heatmap_plot(shift_year, "overlap",
               "Overlap of the densities (year vs other years)"),
  "overlap_year"
)
save_plot(
  heatmap_plot(shift_year, "pct_out_of_range",
               "% of observations outside the range of the other years"),
  "out_of_range_year"
)
save_plot(
  heatmap_plot(shift_zone, "D_ks",
               "Kolmogorov-Smirnov D (zone vs other zones)"),
  "ks_zone"
)
save_plot(
  heatmap_plot(shift_zone, "overlap",
               "Overlap of the densities (zone vs other zones)"),
  "overlap_zone"
)
save_plot(
  heatmap_plot(shift_zone, "pct_out_of_range",
               "% of observations outside the range of the other zones"),
  "out_of_range_zone"
)


# ---- 3) Adversarial validation (multivariate shift) --------------------------

run_adversarial <- function(part, vars) {
  adversarial_validation(dat, part, vars, n_max = n_max, k = n_folds,
                         num_trees = num_trees)
}

adv <- list(
  year_phys = run_adversarial("year", covars_phys),
  year_all  = run_adversarial("year", covars),
  zone_phys = run_adversarial("zone", covars_phys),
  zone_all  = run_adversarial("zone", covars)
)

# AUC per level, and variables carrying the shift
auc_tab <- bind_rows(lapply(adv, `[[`, "auc"), .id = "analysis")
imp_tab <- bind_rows(lapply(adv, `[[`, "imp"), .id = "analysis")
print(auc_tab)
save_table(auc_tab, "adversarial_auc")
save_table(imp_tab, "adversarial_importance")

p_auc <- ggplot(auc_tab, aes(x = level, y = auc)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.5, linetype = 2) +
  facet_wrap(~ analysis, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Adversarial AUC (level vs rest)") +
  theme_bw()
save_plot(p_auc, "adversarial_auc")

save_plot(
  heatmap_plot(adv$year_all$imp, "importance",
               "Permutation importance - years"),
  "adversarial_importance_year"
)
save_plot(
  heatmap_plot(adv$zone_all$imp, "importance",
               "Permutation importance - zones"),
  "adversarial_importance_zone"
)


# ---- Save --------------------------------------------------------------------

saveRDS(
  list(shift_year = shift_year, shift_zone = shift_zone, adv = adv),
  file.path(fig_dir, "distribution_shift_results.rds")
)
cat("Figures and tables saved in:", fig_dir, "\n")
