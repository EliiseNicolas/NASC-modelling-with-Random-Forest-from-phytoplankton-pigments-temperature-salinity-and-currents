# ==============================================================================
# missing_values.R
#
# Data validation: missing values
#
# Input  : learning dataset (NASC, FOD, pigments and FTLE per ESU) at `freq`
#          kHz (01_build_learning_dataset.R)
#          <learning_dataset_dir>/learning_dataset_<years_tag>_<freq>kHz.rds
#          restricted to the ESU with day == `day_code` (config.R)
#
# Steps  : for the variables `vars` (nasc, fod, Chla, ftle)
#          1) map of the number of missing values per cell of a regular grid
#             (km, Lambert azimuthal equal-area projection centred on the
#             study area), one panel per variable
#          2) number of missing values per year
#          3) number of missing values per month and year
#          then, for all the variables of the dataset
#          4) combinations of variables missing together (UpSet plot)
#
# Output : in <fig_dir>
#            na_map_<res_km>km.png
#            na_count_per_year.png
#            na_count_per_month_and_year.png
#            na_upset.png
#
# Requires the packages `naniar` and `UpSetR`, in addition to those loaded
# below.
# ==============================================================================

library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)   # coastline
library(naniar)          # UpSet plot of the missing values

options(scipen = 999)
sf_use_s2(FALSE)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), day_code

freq     <- 38   # kHz
res_km   <- 20   # cell size of the map (km)

in_file <- learning_dataset_file(freq)
fig_dir <- file.path(fig_root_validation, "missing_values")

# Variables shown in the map and in the bar charts
vars <- c("nasc", "fod", "Chla", "ftle")

# One colour per variable (default ggplot colours)
hue_cols <- scales::hue_pal()(4)
var_cols <- c(nasc = hue_cols[1], fod = hue_cols[3],
              Chla = hue_cols[2], ftle = hue_cols[4])

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Save a figure in <fig_dir>, with the same size for all the figures.
save_fig <- function(plot, filename) {
  ggsave(file.path(fig_dir, filename), plot,
         width = 25, height = 22, units = "cm", dpi = 300)
}

# Lambert azimuthal equal-area projection centred on the points `pts` (sf).
laea_crs <- function(pts) {
  lonlat <- st_coordinates(pts)
  sprintf("+proj=laea +lat_0=%.3f +lon_0=%.3f +datum=WGS84 +units=m",
          mean(lonlat[, 2]), mean(lonlat[, 1]))
}

# Land around the points `pts`, projected in `crs`. The land is cropped to the
# area (+ `margin_deg` degrees) before the projection: projecting the whole
# world gives invalid geometries.
projected_land <- function(pts, crs, margin_deg = 20) {
  b    <- st_bbox(pts)
  area <- c(xmin = b[["xmin"]] - margin_deg, ymin = b[["ymin"]] - margin_deg,
            xmax = b[["xmax"]] + margin_deg, ymax = b[["ymax"]] + margin_deg)
  ne_countries(scale = "medium", returnclass = "sf") %>%
    st_crop(area) %>%
    st_union() %>%
    st_transform(crs) %>%
    st_make_valid()
}

# Number of missing values per grid cell and per variable.
#   d      : data frame with the variables `vars`
#   xy     : projected coordinates of the rows of `d` (matrix x, y, in m)
#   res_km : cell size (km)
# Returns a data frame with one row per variable and occupied cell: indices
# ix, iy, number of points n_tot, number of missing values n_na, percentage
# pct_na and cell centre x, y (m).
na_grid <- function(d, xy, vars, res_km) {
  res <- res_km * 1000
  d %>%
    transmute(ix = floor(xy[, 1] / res),
              iy = floor(xy[, 2] / res),
              across(all_of(vars), is.na)) %>%   # is.na() also detects NaN
    pivot_longer(all_of(vars), names_to = "variable", values_to = "missing") %>%
    group_by(variable, ix, iy) %>%
    summarise(n_tot = n(), n_na = sum(missing), .groups = "drop") %>%
    mutate(pct_na = 100 * n_na / n_tot,
           x = (ix + 0.5) * res,
           y = (iy + 0.5) * res)
}

# Map of the number of missing values per cell (log scale), one panel per
# variable. Cells without missing value are light grey. The title of each
# panel gives the overall percentage of missing values of the variable.
#   land, crs : projected land and its projection
#   margin_km : margin added around the grid (km)
na_map <- function(d, xy, vars, res_km, land, crs, margin_km = 0) {
  g <- na_grid(d, xy, vars, res_km)
  res <- res_km * 1000
  margin <- margin_km * 1000

  # Overall % of missing values in the title of each panel
  panel_labels <- g %>%
    group_by(variable) %>%
    summarise(p = 100 * sum(n_na) / sum(n_tot)) %>%
    mutate(lab = sprintf("%s (%.1f %% missing)", variable, p))
  g$variable <- factor(
    g$variable, levels = vars,
    labels = panel_labels$lab[match(vars, panel_labels$variable)]
  )

  ggplot() +
    geom_tile(data = g, aes(x = x, y = y, fill = if_else(n_na == 0, NA, n_na)),
              width = res, height = res) +
    geom_sf(data = land, fill = "grey80", colour = "grey50", linewidth = 0.2) +
    scale_fill_viridis_c(option = "magma", trans = "log10", na.value = "grey90",
                         name = "Number of\nmissing values",
                         labels = scales::label_number(big.mark = " ")) +
    coord_sf(crs = crs,
             xlim = range(g$x) + c(-1, 1) * (res / 2 + margin),
             ylim = range(g$y) + c(-1, 1) * (res / 2 + margin),
             expand = FALSE) +
    facet_wrap(~ variable) +
    labs(title = sprintf("Missing values per cell - %g x %g km grid",
                         res_km, res_km),
         x = NULL, y = NULL) +
    theme_bw()
}


# ---- Read the dataset --------------------------------------------------------

ds <- readRDS(in_file)
str(ds)
ds <- ds[ds$day == day_code, ]

# In the dataset, a missing FOD class is the string "NA"
ds <- ds %>% mutate(fod = na_if(fod, "NA"))


# ---- 1) Map of the missing values --------------------------------------------

# ESU with valid coordinates
d <- ds %>% filter(!is.na(lat_nasc), !is.na(lon_nasc))

pts      <- st_as_sf(d, coords = c("lon_nasc", "lat_nasc"), crs = 4326)
crs_laea <- laea_crs(pts)
xy       <- st_coordinates(st_transform(pts, crs_laea))
land     <- projected_land(pts, crs_laea)

p_map <- na_map(d, xy, vars, res_km, land, crs_laea)
save_fig(p_map, paste0("na_map_", res_km, "km.png"))


# ---- 2) Missing values per year ----------------------------------------------

na_long <- d %>%
  filter(!is.na(time_nasc)) %>%
  transmute(year      = format(time_nasc, "%Y"),
            month_num = as.integer(format(time_nasc, "%m")),
            across(all_of(vars), is.na)) %>%
  pivot_longer(all_of(vars), names_to = "variable", values_to = "missing") %>%
  mutate(variable = factor(variable, levels = vars))

# Months present in the whole dataset (the sampling months)
months_present <- sort(unique(na_long$month_num))
na_long <- na_long %>%
  mutate(month = factor(month.abb[month_num],
                        levels = month.abb[months_present]))

na_year <- na_long %>%
  group_by(year, variable) %>%
  summarise(n_na = sum(missing), .groups = "drop")

p_year <- ggplot(na_year, aes(x = year, y = n_na, fill = variable)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  scale_fill_manual(values = var_cols) +
  labs(x = "Year", y = "Number of missing values", fill = "Variable") +
  theme_bw()
print(p_year)
save_fig(p_year, "na_count_per_year.png")


# ---- 3) Missing values per month and year ------------------------------------

na_year_month <- na_long %>%
  group_by(year, month, variable) %>%
  summarise(n_na = sum(missing), .groups = "drop")

# drop = FALSE: all the sampling months stay on the x axis of every panel
p_year_month <- ggplot(na_year_month,
                       aes(x = month, y = n_na, fill = variable)) +
  geom_col(position = position_dodge(preserve = "single")) +
  facet_wrap(~ year) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  scale_fill_manual(values = var_cols) +
  labs(x = "Month", y = "Number of missing values", fill = "Variable") +
  theme_bw()
print(p_year_month)
save_fig(p_year_month, "na_count_per_month_and_year.png")


# ---- 4) Variables missing together -------------------------------------------

# All the variables with at least one missing value, at most 40 combinations
# On screen
gg_miss_upset(ds, nsets = n_var_miss(ds), nintersects = 40, order.by = "freq")

# In a file
png(file.path(fig_dir, "na_upset.png"),
    width = 30, height = 18, units = "cm", res = 300)
gg_miss_upset(ds, nsets = n_var_miss(ds), nintersects = 40)
dev.off()

cat("Figures saved in:", fig_dir, "\n")
