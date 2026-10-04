# ==============================================================================
# sampling_bias.R
#
# Data validation: sampling bias of the NASC
#
# Input  : learning dataset (NASC, FOD, pigments and FTLE per ESU) at `freq`
#          kHz (01_build_learning_dataset.R)
#          <learning_dataset_dir>/learning_dataset_<years_tag>_<freq>kHz.rds
#          restricted to the ESU with day == `day_code` (config.R)
#
# Steps  : 1) project the NASC points (Lambert azimuthal equal-area projection
#             centred on the study area)
#          2) density map: number of points per cell of a regular grid (km),
#             with the percentage of sea cells without any point
#          3) number of points per year
#          4) number of points per month, one panel per year
#
# Output : in <fig_dir>
#            nasc_density_<res_km>km.png
#            nasc_count_per_year.png
#            nasc_count_per_month_and_year.png
# ==============================================================================

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)   # coastline

options(scipen = 999)
sf_use_s2(FALSE)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, learning_dataset_file(), day_code

freq     <- 38   # kHz
res_km   <- 20   # cell size of the density map (km)

in_file <- learning_dataset_file(freq)
fig_dir <- file.path(fig_root_validation, "sampling_bias")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Save a figure in <fig_dir>, with the same size for all the figures.
save_fig <- function(plot, filename) {
  ggsave(file.path(fig_dir, filename), plot,
         width = 20, height = 18, units = "cm", dpi = 300)
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

# Number of points in each occupied cell of a regular grid.
#   xy     : projected coordinates of the points (matrix x, y, in m)
#   res_km : cell size (km)
# Returns a data frame with one row per occupied cell: cell indices ix, iy and
# number of points n.
count_grid <- function(xy, res_km) {
  res <- res_km * 1000
  data.frame(ix = floor(xy[, 1] / res),
             iy = floor(xy[, 2] / res)) %>%
    count(ix, iy, name = "n")
}

# Full grid covering the points (empty cells included), without the land cells.
#   land, crs : projected land and its projection
# Returns a data frame with one row per sea cell: indices ix, iy, number of
# points n (0 if empty), cell centre x, y (m) and res_km.
full_grid <- function(xy, res_km, land, crs) {
  res <- res_km * 1000
  occupied <- count_grid(xy, res_km)

  grid <- expand.grid(ix = seq(min(occupied$ix), max(occupied$ix)),
                      iy = seq(min(occupied$iy), max(occupied$iy))) %>%
    left_join(occupied, by = c("ix", "iy")) %>%
    mutate(n = coalesce(n, 0L),
           x = (ix + 0.5) * res,   # cell centre
           y = (iy + 0.5) * res,
           res_km = res_km)

  # A cell is on land if its centre is on a continent or an island
  centres <- st_as_sf(grid, coords = c("x", "y"), crs = crs)
  grid[lengths(st_intersects(centres, land)) == 0, ]
}

# Map of the number of NASC points per cell (log scale). Empty cells are white.
#   margin_km : margin added around the grid (km)
density_map <- function(xy, res_km, land, crs, margin_km = 0) {
  g <- full_grid(xy, res_km, land, crs)
  pct_empty <- 100 * mean(g$n == 0)
  res <- res_km * 1000
  margin <- margin_km * 1000

  # Limits = outer edge of the grid cells (+ margin)
  xlim <- range(g$x) + c(-1, 1) * (res / 2 + margin)
  ylim <- range(g$y) + c(-1, 1) * (res / 2 + margin)

  ggplot() +
    geom_tile(data = g, aes(x = x, y = y, fill = if_else(n == 0, NA, n)),
              width = res, height = res,
              colour = if (res_km >= 100) "grey70" else NA, linewidth = 0.1) +
    geom_sf(data = land, fill = "grey80", colour = "grey50", linewidth = 0.2) +
    scale_fill_viridis_c(trans = "log10", na.value = "white",
                         name = "Number of\nNASC points") +
    coord_sf(crs = crs, xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = sprintf("NASC sampling density - %g x %g km grid",
                         res_km, res_km),
         subtitle = sprintf("%.1f %% of empty cells", pct_empty),
         x = NULL, y = NULL) +
    theme_bw()
}


# ---- Read the dataset --------------------------------------------------------

ds <- readRDS(in_file)
str(ds)
ds <- ds[ds$day == day_code, ]


# ---- 1) Projection -----------------------------------------------------------

# Valid NASC points
pts <- ds %>%
  filter(!is.na(nasc), !is.na(lat_nasc), !is.na(lon_nasc)) %>%
  st_as_sf(coords = c("lon_nasc", "lat_nasc"), crs = 4326)

crs_laea <- laea_crs(pts)
xy       <- st_coordinates(st_transform(pts, crs_laea))
land     <- projected_land(pts, crs_laea)


# ---- 2) Density map ----------------------------------------------------------

p <- density_map(xy, res_km, land, crs_laea)
save_fig(p, paste0("nasc_density_", res_km, "km.png"))


# ---- 3) Number of points per year --------------------------------------------

d <- ds %>%
  filter(!is.na(nasc), !is.na(time_nasc)) %>%
  mutate(year  = factor(format(time_nasc, "%Y")),
         month = factor(format(time_nasc, "%m"), levels = sprintf("%02d", 1:12),
                        labels = month.abb))

p <- ggplot(d, aes(x = year)) +
  geom_bar(fill = "steelblue") +
  scale_y_continuous(labels = scales::label_number(big.mark = " ")) +
  labs(x = "Year", y = "Number of NASC points",
       title = "Number of NASC points per year") +
  theme_bw()
save_fig(p, "nasc_count_per_year.png")


# ---- 4) Number of points per month, one panel per year -----------------------

p <- ggplot(d, aes(x = month)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~ year) +
  scale_x_discrete(drop = TRUE) +
  labs(x = "Month", y = "Number of NASC points",
       title = "Number of NASC points per month and year") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p, "nasc_count_per_month_and_year.png")

cat("Figures saved in:", fig_dir, "\n")
