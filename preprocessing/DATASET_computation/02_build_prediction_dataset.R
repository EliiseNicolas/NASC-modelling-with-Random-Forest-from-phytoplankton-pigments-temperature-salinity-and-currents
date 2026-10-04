# ==============================================================================
# 02_build_prediction_dataset.R
#
# Prediction dataset: FTLE, pigments and FOD on one grid, for all the dates
#
# Input  : FTLE at the FOD dates (01_create_ds_ftle_filtered_all_years.R)
#          <ftle_file>, a list with date, lon, lat, ftle (date x lon x lat)
#          pigments at the FOD dates
#          (01_create_ds_pigments_filtered_all_years.R)
#          <pigments_file>, a list with date, lon, lat and one array
#          (date x lon x lat) per variable
#          FOD classes on the grid (04_fod_clustering.R) and FOD grid
#          (02_fod_bspline.R), in <fod_dir>
#            cluster_transition_map_renamed.rds   (lon x lat x time)
#            lon.rds, lat.rds, time.rds
#
# Steps  : 1) reference grid = pigment grid; nearest FTLE and FOD grid point
#             of each of its longitudes and latitudes
#          2) dates common to the three datasets
#          3) extraction on the reference grid and the common dates, in one
#             indexing per variable (no loop over the dates)
#
# Output : <prediction_dataset_file>, a list with
#            date : dates common to the three datasets (Date)
#            lon  : longitudes of the pigment grid
#            lat  : latitudes of the pigment grid
#            ftle : array (date x lon x lat)
#            pig  : list of arrays (date x lon x lat), named like the pigment
#                   columns of the learning dataset: the concentrations
#                   (Chla, Per, ...), Chla_total and the ratios
#                   <pigment>_totpig
#            fod  : array (date x lon x lat), FOD class
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, ftle_file, pigments_file,
                     # prediction_dataset_file

out_file <- prediction_dataset_file


# ---- Functions ---------------------------------------------------------------

# Index of the nearest `source` coordinate of each `target` coordinate (both
# along the same axis, longitude or latitude).
nearest_index <- function(target, source) {
  vapply(target, function(x) which.min(abs(source - x)), integer(1))
}

# Dates as "YYYY-MM-DD" strings, to compare dates of different classes.
date_key <- function(x) format(as.Date(x), "%Y-%m-%d")


# ---- Read the three datasets -------------------------------------------------

ftle     <- readRDS(ftle_file)
pigments <- readRDS(pigments_file)

fod_map  <- readRDS(file.path(fod_dir, "cluster_transition_map_renamed.rds"))
lon_fod  <- as.vector(readRDS(file.path(fod_dir, "lon.rds")))
lat_fod  <- as.vector(readRDS(file.path(fod_dir, "lat.rds")))
time_fod <- readRDS(file.path(fod_dir, "time.rds"))

# Variables of the pigment dataset: concentrations, Chla_total and ratios
pig_vars <- setdiff(names(pigments), c("lon", "lat", "date"))


# ---- 1) Reference grid and nearest grid points -------------------------------

lon <- pigments$lon
lat <- pigments$lat

idx_lon_ftle <- nearest_index(lon, ftle$lon)
idx_lat_ftle <- nearest_index(lat, ftle$lat)
idx_lon_fod  <- nearest_index(lon, lon_fod)
idx_lat_fod  <- nearest_index(lat, lat_fod)


# ---- 2) Common dates ---------------------------------------------------------

dates_ftle <- date_key(ftle$date)
dates_pig  <- date_key(pigments$date)
dates_fod  <- date_key(time_fod)

# ISO strings sort in chronological order
common_dates <- sort(Reduce(intersect, list(dates_ftle, dates_pig, dates_fod)))
cat(length(common_dates), "common dates\n")

idx_date_ftle <- match(common_dates, dates_ftle)
idx_date_pig  <- match(common_dates, dates_pig)
idx_date_fod  <- match(common_dates, dates_fod)


# ---- 3) Extraction on the reference grid -------------------------------------

# FTLE: (date x lon x lat) on its own grid -> reference grid
ftle_all <- ftle$ftle[idx_date_ftle, idx_lon_ftle, idx_lat_ftle]

# FOD: (lon x lat x time) on its own grid -> reference grid, then
# (date x lon x lat) like the other arrays
fod_all <- aperm(fod_map[idx_lon_fod, idx_lat_fod, idx_date_fod], c(3, 1, 2))

# Pigments: already on the reference grid, only the dates are selected.
# "c_cond_Chla" -> "Chla"; Chla_total and the ratios keep their name.
pig_all <- lapply(pig_vars, function(v) pigments[[v]][idx_date_pig, , ])
names(pig_all) <- sub("^c_cond_", "", pig_vars)


# ---- Save --------------------------------------------------------------------

prediction_ds <- list(
  date = as.Date(common_dates),
  lon  = lon,
  lat  = lat,
  ftle = ftle_all,
  pig  = pig_all,
  fod  = fod_all
)

str(prediction_ds, max.level = 2)
cat("Size in memory:", format(object.size(prediction_ds), units = "GB"), "\n")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(prediction_ds, out_file)
cat("File saved:", out_file, "\n")
