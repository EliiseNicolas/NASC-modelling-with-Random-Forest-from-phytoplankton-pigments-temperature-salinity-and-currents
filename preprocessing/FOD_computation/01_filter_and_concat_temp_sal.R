# ==============================================================================
# 01_filter_and_concat_temp_sal.R
#
# Crop and concatenation of the yearly temperature / salinity files
#
# Input  : Copernicus Marine (GLORYS12) daily temperature and salinity
#          (00_download_temp_salinity.R), one NetCDF file per year in
#          <raw_temp_sal_dir>, thetao / so (longitude, latitude, depth, time)
#
# Steps  : 1) crop each year to the study area and depth range
#          2) check that all years share the same longitude / latitude / depth grid
#          3) concatenate the years along time
#          4) save the result in a single NetCDF file
#
# Output : <temp_sal_dir>/thetao_so_crop_<years_tag>.nc
# ==============================================================================

library(ncdf4)
# also requires the `abind` package


# ---- Configuration -----------------------------------------------------------

source("config.R")   # years, study area, directories, temp_sal_file

in_dir <- raw_temp_sal_dir

# One file per year, same area and same period (9 January - 3 March)
in_files <- paste0(
  "cmems_mod_glo_phy_my_0.083deg_P1D-m_thetao-so_40.00E-95.00E_60.00S-20.00S_",
  "0.49-5727.92m_", years, "-01-09-", years, "-03-03.nc"
)
names(in_files) <- years

out_file <- temp_sal_file

# Depth range (bounds included, like the study area)
depth_min <- 20; depth_max <- 500   # m


# ---- Functions ---------------------------------------------------------------

# Read one yearly file, cropped to the study area and depth range.
# Returns a list: thetao and so (lon x lat x depth x time), lon, lat, depth,
# time and time_units.
read_crop_nc <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))

  lon   <- ncvar_get(nc, "longitude")
  lat   <- ncvar_get(nc, "latitude")
  depth <- ncvar_get(nc, "depth")

  i_lon   <- which(lon >= lon_min & lon <= lon_max)
  i_lat   <- which(lat >= lat_min & lat <= lat_max)
  i_depth <- which(depth >= depth_min & depth <= depth_max)

  # Read only the cropped block: variables are (lon, lat, depth, time)
  start <- c(min(i_lon), min(i_lat), min(i_depth), 1)
  count <- c(length(i_lon), length(i_lat), length(i_depth), -1)

  list(
    thetao     = ncvar_get(nc, "thetao", start = start, count = count,
                           collapse_degen = FALSE),
    so         = ncvar_get(nc, "so", start = start, count = count,
                           collapse_degen = FALSE),
    lon        = lon[i_lon],
    lat        = lat[i_lat],
    depth      = depth[i_depth],
    time       = ncvar_get(nc, "time"),
    time_units = nc$dim$time$units
  )
}


# ---- 1) Read and crop --------------------------------------------------------

ds <- lapply(years, function(year) {
  d <- read_crop_nc(file.path(in_dir, in_files[[year]]))

  cat(sprintf(
    "%s: lon %g to %g, lat %g to %g, depth %g to %g m, %d time steps\n",
    year, min(d$lon), max(d$lon), min(d$lat), max(d$lat),
    min(d$depth), max(d$depth), length(d$time)
  ))
  d
})
names(ds) <- years


# ---- 2) Check that all years share the same grid -----------------------------

# Grid of the first year (without the data arrays, to free memory later)
ref <- ds[[1]][c("lon", "lat", "depth", "time_units")]

for (year in years) {
  for (v in c("lon", "lat", "depth")) {
    if (!isTRUE(all.equal(ds[[year]][[v]], ref[[v]]))) {
      stop("`", v, "` of ", year, " differs from ", years[1])
    }
  }
  if (ds[[year]]$time_units != ref$time_units) {
    stop("Time units of ", year, " differ from ", years[1])
  }
}


# ---- 3) Concatenate along time -----------------------------------------------

thetao_all <- do.call(
  abind::abind, c(unname(lapply(ds, `[[`, "thetao")), list(along = 4))
)
so_all <- do.call(
  abind::abind, c(unname(lapply(ds, `[[`, "so")), list(along = 4))
)
time_all <- unlist(lapply(ds, `[[`, "time"), use.names = FALSE)

rm(ds)
gc()

cat("Concatenated dimensions (lon, lat, depth, time):", dim(thetao_all), "\n")
stopifnot(
  identical(dim(thetao_all), dim(so_all)),
  dim(thetao_all)[4] == length(time_all)
)


# ---- 4) Save as NetCDF -------------------------------------------------------

dim_lon   <- ncdim_def("longitude", units = "degrees_east", vals = ref$lon)
dim_lat   <- ncdim_def("latitude", units = "degrees_north", vals = ref$lat)
dim_depth <- ncdim_def("depth", units = "m", vals = ref$depth)
dim_time  <- ncdim_def("time", units = ref$time_units, vals = time_all,
                       unlim = TRUE)

dims <- list(dim_lon, dim_lat, dim_depth, dim_time)

var_thetao <- ncvar_def("thetao", units = "degC", dim = dims,
                        missval = -9999, prec = "float")
var_so     <- ncvar_def("so", units = "psu", dim = dims,
                        missval = -9999, prec = "float")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

nc_out <- nc_create(out_file, vars = list(var_thetao, var_so))
ncvar_put(nc_out, var_thetao, thetao_all)
ncvar_put(nc_out, var_so, so_all)
nc_close(nc_out)

cat("File saved:", out_file, "\n")


# ---- Check the written file --------------------------------------------------

nc <- nc_open(out_file)
dim_written <- nc$var$thetao$varsize
nc_close(nc)

cat("Dimensions in the file (lon, lat, depth, time):", dim_written, "\n")
stopifnot(all(dim_written == dim(thetao_all)))
