# ==============================================================================
# 01_create_ds_ftle_filtered_all_years.R
#
# FTLE dataset at the FOD dates, cropped to the FOD area
#
# Input  : daily FTLE NetCDF files in <raw_ftle_dir> (one file per day, date
#          YYYY-MM-DD in the file name, variable FTLE (lon, lat, time) with a
#          single time step)
#          FOD grid: lon.rds, lat.rds, time.rds in <fod_dir> (02_fod_bspline.R)
#
# Steps  : 1) keep the FTLE files whose date is one of the FOD dates
#          2) crop each file to the longitude / latitude range of the FOD grid
#          3) gather all the dates in a single object
#
# Output : <ftle_file>, a list with
#            lon, lat : cropped FTLE grid
#            date     : date of each file (Date)
#            ftle     : array (date x lon x lat); all NA for a file that was
#                       skipped because its grid differs from the first file
# ==============================================================================

library(ncdf4)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, ftle_file

in_dir   <- raw_ftle_dir
out_file <- ftle_file

# Date in the FTLE file names
date_pattern <- "\\d{4}-\\d{2}-\\d{2}"


# ---- Functions ---------------------------------------------------------------

# Longitude / latitude of an open FTLE file, cropped to the FOD area.
# Returns the cropped coordinates and their position in the file.
crop_grid <- function(ds, lon_range, lat_range) {
  lon <- ncvar_get(ds, "lon")
  lat <- ncvar_get(ds, "lat")

  i_lon <- which(lon >= lon_range[1] & lon <= lon_range[2])
  i_lat <- which(lat >= lat_range[1] & lat <= lat_range[2])

  list(
    lon   = lon[i_lon],
    lat   = lat[i_lat],
    start = c(min(i_lon), min(i_lat), 1),
    count = c(length(i_lon), length(i_lat), 1)
  )
}


# ---- FOD area and dates ------------------------------------------------------

lon_range <- range(readRDS(file.path(fod_dir, "lon.rds")))
lat_range <- range(readRDS(file.path(fod_dir, "lat.rds")))
fod_dates <- unique(format(readRDS(file.path(fod_dir, "time.rds")), "%Y-%m-%d"))


# ---- 1) FTLE files at the FOD dates ------------------------------------------

files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)

# Keep the files whose name contains one of the FOD dates
files <- files[grepl(paste(fod_dates, collapse = "|"), basename(files))]

# Date of each file: the first date found in its name
file_dates <- regmatches(basename(files), regexpr(date_pattern, basename(files)))

if (length(files) == 0) stop("No FTLE file at the FOD dates in ", in_dir)
if (anyDuplicated(file_dates)) {
  warning("Several FTLE files for the same date: ",
          paste(unique(file_dates[duplicated(file_dates)]), collapse = ", "))
}

# FOD dates without FTLE file. Not necessarily a problem: there may be no NASC
# at these dates either.
missing_dates <- setdiff(fod_dates, file_dates)
cat("FOD dates without FTLE file:", length(missing_dates), "/",
    length(fod_dates), "\n")
print(missing_dates)


# ---- 2) and 3) Crop and gather -----------------------------------------------

# Reference grid: the first file
ds  <- nc_open(files[1])
ref <- crop_grid(ds, lon_range, lat_range)
nc_close(ds)

cat("Cropped FTLE grid:", length(ref$lon), "lon x", length(ref$lat), "lat\n")

n_files <- length(files)
ftle <- list(
  lon  = ref$lon,
  lat  = ref$lat,
  date = as.Date(file_dates),
  ftle = array(NA_real_, dim = c(n_files, length(ref$lon), length(ref$lat)))
)
loaded <- logical(n_files)

for (i in seq_len(n_files)) {
  ds   <- nc_open(files[i])
  grid <- crop_grid(ds, lon_range, lat_range)

  same_grid <- isTRUE(all.equal(grid$lon, ref$lon)) &&
    isTRUE(all.equal(grid$lat, ref$lat))

  if (same_grid) {
    # FTLE is (lon, lat, time) with one time step: the result is (lon x lat)
    ftle$ftle[i, , ] <- ncvar_get(ds, "FTLE",
                                  start = grid$start, count = grid$count)
    loaded[i] <- TRUE
    cat(i, "/", n_files, "-", file_dates[i], "\n")
  } else {
    cat(i, "/", n_files, "-", file_dates[i], ": grid differs from",
        file_dates[1], "(", length(grid$lon), "x", length(grid$lat),
        "instead of", length(ref$lon), "x", length(ref$lat),
        ") -> file skipped\n")
  }

  nc_close(ds)
}

str(ftle)
cat(sprintf("FTLE dates loaded: %d / %d\n", sum(loaded), n_files))


# ---- Save --------------------------------------------------------------------

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(ftle, out_file)
cat("File saved:", out_file, "\n")
