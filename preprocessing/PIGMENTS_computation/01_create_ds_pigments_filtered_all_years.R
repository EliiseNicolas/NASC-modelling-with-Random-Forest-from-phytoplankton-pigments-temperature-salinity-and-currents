# ==============================================================================
# 01_create_ds_pigments_filtered_all_years.R
#
# Pigment dataset at the FOD dates, cropped to the FOD area and filtered
#
# Input  : daily PIGMeANN NetCDF files in <raw_pigments_dir> (one file per day,
#          date YYYYMMDD in the file name, variables c_cond_<pigment>,
#          use_<pigment> and in_domain (lon, lat, time) with a single time step)
#          FOD grid: lon.rds, lat.rds, time.rds in <fod_dir> (02_fod_bspline.R)
#
# Steps  : 1) keep the pigment files whose date is one of the FOD dates
#          2) crop each file to the longitude / latitude range of the FOD grid
#          3) filter each pigment, at each date:
#               a. values outside the 1 % - 99 % quantiles of the area -> NA
#               b. use_<pigment> == 0 -> 0
#               c. in_domain == 0 -> NA
#          4) ratio of each pigment to the sum of the pigments other than Chla
#
# Output : <pigments_file>, a list with
#            lon, lat         : cropped pigment grid
#            date             : date of each file (Date)
#            c_cond_<pigment> : array (date x lon x lat), one per pigment
#            Chla_total       : Chla alone
#            <pigment>_totpig : ratio of the pigment to the sum of the pigments
#                               other than Chla (no ratio for Chla)
# ==============================================================================

library(ncdf4)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, pigments_file

in_dir   <- raw_pigments_dir
out_file <- pigments_file

pig_names <- c(
  "Chla", "Per", "But", "Fuco", "Hex",
  "Allo", "Zea", "Chlb", "DvChla"
)

# Values outside these quantiles of the cropped area are set to NA
q_range <- c(0.01, 0.99)

# Date in the pigment file names
date_pattern <- "\\d{8}"


# ---- FOD area and dates ------------------------------------------------------

lon_range <- range(readRDS(file.path(fod_dir, "lon.rds")))
lat_range <- range(readRDS(file.path(fod_dir, "lat.rds")))
fod_dates <- unique(format(readRDS(file.path(fod_dir, "time.rds")), "%Y%m%d"))


# ---- 1) Pigment files at the FOD dates ---------------------------------------

files <- list.files(in_dir, pattern = "\\.nc$", full.names = TRUE)

# Keep the files whose name contains one of the FOD dates
files <- files[grepl(paste(fod_dates, collapse = "|"), basename(files))]

# Date of each file: the first block of 8 digits found in its name
file_names <- basename(files)
file_dates <- regmatches(file_names, regexpr(date_pattern, file_names))

if (length(files) == 0) stop("No pigment file at the FOD dates in ", in_dir)

# FOD dates without pigment file
missing_dates <- setdiff(fod_dates, file_dates)
cat("FOD dates without pigment file:", length(missing_dates), "/",
    length(fod_dates), "\n")
print(missing_dates)


# ---- 2) Cropped grid ---------------------------------------------------------

# Reference grid: the first file. All files must share the same grid.
ds      <- nc_open(files[1])
lon_ref <- ncvar_get(ds, "lon")
lat_ref <- ncvar_get(ds, "lat")
nc_close(ds)

i_lon <- which(lon_ref >= lon_range[1] & lon_ref <= lon_range[2])
i_lat <- which(lat_ref >= lat_range[1] & lat_ref <= lat_range[2])

n_files <- length(files)
n_lon   <- length(i_lon)
n_lat   <- length(i_lat)

cat("Cropped pigment grid:", n_lon, "lon x", n_lat, "lat\n")

# Block read in each file: variables are (lon, lat, time) with one time step
start <- c(min(i_lon), min(i_lat), 1)
count <- c(n_lon, n_lat, 1)

pigments <- list(
  lon  = lon_ref[i_lon],
  lat  = lat_ref[i_lat],
  date = as.Date(file_dates, format = "%Y%m%d")
)
for (pig in pig_names) {
  pigments[[paste0("c_cond_", pig)]] <- array(
    NA_real_, dim = c(n_files, n_lon, n_lat)
  )
}


# ---- 3) Read, crop and filter each file --------------------------------------

for (i in seq_len(n_files)) {
  ds <- nc_open(files[i])

  same_grid <- isTRUE(all.equal(ncvar_get(ds, "lon"), lon_ref)) &&
    isTRUE(all.equal(ncvar_get(ds, "lat"), lat_ref))
  if (!same_grid) {
    nc_close(ds)
    stop("Grid of ", basename(files[i]), " differs from ", basename(files[1]))
  }

  in_domain <- ncvar_get(ds, "in_domain", start = start, count = count)

  for (pig in pig_names) {
    c_cond <- ncvar_get(ds, paste0("c_cond_", pig),
                        start = start, count = count)
    use    <- ncvar_get(ds, paste0("use_", pig),
                        start = start, count = count)

    # a. values outside the quantiles of the area
    q <- quantile(c_cond, q_range, na.rm = TRUE)
    c_cond[c_cond < q[1] | c_cond > q[2]] <- NA

    # b. pigment flagged as not usable
    c_cond[use == 0] <- 0

    # c. pixels outside the domain of validity
    c_cond[in_domain == 0] <- NA

    pigments[[paste0("c_cond_", pig)]][i, , ] <- c_cond
  }

  nc_close(ds)
  cat(i, "/", n_files, "-", format(pigments$date[i], "%Y-%m-%d"), "\n")
}


# ---- 4) Ratios ---------------------------------------------------------------

# Chla_total = Chla alone
pigments$Chla_total <- pigments$c_cond_Chla

# Sum of the pigments other than Chla
other_pigs <- setdiff(pig_names, "Chla")
sum_others <- Reduce(`+`, pigments[paste0("c_cond_", other_pigs)])

# Ratio of each pigment to this sum. No ratio for Chla: it is not part of
# the sum.
for (pig in other_pigs) {
  pigments[[paste0(pig, "_totpig")]] <-
    pigments[[paste0("c_cond_", pig)]] / sum_others
}

str(pigments)


# ---- Save --------------------------------------------------------------------

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(pigments, out_file)
cat("File saved:", out_file, "\n")
