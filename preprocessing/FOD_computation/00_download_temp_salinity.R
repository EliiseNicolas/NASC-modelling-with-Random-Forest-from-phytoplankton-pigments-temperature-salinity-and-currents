# ==============================================================================
# 00_download_temp_salinity.R
#
# Download of the daily temperature and salinity (Copernicus Marine)
#
# Input  : Copernicus Marine credentials, in a `.env` file at the project root
#          (see .env.example)
#            COPERNICUS_USERNAME=...
#            COPERNICUS_PASSWORD=...
#
# Steps  : call the Copernicus Marine Toolbox (`copernicusmarine subset`) to
#          download thetao and so from the GLORYS12V1 reanalysis
#          (cmems_mod_glo_phy_my_0.083deg_P1D-m: daily, 1/12 degree, 50 levels)
#          over the requested period and area
#
# Output : one NetCDF file in <raw_temp_sal_dir>, named by the Toolbox
#
# Requires the Copernicus Marine Toolbox (`copernicusmarine` executable; it is
# not an R package) and the R package `dotenv`.
#
# Usage  : from the project root, in R
#            source("FOD_computation/00_download_temp_salinity.R")
#            download_temp_salinity("2022-01-09", "2022-03-03",
#                                   -60, -20, 40, 95)
#          or on the command line
#            Rscript FOD_computation/00_download_temp_salinity.R \
#              --date-start 2022-01-09 --date-end 2022-03-03 \
#              --lat-min -60 --lat-max -20 --lon-min 40 --lon-max 95 \
#              [--output-dir DIR]
#          01_filter_and_concat_temp_sal.R expects one file per year, for the
#          area and the period (9 January - 3 March) of this example.
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories

out_dir    <- raw_temp_sal_dir
dataset_id <- "cmems_mod_glo_phy_my_0.083deg_P1D-m"

# Credentials: read the .env file of the project root, if any
if (file.exists(".env")) dotenv::load_dot_env(".env")


# ---- Functions ---------------------------------------------------------------

# Download the daily temperature and salinity from Copernicus Marine.
#   date_start, date_end : first and last day, "YYYY-MM-DD"
#   lat_min, lat_max     : latitude range (degrees North)
#   lon_min, lon_max     : longitude range (degrees East)
#   output_dir           : destination folder
#   toolbox              : path to the `copernicusmarine` executable
#                          (default: the one found in the PATH)
# Returns the path of the destination folder (invisibly).
download_temp_salinity <- function(date_start,
                                   date_end,
                                   lat_min,
                                   lat_max,
                                   lon_min,
                                   lon_max,
                                   output_dir = out_dir,
                                   toolbox = Sys.which("copernicusmarine")) {
  username <- Sys.getenv("COPERNICUS_USERNAME")
  password <- Sys.getenv("COPERNICUS_PASSWORD")

  if (username == "" || password == "") {
    stop(
      "Copernicus credentials not found.\n",
      "1. Copy .env.example to .env\n",
      "2. Fill in COPERNICUS_USERNAME and COPERNICUS_PASSWORD"
    )
  }

  if (toolbox == "" || !file.exists(toolbox)) {
    stop(
      "`copernicusmarine` executable not found.\n",
      "Install the Copernicus Marine Toolbox, or give its path with the ",
      "`toolbox` argument."
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_dir)

  cat("Download parameters:\n")
  cat("  Period    :", date_start, "to", date_end, "\n")
  cat("  Latitude  :", lat_min, "to", lat_max, "degrees N\n")
  cat("  Longitude :", lon_min, "to", lon_max, "degrees E\n")
  cat("  Folder    :", output_path, "\n\n")

  # The Toolbox reads the credentials from these environment variables, so
  # that they do not appear on the command line.
  Sys.setenv(
    COPERNICUSMARINE_SERVICE_USERNAME = username,
    COPERNICUSMARINE_SERVICE_PASSWORD = password
  )

  args <- c(
    "subset",
    "--dataset-id", dataset_id,
    "--variable", "thetao",
    "--variable", "so",
    "--start-datetime", paste0(date_start, "T00:00:00"),
    "--end-datetime", paste0(date_end, "T23:59:59"),
    "--minimum-latitude", lat_min,
    "--maximum-latitude", lat_max,
    "--minimum-longitude", lon_min,
    "--maximum-longitude", lon_max,
    "--output-directory", shQuote(output_path)
  )

  status <- system2(toolbox, args)
  if (status != 0) {
    stop("Download failed (exit status ", status, ").")
  }

  cat("\nDownload finished.\n")
  cat("Files saved in:", output_path, "\n")
  invisible(output_path)
}


# ---- Command line (Rscript) --------------------------------------------------
# This block is not run when the file is loaded with source().

if (sys.nframe() == 0 && !interactive()) {
  usage <- paste(
    "Download the daily temperature and salinity from Copernicus Marine.",
    "",
    "Arguments:",
    "  --date-start  First day (YYYY-MM-DD)",
    "  --date-end    Last day  (YYYY-MM-DD)",
    "  --lat-min     Minimum latitude  (degrees N)",
    "  --lat-max     Maximum latitude  (degrees N)",
    "  --lon-min     Minimum longitude (degrees E)",
    "  --lon-max     Maximum longitude (degrees E)",
    "  --output-dir  Destination folder",
    paste0("                (default: ", out_dir, ")"),
    sep = "\n"
  )

  # Arguments are read as `--name value` pairs
  cli <- commandArgs(trailingOnly = TRUE)
  if (any(cli %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(status = 0)
  }
  if (length(cli) %% 2 != 0) stop("Invalid arguments.\n\n", usage)

  opts <- as.list(cli[c(FALSE, TRUE)])
  names(opts) <- cli[c(TRUE, FALSE)]

  required <- c("--date-start", "--date-end",
                "--lat-min", "--lat-max", "--lon-min", "--lon-max")
  allowed  <- c(required, "--output-dir")

  if (!all(names(opts) %in% allowed)) {
    stop("Unknown argument: ",
         paste(setdiff(names(opts), allowed), collapse = ", "), "\n\n", usage)
  }
  if (!all(required %in% names(opts))) {
    stop("Missing argument: ",
         paste(setdiff(required, names(opts)), collapse = ", "), "\n\n", usage)
  }

  bounds <- sapply(opts[c("--lat-min", "--lat-max", "--lon-min", "--lon-max")],
                   as.numeric)
  if (anyNA(bounds)) stop("Latitudes and longitudes must be numbers.")

  if (is.null(opts[["--output-dir"]])) opts[["--output-dir"]] <- out_dir

  download_temp_salinity(
    date_start = opts[["--date-start"]],
    date_end   = opts[["--date-end"]],
    lat_min    = bounds[["--lat-min"]],
    lat_max    = bounds[["--lat-max"]],
    lon_min    = bounds[["--lon-min"]],
    lon_max    = bounds[["--lon-max"]],
    output_dir = opts[["--output-dir"]]
  )
}
