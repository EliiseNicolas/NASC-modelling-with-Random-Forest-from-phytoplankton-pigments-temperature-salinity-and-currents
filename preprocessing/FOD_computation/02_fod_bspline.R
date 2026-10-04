# ==============================================================================
# 02_fod_bspline.R
#
# Functional Oceanographic Domains (FOD), step 1/3: B-spline estimation
#
# Input  : cropped and concatenated temperature / salinity
#          (01_filter_and_concat_temp_sal.R), thetao and so
#          (longitude, latitude, depth, time)
#
# Steps  : 1) flatten the data to one row per profile (lon x lat x time)
#          2) keep the profiles with no missing value in temperature and salinity
#          3) project each profile on a penalised B-spline basis
#          4) diagnostic figures
#
# Output : in <fod_dir>
#            lon.rds, lat.rds, depth.rds, time.rds   grid (time in POSIXct UTC)
#            mask_common.rds, valid_idx.rds          profiles kept (flattened grid)
#            phi.rds                                 B-spline basis (n_depth x K)
#            coef_thetao.rds, coef_so.rds            coefficients (n_profiles x K)
#            temperature_profiles.rds,
#            salinity_profiles.rds                   smoothed profiles
#                                                    (n_profiles x n_depth)
#          in <fig_dir>
#            bspline_basis.png, bspline_fit_examples.png, bspline_rmse_by_depth.png
# ==============================================================================

library(ncdf4)
library(splines)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories, temp_sal_file

in_file <- temp_sal_file
fig_dir <- file.path(fig_root, "FOD_computation", "02_fod_bspline")

# B-spline basis
K             <- 25     # number of basis functions
spline_order  <- 4      # order of the basis functions (4 = cubic)
lambda_spline <- 0.25   # weight of the roughness penalty

# Diagnostic figures (computed on a random sample of profiles)
n_check    <- 100000    # profiles used for the RMSE by depth
n_examples <- 6         # profiles shown in the example figure
seed       <- 1

dir.create(fod_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Functions ---------------------------------------------------------------

# Read a variable (lon, lat, depth, time) and flatten it to a matrix with one
# row per profile (lon varies fastest, then lat, then time) and one column per
# depth.
read_profiles <- function(ds, varname) {
  x <- ncvar_get(ds, varname)
  matrix(aperm(x, c(1, 2, 4, 3)), ncol = dim(x)[3])
}

# Project the profiles `x` on the basis, save the coefficients and the smoothed
# profiles, and return the raw / smoothed values of the profiles `idx_check`.
fit_bspline <- function(x, phi, B, idx_check, coef_file, profiles_file) {
  coef <- x %*% t(B)            # (n_profiles x K)
  saveRDS(coef, coef_file)

  rec <- coef %*% t(phi)        # (n_profiles x n_depth)
  saveRDS(rec, profiles_file)

  list(
    raw = x[idx_check, , drop = FALSE],
    rec = rec[idx_check, , drop = FALSE]
  )
}

# Raw values (points) and B-spline fit (lines) of a few profiles
plot_fit_examples <- function(raw, rec, depth, xlab) {
  plot(
    NULL,
    xlim = range(raw, rec), ylim = rev(range(depth)),
    xlab = xlab, ylab = "Depth (m)"
  )
  for (i in seq_len(nrow(raw))) {
    points(raw[i, ], depth, col = i, pch = 16, cex = 0.7)
    lines(rec[i, ], depth, col = i, lwd = 2)
  }
  legend(
    "bottomright", legend = c("Data", "B-spline fit"),
    pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2), bty = "n", cex = 0.8
  )
}


# ---- 1) Read and flatten -----------------------------------------------------

ds <- nc_open(in_file)

lon        <- ncvar_get(ds, "longitude")
lat        <- ncvar_get(ds, "latitude")
depth      <- ncvar_get(ds, "depth")
time_raw   <- ncvar_get(ds, "time")
time_units <- ds$dim$time$units

thetao <- read_profiles(ds, "thetao")   # (n_grid x n_depth)
so     <- read_profiles(ds, "so")       # (n_grid x n_depth)

nc_close(ds)

if (!grepl("^hours since 1950-01-01", time_units)) {
  stop("Unexpected time units: '", time_units,
       "' (expected 'hours since 1950-01-01')")
}
time <- as.POSIXct(time_raw * 3600, origin = "1950-01-01", tz = "UTC")

stopifnot(
  nrow(thetao) == length(lon) * length(lat) * length(time),
  ncol(thetao) == length(depth),
  identical(dim(thetao), dim(so))
)
cat("Grid:", length(lon), "lon x", length(lat), "lat x", length(depth),
    "depths x", length(time), "time steps\n")


# ---- 2) Keep the profiles complete in temperature and salinity ---------------

# A row sum is NA as soon as the profile has one missing value
na_thetao <- is.na(rowSums(thetao))
na_so     <- is.na(rowSums(so))

mask_common <- !na_thetao & !na_so   # one value per (lon, lat, time)

cat(sprintf("Profiles with NA: temperature %.2f %%, salinity %.2f %%\n",
            100 * mean(na_thetao), 100 * mean(na_so)))
cat(sprintf("Profiles removed: %.2f %% (%d kept)\n",
            100 * mean(!mask_common), sum(mask_common)))

thetao <- thetao[mask_common, , drop = FALSE]
so     <- so[mask_common, , drop = FALSE]
rm(na_thetao, na_so)
gc()

saveRDS(lon, file.path(fod_dir, "lon.rds"))
saveRDS(lat, file.path(fod_dir, "lat.rds"))
saveRDS(depth, file.path(fod_dir, "depth.rds"))
saveRDS(time, file.path(fod_dir, "time.rds"))
saveRDS(mask_common, file.path(fod_dir, "mask_common.rds"))
saveRDS(which(mask_common), file.path(fod_dir, "valid_idx.rds"))


# ---- 3) Penalised B-spline estimation ----------------------------------------

# Basis evaluated at the data depths (n_depth x K)
phi <- as.matrix(bs(depth, df = K, degree = spline_order - 1, intercept = TRUE))

# Roughness penalty: second-order differences of the coefficients
D2 <- diff(diag(K), differences = 2)
R  <- crossprod(D2)

# Penalised least squares: coef = (phi'phi + lambda R)^-1 phi' x
A <- crossprod(phi) + lambda_spline * R
cat("Depth levels:", length(depth), "- basis functions:", K,
    "- rank of the penalised system:", qr(A)$rank, "\n")
B <- solve(A, t(phi))   # (K x n_depth)

saveRDS(phi, file.path(fod_dir, "phi.rds"))

# Profiles used for the diagnostic figures
set.seed(seed)
idx_check <- sample(nrow(thetao), min(n_check, nrow(thetao)))

check_temp <- fit_bspline(
  thetao, phi, B, idx_check,
  coef_file     = file.path(fod_dir, "coef_thetao.rds"),
  profiles_file = file.path(fod_dir, "temperature_profiles.rds")
)
rm(thetao)
gc()

check_sal <- fit_bspline(
  so, phi, B, idx_check,
  coef_file     = file.path(fod_dir, "coef_so.rds"),
  profiles_file = file.path(fod_dir, "salinity_profiles.rds")
)
rm(so)
gc()


# ---- 4) Figures --------------------------------------------------------------

# B-spline basis
png(file.path(fig_dir, "bspline_basis.png"),
    width = 2400, height = 1500, res = 300)
matplot(
  depth, unclass(phi), type = "l", lty = 1,
  xlab = "Depth (m)", ylab = "Basis function",
  main = paste0("B-spline basis (K = ", K, ", order ", spline_order, ")")
)
rug(depth)   # depth levels of the data
dev.off()

# Raw vs smoothed profiles, for a few profiles
i_ex <- seq_len(min(n_examples, length(idx_check)))

png(file.path(fig_dir, "bspline_fit_examples.png"),
    width = 2800, height = 2000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plot_fit_examples(check_temp$raw[i_ex, , drop = FALSE],
                  check_temp$rec[i_ex, , drop = FALSE],
                  depth, xlab = "Temperature (°C)")
plot_fit_examples(check_sal$raw[i_ex, , drop = FALSE],
                  check_sal$rec[i_ex, , drop = FALSE],
                  depth, xlab = "Salinity (PSU)")
mtext("B-spline fit of randomly chosen profiles", outer = TRUE, cex = 1.2)
dev.off()

# RMSE of the fit by depth
rmse_temp <- sqrt(colMeans((check_temp$raw - check_temp$rec)^2))
rmse_sal  <- sqrt(colMeans((check_sal$raw - check_sal$rec)^2))

png(file.path(fig_dir, "bspline_rmse_by_depth.png"),
    width = 2800, height = 2000, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
plot(rmse_temp, depth, type = "b", pch = 16, ylim = rev(range(depth)),
     xlab = "RMSE temperature (°C)", ylab = "Depth (m)")
plot(rmse_sal, depth, type = "b", pch = 16, ylim = rev(range(depth)),
     xlab = "RMSE salinity (PSU)", ylab = "Depth (m)")
mtext(paste0("B-spline fit error by depth (", length(idx_check), " profiles)"),
      outer = TRUE, cex = 1.2)
dev.off()

cat("Results saved in:", fod_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
