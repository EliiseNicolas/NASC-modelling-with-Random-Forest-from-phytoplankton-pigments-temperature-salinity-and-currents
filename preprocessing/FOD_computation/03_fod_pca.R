# ==============================================================================
# 03_fod_pca.R
#
# Functional Oceanographic Domains (FOD), step 2/3: PCA
#
# Input  : B-spline coefficients of temperature and salinity (02_fod_bspline.R)
#          coef_thetao.rds, coef_so.rds, phi.rds, depth.rds in <fod_dir>
#
# Steps  : 1) PCA of the joint coefficients [temperature, salinity]
#          2) table of explained variance
#          3) scores of the profiles on the first `nharm` components
#          4) figures
#
# Output : in <fod_dir>
#            pca_scores.rds          scores (n_profiles x nharm)
#            pca_results.rds         mean coefficients, eigenvectors, eigenvalues,
#                                    explained variance
#            pca_variance_table.csv  explained variance of every component
#          in <fig_dir>
#            pca_explained_variance.png, pca_eigenfunctions.png
# ==============================================================================


# ---- Configuration -----------------------------------------------------------

source("config.R")   # directories

fig_dir <- file.path(fig_root, "FOD_computation", "03_fod_pca")

nharm  <- 6    # number of principal components kept
n_show <- 15   # number of components shown in the explained-variance figure

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- 1) PCA of the B-spline coefficients -------------------------------------

coef_thetao <- readRDS(file.path(fod_dir, "coef_thetao.rds"))   # (n x K)
coef_so     <- readRDS(file.path(fod_dir, "coef_so.rds"))       # (n x K)
K <- ncol(coef_thetao)
stopifnot(identical(dim(coef_thetao), dim(coef_so)))

# Joint coefficients: temperature in columns 1..K, salinity in K+1..2K
coef_biv <- cbind(coef_thetao, coef_so)
rm(coef_thetao, coef_so)
gc()

# Centring
alpha_mean <- colMeans(coef_biv)
Xc <- scale(coef_biv, center = alpha_mean, scale = FALSE)
rm(coef_biv)
gc()

# Covariance of the coefficients and spectral decomposition
# (eigen() returns the eigenvalues in decreasing order)
V   <- crossprod(Xc) / (nrow(Xc) - 1)
eig <- eigen(V, symmetric = TRUE)

lambda   <- eig$values
prop_var <- lambda / sum(lambda)
cum_var  <- cumsum(prop_var)


# ---- 2) Table of explained variance ------------------------------------------

variance_table <- data.frame(
  PC         = seq_along(lambda),
  eigenvalue = lambda,
  prop_var   = prop_var,
  cum_var    = cum_var
)
print(head(variance_table, 10))

write.csv(variance_table, file.path(fod_dir, "pca_variance_table.csv"),
          row.names = FALSE)


# ---- 3) Scores on the first components ---------------------------------------

U      <- eig$vectors[, seq_len(nharm), drop = FALSE]   # (2K x nharm)
scores <- Xc %*% U                                      # (n x nharm)
rm(Xc)
gc()

cat(sprintf("%d components kept: %.2f %% of the variance\n",
            nharm, 100 * cum_var[nharm]))

saveRDS(scores, file.path(fod_dir, "pca_scores.rds"))
saveRDS(
  list(
    alpha_mean     = alpha_mean,      # mean coefficients (2K)
    eigenvectors   = U,               # (2K x nharm)
    eigenvalues    = lambda,          # all eigenvalues
    variance_table = variance_table,
    nharm          = nharm,
    K              = K
  ),
  file.path(fod_dir, "pca_results.rds")
)


# ---- 4) Figures --------------------------------------------------------------

# Explained variance
i_show <- seq_len(min(n_show, length(lambda)))

png(file.path(fig_dir, "pca_explained_variance.png"),
    width = 2400, height = 1800, res = 300)
plot(
  i_show, 100 * cum_var[i_show], type = "b", pch = 16, ylim = c(0, 100),
  xlab = "Principal component", ylab = "Explained variance (%)",
  main = "PCA of the temperature-salinity B-spline coefficients"
)
lines(i_show, 100 * prop_var[i_show], type = "b", pch = 1, lty = 2)
abline(v = nharm + 0.5, lty = 3)
legend(
  "right", legend = c("Cumulative", "Per component", "Components kept"),
  pch = c(16, 1, NA), lty = c(1, 2, 3), bty = "n"
)
dev.off()

# Eigenfunctions: eigenvectors brought back to the depth space
phi   <- unclass(readRDS(file.path(fod_dir, "phi.rds")))   # (n_depth x K)
depth <- readRDS(file.path(fod_dir, "depth.rds"))
stopifnot(ncol(phi) == K, nrow(phi) == length(depth))

eigfun_temp <- phi %*% U[1:K, , drop = FALSE]                # (n_depth x nharm)
eigfun_sal  <- phi %*% U[(K + 1):(2 * K), , drop = FALSE]    # (n_depth x nharm)

pc_cols   <- seq_len(nharm)
pc_labels <- sprintf("PC%d (%.1f %%)", pc_cols, 100 * prop_var[pc_cols])

png(file.path(fig_dir, "pca_eigenfunctions.png"),
    width = 2800, height = 2200, res = 300)
par(mfrow = c(1, 2), oma = c(0, 0, 2, 0))
matplot(
  eigfun_temp, depth, type = "l", lty = 1, lwd = 2, col = pc_cols,
  ylim = rev(range(depth)), xlab = "Temperature component", ylab = "Depth (m)"
)
abline(v = 0, lty = 3)
matplot(
  eigfun_sal, depth, type = "l", lty = 1, lwd = 2, col = pc_cols,
  ylim = rev(range(depth)), xlab = "Salinity component", ylab = "Depth (m)"
)
abline(v = 0, lty = 3)
legend("bottomright", legend = pc_labels, col = pc_cols, lty = 1, lwd = 2,
       bty = "n", cex = 0.8)
mtext("Eigenfunctions of the first principal components", outer = TRUE,
      cex = 1.2)
dev.off()

cat("Results saved in:", fod_dir, "\n")
cat("Figures saved in:", fig_dir, "\n")
