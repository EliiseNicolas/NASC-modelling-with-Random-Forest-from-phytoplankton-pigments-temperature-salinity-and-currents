# ==============================================================================
# 03_compute_nasc.R
#
# NASC computation, step 3/3: NASC per ESU
#
# Input  : Sv profiles of all years, per frequency
#          (02_concat_sv_profiles_all_years.R)
#          <sv_all_years_dir>/Sv_<years_tag>_<freq>kHz.rds (one profile per ESU)
#
# Steps  : for each frequency
#            1) remove the profiles with too many NA (> max_na_fraction)
#            2) diagnostic figure: remaining NA by depth and year
#            3) fill the remaining NA along depth (linear interpolation,
#               nearest value beyond the first / last valid depth)
#            4) integrate over depth and compute the NASC of each profile
#            5) diagnostic figures: NASC over time, per year
#
# Output : <nasc_dir>/NASC_per_esu_<years_tag>_<freq>kHz.rds,
#          a data frame with time, lat, lon, day, NASC (one row per ESU)
#          in <fig_dir>
#            diag_NA_by_depth_per_esu_<freq>kHz.png
#            diag_NASC_time_per_esu_<freq>kHz_<year>.png
# ==============================================================================

library(ggplot2)


# ---- Configuration -----------------------------------------------------------

source("config.R")   # years, freqs, sv_all_years_file(), nasc_file()

fig_dir <- file.path(fig_root, "NASC_computation", "03_diag_NASC")

# Profiles with more than this fraction of NA are removed
max_na_fraction <- 0.5

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- NASC per frequency ------------------------------------------------------

for (freq in freqs) {
  cat("\n---", freq, "kHz ---\n")

  sv_data  <- readRDS(sv_all_years_file(freq))
  profiles <- sv_data$profiles   # (n_profiles x n_depth), Sv in dB
  depth    <- sv_data$depth

  # ---- 1) Remove profiles with too many NA -----------------------------------

  n_depth <- ncol(profiles)
  n_na    <- rowSums(is.na(profiles))   # number of NA in each profile

  cat(sprintf("Profiles entirely NA   : %.2f %%\n", 100 * mean(n_na == n_depth)))
  cat(sprintf("Profiles with >20%% NA  : %.2f %%\n", 100 * mean(n_na > 0.2 * n_depth)))

  keep <- n_na <= max_na_fraction * n_depth
  cat("Profiles removed       :", sum(!keep), "out of", length(keep), "\n")

  profiles <- profiles[keep, , drop = FALSE]
  time     <- sv_data$time[keep]
  lat      <- sv_data$lat[keep]
  lon      <- sv_data$lon[keep]
  day      <- sv_data$day[keep]
  year     <- format(as.Date(time), "%Y")

  # ---- 2) Where are the remaining NA? ----------------------------------------

  na_df <- do.call(rbind, lapply(years, function(y) {
    profiles_year <- profiles[year == y, , drop = FALSE]
    if (nrow(profiles_year) == 0) return(NULL)

    n_na_depth <- colSums(is.na(profiles_year))
    data.frame(
      depth      = depth,
      n_NA       = n_na_depth,
      pct_NA     = 100 * n_na_depth / nrow(profiles_year),
      n_profiles = nrow(profiles_year),
      year       = y
    )
  }))
  rownames(na_df) <- NULL

  p_na <- ggplot(na_df, aes(x = depth, y = pct_NA, color = year)) +
    geom_line() +
    geom_point(size = 0.5) +
    labs(
      x = "Depth (m)", y = "Profiles with NA (%)", color = "Year",
      title = "Remaining NA by depth",
      subtitle = paste0(
        "After removing profiles with >", 100 * max_na_fraction, "% NA - ",
        freq, " kHz"
      )
    ) +
    theme_minimal()
  print(p_na)

  ggsave(
    filename = file.path(
      fig_dir, paste0("diag_NA_by_depth_per_esu_", freq, "kHz.png")
    ),
    plot = p_na, width = 10, height = 6, dpi = 300, units = "in"
  )

  # ---- 3) Fill the remaining NA along depth ----------------------------------

  cat("NA before interpolation:", sum(is.na(profiles)), "\n")

  for (i in seq_len(nrow(profiles))) {
    sv_i <- profiles[i, ]
    ok   <- !is.na(sv_i)

    # rule = 2: beyond the valid range, use the nearest valid value
    if (sum(ok) >= 2 && !all(ok)) {
      profiles[i, !ok] <- approx(
        x = depth[ok], y = sv_i[ok], xout = depth[!ok],
        method = "linear", rule = 2
      )$y
    }
  }

  cat("NA after interpolation :", sum(is.na(profiles)), "\n")

  # ---- 4) NASC ---------------------------------------------------------------

  sv_linear  <- 10^(profiles / 10)                            # dB -> linear
  depth_step <- mean(diff(depth))                             # m
  sa         <- rowSums(sv_linear, na.rm = TRUE) * depth_step # depth integral
  NASC       <- 4 * pi * 1852^2 * sa

  cat("Depth step (m):", depth_step, "\n")
  print(summary(NASC))

  nasc_df <- data.frame(
    time = time,
    lat  = lat,
    lon  = lon,
    day  = day,
    NASC = NASC
  )
  str(nasc_df)

  out_file <- nasc_file(freq)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(nasc_df, out_file)

  # ---- 5) NASC over time, per year -------------------------------------------

  for (y in years) {
    df_year <- nasc_df[year == y, ]
    if (nrow(df_year) == 0) next

    p_nasc <- ggplot(df_year, aes(x = time, y = log(NASC))) +
      geom_point(size = 1) +
      labs(
        x = "Time", y = "log(NASC)",
        title = paste("NASC over time in", y),
        subtitle = paste0("Transect data, ", freq, " kHz"),
        caption = paste0(
          "Profiles with >", 100 * max_na_fraction,
          "% NA removed, linear interpolation + extrapolation"
        )
      ) +
      theme_minimal() +
      theme(plot.caption = element_text(hjust = 0))
    print(p_nasc)

    ggsave(
      filename = file.path(
        fig_dir, paste0("diag_NASC_time_per_esu_", freq, "kHz_", y, ".png")
      ),
      plot = p_nasc, width = 10, height = 6, dpi = 300, units = "in"
    )
  }
}
