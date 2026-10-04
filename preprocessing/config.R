# ==============================================================================
# config.R
#
# Shared configuration of the data-preparation scripts
#
# Sourced at the top of every script with `source("config.R")`: the scripts
# must be run with the project root (the folder of this file) as working
# directory.
#
# Defines : roots of the data and figure folders (machine-specific)
#           years, acoustic frequencies and study area
#           raw and processed data directories
#           files exchanged between scripts (written by one, read by another)
# ==============================================================================


# ---- Roots (machine-specific: edit these two lines) --------------------------

data_root <- "/run/media/mmolinet/DATA1/data_elise_stage"
fig_root  <- "~/Elisou/part_V_final_model_clean_rep/figures/data_analysis"


# ---- Campaigns ---------------------------------------------------------------

years     <- c("2018", "2021", "2022", "2023")
years_tag <- paste(years, collapse = "_")   # used in file and folder names
freqs     <- c(18, 38, 70, 120, 200)        # kHz


# ---- Study area --------------------------------------------------------------

lat_min <- -60; lat_max <- -30
lon_min <-  45; lon_max <-  90


# ---- Raw data directories ----------------------------------------------------

raw_acoustic_dir <- file.path(data_root, "raw", "acoustic")
raw_temp_sal_dir <- file.path(data_root, "raw", "temperature_salinity")
raw_ftle_dir     <- file.path(data_root, "raw", "ftle")
raw_pigments_dir <- file.path(data_root, "raw", "PIGMeANN", "daily")


# ---- Processed data directories ----------------------------------------------
# One folder per pipeline, one sub-folder per script (same number as the script)

nasc_root     <- file.path(data_root, "processed", "NASC_computation")
fod_root      <- file.path(data_root, "processed", "FOD_computation")
ftle_root     <- file.path(data_root, "processed", "FTLE_computation")
pigments_root <- file.path(data_root, "processed", "PIGMENTS_computation")

sv_per_year_dir  <- file.path(nasc_root, "01_sv_cropped_per_year")
sv_all_years_dir <- file.path(nasc_root, "02_sv_cropped_all_years")
nasc_dir         <- file.path(nasc_root, "03_NASC_per_esu")

temp_sal_dir <- file.path(fod_root, "01_concat_temp_sal")
fod_dir      <- file.path(
  fod_root, paste0("02_03_04_mclust_results_", years_tag)
)
fod_nasc_dir <- file.path(fod_root, "05_fod_colocated_nasc")

ftle_dir      <- file.path(ftle_root, "01_ftle_cropped_all_years")
ftle_nasc_dir <- file.path(ftle_root, "02_ftle_colocated_nasc")

pigments_dir      <- file.path(pigments_root, "01_pigments_cropped_all_years")
pigments_nasc_dir <- file.path(pigments_root, "02_pigments_colocated_nasc")


# ---- Files exchanged between scripts -----------------------------------------

# Cropped Sv profiles of one year and one frequency (NASC 01 -> NASC 02)
sv_year_file <- function(year, freq) {
  file.path(sv_per_year_dir, paste0(freq, "kHz"),
            paste0("Sv_", year, "_", freq, "kHz.rds"))
}

# Cropped Sv profiles of all years, one frequency (NASC 02 -> NASC 03)
sv_all_years_file <- function(freq) {
  file.path(sv_all_years_dir, paste0("Sv_", years_tag, "_", freq, "kHz.rds"))
}

# NASC per ESU, one frequency (NASC 03 -> the three colocation scripts)
nasc_file <- function(freq) {
  file.path(nasc_dir, paste0("NASC_per_esu_", years_tag, "_", freq, "kHz.rds"))
}

# Cropped and concatenated temperature / salinity (FOD 01 -> FOD 02)
temp_sal_file <- file.path(
  temp_sal_dir, paste0("thetao_so_crop_", years_tag, ".nc")
)

# FTLE at the FOD dates, cropped to the FOD area (FTLE 01)
ftle_file <- file.path(ftle_dir, paste0("ftle_", years_tag, "_cropped.rds"))

# Pigments at the FOD dates, cropped to the FOD area (PIGMENTS 01 -> 02)
pigments_file <- file.path(
  pigments_dir, paste0("pigments_", years_tag, "_cropped.rds")
)
