# Acoustic NASC and environmental covariates: data preparation and validation

R scripts that prepare the datasets used to relate acoustic backscatter to the
environment along the OBSAUSTRAL transects of the *Marion Dufresne*
(2018, 2021, 2022 and 2023; 45-90°E, 60-30°S).

They compute the NASC of each elementary sampling unit (ESU) at five
frequencies (18, 38, 70, 120 and 200 kHz), then attach three environmental
descriptors to each ESU:

- **FOD**: Functional Oceanographic Domains, a classification of the
  temperature / salinity profiles (B-splines, PCA, Gaussian mixture model);
- **FTLE**: finite-time Lyapunov exponents;
- **Pigments**: PIGMeANN phytoplankton pigment concentrations and ratios.

A second set of scripts checks the resulting data before modelling: sampling
bias, missing values, distribution shift and spatio-temporal autocorrelation.

## Repository layout

```
config.R                    shared configuration, sourced by every script
.env.example                template for the Copernicus Marine credentials
preprocessing/
  NASC_computation/
    01_filter_crop_sv_profiles.R
    02_concat_sv_profiles_all_years.R
    03_compute_nasc.R
  FOD_computation/
    00_download_temp_salinity.R
    01_filter_and_concat_temp_sal.R
    02_fod_bspline.R
    03_fod_pca.R
    04_fod_clustering.R
    05_colocate_fod_nasc.R
  FTLE_computation/
    01_create_ds_ftle_filtered_all_years.R
    02_colocate_ftle_nasc.R
  PIGMENTS_computation/
    01_create_ds_pigments_filtered_all_years.R
    02_colocate_pigments_nasc.R
data_validation/
  sampling_bias.R
  missing_values.R
  distribution_shift.R
  spatio_temporal_autocorrelation_sv.R
  spatio_temporal_autocorrelation_pigments.R
```

## Setup

1. **R packages**

   ```r
   # preprocessing
   install.packages(c("ncdf4", "abind", "mclust", "ggplot2", "sf",
                      "rnaturalearth", "rnaturalearthdata", "dotenv"))
   # data validation
   install.packages(c("dplyr", "tidyr", "patchwork", "gstat", "ranger",
                      "naniar", "UpSetR", "scales"))
   ```

2. **Paths.** Edit the two roots at the top of `config.R`:

   ```r
   data_root   <- "/path/to/data"      # contains raw/ and processed/
   figures_dir <- "/path/to/figures"
   ```

   Every other path is built from these two. Figures go to
   `<figures_dir>/preprocessing/` and `<figures_dir>/data_validation/`.

3. **Working directory.** Every script starts with `source("config.R")`, so
   it must be run from the project root:

   ```sh
   Rscript preprocessing/NASC_computation/01_filter_crop_sv_profiles.R
   ```

   In RStudio, open the project at the root (or `setwd()` to it) before
   sourcing a script.

4. **Copernicus Marine (only for `00_download_temp_salinity.R`).** Install the
   [Copernicus Marine Toolbox](https://help.marine.copernicus.eu/en/collections/9080063-copernicus-marine-toolbox)
   (the `copernicusmarine` executable), copy `.env.example` to `.env` and fill
   in your credentials. `.env` is ignored by git.

## Preprocessing

### Running order

The NASC and FOD pipelines are independent. The three colocation scripts need
the NASC (`NASC 03`); the FTLE and pigment datasets need the FOD grid
(`FOD 02`), and the FOD colocation needs the FOD classes (`FOD 04`).

```
NASC 01 -> NASC 02 -> NASC 03 ---------------------+
                                                   |
FOD 00 -> FOD 01 -> FOD 02 -> FOD 03 -> FOD 04 ----+-> FOD 05
                      |                            |
                      +-> FTLE 01                  +-> FTLE 02
                      |                            |
                      +-> PIGMENTS 01 -------------+-> PIGMENTS 02
```

### NASC (`preprocessing/NASC_computation/`)

| Script | Reads | Writes |
|---|---|---|
| `01_filter_crop_sv_profiles.R` | `raw/acoustic/*.nc` | `01_sv_cropped_per_year/<freq>kHz/Sv_<year>_<freq>kHz.rds`, NA diagnostic figures |
| `02_concat_sv_profiles_all_years.R` | output of 01 | `02_sv_cropped_all_years/Sv_<years>_<freq>kHz.rds` |
| `03_compute_nasc.R` | output of 02 | `03_NASC_per_esu/NASC_per_esu_<years>_<freq>kHz.rds`, diagnostic figures |

### FOD (`preprocessing/FOD_computation/`)

| Script | Reads | Writes |
|---|---|---|
| `00_download_temp_salinity.R` | Copernicus Marine (GLORYS12V1) | `raw/temperature_salinity/*.nc` (run once per year) |
| `01_filter_and_concat_temp_sal.R` | output of 00 | `01_concat_temp_sal/thetao_so_crop_<years>.nc` |
| `02_fod_bspline.R` | output of 01 | grid, B-spline basis, coefficients and smoothed profiles in `02_03_04_mclust_results_<years>/` |
| `03_fod_pca.R` | output of 02 | PCA scores and results, same folder |
| `04_fod_clustering.R` | outputs of 02 and 03 | model, classes per profile and on the grid, same folder; daily maps |
| `05_colocate_fod_nasc.R` | output of 04, NASC | `05_fod_colocated_nasc/fod_colocated_NASC_per_esu_<years>_<freq>kHz.rds` |

`04_fod_clustering.R` reuses the saved model unless `refit_gmm <- TRUE`. After
a new fit, check `cluster_rename` (numbering of the clusters from South to
North), which is specific to the fitted model.

### FTLE (`preprocessing/FTLE_computation/`)

| Script | Reads | Writes |
|---|---|---|
| `01_create_ds_ftle_filtered_all_years.R` | `raw/ftle/*.nc`, FOD grid | `01_ftle_cropped_all_years/ftle_<years>_cropped.rds` |
| `02_colocate_ftle_nasc.R` | `raw/ftle/*.nc`, NASC | `02_ftle_colocated_nasc/ftle_colocated_NASC_per_esu_<years>_<freq>kHz.rds` |

`02_colocate_ftle_nasc.R` reads the raw daily maps, not the dataset written by
`01`.

### Pigments (`preprocessing/PIGMENTS_computation/`)

| Script | Reads | Writes |
|---|---|---|
| `01_create_ds_pigments_filtered_all_years.R` | `raw/PIGMeANN/daily/*.nc`, FOD grid | `01_pigments_cropped_all_years/pigments_<years>_cropped.rds` |
| `02_colocate_pigments_nasc.R` | output of 01, NASC | `02_pigments_colocated_nasc/pigments_colocated_NASC_per_esu_3x3_<years>_<freq>kHz.rds` |

`<years>` stands for `2018_2021_2022_2023` (`years_tag` in `config.R`).

## Data validation

The scripts of `data_validation/` are independent of each other and can be
run in any order, once the data they read exist. They only write figures (and
a few tables), in `<figures_dir>/data_validation/`.

| Script | Reads | Shows |
|---|---|---|
| `sampling_bias.R` | merged dataset, 38 kHz | number of NASC points per 20 km cell, per year, per month |
| `missing_values.R` | merged dataset, 38 kHz | where and when `nasc`, `fod`, `Chla` and `ftle` are missing; variables missing together |
| `distribution_shift.R` | merged dataset, 38 kHz | shift of the covariates and of the NASC between years and between FOD zones (KS statistic, density overlap, adversarial validation) |
| `spatio_temporal_autocorrelation_sv.R` | Sv of all years (`NASC 02`) | variogram, spatial and temporal correlograms of the mean Sv per ping, per frequency, for all / day / night pings |
| `spatio_temporal_autocorrelation_pigments.R` | pigment dataset (`PIGMENTS 01`) | directional variogram and correlogram, temporal ACF of each pigment |

The **merged dataset** gathers the NASC and the three colocated covariates in
one table per frequency:
`processed/ds_NASC_pig_ftle_fod/NASC_per_esu_pig_ftle_fod_<years>_transect_<freq>kHz.rds`
(`dataset_file()` in `config.R`). The script that builds it is not in this
repository yet. The three scripts that read it keep the
daytime ESU (`day == 3`) and expect the columns `nasc`, `time_nasc`,
`lat_nasc`, `lon_nasc`, `day`, `fod`, `ftle`, `Chla`, `total_pig` and
`<pigment>_totpig`.

## Data folders

```
<data_root>/
  raw/
    acoustic/                 OBSAUSTRAL echo-integration NetCDF, one per year
    temperature_salinity/     GLORYS12V1 thetao / so, one NetCDF per year
    ftle/                     daily FTLE maps, map_<YYYY-MM-DD>*.nc
    PIGMeANN/daily/           daily pigment maps, date YYYYMMDD in the name
  processed/
    NASC_computation/         01_..., 02_..., 03_...
    FOD_computation/          01_..., 02_03_04_..., 05_...
    FTLE_computation/         01_..., 02_...
    PIGMENTS_computation/     01_..., 02_...
    ds_NASC_pig_ftle_fod/     merged dataset, one file per frequency
<figures_dir>/
  preprocessing/
    NASC_computation/         01_diag_NA/, 03_diag_NASC/
    FOD_computation/          02_fod_bspline/, 03_fod_pca/, 04_fod_clustering/
  data_validation/
    sampling_bias/
    missing_values/
    distribution_shift/       figures, CSV tables and results (.rds)
    spatio_temporal_autocorrelation/
```

Each processed sub-folder carries the number of the script that writes it.
Folders are created by the scripts.

## Colocated outputs

The three colocation scripts write one file per frequency, with one row per
ESU, in the same order as the NASC file. They share their first columns, so
that they can be bound to the NASC and to each other:

| Column | Content |
|---|---|
| `time` | time of the ESU (POSIXct, UTC) |
| `lat_sv`, `lon_sv` | position of the ESU |
| `lat_<x>`, `lon_<x>` | matched pixel of the covariate (`fod`, `ftle`, `pig`) |
| ... | covariate: `fod_cluster`, `ftle`, or pigment concentrations and ratios |

FOD and FTLE take the nearest pixel on the same day; pigments are averaged
over a 3 x 3 pixel window centred on the nearest pixel.

FOD classes: 1-6 are the clusters C1 to C6 (South to North), 7-13 the
transitions between two clusters, 0 a transition that is not in the list
(see the header of `04_fod_clustering.R`).

## Conventions

- **Header.** Every script starts with the same block: file name, one-line
  title, `Input`, `Steps`, `Output`.
- **Sections.** `# ---- Title ----`, in this order: libraries, configuration,
  functions, numbered steps matching the header.
- **Paths.** No absolute path in the scripts. Directories and the files
  exchanged between scripts are defined once in `config.R`; scripts use
  `in_dir` / `in_file`, `out_dir` / `out_file` and `fig_dir` for their own
  input, output and figures.
- **File names.** Preprocessing scripts are numbered in running order within
  their folder; data-validation scripts are independent and not numbered.
  Output files of the preprocessing end with `<years>_<freq>kHz`.
- **Language.** Comments, console messages, figure labels and output file
  names are in English.
