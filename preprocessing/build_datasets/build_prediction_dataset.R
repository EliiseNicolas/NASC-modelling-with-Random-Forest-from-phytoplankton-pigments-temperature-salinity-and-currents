# ============================================================
# Extraction FTLE, pigments et FOD pour TOUTES les dates
# -> un seul objet / un seul fichier RDS en sortie
# Version entièrement vectorisée : aucune boucle sur les dates.
#
# Les trois jeux de données sont ramenés sur la grille des pigments (plus
# proche voisin pour FTLE et FOD) et restreints à leurs dates communes.
#
# Entrées : ftle_2018_2021_2022_2023_cropped.rds (date, lon, lat, ftle),
#           pigments_2018_2021_2022_2023_crop.rds (crop_pigments_to_fod.R),
#           cluster_transition_map_renamed.rds, lon.rds, lat.rds, time.rds (FOD)
# Sortie  : out_file, une liste contenant
#             date  dates communes aux trois jeux de données
#             lon   longitudes de la grille des pigments
#             lat   latitudes de la grille des pigments
#             ftle  array [date, lon, lat]
#             pig   liste d'arrays [date, lon, lat] : les 9 pigments,
#                   chla_total et les ratios <pigment>_totpig
#             fod   array [lon, lat, date]  (ordre des dimensions différent
#                   de ftle et pig)
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------

out_file <- "F:/data_elise/prediction_ds/ds_ftle_pig_fod_ALL_DATES.rds"

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

path_ftle <- "F:/data_elise/ftle/ftle_2018_2021_2022_2023_cropped.rds"
path_pig  <- "F:/data_elise/pigmeann/pigments_2018_2021_2022_2023_crop.rds"

path_fod_clusters <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/transitions_upgraded/cluster_transition_map_renamed.rds"
path_fod_lon  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lon.rds"
path_fod_lat  <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/lat.rds"
path_fod_time <- "F:/data_elise/fod_elise_2018_2021_2022_2023/mclust_results_2018_2021_2022_2023/time.rds"

# ------------------------------------------------------------
# Chargement
# ------------------------------------------------------------

ftle_grid <- readRDS(path_ftle)
pigs_grid <- readRDS(path_pig)

fod_clusters <- readRDS(path_fod_clusters)
fod_lon  <- readRDS(path_fod_lon)
fod_lat  <- readRDS(path_fod_lat)
fod_time <- readRDS(path_fod_time)
fod_date <- as.Date(fod_time)

list_pigs <- c("Chla", "Per", "But", "Fuco", "Hex", "Allo", "Zea", "Chlb", "DvChla")

# ------------------------------------------------------------
# Grille de référence = pigments (1080 x 720)
# ------------------------------------------------------------

grid_lon <- pigs_grid$lon
grid_lat <- pigs_grid$lat

# ------------------------------------------------------------
# Indices de reprojection spatiale (nearest neighbor), 1 seule fois
# ------------------------------------------------------------

#' Indice du plus proche voisin de chaque coordonnée cible.
#'
#' @param target coordonnées de la grille de référence (lon ou lat)
#' @param source coordonnées de la grille d'origine, sur le même axe
#' @return pour chaque valeur de target, l'indice de la valeur de source la
#'   plus proche
nearest_index <- function(target, source) {
  vapply(target, function(x) which.min(abs(source - x)), integer(1))
}

idx_lon_ftle <- nearest_index(grid_lon, ftle_grid$lon)
idx_lat_ftle <- nearest_index(grid_lat, ftle_grid$lat)

idx_lon_fod <- nearest_index(grid_lon, fod_lon)
idx_lat_fod <- nearest_index(grid_lat, fod_lat)

# ------------------------------------------------------------
# Dates communes aux trois jeux de données
# ------------------------------------------------------------

dates_communes <- sort(as.Date(Reduce(intersect, list(ftle_grid$date, pigs_grid$date, fod_date))))
n_date <- length(dates_communes)
cat(n_date, "dates communes.\n")

idx_date_ftle <- match(dates_communes, ftle_grid$date)
idx_date_pig  <- match(dates_communes, pigs_grid$date)
idx_date_fod  <- match(dates_communes, fod_date)

# ------------------------------------------------------------
# Extraction + reprojection en UNE SEULE opération d'indexation
# par variable (produit cartésien des indices sur les 3 dimensions)
# ------------------------------------------------------------

# --- FTLE : natif [ndate_tot, 901, 600] -> [n_date, 1080, 720] ---
ftle_all <- ftle_grid$ftle[idx_date_ftle, idx_lon_ftle, idx_lat_ftle]

# --- FOD : natif [541, 361, ndate_tot] -> [1080, 720, n_date] ---
fod_all <- fod_clusters[idx_lon_fod, idx_lat_fod, idx_date_fod]

# --- Pigments : déjà sur la grille de référence, juste sélection des dates ---
pig_all <- lapply(list_pigs, function(p) {
  pigs_grid[[paste0("c_cond_", p)]][idx_date_pig, , ]   # [n_date, 1080, 720]
})
names(pig_all) <- list_pigs

# ------------------------------------------------------------
# chla_total = Chla seule (DvChla volontairement exclue)
# ------------------------------------------------------------

pig_all[["chla_total"]] <- pig_all[["Chla"]]

# ------------------------------------------------------------
# Somme des concentrations de tous les pigments SAUF Chla
# (calculee une seule fois, reutilisee comme denominateur commun)
# ------------------------------------------------------------

pigs_sans_chla <- setdiff(list_pigs, "Chla")

pig_sum_others <- Reduce(`+`, pig_all[pigs_sans_chla])   # somme elementwise, meme dims [n_date, 1080, 720]

# ------------------------------------------------------------
# Ratio de chaque pigment (les 9 pigments bruts, Chla incluse) sur
# la somme des pigments hors Chla, ajoutes dans pig_all avec un
# suffixe "_totpig" (ex : chla_totpig, but_totpig, fuco_totpig, ...).
# Vectorise via lapply : une division element par element par
# pigment, aucune boucle sur les dates/pixels.
# ------------------------------------------------------------

ratio_list <- lapply(list_pigs, function(p) pig_all[[p]] / pig_sum_others)
names(ratio_list) <- paste0(tolower(list_pigs), "_totpig")

pig_all <- c(pig_all, ratio_list)

# ------------------------------------------------------------
# Assemblage final : un seul objet pour toutes les dates
# ------------------------------------------------------------

all_ds <- list(
  date = dates_communes,   # vecteur de n_date dates
  lon  = grid_lon,
  lat  = grid_lat,
  ftle = ftle_all,         # array [n_date, 1080, 720]
  pig  = pig_all,          # liste : 9 pigments bruts + chla_total + ratios *_totpig
  fod  = fod_all           # array [1080, 720, n_date]
)

str(all_ds, max.level = 2)
cat("Taille en mémoire :", format(object.size(all_ds), units = "GB"), "\n")

# ------------------------------------------------------------
# Sauvegarde (un seul fichier)
# ------------------------------------------------------------

saveRDS(all_ds, out_file)
cat("Sauvegardé dans :", out_file, "\n")