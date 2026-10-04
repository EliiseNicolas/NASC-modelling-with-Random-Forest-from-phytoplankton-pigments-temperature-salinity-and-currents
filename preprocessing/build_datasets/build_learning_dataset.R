# ============================================================
# Assemblage du dataset d'entrée des arbres de régression
#
# Pour chaque fréquence acoustique, réunit dans un seul data.frame (une ligne
# par ESU) le NASC et les variables colocalisées par les scripts
# colocate_*_with_nasc.R : classe FOD, pigments et FTLE. Les quatre fichiers
# ont les mêmes lignes dans le même ordre ; le script le vérifie et s'arrête
# sinon.
#
# Entrées (une par fréquence) :
#   NASC_per_ESU_2018_2021_2022_2023_<freq>kHz.rds
#   NASC_per_esu_FOD_with_transitions_cluster_match_2018_2021_2022_2023_<freq>kHz.rds
#   NASC_per_esu_pig_conc_ratio_9_1d_2018_2021_2022_2023_<freq>kHz.rds
#   ftle_colocated_with_NASC_per_esu_2018_2021_2022_2023_<freq>kHz.rds
# Sortie :
#   ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_<freq>kHz_mask9.rds
#   avec les colonnes
#     time_nasc, lat_nasc, lon_nasc  date et position de l'ESU
#     day                            3 = jour, 1 = nuit
#     nasc                           NASC
#     lat_fod, lon_fod, fod          point de grille FOD associé et classe FOD
#                                    (texte ; "NA" si absente)
#     lat_pig, lon_pig, ...          pixel de pigments associé, puis toutes
#                                    les variables du fichier de pigments
#     lat_ftle, lon_ftle, ftle       pixel FTLE associé et valeur de FTLE
# ============================================================

# Libraries
library(rpart)

# Global variables
rm(list = ls())

freqs <- c(120, 200)# ,
for (freq in freqs) {
  # ds NASC per ESU
  path_pig <- paste0("F:/data_elise/pigmeann/pigs_colocated_NASC_per_esu/NASC_per_esu_pig_conc_ratio_9_1d_2018_2021_2022_2023_", freq, "kHz.rds")
  path_fod <- paste0("F:/data_elise/fod_elise_2018_2021_2022_2023/fod_colocated_nasc_2018_2021_2022_2023_transect/fod_colocated_NASC_per_esu/NASC_per_esu_FOD_with_transitions_cluster_match_2018_2021_2022_2023_", freq, "kHz.rds")
  path_nasc <- paste0("F:/data_elise/NASC/NASC_all_ESU/NASC_per_ESU_2018_2021_2022_2023_", freq, "kHz.rds")
  path_ftle <- paste0("F:/data_elise/ftle/ftle_colocated_transect/ftle_colocated_NASC_per_esu/ftle_colocated_with_NASC_per_esu_2018_2021_2022_2023_", freq, "kHz.rds")
  print(path_nasc)
  
  pig <- readRDS(path_pig)
  fod <- readRDS(path_fod)
  nasc <- readRDS(path_nasc)
  ftle <- readRDS(path_ftle)
  
  str(pig)
  str(fod)
  str(nasc)
  str(ftle)
  print(length(unique(as.Date(nasc$time))))
  
  # --------------- Creation du dataset final
  regression_ds <- data.frame(
    time_nasc = nasc$time,
    lat_nasc = nasc$lat,
    lon_nasc = nasc$lon,
    day = nasc$day,
    nasc = nasc$NASC
  )
  
  #----------------  Vérification alignement
  # (arrêt du script si une ESU ne correspond pas entre les fichiers)
  stopifnot(
    # verif fod
    regression_ds$time_nasc == fod$time_nasc,
    regression_ds$lat_nasc == fod$lat_nasc,
    regression_ds$lon_nasc == fod$lon_nasc,
    # verif pig
    regression_ds$lat_nasc == pig$lat_sv,
    regression_ds$lon_nasc == pig$lon_sv,
    # verif ftle
    regression_ds$time_nasc == ftle$time,
    regression_ds$lat_nasc == ftle$lat_sv,
    regression_ds$lon_nasc == ftle$lon_sv
  )
  
  # ------------------- Ajout des variables ftle, pig et fod dans le dataset final
  # ---- FOD
  # Ajout coordonnées FOD + cluster
  print(sum(!is.na(fod$lat_fod)))
  regression_ds$lat_fod <- fod$lat_fod
  regression_ds$lon_fod <- fod$lon_fod
  regression_ds$fod <- format(fod$fod_cluster)
  
  # ----- Pigments
  vars_pig <- setdiff(
    names(pig),
    c("time", "lat_sv", "lon_sv")
  )
  print(vars_pig)
  regression_ds[vars_pig] <- pig[vars_pig]
  
  # ---- FTLE
  regression_ds$lat_ftle <- ftle$lat_ftle
  regression_ds$lon_ftle <- ftle$lon_ftle
  regression_ds$ftle <- ftle$ftle
  
  # ---------------------- VERIFS
  str(regression_ds)
  
  # verif années contenues
  print(unique(format(regression_ds$time_nasc, "%Y"))) # 2018 2021 2022 2023
  
  # différence de lat/lon entre nasc fod et pig
  dlat_fod <- regression_ds$lat_fod - regression_ds$lat_nasc
  dlon_fod <- regression_ds$lon_fod - regression_ds$lon_nasc
  print(summary(dlat_fod))
  print(summary(dlon_fod)) # OK
  
  dlat_pig <- regression_ds$lat_pig - regression_ds$lat_nasc
  dlon_pig <- regression_ds$lon_pig - regression_ds$lon_nasc
  print(summary(dlat_pig))
  print(summary(dlon_pig)) # OK
  
  dlat_ftle <- regression_ds$lat_ftle - regression_ds$lat_nasc
  dlon_ftle <- regression_ds$lon_ftle - regression_ds$lon_nasc
  print(summary(dlat_ftle))
  print(summary(dlon_ftle)) # OK
  
  # ------------ SAVE
  saveRDS(regression_ds, paste0("F:/data_elise/ds_NASC_pig_ftle_fod/ds_NASC_per_esu_all/ds_NASC_per_esu_pig_ftle_fod_2018_2021_2022_2023_transect_", freq, "kHz_mask9.rds"))
}