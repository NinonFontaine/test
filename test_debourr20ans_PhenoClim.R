library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)
library(terra)

# OBJECTIFS ----
# Concernant les données PhenoClim et leur analyse, certains éléments ont déjà été explorés (notamment par Marjorie Bison).
# Les analyses se sont principalement concentrées sur les Alpes, la phase de débourrement, et sur une première série de données, de 
# 2004 à 2016. Les phases de floraison et feuillaison restent à explorer, tout comme le cas des autres massifs (Pyrénées, Massif Central,
# Vosges...).
# 
# Par rapport aux analyses de débourrement dans les Alpes, il y a désormais 8 années supplémentaires de données, mais les premières 
# explorations de ces données supplémentaires tendent à montrer que les tendances d'avancement du débourrement ne se retrouvent plus 
# avec les 8 ans supplémentaires, peut-être à cause d'une grande variabilité interzone avec des outliers.
# Cela reste à explorer, avec plusieurs pistes : 
# - Piste 1 : éliminer les outliers (Vitasse et al 2017)
# - Piste 2 : contraindre la distribution en jouant sur les priors dans les méthodes bayésiennes
# - Piste 3 : intégrer des données de température comme variable explicative, plus directement que le paramètre 'altitude'


# DATA ----
#---- Fonds carto
cle_ggmap ="" # à récupérer
ggmap::register_google(key=cle_ggmap)
# Alpes
map_base_Alp <- get_map(location = c(lon = 6.5, lat = 45.3), zoom = 7,
                     maptype = "terrain", scale = 2)
# Alpes du Nord
map_base_AlpN <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
                         maptype = "terrain", scale = 2)
# Massif du Mont Blanc
map_base_MtBlc <- get_map(location = c(lon = 6.8, lat = 45.95), zoom = 11,
                    maptype = "terrain", scale = 2)

#---- Données Phénoclim
phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
# Données obtenues avec le script 1_mise_en_forme_BDD.R



#
# Pour la piste 3, pour corriger la variabilité interzone par la température... ----
#*---- Données climatiques potentielles : comparaison des couvertures spatiales et temporelles ----
#
#         DONNEES REANALYSEES
# - Copernicus CERRA : données météo Europe, jusqu'à 2022, à une résoluton de 5.5 km 
#                     (https://climate.copernicus.eu/copernicus-regional-reanalysis-europe-cerra)
# - ERA 5-Land : données météo Monde, en temps réel j-5, à une résolution de 9 km
#                     (https://climate.copernicus.eu/climate-reanalysis)
# - CHELSA : données météo Europe, 1979-2013, à une résolution de 1 km
#                     (https://chelsa-climate.org/timeseries/)
# - SAFRAN - MeteoFrance : données météo France, 1958-2024, à une résolution de 8 km (calcul par massif x altitude : il faudrait 
#                     récupérer la base de ce calcul... IsaB ?)
#                     (https://meteo.data.gouv.fr/datasets/donnees-changement-climatique-sim-quotidienne/)
#
#         DONNEES MESUREES
# - Stations Phénoclim : données calées sur la période d'étude, avec une méthode pour boucher les 'manques' liés à des bugs techniques
#                     mais couverture spatiale relativement limitée (même si il y a qq stations en altitude)...



#*---- Couverture spatiale des données dont la couverture temporelle est suffisante : SAFRAN + Phénoclim ----
#
# grille_SAFRAN = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/coordonnees_grille_safran_lambert-2-etendu.csv",
#                          sep=";", dec=",")
# # ggmap(map_base_Alp) +
# #   geom_point(data=grille_SAFRAN[grille_SAFRAN$LAMBX..hm. <1200 & grille_SAFRAN$LAMBX..hm.>1000,], aes(x=LON_DG, y=LAT_DG))
# ggmap(map_base_MtBlc) +
#   geom_point(data=grille_SAFRAN, aes(x=LON_DG, y=LAT_DG))
# 
# stt_Phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Stations-phenoclim/analyse/coordonnees_crea.csv",
#                          sep=";")
# stt_Phenoclim_sp = vect(stt_Phenoclim, geom=c("x_crea", "y_crea"), crs="+init=epsg:4326")
# stt_Phenoclim_buff5km = project(aggregate(buffer(project(stt_Phenoclim_sp, "epsg:2154"), 5000)), "epsg:4326")
# 
# ggmap(map_base_Alp) +
#   # geom_polygon(data=geom(stt_Phenoclim_buff5km), aes(x=x, y=y, group = part), fill="red", alpha=0.7)+
#   geom_point(data=phenoclim, aes(x=coord_x_4326, y=coord_y_4326), col="orange")+
#   geom_point(data=stt_Phenoclim, aes(x=x_crea, y=y_crea))
# 
# # La résolution spatiale des stations phénoclim et des points SAFRAN ne semble pas suffisante pour couvrir l'ensemble des observations
# # phénologiques, pour pouvoir inclure le paramètre température dans les modèles phénologiques...
# # Deux solutions sont envisagées : 
# # 1) interpoler les mesures de température des stations Phénoclim en fonction de l'altitude (cf Asse et al. 2018)
# # 2) utiliser les températures du modèle S2M (MétéoFrance - SAFRAN) par massif x altitude x pente x exposition 

#*---- Interpolation des données existantes pour améliorer la couverture spatiale des informations 'température' ----
#*-------- 1) Utilisation des données des stations PhénoClim interpolées pour tous les sites d'observation ----

# En se basant sur le machine learning, sur le même principe que son utilisation pour boucher les 'gros trous' dans les données des 
# stations météo (bugs techniques), on peut obtenir des valeurs de température pour tous les sites en se basant sur leurs coordonnées 
# x, y, z (cf scripts Geoffrey Klein dans le projet R "Stations-phenoclim.Rproj").
# NOTE : ce qui m'étonne, c'est de pouvoir reconstruire ces données de température sans prendre en compte l'exposition des sites, 
#         seulement leur altitude et leur localisation (parce qu'on se base sur les stations PhénoClim les plus proches).
#*         => À VOIR AVEC GEOFFREY ?? ----
#
# Ce qui nous intéresse dans ces données c'est :
# - le sensor 3 (T_air à 2m de haut)
# - les températures moyennes, minimales et maximales quotidiennes (qu'on pourra utiliser ensuite pour calculer une température moyenne
#   sur une période donnée, comme le printemps par exemple, ou des GDD)

# stt_Phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Stations-phenoclim/analyse/coordonnees_crea.csv",
#                          sep=";")
donnees_completes_2023 = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Stations-phenoclim/data/fichiers_reconstruits_annee/2023_reconstruit.csv")

test = donnees_completes_2023 %>% group_by(station_name, month, day, x_crea, y_crea, z_crea) %>%
  summarise(Tair_moy = mean(sensor_value4),
            Tair_min = min(sensor_value4),
            Tair_max = max(sensor_value4))
test$date = as.Date(paste(test$day, test$month, "2023", sep="/"), format = "%d/%m/%Y")

ggplot(test) + 
  geom_line(aes(x=date, y=Tair_moy, group=station_name), col="black") +
  geom_line(aes(x=date, y=Tair_min, group=station_name), col="blue", linetype=2)+
  geom_line(aes(x=date, y=Tair_max, group=station_name), col="red", linetype=2) + 
  xlim(as.Date(c("15/02/2023","15/04/2023"), format = "%d/%m/%Y"))

#*-------- 2) Utilisation des données SAFRAN par massif x topographie ----
# PISTE 2 : DONNÉES SAFRAN 
# # Si on veut les données SAFRAN SIM (/!\ la résolution à 8km n'est pas optimale) :
# 
# sim_2020a24 = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/QUOT_SIM2_previous-2020-202412.csv", sep=";")
# # https://www.aeris-data.fr/catalogue/?keywords=%5B%22safran%22%5D&currentSelection=865730e8-edeb-4c6b-ae58-80f95166509b
# # https://meteo.data.gouv.fr/datasets/donnees-changement-climatique-sim-quotidienne/ 
# 
# # Sélection des données pour les Alpes (8000 < X < 10500, 16000 < Y < 22000)
# sim_2020a24_Alp = sim_2020a24[sim_2020a24$LAMBX <10500 & sim_2020a24$LAMBX>8000 & 
#                                 sim_2020a24$LAMBY<22000 & sim_2020a24$LAMBY>16000, ]
# pts = grille_SAFRAN[grille_SAFRAN$LAMBX..hm. <10500 & grille_SAFRAN$LAMBX..hm.>8000 & 
#                       grille_SAFRAN$LAMBY..hm. <22000 & grille_SAFRAN$LAMBY..hm.>16000,]
# pts$id = 1:nrow(pts)
# write.csv(pts, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/pts_SAFRAN_Alps.csv", row.names = F)
# 
# sim_2020a24_Alp$id = apply(sim_2020a24_Alp[,c("LAMBX", "LAMBY")], 1, 
#                            FUN = function(X){output = pts$id[pts$LAMBX..hm. == X[1] & pts$LAMBY..hm.==X[2]]})
# write.csv(sim_2020a24_Alp, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/sim2020a24_pts_SAFRAN_Alps.csv", row.names = F)


# A la place, on peut utiliser les données SAFRAN S2M, avec une simulation par massif et classe d'altitude
# Dans un premier temps, on attribue un massif SAFRAN à chaque observation Phénoclim, et sa classe d'altitude x exposition
# (en gardant en tête le fait qu'il va falloir interpoler les températures pour éviter les valeurs 'en escalier' liées à ce système
# de classe)

massifs_safran = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/alp_allslopes/massifs_alpes_2154.shp")
MNT = rast("/Volumes/RESSOURCES_SIG/MNT/MNT_Alpes/MNT_Alpes_25m_L93.asc")
pente = terrain(MNT, v="slope", unit="degrees")
aspect = terrain(MNT, v="aspect", unit="degrees")
stacktopo = c(MNT, pente, aspect)


sites_cl_safran = phenoclim[phenoclim$nom_massif=="Alpes",] %>% distinct(id_base_site, coord_x_2154, coord_y_2154, altitude)

plot(stt_Phenoclim_buff5km)
points(phenoclim$coord_x_4326[phenoclim$nom_massif=="Alpes"], phenoclim$coord_y_4326[phenoclim$nom_massif=="Alpes"], cex=0.3, col="red")
lines(project(massifs_safran, "epsg:4326"), col="darkgreen")

# En fait même les données SAFRAN ne couvriront pas toute la zone, puisque de nombreuses observations Phénoclim sont en dehors de la zone
# 'Alpes' où sont définis les massifs SAFRAN... Soit on leur associe le massif le plus proche, soit on trouve d'autres sources de données T

sites_cl_safran[,c("massif_num","massif_nom")] = extract(massifs_safran, sites_cl_safran[,c("coord_x_2154","coord_y_2154")])[,c("massif_num","nom")]
sites_cl_safran[,c("pente","expo")] = extract(stacktopo[[2:3]], sites_cl_safran[,c("coord_x_2154","coord_y_2154")])[,2:4] #on ne prend pas l'altitude qu'on a déjà avec les infos des sites
write.csv(sites_cl_safran, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_sites_SAFRAN.csv")


# ggplot(sites_cl_safran) +
#   # geom_polygon(data=geom(massifs_safran), aes(x=x, y=y, group=part), fill="lightgray", alpha=0.5)+
#   geom_point(aes(x=coord_x_2154, y=coord_y_2154, color=massif_nom))+
#   coord_equal() + xlim(760000, 1100000)



