library(ggplot2)
library(tidyr)
library(statip)
library(terra)
library(sf)
library(RJSONIO)
library(rinat)
library(leaflet)
library(ggmap)
library(lubridate)
library(dplyr)
library(readxl)
library(readODS)
# library(RStoolbox)
# library(ade4)
library(cluster)
library(StatMatch)
library(FactoMineR)
library(factoextra)

# ----------------------------------------------------------------------------------------*

# Le protocole Floraison d'Altitude évolue : deux développements envisagés sont :
#
# 1) de proposer à d'autres citoyen.nes scientifiques de rejoindre le programme, notamment les observateur.rices iNaturalist
#       => Pour cela, on souhaite créer un 'projet' iNaturalist, basé sur une zone donnée, pour que les personnes ayant déjà fait / faisant
#           des observations dans cette zone soient contactées pour participer à SPOT. Il faut donc créer un shapefile sous forme de buffer
#           autour des points d'observation SPOT.
#
# 2) d'ajouter des sites de quantification de la production de fruits de myrtille, pour évaluer les effets de l'environnement abiotique
#     sur la production de fruits (les effets biotiques, et notamment la pollinisation, ne seront suivis que sur certains sites avec 
#     l'écoacoustique)
#       => Pour cela, on souhaite échantillonner une diversité de contextes abiotiques où est présente la myrtille, afin d'avoir un plan
#           d'échantillonnage incluant ces différentes conditions.


# ----------------------------------------------------------------------------------------*


# PROJET iNATURALIST (développement n°1) ----

# Localisation des points SPOT Floraison d'Altitude :
FlodAlti = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/FloraisonAltitude/data_long.csv") # Fichier issu de GeoNature (2025_01_13_19h23m43_Floraison_daltitude_private.csv) mis en forme
sites_FlodAlti = FlodAlti %>% distinct(altitude, base_site_name, coord_x_2154, coord_y_2154, id_site, cd_nom, nom_cite, .keep_all = F)
# NB : on prend les coordonnées en lambert 93, pour pouvoir faire des buffers kilométriques

sites_FlodAlti_sp = vect(sites_FlodAlti, geom=c("coord_x_2154", "coord_y_2154"), crs="+init=epsg:2154")
sites_FlodAlti_buff2km = aggregate(buffer(sites_FlodAlti_sp, 2000))
# -----------------  /!\ décider la taille du buffer pour le projet iNat -----------------*

writeVector(sites_FlodAlti_buff2km, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/SPOT_FlodAlti_buff2km_l93.shp")



# plot(sites_FlodAlti_buff2km, ylim=c(6430000, 6450000))
# points(sites_FlodAlti_sp)

# # Aperçu de la localisation des points
# ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
# map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
#                     maptype = "terrain", scale = 2)
# ggmap(map_base) +
#   geom_polygon(data=geom(project(sites_FlodAlti_buff2km, "epsg:4326")), aes(x=x, y=y, group = part))+
#   geom_point(data=geom(project(sites_FlodAlti_sp, "epsg:4326")), aes(x=x, y=y), col="red", size=0.5)



# Plutôt que de se concentrer sur les points SPOT déjà identifiés, on travaille sur des zones où il y a animation possible :
# - Massif du Mont Blanc (délimitation existante sur iNaturalist)
# - Belledonne (délimitation en croisant les limites du massif tel que défini dans les modèles météo SAFRAN, et les lignes de niveau 1000m)
# - Mont Aiguille (pelouse sommitale - délimitation 'à la main')
# - Lautaret et alentours (zone LTSER de la ZAA élargie)
# - coeur du PN Ecrins (délimitation existante sur iNaturalist)

# Récupération de Belledonne dans les massifs SAFRAN :
massifs_SAFRAN = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/alp_allslopes/massifs_alpes_2154.shp")
Belledonne = massifs_SAFRAN[massifs_SAFRAN$nom == "Belledonne"]
Belledonne_WGS84 = project(Belledonne, "epsg:4326")
writeVector(Belledonne_WGS84, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/Belledonne_WGS84.kml")

# Récupération du Mt Aiguille d'après la délimitation manuelle sur mymaps :
MtAiguille = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/MtAiguille_mymap_WGS84.kml")

# Filtre altitudinal : on évite les plaines avec un filtre à 1000m d'altitude (pour Belledonne)
MNT_38 = rast("/Volumes/RESSOURCES_SIG/MNT/MNT_Alpes/MNT_FRANCE_25m_L93/Dpt_38_asc.asc")
MNT_Belledonne = mask(MNT_38, Belledonne)
MNT_Belledonne_1000m = ifel(MNT_Belledonne >= 1000,1,NA)
vect_altiBelledonne = as.polygons(MNT_Belledonne_1000m)
crs(vect_altiBelledonne) = "epsg:2154"
Belledonne_alti1000_WGS84 = project(vect_altiBelledonne, "epsg:4326")

writeVector(Belledonne_alti1000_WGS84, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/Belledonne_alti1000_WGS84.kml")

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
                    maptype = "terrain", scale = 2)
ggmap(map_base) +
  geom_polygon(data=geom(Belledonne_alti1000_WGS84), aes(x=x, y=y, group = part), col="red", fill=NA)+
  geom_polygon(data=geom(MtAiguille), aes(x=x, y=y, group = part), col="red", fill=NA)


# ----------------------------------------------------------------------------------------*

# Un projet iNaturalist a été créé pour regrouper toutes les observations des 7 espèces d'intérêt (myrtille, airelle bleue, dryade à huit pétales, marguerite des 
# Alpes, lis de St-Bruno, soldanelle, rhododendron ferrugineux) avec une information phénologique dans les zones où il peut y avoir une animation (Mont-Blanc,
# Ecrins, Belledonne, Mont Aiguille). Les données récupérées via ce projet peuvent être téléchargées.

projFA_iNat = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/data/FlodAlti/iNat_projetFlodAlti__export.csv")
paste("@", unique(projFA_iNat$user_login), sep="", collapse = " ")


# ----------------------------------------------------------------------------------------*
# ----------------------------------------------------------------------------------------*
# ----------------------------------------------------------------------------------------*


# PLAN D'ÉCHANTILLONNAGE STRUCTURÉ POUR LES SITES SUPPLEMENTAIRES DE SUIVI DE LA PRODUCTION DE FRUITS (développement n°2) ----

# On souhaite échantillonner différents environnements où la myrtille est présente, en supposant que la production de fruits dépend de
# ces conditions environnementales. On croise donc :
# - la topographie
#       => altitude, exposition, rayonnement incident, relief en creux ou en bosse, peuvent jouer sur T et lumière dispo pour le 
#           développement des fruits
# - la structure de la végétation 
#       => en sous-bois ou en lande, humidité, accès à la lumière, disponibilité de ressources varient et peuvent jouer sur le 
#           développement des fruits
#
# (à ces paramètres qui changent peu d'une année à l'autre s'ajoutent des effets annuels dont l'effet du climat sur la phénologie, la 
# pollinisation, qu'on étudie par ailleurs sur davantage de sites)


#* 0. Données utiles ----

# /!\ Pour les fichiers dispo sur le NAS (topo IGN + habitats ORION + lim CC), à voir pour les accès distants !!

# zone_etude = rbind(vect("/Volumes/RESSOURCES_SIG/Perim_administratif/Communes_EMB/CC_Pays_Mont_Blanc.shp"),
#                    vect("/Volumes/RESSOURCES_SIG/Perim_administratif/Communes_EMB/CC_Vallee_CHAMONIX.shp"))
zone_etude = rbind(vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_lim_admin/CC_Pays_Mont_Blanc.shp"),
                   vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_lim_admin/CC_Vallee_CHAMONIX.shp"))
zone_etude = project(zone_etude, "epsg:2154")


# Fonds carto
cle_ggmap ="" # à récupérer
ggmap::register_google(key=cle_ggmap)
# # Alpes du Nord
map_base_AlpN <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
                    maptype = "terrain", scale = 2)
# Mont-Blanc
map_base_MtBlc <- get_map(location = c(lon = 6.76, lat = 45.90), zoom = 11,
                    maptype = "terrain", scale = 2)


#*---- 0.1. Topo (IGN) ----

# MNT_MB = rast("/Volumes/RESSOURCES_SIG/ATLAS_MB/Topo/MB_MNT_25m.tif")
MNT_MB = rast("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_topo_IGN/MB_MNT_25m.tif")
MNT_MB = project(MNT_MB, "epsg:2154")
MNT_MB = crop(MNT_MB, zone_etude, mask=T)

# Calcul de variables topographiques
# PENTE
slope_MB = terrain(MNT_MB, v='slope', unit='degrees')
# EXPOSITION
aspect_MB = terrain(MNT_MB, v='aspect', unit='radians')
northness_MB = cos(aspect_MB)
eastness_MB = sin(aspect_MB)

rast_topo = c(MNT_MB, slope_MB, aspect_MB, northness_MB, eastness_MB)
names(rast_topo) = c("altitude","slope", "aspect", "northness","eastness")


#*---- 0.2. Type d'habitat / lande vs forêt (ORION) ----

# classes_ORION = read_xlsx("/Volumes/RESSOURCES_SIG/ATLAS_MB/Habitat/ORION/ORION_classes.xlsx")
# hab_ORION_CCPMB = rast("/Volumes/RESSOURCES_SIG/ATLAS_MB/Habitat/ORION/LC_pred_CCPMB_2806023_FINAL.tif")
# hab_ORION_CCVCMB = rast("/Volumes/RESSOURCES_SIG/ATLAS_MB/Habitat/ORION/LC_pred_CCVCMB_2806023_FINAL.tif")
# hab_ORION = merge(extend(hab_ORION_CCPMB, hab_ORION_CCVCMB),extend(hab_ORION_CCVCMB, hab_ORION_CCPMB))
# # PB : l'extent n'est pas le même, et pour concaténer les 2 rasters il faut du resampling, qui pose problème pour les classes, puisqu'on
# #      n'a plus les valeurs des classes mais des "moyennes"... 
# # => on utilise plutôt le fichier plus large, et on crop à la zone d'étude envisagée
# hab_ORION = rast("/Volumes/RESSOURCES_SIG/ATLAS_MB/Habitat/ORION/LC_pred_EMB_2806023_FINAL.tif")

# hab_ORION = rast("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_hab_flore/ORION/LC_pred_EMB_2806023_FINAL.tif")
# hab_ORION_l93 = crop(project(hab_ORION, "epsg:2154", method="near"), zone_etude, mask=T)
# writeRaster(hab_ORION_l93, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/habORION_MB_l93.tif")

classes_ORION = read_xlsx("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_hab_flore/ORION/ORION_classes.xlsx")
hab_ORION_l93 = rast("/Users/ninonfontaine/Desktop/projetsR/TEST/output/habORION_MB_l93.tif")

# On souhaite une résolution et un extent compatibles avec les rasters topo (25m)
hab_ORION = resample(hab_ORION_l93, MNT_MB, method="near") 
# /!\ la meilleure méthode serait "modal", pour avoir la classe majoritaire et pas la plus proche du point de resampling (comme avec
#     "near")... mais la fonction "modal" n'est pas implémentée dans resample... aggregate(fun="modal") avant resample ?

levels(hab_ORION) = classes_ORION # pour avoir une version en classe, pas numérique
names(hab_ORION) = "habitat"
# hist(MNT_MB[hab_ORION %in% 2:7])


#*---- 0.3. Présence de chemins pour l'accessibilité des points sélectionnés (BD TOPO IGN) ----

chemins = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_chemins_IGN/ITINERAIRE_AUTRE.shp")
routes = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_chemins_IGN/TRONCON_DE_ROUTE.shp")

rast_dist_acces = sum(rasterize(chemins, MNT_MB, background=NA), rasterize(routes, MNT_MB, background=NA), na.rm=T)
rast_dist_acces = distance(rast_dist_acces)
rast_dist_acces = mask(rast_dist_acces, MNT_MB)

rast_dist_acces1500 = mask(rast_dist_acces, ifel(rast_dist_acces > 1500, NA,1))

#------------------ En plus de ça, ça pourrait être pertinent de se concentrer sur les zones où le CREA a déjà concentré ses efforts (optimisation 
#                   déplacements & co) SSI ça couvre les différentes conditions habitat x topo.
#                   => Blaitière, Péclerey, Loriaz + Servoz


#*---- 0.4. Présences de myrtilles (CBNA + iNaturalist + relevés landes CREA) ----

# On ne regarde les données que sur la zone d'étude
zone = project(ext(zone_etude), from="epsg:2154",to="epsg:4326")

# # data_CBNA = read.csv("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/LANDE/Data/export_cbna_20231116/export_cbna_crea.csv")
data_CBNA = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_hab_flore/bota/export_cbna_crea.csv")
data_CBNA = data_CBNA[data_CBNA$cd_nom %in% c(128345, 128354),]

data_CBNA = data_CBNA[data_CBNA$lon_wgs84 <= zone$xmax & data_CBNA$lon_wgs84 >= zone$xmin &
                        data_CBNA$lat_wgs84 <= zone$ymax & data_CBNA$lat_wgs84 >= zone$ymin,]
data_CBNA = data_CBNA %>% rename("longitude" = lon_wgs84, "latitude" = lat_wgs84)

# 
# # Aperçu de la localisation des points
# cle_ggmap ="" # à récupérer
# ggmap::register_google(key=cle_ggmap)
# # # Alpes du Nord
# # map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
# #                     maptype = "terrain", scale = 2)
# # ggmap(map_base) +
# #   geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84))
# # Mont-Blanc
# map_base_MtBlc <- get_map(location = c(lon = 6.8, lat = 45.95), zoom = 11,
#                     maptype = "terrain", scale = 2)
# ggmap(map_base_MtBlc) +
#   geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84))
# 
# # Il y a trop peu de points pour se baser là-dessus pour sélectionner les points supplémentaires !
# #   => se baser plutôt sur des modèles de distribution ? (voir avec Isa)
# #   => retour Isa : les modèles de distribution myrtille ne marchent pas très bien, ce n'est pas optimal de se baser là-dessus


data_iNat = rbind(get_inat_obs(taxon_name = "Vaccinium myrtillus",
                         bounds = c(zone$ymin, zone$xmin, zone$ymax, zone$xmax)),
                  get_inat_obs(taxon_name = "Vaccinium uliginosum",
                               bounds = c(zone$ymin, zone$xmin, zone$ymax, zone$xmax)))
data_iNat[,c("x_l93","y_l93")] = crds(project(vect(data_iNat, geom=c("longitude","latitude"), "epsg:4326"),"epsg:2154"))


data_landesCREA = RJSONIO::fromJSON("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/05. DONNEES/Floraison_altitude_et_myrtille/Relevé_végétation/releve_cleaned_301224.json")
codes_esp = names(data_landesCREA$metadata_global$vegetation)
data_landesCREA = as.data.frame(do.call("rbind", lapply(data_landesCREA$area, function(x){c(unlist(x$metadata[c("date","ref_project","elevation","latitude","longitude")]), unlist(x$PROTOCOL_PPFG_SUMMARY))})))
colnames(data_landesCREA) = c("date","ref_project","elevation","latitude","longitude",
                              codes_esp)
data_landesCREA$site = rownames(data_landesCREA)
data_landesCREA[,3:24] = apply(data_landesCREA[,3:24], 2, as.numeric)
data_landesCREA[,c("x_l93","y_l93")] = crds(project(vect(data_landesCREA, geom=c("longitude","latitude"), "epsg:4326"),"epsg:2154"))



pts_myrtille = rbind(data_iNat[grep("myrtillus",data_iNat$scientific_name),c("x_l93","y_l93", "longitude", "latitude")],
                     data_CBNA[data_CBNA$cd_nom==128345, c("x_l93", "y_l93", "longitude", "latitude")],
                     data_landesCREA[data_landesCREA$LVM != 0, c("x_l93", "y_l93", "longitude", "latitude")])
pts_airelle = rbind(data_iNat[grep("uliginosum",data_iNat$scientific_name),c("x_l93","y_l93", "longitude", "latitude")],
                     data_CBNA[data_CBNA$cd_nom==128354, c("x_l93", "y_l93", "longitude", "latitude")],
                     data_landesCREA[data_landesCREA$LVU != 0, c("x_l93", "y_l93", "longitude", "latitude")])

#*-------- PISTE : modèle de distribution de myrtille (à utiliser pour filtrer les sites potentiels) ----
library(biomod2)

# FORMATTING DATA
myResp <- rep(1, nrow(pts_myrtille))
myRespXY <- pts_myrtille[, c('x_l93', 'y_l93')]
# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12) 
# data(bioclim_current)
# myExpl <- terra::rast(bioclim_current)
myExpl = c(rast_topo, hab_ORION_l93)
# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = as.vector(myResp),expl.var = myExpl,resp.xy = myRespXY,resp.name = "Myrtille",
                                     PA.nb.rep = 4, PA.strategy = 'disk', PA.dist.min = 50, PA.dist.max = 5000, PA.nb.absences = 400)#, filter.raster = TRUE)
# biomod2::plot(myBiomodData)
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    models = c('RF', 'GLM'),
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3,
                                    seed.val = 42)
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'PA+run',
                                      em.algo = c('EMmean', 'EMca'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.5),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      seed.val = 42)
# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = 'CurrentEM',
                                             new.env = myExpl,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
biomod2::plot(myBiomodEMProj)
# models.proj <- get_built_models(myBiomodModelOut, algo = "RF")
# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'CurrentRangeSize',
                                  new.env = myExpl,
                                  models.chosen = 'all',
                                  metric.binary = 'all')
biomod2::plot(myBiomodProj)


test = rast("/Users/ninonfontaine/Desktop/projetsR/TEST/Myrtille/proj_CurrentEM/proj_CurrentEM_Myrtille_ensemble.tif")
test1 = test[["Myrtille_EMmeanByTSS_PA1_RUN2_mergedAlgo"]]



#*---- 0.5. Points pré-envisagés, pour l'écoacoustique notamment (ORCHAMP + Colin + CamTrap ?) ----
pts_envisages = read_xlsx("/Users/ninonfontaine/Library/CloudStorage/GoogleDrive-nfontaine@creamontblanc.org/Drive partagés/SoPheno/Sites_spots.xlsx",
                          sheet = "preselection", col_types = c("text",rep("numeric",3), rep("text",10)))
pts_envisages = pts_envisages[!is.na(pts_envisages$Lat),]
pts_envisages$selec = paste0("selection initiale",ifelse(!is.na(pts_envisages$Suivi_acoustique)," - ecoacoustique",""))
pts_envisages = pts_envisages %>% rename("Name"="site")


# On inclut aussi les gradients ORCHAMP 'doubles', où il y a des camtrap et équipements déjà installés (et qui font un doublon de conditions similaires). 
# Sur ces points ORCHAMP, il y a des CamTrap : on les récupère par ce biais (les noms des CamTrap doublons ont une indication d'orientation).

# # ORCHAMP = st_read("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/CameraTrap/GPX:KML/camerainfo_running_20221019.kml")
# ORCHAMP = st_read("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_camtrap/camerainfo_running_20221019.kml")
# ORCHAMP[,c("Long","Lat")] = st_coordinates(ORCHAMP)
# ORCHAMP_df = as.data.frame(ORCHAMP)
# ORCHAMP_df$selec = "doublure ORCHAMP"

CamTrap = read_ods("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/CameraTrap/DATA/camerainfo_20240113.ods", sheet = "Tous")
CamTrap = CamTrap[!is.na(CamTrap$site) & !grepl("Campagnol", CamTrap$site) & CamTrap$running == "Y",]
CamTrap = CamTrap %>% rename("Long"="long","Lat"="lat")
CamTrap = CamTrap %>% rename("Name"="Station")
CamTrap$selec = "reseau CamTrap"
CamTrap$selec[CamTrap$site %in% c("Blaitiere","Loriaz","Peclerey")] = "ORCHAMP"
CamTrap$selec[CamTrap$selec == "ORCHAMP" & grepl("Ouest|Est|Sud|Nord", CamTrap$Name)] = "doublure ORCHAMP"
CamTrap$selec[CamTrap$site == "Para"] = "doublure ORCHAMP"

# ggmap(map_base_MtBlc) +
#   # geom_point(data=pts_myrtille, aes(x=lon_wgs84, y=lat_wgs84), col="purple")+
#   geom_point(data=pts_envisages, aes(x=Long, y=Lat))+
#   geom_point(data=ORCHAMP, aes(x=Long, y=Lat), col="orange")+
#   geom_point(data=CamTrap, aes(x=long, y=lat), col="red")


# /!\ on n'a pas les infos des identifiants des points ORCHAMP dans le kml... on a les noms dans le fichier des CamTrap, mais pas 
#     tous les points !
# /!\ il y a des doublons dans les points entre 'pts_envisages' et 'ORCHAMP' ! (ça ne devrait plus apparaître une fois rasterisé)

pts_pot = rbind(pts_envisages[!is.na(pts_envisages$Lat), c("Name","Long","Lat","selec")],
                # ORCHAMP_df[!is.na(ORCHAMP_df$Lat), c("Name","Long","Lat","selec")],
                CamTrap[!is.na(CamTrap$Lat), c("Name","Long","Lat","selec")])
pts_pot = pts_pot %>% distinct(Long, Lat, .keep_all = T)
pts_pot = project(vect(pts_pot, geom=c("Long","Lat"), crs="epsg:4326"), "epsg:2154")              


rast_pts_pot = rasterize(pts_pot, MNT_MB, background=0)
names(rast_pts_pot) = "potential_pts"



#* 1. Diversité des conditions dans le Mont-Blanc ----

#*---- 1.1. Typologie 'basique' en créant des catégories exposition x habitat ----

#---- Aperçu des conditions au niveau des points envisagés :
pts_pot_cond = cbind(crds(pts_pot),as.data.frame(pts_pot), extract(c(rast_topo, hab_ORION), pts_pot)[,-1])

par(mfrow=c(2,3))
for (i in 4:ncol(pts_pot_cond)){
  if(is.numeric(pts_pot_cond[,i])){par(mar=c(2,2,4,0))
    hist(pts_pot_cond[,i], main=colnames(pts_pot_cond)[i], xlab="")}
  else if(is.factor(pts_pot_cond[,i])){par(mar=c(2,9,4,0))
    barplot(table(pts_pot_cond[,i]), main=colnames(pts_pot_cond)[i], las=2, horiz=T)}
}

#---- Typologie 'basique' en créant des catégories exposition x habitat

# Simplification des habitats en 3 classes : 
# - foret pour 2 Forêt, 3 Limite de la forêt
# - lande pour 5 Lande, 6 Ecotone lande-prairie
# - prairie pour 4 Prairie subalpine, 7 Prairie alpine
hab_ORION_class = classify(hab_ORION, cbind(0:10, 
                                            c(NA,NA,1,1,3,2,2,3,NA,NA,NA)))
#c(NA, NA,"foret","foret","prairie","lande","lande","prairie",NA, NA, NA)))

# Simplification des expositions en 4 classes : 
# - N=1 si aspect entre 0 - pi/4 rad, ou 7pi/4 - 2pi rad
# - E=2 si aspect entre pi/4 - 3*pi/4 rad
# - S=3 si aspect entre 3*pi/4 - 5*pi/4 rad
# - W=4 si aspect entre 5*pi/4 - 7*pi/4 rad
asp_class = classify(aspect_MB, matrix(c(0, pi/4, 1,
                                         pi/4, 3*pi/4, 2,
                                         3*pi/4, 5*pi/4, 3,
                                         5*pi/4, 7*pi/4, 4,
                                         7*pi/4, 2*pi, 1), ncol=3, byrow=TRUE), include.lowest=T)

# => on obtient ainsi 12 classes habitat x exposition
typo = as.factor(hab_ORION_class*10 + asp_class)
names(typo) = "typo_Hab_x_Asp"
plot(typo)
table(values(typo))
# # Filtre par l'accessibilité ?
# table(values(mask(typo, ifel(rast_dist_acces > 1500, NA,1))))


#*---- 1.2. Clustering des types de milieux ----

# # V1
# # Les habitats sur lesquels on se focalise sont les habitats 2 à 7 
# # (2 Forêt, 3 Limite de la forêt, 4 Prairie subalpine, 5 Lande, 6 Ecotone lande-prairie, 7 Prairie alpine)
# # V2 
# # Au vu des points à myrtille identifiés, les habitats sur lesquels on se focalise sont les habitats 2, 3, 5 et 6
# # (2 Forêt, 3 Limite de la forêt, 5 Lande, 6 Ecotone lande-prairie)
# 
# # par(mfrow=c(3,2), mar=c(4,2,4,1))
# # for (i in 2:7){
# #   hist(MNT_MB[hab_ORION == i], main=classes_ORION$NOM[classes_ORION$CODE ==i], xlab="Altitude", xlim=c(minmax(MNT_MB)))
# # }
# 
# stackPCA = c(rast_topo, hab_ORION, typo, rast_pts_pot)
# # stackPCA = mask(stackPCA, hab_ORION, maskvalues=2:7, inverse=T) # on ne garde que les habitats d'intérêt
# stackPCA = mask(stackPCA, hab_ORION, maskvalues=c(2,3,5,6), inverse=T) # on ne garde que les habitats d'intérêt
# stackPCA = mask(stackPCA, app(stackPCA, fun = sum)) # on ne garde que les pixels où on a toutes les infos
# stackPCA[[5]] = droplevels(stackPCA[[5]])
# stackPCA[[6]] = droplevels(stackPCA[[6]])
# 
# # plot(stackPCA)
# 
# 
# # Différentes méthodes testées pour le clustering : 
# # 1) ACP 'classique' MAIS on ne peut pas mettre de variables catégorielles comme les habitats !
# 
# # # PCA on topography (only on focus habitats)
# # PCA = prcomp(stackPCA[[1:4]], scale. = T)
# # 
# # # Eigenvalues graph
# # par(mfrow=c(1,1))
# # barplot(summary(PCA)[[6]][2,], main="Proportion of variance")
# # # Arrow graph
# # plot(c(-2,2), c(-2,2), type="n")
# # arrows(x0=0, y0=0, x1=PCA$rotation[,1], y1=PCA$rotation[,2])
# # text(x=PCA$rotation[,1]+0.2, y=PCA$rotation[,2], rownames(PCA$rotation))
# # # points(PCA$x[na.omit(values(stackPCA[[5]]))==4,1:2], cex=0.05)
# # # points(PCA$x[na.omit(values(stackPCA[[5]]))==4,1:2], cex=0.05)
# # points(PCA$x[na.omit(values(stackPCA[[6]]))==1,1:2], cex=0.5, col="red", pch=20) # points envisagés
# # 
# # # Il y a un trou / il manquerait des points à faibles pentes et altitudes => à voir là où ce sont des zones favorables aux myrtilles
# # 
# # # Problème : on n'a pas intégré les différents habitats, puisque l'ACP ne prend pas de variables catégorielles...
# # #       => choix de méthodes mixtes ?
# 
# 
# # 2) Analyse multivariée mixte + clustering
# 
# # # Analyse multivariée ade4
# # AnMult = dudi.mix(na.omit(values(stackPCA[[1:5]])), nf=4)
# # 
# # # Clustering kmeans du package terra
# # # Normalisation des données (sinon l'altitude a beaucoup plus de poids)
# # nx <- minmax(stackPCA[[1:4]])    
# # rn <- (stackPCA[[1:4]] - nx[1,]) / (nx[2,] - nx[1,])
# # stackNorm = c(rn, stackPCA[[5:6]])
# # # Calcul des kmeans /!\ en théorie ça ne marche pas avec des variables catégorielles !! /!\
# # kmeans = k_means(na.omit(values(stackNorm[[1:5]])), centers=6)
# # clusters = stackNorm[[1]]
# # values(clusters)[!is.na(values(clusters))] = kmeans$cluster
# # plot(clusters)
# # 
# # # Matrice de distance puis hclust
# # matdist = daisy(na.omit(values(stackPCA[[1:5]])), metric = "gower")
# # clust = hclust(matdist)
# # # Trop de pixels = trop lourd !
# # 
# # 
# # # Analyses multivariée de données mixtes (FAMD)
# # tab_FAMD = as.data.frame(na.omit(values(stackPCA[[1:5]])))
# # tab_FAMD$habitat = as.character(tab_FAMD$habitat)
# # famd = FAMD(tab_FAMD, graph=F, ncp=9)
# # 
# # # fviz_famd_var(famd)
# # fviz_eig(famd) # 4 ou 5 dimensions
# # # fviz_famd_ind(famd, habillage=5)
# # 
# # clustering = HCPC(famd, graph=F) # TROP LOURD !!!
# 
# 
# # Utilisation de hclust avec la distance de Gower (fonctionne avec les variables catégorielles)
# h_clust <- function(x, database=as.data.frame(na.omit(values(x))), 
#                     ngroups,  clust_method="complete", #dist_metric="euclidean",
#                     types=list(numeric=1:dim(x)[3]), agfuns=c(mean,mfv1), 
#                     matchfun=gower.dist, ..., #matchfun="squared"
#                     maxcell=10000, filename="", overwrite=FALSE, wopt=list()) {
#   stopifnot(maxcell > 0)
#   stopifnot(ngroups > 0)
#   stopifnot(ngroups < maxcell)
#   d <- na.omit(spatSample(x, maxcell, "regular"))
#   # print(dim(d))
#   # dd <- stats::dist(d, dist_metric) # adaptation par rapport au script initial
#   # hc <- stats::hclust(dd, clust_method) # adaptation par rapport au script initial
#   dd <- StatMatch::gower.dist(d) # adaptation par rapport au script initial
#   hc <- stats::hclust(as.dist(dd), clust_method) # adaptation par rapport au script initial
#   th <- sort(hc$height, TRUE)[ngroups]
#   cls <- stats::cutree(hc, h = th)
#   hc <- cut(stats::as.dendrogram(hc), h=th)$upper
#   # d <- aggregate(d, list(cls=cls), agfun)
#   d <- merge(as.data.frame(aggregate(d[,types[["numeric"]]], list(cls=cls), agfuns[[1]])),
#              as.data.frame(aggregate(d[,types[["factor"]]], list(cls=cls), agfuns[[2]])), by="cls")# adaptation par rapport au script initial
#   # print(d)
#   cls <- d$cls 
#   d$cls <- NULL
#   colnames(d)=names(x)
#   # b <- bestMatch(x, d, fun=matchfun, ..., filename=filename, overwrite=overwrite, wopt=wopt)	
#   # print(head(d))
#   bdist = gower.dist(database,d)
#   b=apply(bdist,1,FUN = function(X){c(1:ncol(bdist))[X==min(X)]})
#   brast=x[[1]]
#   brast[!is.na(values(brast))]=b
#   # b <- bestMatch(x, d, fun=matchfun, ..., filename=filename, overwrite=overwrite, wopt=wopt)
#   return(list(barycentres=d,  dendrogram=hc, clusterrast=brast))
# } # https://rdrr.io/github/rspatial/terra/src/R/k_means.R 
# 
# database = as.data.frame(na.omit(values(stackPCA)))
# database$habitat = droplevels(factor(database$habitat, levels=classes_ORION$CODE, labels=classes_ORION$NOM))
# 
# clustering = h_clust(stackPCA[[c(1,2,4,5,6)]], database = database[c(1,2,4,5,6)],
#                      ngroups=8, types=list(numeric=1:4, factor=5), maxcell=50000) 
# # clustering = h_clust(stackPCA[[1:5]], database = database,
# #                      ngroups=6, types=list(numeric=1:4, factor=5), maxcell=50000) 
# 
# clust_rast = clustering$clusterrast
# terra::plot(clust_rast)
# 
# database$clust8 = na.omit(values(clust_rast))
# 
# # Au final le clustering ressemble aux différents habitats, MAIS 
# # - l'habitat "forêt" est subdivisé en limite de la forêt (cluster 5)  VS   forêt en fonction de l'exposition (4 expositions - clusters 1 6 7 8)
# # - l'habitat "lande" est subdivisé en écotone lande-prairie (un peu plus haut en altitude - cluster 2)  VS  lande NW (cluster 3)   VS   lande SE (4)



#---------------- => ce système de clustering ne semble donc pas apporter grand-chose relativement au système en classe habitat x exposition, 
#                    qui a le mérite d'être plus lisible => on abandonne la démarche clustering.



#* 2. Représentativité des points présélectionnés ----
# Points d'écoacoustique ORCHAMP + 4 sites + doublons ORCHAMP

pts_pot_cond$typo = extract(typo, pts_pot_cond[,c("x","y")])[,-1]
table(pts_pot_cond$selec, pts_pot_cond$typo)

plot(typo, alpha=0.8)
points(pts_pot_cond$x, pts_pot_cond$y, pch=21, 
       cex=ifelse(pts_pot_cond$selec=="selection initiale - ecoacoustique",1.2, 0.8),
       bg=as.character(factor(pts_pot_cond$selec,
                              levels=c("selection initiale - ecoacoustique","selection initiale","doublure ORCHAMP","reseau CamTrap"),
                              labels=c("darkgreen","yellowgreen","lightgray","darkgrey"))))


par(mfrow=c(1,1), mar=c(2,2,2,2))
barplot(table(pts_pot_cond$selec, pts_pot_cond$typo), las=2, horiz=T, col=c("lightgray","darkgrey","yellowgreen","darkgreen"))
# Dans la sélection initiale, il n'y a pas du tout de point dans l'habitat "prairie" (typo 3X)
# --> est-ce qu'on n'abandonnerait pas cette catégorie ?
# La plupart des points présélectionnés pour l'écoacoustique (vert foncé) sont en forêt (n=7) ou en lande (n=3), majoritairement en exposition W


# Diversité au niveau des points à myrtille du CBNA + iNat
pts_myrtille$typo = extract(typo, vect(pts_myrtille, geom=c("x_l93","y_l93"), crs="epsg:2154"))[,-1]
par(mfrow=c(1,1), mar=c(2,2,2,2))
barplot(table(pts_myrtille$typo), las=2, horiz=T)
# Dans les points à myrtille, il y a 1 seul point dans l'habitat "prairie" (typo 3X)... En plus dans le relevé CBNA, l'habitat indiqué est 
# "Landine rase à Loiseleuria procumbens située dans une concavité acide temporairement humide" !
# --> ça va dans le sens de l'abandon de la catégorie "prairie"...


#---------------- => Au vu des données qu'on a, on se limite à 8 des catégories initialement envisagées : on exclut les catégories "prairie" et on se
#                    concentre sur les landes et forêts (en sachant qu'il y aura peut-être malgré tout des milieux intermédiaires lande/prairie)

typo_sel = mask(typo, typo , maskvalues=c(31,32,33,34))
colors = data.frame(ID = levels(typo_sel)[[1]]$habitat[1:8], col = palette("Paired")[1:8])
# on remplace juste la dernière couleur, qui bugue
colors[8,2] = "darkorange"
plot(typo_sel, alpha=0.8, col=colors$col)
points(pts_pot_cond$x, pts_pot_cond$y, pch=21, 
       cex=ifelse(pts_pot_cond$selec=="selection initiale - ecoacoustique",1.2, 0.8),
       bg=as.character(factor(pts_pot_cond$selec,
                              levels=c("selection initiale - ecoacoustique","selection initiale","doublure ORCHAMP","reseau CamTrap"),
                              labels=c("darkgreen","yellowgreen","lightgray","darkgrey"))))
points(pts_myrtille$x_l93, pts_myrtille$y_l93, cex=0.5, pch=16,col="black")


typo_sel_leaflet = project(typo_sel, "epsg:4326", method="near")
pts_pot_cond[,c("longitude","latitude")] = crds(project(vect(pts_pot_cond, geom=c("x","y"), "epsg:2154"),"epsg:4326"))



carto = leaflet() %>% #addCircleMarkers(color = ~c("white","purple")[points_couchades_potentiels_proj$pasto_couchades +1]) %>% 
  addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>% 
  addRasterImage(typo_sel_leaflet, colors="Paired", opacity=0.7) %>%
  addCircleMarkers(lng=pts_myrtille$longitude, lat=pts_myrtille$latitude, col="black", radius=1, opacity=1)%>%
  addCircleMarkers(lng=pts_pot_cond$longitude, lat=pts_pot_cond$latitude, 
                   fillColor =as.character(factor(pts_pot_cond$selec,
                              levels=c("selection initiale - ecoacoustique","selection initiale","doublure ORCHAMP","reseau CamTrap"),
                              labels=c("darkgreen","yellowgreen","lightgray","darkgrey"))), 
                   radius=ifelse(pts_pot_cond$selec=="selection initiale - ecoacoustique",6,4), 
                   color="black",opacity=1,fillOpacity = 1, weight=0.5)
carto




#* 3. Échantillonnage aléatoire de points supplémentaires où faire du suivi de production de fruits ----

# # Pour optimiser la logistique terrain, on se limite à quelques zones où il y a de la myrtille d'après les relevés CBNA ou les pointages iNaturalist
# # ou d'autres sources d'information (ex. CamTrap, AltitudeRando !).
# zone_myrtille = aggregate(buffer(vect(pts_myrtille, geom=c("x_l93","y_l93"), "epsg:2154"), 500))
# rast_myrtille = rasterize(zone_myrtille, typo)
# # plot(mask(typo, rast_myrtille), col=palette("Paired"))
# # /!\ BIFBOF, très peu de points quand même ! Se baser plutôt sur des zones préidentifiées ?

# L'idée est de parcourir des gradients altitudinaux pour faire à la fois des sites en forêt et en lande.
# Les zones pré-identifiées d'après la présence avérée de myrtille ET pour prendre en compte les différentes situations habitat x exposition sont :
# - le fond du vallon des Contamines, en exposition W sous le refuge de Tré la Tête
# - le fond du vallon des Contamines, en exposition E sous le refuge des Prés
# - du côté des réserves de Passy et des Aiguilles rouges, en montant au refuge de Moede par Plaine Joux puis en continuant par le chemin sous le
#   lac Cornu (i.e. au nord de Pormenaz) --> expo SW 
# - du col des Montets au Lac Blanc --> expo E
# - à Loriaz, avec le double gradient ORCHAMP --> expo SE
# - à Péclerey, avec le double gradient ORCHAMP --> expo W
# - vers le Mont Lachat (ou vers le Prarion) --> expo N
# - du côté de Blaitière (déjà des équipements CamTrap + ORCHAMP avec les infos associées, tous les points sauf le plus haut où il n'y a pas de 
#   myrtille) --> expo W
#   
# - À DISCUTER en plus, si on veut aussi avoir des sites du côté de l'autre Comcom : on peut viser la zone au-dessus de Combloux, vers le Planay ?

zones_preid = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/FlodAlti/zonesMyrtillePreSelec_mymaps.kml")

# ggmap(map_base_MtBlc) +
#   geom_polygon(data=geom(zones_preid), aes(x=x, y=y, group=geom)) + 
#   geom_point(data=pts_myrtille, aes(x=longitude, y=latitude), col="darkviolet")

zones_preid = project(zones_preid, "epsg:2154")
rast_preid = rasterize(zones_preid, typo)
# plot(typo_sel, col=palette("Paired"))
# lines(zones_preid, col="black")


# On se concentre aussi sur les zones pas trop éloignées des chemins (200 m max ?).
typo_preid = mask(typo_sel, mask(rast_preid, ifel(rast_dist_acces<200,1,NA)))
typo_preid = droplevels(typo_preid)
plot(typo_preid, col=colors$col)

# Dans ces zones, on pioche aléatoirement 15 points par catégorie. Dans un second temps, on sélectionnera parmi ces 15 points ceux (6) qui sont les 
# plus pratiques niveau logistique, et ça permet d'avoir des back-ups s'il n'y a pas de myrtille.
# Cette resélection peut aussi permettre de prendre en compte les CamTrap et points ORCHAMP préexistants --> ex. si un point a été sélectionné et se 
# trouve à proximité immédiate (moins de 200m (?)) d'une CamTrap ou d'un point ORCHAMP, on le remplace par le site CamTrap s'il y a de la myrtille.

# ech = spatSample(typo_preid, size=15, method="stratified", xy=T)
# # ça ne marche pas !! Peut-être parce que trop de NA ? Alternative : faire par catégorie
ech=data.frame("x"=NULL,"y"=NULL, "typo_Hab_x_Asp"=NULL)
for (type in levels(typo_preid)[[1]]$ID){
  # print(type)
  ech = rbind(ech,spatSample(mask(typo_preid, ifel(typo_preid==type,type,NA)), 
                             size=15, method="random", xy=T, exhaustive=T, na.rm=T))
}
dist = distance(vect(ech, geom=c("x","y"), "epsg:2154"), 
                vect(pts_pot_cond, geom=c("x","y"), "epsg:2154"))
distproxim = 200
ech$proxim = unlist(apply(dist, 1, function(x){prox = pts_pot_cond[x<distproxim,c("Name","typo_Hab_x_Asp")]
return(ifelse(nrow(prox)==0,NA,
              paste0(prox[,1]," (type ", prox[,2],") à ",round(x[x<distproxim],0)," m")))}))

dist_airelle = distance(vect(ech, geom=c("x","y"), "epsg:2154"), 
                        vect(pts_airelle, geom=c("x_l93","y_l93"), "epsg:2154"))
ech$airelle_proxim = apply(dist_airelle,1, min)

ech$zone = extract(zones_preid, ech[,c("x","y")])$Name
ech[,c("longitude","latitude")] = crds(project(vect(ech, geom=c("x","y"), "epsg:2154"),"epsg:4326"))
ech$altitude = extract(MNT_MB, ech[,c("x","y")])$MB_MNT_25m
ech$pente = extract(slope_MB, ech[,c("x","y")])$slope


write.csv(ech, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/PlEchantillonnage_FruitsMyrt.csv" ,row.names = F)  



# Visualisation de la sélection

ech = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/PlEchantillonnage_FruitsMyrt.csv")  

carto_ech = leaflet() %>%
  addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>% 
  addRasterImage(typo_sel_leaflet, colors="Paired", opacity=0.7) %>%
  addCircleMarkers(lng=pts_myrtille$longitude, lat=pts_myrtille$latitude, col="darkviolet", radius=0.5, opacity=1)%>%
  addCircleMarkers(lng=pts_airelle$longitude, lat=pts_airelle$latitude, col="darkblue", radius=0.5, opacity=1)%>%
  addCircleMarkers(lng=pts_pot_cond$longitude, lat=pts_pot_cond$latitude, 
                   popup = paste(pts_pot_cond$Name,"- Typo :", pts_pot_cond$typo),
                   fillColor =as.character(factor(pts_pot_cond$selec,
                                                  levels=c("selection initiale - ecoacoustique","selection initiale","ORCHAMP","doublure ORCHAMP","reseau CamTrap"),
                                                  labels=c("darkgreen","yellowgreen","yellowgreen","darkgrey","lightgray"))), 
                   radius=ifelse(pts_pot_cond$selec=="selection initiale - ecoacoustique",8,4), 
                   color="black",opacity=1,fillOpacity = 1, weight=1)%>%
  addCircleMarkers(lng=ech$longitude, lat=ech$latitude, 
                   popup = paste(ech$zone,"- Alt :",ech$altitude,", Typo :", ech$typo_Hab_x_Asp),
                   fillColor = as.character(factor(ech$typo_Hab_x_Asp,
                                                   levels=as.numeric(colors$ID),
                                                   labels=colors$col)),
                   radius=4, color="black",opacity=1,fillOpacity = 1, weight=3) %>%
  addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Classes hab. x expo.")
carto_ech
  
  

# Concaténation sélection + points acoustique + ORCHAMP

ech$selec = "selection aleatoire"
ech = ech  %>%
  group_by(typo_Hab_x_Asp) %>%
  mutate(Name = paste(typo_Hab_x_Asp, 1:length(typo_Hab_x_Asp), sep="_") )
ech$typo_Hab_x_Asp = as.factor(ech$typo_Hab_x_Asp)

pts_pot_cond$zone = extract(zones_preid, pts_pot_cond[,c("x","y")])$Name
pts_pot_cond$proxim = NA

dist_airelle = distance(vect(pts_pot_cond, geom=c("x","y"), "epsg:2154"), 
                        vect(pts_airelle, geom=c("x_l93","y_l93"), "epsg:2154"))
pts_pot_cond$airelle_proxim = apply(dist_airelle,1, min)

pts_pot_cond = pts_pot_cond %>% rename("pente" = slope, "typo_Hab_x_Asp" = typo)

recap_sites = rbind(ech,
                    pts_pot_cond[,colnames(ech)])


write.csv(recap_sites, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/PlEchantillonnage_FruitsMyrt_RECAPtotal.csv" ,row.names = F) 

rownames(recap_sites) = recap_sites$Name
recap_sites_vect = vect(recap_sites[,c("Name","altitude","x","y")], geom=c("x","y"), "epsg:2154")
writeVector(recap_sites_vect, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/PlEchantillonnage_FruitsMyrt.kml", overwrite=T)
# writeVector(recap_sites_vect, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/PlEchantillonnage_FruitsMyrt.gpx", filetype="GPX", overwrite=T)

table(recap_sites$zone, recap_sites$selec, recap_sites$typo_Hab_x_Asp)

# carto_acoustique = leaflet(recap_sites[recap_sites$selec == "selection initiale - ecoacoustique",]) %>%
#   addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>% 
#   addCircleMarkers(lng=recap_sites[recap_sites$selec == "selection initiale - ecoacoustique","longitude"], lat=recap_sites[recap_sites$selec == "selection initiale - ecoacoustique","latitude"], 
#                    fillColor ="darkgreen", 
#                    radius=8, 
#                    color="black",opacity=1,fillOpacity = 1, weight=1)
# carto_acoustique



# TEST avec modele distribution myrtille
recap_sites = read_xlsx("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/05. DONNEES/Floraison_altitude_et_myrtille/Sites_spots.xlsx",sheet="recap_toutesselections") 
recap_sites$modTEST = terra::extract(test1, recap_sites[,c("x","y")])[,2]

write.csv(recap_sites, "/Users/ninonfontaine/Google Drive/Mon Drive/PlEchantillonnage_FruitsMyrt_RECAPtotalMOD.csv", row.names=F)


