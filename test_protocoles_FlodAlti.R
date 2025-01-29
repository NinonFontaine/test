library(ggplot2)
library(tidyr)
library(terra)
library(sf)
library(ggmap)
library(lubridate)
library(dplyr)
library(readxl)

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

writeVector(sites_FlodAlti_buff2km, "/Users/ninonfontaine/Desktop/projetsR/FloraisonAltitude/SPOT_FlodAlti_buff2km_l93.shp")



# plot(sites_FlodAlti_buff2km, ylim=c(6430000, 6450000))
# points(sites_FlodAlti_sp)

# # Aperçu de la localisation des points
# ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
# map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
#                     maptype = "terrain", scale = 2)
# ggmap(map_base) +
#   geom_polygon(data=geom(project(sites_FlodAlti_buff2km, "epsg:4326")), aes(x=x, y=y, group = part))+
#   geom_point(data=geom(project(sites_FlodAlti_sp, "epsg:4326")), aes(x=x, y=y), col="red", size=0.5)
# 


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


#*---- 0. Données utiles ----
#*-------- 0.1. Topo (IGN) ----



#*-------- 0.2. Type d'habitat / lande vs forêt (ORION) ----



#*-------- 0.3. Présences de myrtilles (CBNA) ----
data_CBNA = read.csv("/Users/ninonfontaine/Desktop/projetsR/FloraisonAltitude/data/export_cbna_crea.csv")
data_CBNA = data_CBNA[data_CBNA$cd_nom == 128345,]

# Aperçu de la localisation des points
ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
# # Alpes du Nord
# map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
#                     maptype = "terrain", scale = 2)
# ggmap(map_base) +
#   geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84))
# Mont-Blanc
map_base_MtBlc <- get_map(location = c(lon = 6.8, lat = 45.95), zoom = 11,
                    maptype = "terrain", scale = 2)
ggmap(map_base_MtBlc) +
  geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84))

# Il y a trop peu de points pour se baser là-dessus pour sélectionner les points supplémentaires !
#   => se baser plutôt sur des modèles de distribution ? (voir avec Isa)


#*-------- 0.4. Points pré-envisagés, pour l'écoacoustique notamment (ORCHAMP + Colin) ----
pts_envisages = read_xlsx("/Users/ninonfontaine/Library/CloudStorage/GoogleDrive-nfontaine@creamontblanc.org/Drive partagés/SoPheno/Sites_spots.xlsx",
                          sheet = "Feuil1", col_types = c("text",rep("numeric",3), rep("text",9)))
ggmap(map_base_MtBlc) +
  geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84), col="purple")+
  geom_point(data=pts_envisages, aes(x=Long, y=Lat))


#*---- 1. Diversité des conditions dans le Mont-Blanc où se rencontre la myrtille ----



#*---- 2. Représentativité des points d'écoacoustique présélectionnés ----
# Points d'écoacoustique = ORCHAMP + 4 sites



#*---- 3. Échantillonnage aléatoire de points supplémentaires où faire du suivi de production de fruits ----




