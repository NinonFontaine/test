library(ggplot2)
library(tidyr)
library(statip)
library(terra)
library(sf)
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
map_base_MtBlc <- get_map(location = c(lon = 6.8, lat = 45.95), zoom = 11,
                    maptype = "terrain", scale = 2)


#*-------- 0.1. Topo (IGN) ----

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

rast_topo = c(MNT_MB, slope_MB, northness_MB, eastness_MB)
names(rast_topo) = c("altitude","slope", "northness","eastness")


#*-------- 0.2. Type d'habitat / lande vs forêt (ORION) ----

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




#*-------- 0.3. Présences de myrtilles (CBNA) ----

# # data_CBNA = read.csv("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/LANDE/Data/export_cbna_20231116/export_cbna_crea.csv")
# data_CBNA = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_hab_flore/bota/export_cbna_crea.csv")
# data_CBNA = data_CBNA[data_CBNA$cd_nom == 128345,]
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




#*-------- 0.4. Points pré-envisagés, pour l'écoacoustique notamment (ORCHAMP + Colin + CamTrap ?) ----
pts_envisages = read_xlsx("/Users/ninonfontaine/Library/CloudStorage/GoogleDrive-nfontaine@creamontblanc.org/Drive partagés/SoPheno/Sites_spots.xlsx",
                          sheet = "Feuil1", col_types = c("text",rep("numeric",3), rep("text",10)))

# On inclut aussi les gradients ORCHAMP 'doubles', où il y a des camtrap et équipements déjà installés (et qui font un doublon de 
# conditions similaires)
# ORCHAMP = st_read("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/CameraTrap/GPX:KML/camerainfo_running_20221019.kml")
ORCHAMP = st_read("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_camtrap/camerainfo_running_20221019.kml")
ORCHAMP[,c("Long","Lat")] = st_coordinates(ORCHAMP)
ORCHAMP_df = as.data.frame(ORCHAMP)

CamTrap = read_ods("/Volumes/EN_COURS/2.Recherche/3.Mont_Blanc/CameraTrap/DATA/camerainfo_20240113.ods", sheet = "Tous")
CamTrap = CamTrap[!is.na(CamTrap$site) & !grepl("Campagnol", CamTrap$site) & CamTrap$running == "Y",]
CamTrap = CamTrap %>% rename("Long"="long","Lat"="lat")


# ggmap(map_base_MtBlc) +
#   # geom_point(data=data_CBNA, aes(x=lon_wgs84, y=lat_wgs84), col="purple")+
#   geom_point(data=pts_envisages, aes(x=Long, y=Lat))+
#   geom_point(data=ORCHAMP, aes(x=Long, y=Lat), col="orange")+
#   geom_point(data=CamTrap, aes(x=long, y=lat), col="red")


# /!\ on n'a pas les infos des identifiants des points ORCHAMP dans le kml... on a les noms dans le fichier des CamTrap, mais pas 
#     tous les points !
# /!\ il y a des doublons dans les points entre 'pts_envisages' et 'ORCHAMP' ! (ça ne devrait plus apparaître une fois rasterisé)

pts_pot = rbind(pts_envisages[!is.na(pts_envisages$Lat), c("Long","Lat")],
                CamTrap[!is.na(CamTrap$Lat), c("Long","Lat")],
                ORCHAMP_df[!is.na(ORCHAMP_df$Lat), c("Long","Lat")])
pts_pot = pts_pot %>% distinct()
pts_pot = project(vect(pts_pot, geom=c("Long","Lat"), crs="epsg:4326"), "epsg:2154")              


rast_pts_pot = rasterize(pts_pot, MNT_MB, background=0)
names(rast_pts_pot) = "potential_pts"



#*---- 1. Diversité des conditions dans le Mont-Blanc ----

#*-------- 1.1. ... dans les habitats type forêt, lande ou prairie d'altitude ----

#---- Aperçu des conditions au niveau des points envisagés :
pts_pot_cond = cbind(crds(pts_pot),extract(c(rast_topo, hab_ORION), pts_pot)[,-1])

par(mfrow=c(2,3))
for (i in 3:ncol(pts_pot_cond)){
  if(is.numeric(pts_pot_cond[,i])){par(mar=c(2,2,4,0))
    hist(pts_pot_cond[,i], main=colnames(pts_pot_cond)[i], xlab="")}
  else if(is.factor(pts_pot_cond[,i])){par(mar=c(2,9,4,0))
    barplot(table(pts_pot_cond[,i]), main=colnames(pts_pot_cond)[i], las=2, horiz=T)}
}

#---- Clustering des types de milieux

# Les habitats sur lesquels on se focalise sont les habitats 2 à 7 
# (2 Prairie montagnarde, 3 Forêt, 4 Limite de la forêt, 5 Prairie subalpine, 6 Lande, 7 Ecotone lande-prairie)

# par(mfrow=c(3,2), mar=c(4,2,4,1))
# for (i in 2:7){
#   hist(MNT_MB[hab_ORION == i], main=classes_ORION$NOM[classes_ORION$CODE ==i], xlab="Altitude", xlim=c(minmax(MNT_MB)))
# }

stackPCA = c(rast_topo, hab_ORION, rast_pts_pot)
stackPCA = mask(stackPCA, hab_ORION, maskvalues=2:7, inverse=T) # on ne garde que les habitats d'intérêt
stackPCA = mask(stackPCA, app(stackPCA, fun = sum)) # on ne garde que les pixels où on a toutes les infos
stackPCA[[5]] = droplevels(stackPCA[[5]])

# plot(stackPCA)


# Différentes méthodes testées pour le clustering : 
# 1) ACP 'classique' MAIS on ne peut pas mettre de variables catégorielles comme les habitats !

# PCA on topography (only on focus habitats)
PCA = prcomp(stackPCA[[1:4]], scale. = T)

# Eigenvalues graph
par(mfrow=c(1,1))
barplot(summary(PCA)[[6]][2,], main="Proportion of variance")
# Arrow graph
plot(c(-2,2), c(-2,2), type="n")
arrows(x0=0, y0=0, x1=PCA$rotation[,1], y1=PCA$rotation[,2])
text(x=PCA$rotation[,1]+0.2, y=PCA$rotation[,2], rownames(PCA$rotation))
# points(PCA$x[na.omit(values(stackPCA[[5]]))==4,1:2], cex=0.05)
# points(PCA$x[na.omit(values(stackPCA[[5]]))==4,1:2], cex=0.05)
points(PCA$x[na.omit(values(stackPCA[[6]]))==1,1:2], cex=0.5, col="red", pch=20) # points envisagés

# Il y a un trou / il manquerait des points à faibles pentes et altitudes => à voir là où ce sont des zones favorables aux myrtilles

# Problème : on n'a pas intégré les différents habitats, puisque l'ACP ne prend pas de variables catégorielles...
#       => choix de méthodes mixtes ?


# 2) Analyse multivariée mixte + clustering

# # Analyse multivariée ade4
# AnMult = dudi.mix(na.omit(values(stackPCA[[1:5]])), nf=4)
# 
# # Clustering kmeans du package terra
# # Normalisation des données (sinon l'altitude a beaucoup plus de poids)
# nx <- minmax(stackPCA[[1:4]])    
# rn <- (stackPCA[[1:4]] - nx[1,]) / (nx[2,] - nx[1,])
# stackNorm = c(rn, stackPCA[[5:6]])
# # Calcul des kmeans /!\ en théorie ça ne marche pas avec des variables catégorielles !! /!\
# kmeans = k_means(na.omit(values(stackNorm[[1:5]])), centers=6)
# clusters = stackNorm[[1]]
# values(clusters)[!is.na(values(clusters))] = kmeans$cluster
# plot(clusters)
# 
# # Matrice de distance puis hclust
# matdist = daisy(na.omit(values(stackPCA[[1:5]])), metric = "gower")
# clust = hclust(matdist)
# # Trop de pixels = trop lourd !
# 
# 
# # Analyses multivariée de données mixtes (FAMD)
# tab_FAMD = as.data.frame(na.omit(values(stackPCA[[1:5]])))
# tab_FAMD$habitat = as.character(tab_FAMD$habitat)
# famd = FAMD(tab_FAMD, graph=F, ncp=9)
# 
# # fviz_famd_var(famd)
# fviz_eig(famd) # 4 ou 5 dimensions
# # fviz_famd_ind(famd, habillage=5)
# 
# clustering = HCPC(famd, graph=F) # TROP LOURD !!!


# Utilisation de hclust avec la distance de Gower (fonctionne avec les variables catégorielles)
h_clust <- function(x, database=as.data.frame(na.omit(values(x))), 
                    ngroups,  clust_method="complete", #dist_metric="euclidean",
                    types=list(numeric=1:dim(x)[3]), agfuns=c(mean,mfv1), 
                    matchfun=gower.dist, ..., #matchfun="squared"
                    maxcell=10000, filename="", overwrite=FALSE, wopt=list()) {
  stopifnot(maxcell > 0)
  stopifnot(ngroups > 0)
  stopifnot(ngroups < maxcell)
  d <- na.omit(spatSample(x, maxcell, "regular"))
  # print(dim(d))
  # dd <- stats::dist(d, dist_metric) # adaptation par rapport au script initial
  # hc <- stats::hclust(dd, clust_method) # adaptation par rapport au script initial
  dd <- StatMatch::gower.dist(d) # adaptation par rapport au script initial
  hc <- stats::hclust(as.dist(dd), clust_method) # adaptation par rapport au script initial
  th <- sort(hc$height, TRUE)[ngroups]
  cls <- stats::cutree(hc, h = th)
  hc <- cut(stats::as.dendrogram(hc), h=th)$upper
  # d <- aggregate(d, list(cls=cls), agfun)
  d <- merge(as.data.frame(aggregate(d[,types[["numeric"]]], list(cls=cls), agfuns[[1]])),
             as.data.frame(aggregate(d[,types[["factor"]]], list(cls=cls), agfuns[[2]])), by="cls")# adaptation par rapport au script initial
  # print(d)
  cls <- d$cls 
  d$cls <- NULL
  colnames(d)=names(x)
  # b <- bestMatch(x, d, fun=matchfun, ..., filename=filename, overwrite=overwrite, wopt=wopt)	
  # print(head(d))
  bdist = gower.dist(database,d)
  b=apply(bdist,1,FUN = function(X){c(1:ncol(bdist))[X==min(X)]})
  brast=x[[1]]
  brast[!is.na(values(brast))]=b
  # b <- bestMatch(x, d, fun=matchfun, ..., filename=filename, overwrite=overwrite, wopt=wopt)
  return(list(barycentres=d,  dendrogram=hc, clusterrast=brast))
} # https://rdrr.io/github/rspatial/terra/src/R/k_means.R 

database = as.data.frame(na.omit(values(x)))
database$habitat = droplevels(factor(database$habitat, levels=classes_ORION$CODE, labels=classes_ORION$NOM))

clustering = h_clust(stackPCA[[1:5]], database = database,
                     ngroups=6, types=list(numeric=1:4, factor=5), maxcell=50000) 

clust_rast = clustering$clusterrast
terra::plot(clust_rast)

# Au final le clustering ressemble ++ aux différents habitats... --> essayer de voir pourquoi ?


#*---- 2. Représentativité des points présélectionnés ----
# Points d'écoacoustique ORCHAMP + 4 sites + doublons ORCHAMP

table(clust_rast[stackPCA[["potential_pts"]]])



#*---- 3. Échantillonnage aléatoire de points supplémentaires où faire du suivi de production de fruits ----




