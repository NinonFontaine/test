##############################################-
# ANALYSE DES 20 ANS DE DONNÉES PHÉNOCLIM  ----
#  Ninon Fontaine - printemps 2025            -
##############################################-


library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(data.table)
library(reshape2)
library(dplyr)
library(terra)
library(lme4)
library(lmerTest) # permet d'avoir les pvalues qui s'affichent pour les lmer
library(visreg)
library(MuMIn) # fonction r.squaredGLMM notamment
library(variancePartition)


##############################################-
# OBJECTIFS ----
##############################################-

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
# 
# Avec l'exploration des effets 'température', d'autres questions peuvent être posées : 
#   (1) l'ajout d'une variable de température améliore-t-elle les modèles phénologiques ?
#   (2) y a-t-il une tendance à l'avancement de la phénologie des arbres ?
#   (3) cette tendance est-elle différente selon les périodes étudiées (2006-2016 vs 2006-2024 vs 2016-2024) ?
#   (4) les tendances dans les décalages phénologiques sont-elles différentes selon les altitudes considérées ?
#   (5) les tendances sont-elles différentes selon les catégories d'observateur·rices (scolaires, professionnels, particuliers) ?



##############################################-
# DATA ----
##############################################-

#---- Fonds carto
cle_ggmap ="" # à récupérer "AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes"
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
# phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/data/_CLEANED_data_pheno.csv")
# Données obtenues avec le script 1_mise_en_forme_BDD.R


############################################################################################-
# VISUALISATION DES DONNEES PAR ESPÈCE ----
############################################################################################-

stades = grep("Ok 10%",unique(phenoclim$pheno_etape_value), value = T)

# En prenant toutes les données
resume_stadeesp = phenoclim %>% filter(pheno_etape_value %in% stades) %>% group_by(species, year, pheno_stade_value, pheno_etape_value) %>%
  summarise(nb_obs = length(julian_day),
            julian_day = mean(julian_day))
resume_stadeesp$date_moyenne = as.Date(resume_stadeesp$julian_day, origin = paste0(resume_stadeesp$year,"-01-01"))
write.csv(resume_stadeesp, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resume_datepheno_paresp.csv")

ggplot(resume_stadeesp, aes(x=julian_day, y=year, col=pheno_stade_value)) + geom_point() + facet_wrap(~ species) + 
  scale_color_manual(values = list("Débourrement" = "royalblue3", "Floraison" = "skyblue2", "Feuillaison" = "olivedrab3", "Changement couleur"="chocolate2")) +
  scale_x_continuous(breaks= c(79,171,263), labels=c("mars","juin","sept.")) + 
  labs(x="", y="")+
  theme(legend.position = "none", 
        panel.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major.y = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

# Pour une zone en particulier
zone_ENS = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_lim_admin/lim_zone_ENS_exemple/limite_ens.shp")
zone_ENS_buffer5km = buffer(zone_ENS, 5000)
phenoclim_ENS = phenoclim[!is.na(terra::extract(zone_ENS_buffer5km, vect(phenoclim, geom=c("coord_x_4326","coord_y_4326"), "epsg:4326"))$NOM.SITE),]

resume_stadeesp = phenoclim_ENS %>% filter(pheno_etape_value %in% stades) %>% filter(pheno_etape_value %in% stades) %>% group_by(species, year, pheno_stade_value, pheno_etape_value) %>%
  summarise(nb_obs = length(julian_day),
            julian_day = mean(julian_day))
resume_stadeesp$date_moyenne = as.Date(resume_stadeesp$julian_day, origin = paste0(resume_stadeesp$year,"-01-01"))
write.csv(resume_stadeesp, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resume_datepheno_paresp_zoneENSGresivaudan.csv")

ggplot(resume_stadeesp, aes(x=julian_day, y=year, col=pheno_stade_value)) + geom_point() + facet_wrap(~ species) + 
  scale_color_manual(values = list("Débourrement" = "royalblue3", "Floraison" = "skyblue2", "Feuillaison" = "olivedrab3", "Changement couleur"="chocolate2")) +
  scale_x_continuous(breaks= c(79,171,263), labels=c("mars","juin","sept.")) + 
  labs(x="", y="")+
  theme(legend.position = "none", 
        panel.background = element_blank(),
        panel.border = element_blank(),  
        panel.grid.major.y = element_line(linewidth = 0.5, linetype = 'solid', colour = "grey"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

############################################################################################-
# INTÉGRATION DE VARIABLES EXPLICATIVES ASSOCIÉES AUX TEMPÉRATURES ----
############################################################################################-

# Pour la piste 3, on veut corriger la variabilité interzone par la température... 
#*---- Données climatiques potentielles : comparaison des couvertures spatiales et temporelles ----
#
#         DONNEES REANALYSEES
# - Copernicus CERRA : données météo Europe, jusqu'à 2022, à une résoluton de 5.5 km 
#                     (https://climate.copernicus.eu/copernicus-regional-reanalysis-europe-cerra)
# - ERA 5-Land : données météo Monde, en temps réel j-5, à une résolution de 9 km
#                     (https://climate.copernicus.eu/climate-reanalysis)
# - CHELSA : données météo Europe, 1979-2013, à une résolution de 1 km
#                     (https://chelsa-climate.org/timeseries/)
# - SAFRAN - MeteoFrance : données météo France, 1958-2023, à une résolution de 8 km (calcul par massif x altitude : il faudrait 
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
# # 1) utiliser les températures du modèle S2M (MétéoFrance - SAFRAN) par massif x altitude x pente x exposition 
# # 2) interpoler les mesures de température des stations Phénoclim en fonction de l'altitude (cf Asse et al. 2018)

#*---- Interpolation des données existantes pour améliorer la couverture spatiale des informations 'température' ----
#*-------- 1) Utilisation des données SAFRAN par massif x topographie ----


# # PISTE 1 : DONNÉES SAFRAN SIM
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
#
# # PISTE 2 : DONNÉES SAFRAN S2M
#
# # A la place, on peut utiliser les données SAFRAN S2M, avec une simulation par massif et classe d'altitude
# # Dans un premier temps, on attribue un massif SAFRAN à chaque observation Phénoclim, et sa classe d'altitude x exposition
# # (en gardant en tête le fait qu'il va falloir interpoler les températures pour éviter les valeurs 'en escalier' liées à ce système
# # de classe)
# 
# massifs_safran = vect("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/alp_allslopes/massifs_alpes_2154.shp")
# MNT = rast("/Volumes/RESSOURCES_SIG/MNT/MNT_Alpes/MNT_Alpes_25m_L93.asc")
# pente = terrain(MNT, v="slope", unit="degrees")
# aspect = terrain(MNT, v="aspect", unit="degrees")
# stacktopo = c(MNT, pente, aspect)
# 
# 
# sites_cl_safran = phenoclim[phenoclim$nom_massif=="Alpes",] %>% distinct(id_base_site, coord_x_2154, coord_y_2154, altitude)
# 
# plot(stt_Phenoclim_buff5km)
# points(phenoclim$coord_x_4326[phenoclim$nom_massif=="Alpes"], phenoclim$coord_y_4326[phenoclim$nom_massif=="Alpes"], cex=0.3, col="red")
# lines(project(massifs_safran, "epsg:4326"), col="darkgreen")
# 
# # En fait même les données SAFRAN ne couvriront pas toute la zone, puisque de nombreuses observations Phénoclim sont en dehors de la zone
# # 'Alpes' où sont définis les massifs SAFRAN... Soit on leur associe le massif le plus proche, soit on trouve d'autres sources de données T
# 
# sites_cl_safran[,c("massif_num","massif_nom")] = extract(massifs_safran, sites_cl_safran[,c("coord_x_2154","coord_y_2154")])[,c("massif_num","nom")]
# sites_cl_safran[,c("pente","expo")] = extract(stacktopo[[2:3]], sites_cl_safran[,c("coord_x_2154","coord_y_2154")])[,2:4] #on ne prend pas l'altitude qu'on a déjà avec les infos des sites
# write.csv(sites_cl_safran, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_sites_SAFRAN.csv")
# 
# 
# # ggplot(sites_cl_safran) +
# #   # geom_polygon(data=geom(massifs_safran), aes(x=x, y=y, group=part), fill="lightgray", alpha=0.5)+
# #   geom_point(aes(x=coord_x_2154, y=coord_y_2154, color=massif_nom))+
# #   coord_equal() + xlim(760000, 1100000)



#*-------- 2) Utilisation des données des stations PhénoClim interpolées pour tous les sites d'observation ----

# En se basant sur le machine learning, sur le même principe que son utilisation pour boucher les 'gros trous' dans les données des 
# stations météo (bugs techniques), on peut obtenir des valeurs de température pour tous les sites en se basant sur leurs coordonnées 
# x, y, z (cf scripts Geoffrey Klein dans le projet R "Stations-phenoclim.Rproj").
# NOTE : ce qui m'étonne, c'est de pouvoir reconstruire ces données de température sans prendre en compte l'exposition des sites, 
#         seulement leur altitude et leur localisation (parce qu'on se base sur les stations PhénoClim les plus proches).
#   => MISE À JOUR DE LA NOTE : au final, il semble qu'ajouter l'exposition n'améliore pas les prédictions (RMSE plus élevée, R2 plus faible)
#         (cf script "/Users/ninonfontaine/Desktop/projetsR/TEST/_scripts_meteo/adaptNF__(test_3 )_test_qualite_reconstruction_ML_xgboost.R")
# 
#
# Ce qui nous intéresse dans ces données c'est :
# - le sensor 4 (T_air à 2m de haut)
# - les températures moyennes, minimales et maximales quotidiennes (qu'on pourra utiliser ensuite pour calculer une température moyenne
#   sur une période donnée, comme le printemps par exemple, ou des GDD)

# #*------------ Aperçu des données par station et de leur localisation : ----
# stt_Phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Stations-phenoclim/analyse/coordonnees_crea.csv",
#                          sep=";")
# donnees_completes_2023 = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Stations-phenoclim/data/fichiers_reconstruits_annee/2023_reconstruit.csv")
# 
# test = donnees_completes_2023 %>% group_by(station_name, month, day, x_crea, y_crea, z_crea) %>%
#   summarise(Tair_moy = mean(sensor_value4),
#             Tair_min = min(sensor_value4),
#             Tair_max = max(sensor_value4))
# test$date = as.Date(paste(test$day, test$month, "2023", sep="/"), format = "%d/%m/%Y")
# 
# ggplot(test) + 
#   geom_line(aes(x=date, y=Tair_moy, group=station_name), col="black") +
#   geom_line(aes(x=date, y=Tair_min, group=station_name), col="blue", linetype=2)+
#   geom_line(aes(x=date, y=Tair_max, group=station_name), col="red", linetype=2) + 
#   xlim(as.Date(c("15/02/2023","15/04/2023"), format = "%d/%m/%Y"))

#*------------ Récupération des données météo reconstruites et calcul de variables pertinentes pour la phénologie ----

# # Ces données reconstruites sont résumées en Tmin, Tmoy, Tmax par jour et par site, pour chaque année
# # (Cela représente 3503 sites x 365 jours...)
# reconstruc_sites_2023 = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/2023_Tairday_SitesPheno_reconstruit.csv")

# Pour les études phénologiques, on ne veut pas des températures quotidiennes mais plutôt une valeur de température informative,
# par exemple la température moyenne pendant la période avant la phase phénologique d'intérêt.
# Comment choisir cette période ? 
# - Se baser sur une période unique par année ? ex : on choisit une méthode de calcul commune à tous les sites, basée sur les 60j
#   avant la période globale de débourrement (cf quantiles de la période de débourrement) ?
# - Se baser sur un calcul par site x année ?

# debourrBpen_Alps = phenoclim[phenoclim$pheno_etape_value == "Debourrement - Ok 10%" & 
#                             phenoclim$nom_massif == "Alpes" &
#                               phenoclim$nom_cite == "Betula pendula Roth, 1788",]
# # ggplot(debourrBpen_Alps) + geom_histogram(aes(x=julian_day))
# ggplot(debourrBpen_Alps, aes(y=julian_day, x=altitude)) + geom_point(aes(col=as.factor(year))) + geom_smooth(method=lm, col="black")

# # APERÇU
# quantile_df <- function(x, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)) {
#   tibble(
#     val = quantile(x, probs, na.rm = TRUE),
#     quant = probs
#   )
# }
# test = debourrBpen_Alps %>% group_by(year) %>%
#   reframe(quantile_df(julian_day))
# test = spread(test, key=quant, value = val, fill=NA)
# test$duree = test$`1` - test$`0`
# test$duree80pourc = test$`0.9` - test$`0.1`
# test$duree50pourc = test$`0.75` - test$`0.25`
# #     => Le débourrement du bouleau s'étale sur 120 jours certaines années ! Quelle période sélectionner ?? 60 jours avant la médiane, 30 ?
# #         Fonctionner par site, avec 30 jours avant le débourrement de chaque site, semblerait plus cohérent
# ggplot(debourrBpen_Alps) +
#   geom_boxplot(aes(x=julian_day, y=year, group=year)) +
#   geom_point(data=test, aes(x=`0.5`-60, y=year), col="red")



# #=============================================================================================================================*
# # VERSION 1 - ON CALCULE DES TEMPÉRATURES BASÉES SUR LA DATE DE DÉBOURREMENT PAR SITE ET LA PÉRIODE PRÉCÉDENT LE DÉBOURREMENT
# Tperiod <- function(jourcible, site,  tabT, 
#                     vartabT = c("julian_day","station_name","Tair_moy"),
#                     dureej, tempcalc="Tmoy_Xj"){
#   # NB : le tableau de températures en entrée est un tableau avec une seule année, donc on peut fonctionner en "julian day"
#   tabT = rename(tabT, c(jul_day=vartabT[1], sitePheno=vartabT[2], T=vartabT[3]))
#   output = matrix(ncol=length(tempcalc))
#   colnames(output)=tempcalc
#   # tabT$date = as.Date(tabT$date)
#   # jourcible = as.Date(jourcible)
#   if ("Tmoy_Xj" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#     output[,"Tmoy_Xj"] = mean(Tcalc$T, na.rm=T)
#   }
#   if ("Tmin" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#     output[,"Tmin"]=min(Tcalc$T, na.rm=T)
#   }
#   if ("Tmax" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#     output[,"Tmax"]=max(Tcalc$T, na.rm=T)
#   }
#   if ("GDD0" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
#     output[,"GDD0"]=sum(Tcalc$T[Tcalc$T > 0], na.rm=T)
#   }
#   if ("GDD5" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
#     output[,"GDD5"]=sum(Tcalc$T[Tcalc$T >=5]-5, na.rm=T)
#   }
#   if ("Tmoyhiv" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= 31 & jul_day >= -60 & sitePheno ==site) #on définit l'hiver comme la période novembre-janvier (Vitasse)
#     output[,"Tmoyhiv"]=mean(Tcalc$T, na.rm=T)
#   }
#   if ("dChill" %in% tempcalc){
#     Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= -60 & sitePheno ==site) #on calcule les jours de chilling du 1/11 au débourrement (Vitasse) - en proportion
#     output[,"dChill"]=nrow(Tcalc[Tcalc$T <=8 & Tcalc$T>=0,]) / nrow(Tcalc)
#   }
#   return(output)
# }
# 
# 
# # Construction d'un tableau sur le débourrement dans les Alpes, où chaque donnée (site x année) est associée à une valeur de température supposée
# # pertinente pour la phénologie : la température moyenne pendant les 30 jours précédant la date de débourrement.
# 
# debourr_Alps = phenoclim[phenoclim$pheno_etape_value == "Debourrement - Ok 10%" & 
#                            phenoclim$nom_massif == "Alpes",]
# # /!\ il y a beaucoup de lignes remplies de NA -> à supprimer :
# debourr_Alps = debourr_Alps[!is.na(debourr_Alps$ids_observers),]
# 
# varTselec = c("Tmoy_Xj","GDD0","GDD5","Tmoyhiv","dChill")
# debourr_Alps[,varTselec] = NA
# 
# for(annee in c(2006:2024)){#unique(debourr_Alps$year)){
#   print(annee)
#   # reconstruc_Tday_sites = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
#   # reconstruc_Tday_sites$julian_day = yday(as.Date(reconstruc_Tday_sites$date))
#   # Tcalcs = apply(debourr_Alps[debourr_Alps$year == annee, c("julian_day", "id_base_site")],1, 
#   #                FUN=function(x){Tperiod(jourcible = x[1], site=x[2], dureej = 30, tabT=reconstruc_Tday_sites)})
#   # # ---------------------------------------------------------------------------------------------------------------------------*
#   # # /!\ pour certaines espèces précoces, pour 30 jours avant la date de débourrement on passe à l'année précédente !!! 
#   # #     => corriger le calcul pour le cas de ces espèces précoces !!
#   # # ---------------------------------------------------------------------------------------------------------------------------*
#   if (annee != 2006){ # si on est en 2006, c'est la première année où on a des infos météo...
#     reconstruc_Tday_sites_an0 = reconstruc_Tday_sites_an1
#     # on renumérote les jours de l'année précédente en négatif
#     reconstruc_Tday_sites_an0$julian_day = reconstruc_Tday_sites_an0$julian_day - yday(as.Date(paste0(annee-1,"/12/31"))) # pour prendre en compte les années bisextiles
#     }
#   reconstruc_Tday_sites_an1 = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
#   reconstruc_Tday_sites_an1$julian_day = yday(as.Date(reconstruc_Tday_sites_an1$date))
#   if (annee != 2006){
#     reconstruc_Tday_sites = rbind(reconstruc_Tday_sites_an0, reconstruc_Tday_sites_an1)}
#   else {
#     reconstruc_Tday_sites = reconstruc_Tday_sites_an1
#   }
#   
#   Tcalcs = apply(debourr_Alps[debourr_Alps$year == annee, c("julian_day", "id_base_site")],1, 
#                  FUN=function(x){Tperiod(jourcible = x[1], site=x[2], dureej = 30, 
#                                          tabT=reconstruc_Tday_sites,
#                                          tempcalc = varTselec)})
#   
#   debourr_Alps[debourr_Alps$year == annee, varTselec] = matrix(unlist(Tcalcs),ncol=length(varTselec),byrow=T)
# }
# 
# debourr_Alps = debourr_Alps %>% rename("Tmoy30j"="Tmoy_Xj")
# 
# write.csv(debourr_Alps, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Debourr_Alps_T.csv")



#=============================================================================================================================*
# VERSION 2 - ON CALCULE DES TEMPÉRATURES BASÉES SUR LA DATE *MÉDIANE* DE DÉBOURREMENT ET LA PÉRIODE PRÉCÉDENT CETTE MÉDIANE

debourr_Alps = phenoclim[phenoclim$pheno_etape_value == "Debourrement - Ok 10%" &
                           phenoclim$nom_massif == "Alpes",]
# /!\ il y a beaucoup de lignes remplies de NA -> à supprimer :
debourr_Alps = debourr_Alps[!is.na(debourr_Alps$ids_observers),]


# 1) Calcul de la date de débourrement médiane par espèce, sur différentes périodes
dates_deb_med = debourr_Alps %>% group_by(species) %>% 
  summarise(med_0616 = median(julian_day[year %in% 2006:2016], na.rm=T),
            med_1624 = median(julian_day[year %in% 2016:2024], na.rm=T),
            med_0624 = median(julian_day[year %in% 2006:2024], na.rm=T))

# Tperiod_med__ <- function(site, esp, dates_medianes, tabT,
#                     vartabT = c("julian_day","station_name","Tair_moy"),
#                     dureej, tempcalc="Tmoy_Xj"){
#   # NB : le tableau de températures en entrée est un tableau avec 2 années, on peut fonctionner en "julian day" parce qu'on a mis des jours de -365 à +365
#   tabT = rename(tabT, c(jul_day=vartabT[1], sitePheno=vartabT[2], T=vartabT[3]))
#   output = matrix(ncol=length(tempcalc) * (ncol(dates_medianes)-1))
#   dimnames(output)[[1]]=paste0(tempcalc,"_",rep(colnames(dates_medianes)[-1],each=length(tempcalc)))
#   
#   jourscibles = dates_medianes[dates_medianes$species==esp,-1]
# 
#   if ("Tmoy_Xj" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#         output[,paste0("Tmoy_Xj_",colnames(jourscibles)[i])] = mean(Tcalc$T, na.rm=T) 
#       }
#   }}
#   if ("Tmin" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#         output[,paste0("Tmin_",colnames(jourscibles)[i])] = min(Tcalc$T, na.rm=T)
#       }
#   }}
#   if ("Tmax" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
#         output[,paste0("Tmax_",colnames(jourscibles)[i])] = max(Tcalc$T, na.rm=T)
#       }
#   }}
#   if ("GDD0" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
#         output[,paste0("GDD0_",colnames(jourscibles)[i])] = sum(Tcalc$T[Tcalc$T > 0], na.rm=T)
#       }
#   }}
#   if ("GDD5" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
#         output[,paste0("GDD5_",colnames(jourscibles)[i])] = sum(Tcalc$T[Tcalc$T >=5]-5, na.rm=T)
#       }
#   }}
#   if ("Tmoyhiv" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       Tcalc = tabT %>% filter(jul_day <= 31 & jul_day >= -60 & sitePheno ==site) #on définit l'hiver comme la période novembre-janvier (Vitasse)
#       output[,paste0("Tmoyhiv_",colnames(jourscibles)[i])] = mean(Tcalc$T, na.rm=T)
#       }
#   }
#   if ("dChill" %in% tempcalc){
#     for (i in 1:length(jourscibles)){
#       jourcible = as.numeric(jourscibles[i])
#       if (!is.na(jourcible)){
#         Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= -60 & sitePheno ==site) #on calcule les jours de chilling du 1/11 au débourrement (Vitasse) - en proportion
#         output[,paste0("dChill_",colnames(jourscibles)[i])] = nrow(Tcalc[Tcalc$T <=8 & Tcalc$T>=0,]) / nrow(Tcalc)
#       }
#   }}
#   return(output)
# }
# 


# Construction d'un tableau sur le débourrement dans les Alpes, où chaque donnée (site x année) est associée à une valeur de température supposée
# pertinente pour la phénologie : la température moyenne pendant les 30 jours précédant la date de débourrement.

varTselec = c("Tmoy30j","Tmoy40j","Tmoy50j","gel30j","gel40j","gel50j","GDD0","GDD5","Tmoyhiv","dChill")
output = data.frame(matrix(ncol=4+length(varTselec)))
colnames(output) = c("sitePheno",varTselec, "periode","esp","annee")
# debourr_Alps[,paste0(varTselec,"_",rep(colnames(dates_deb_med)[-1],each=length(varTselec)))] = NA

vartabT = c("julian_day","station_name","Tair_moy", "Tair_min")

for(annee in c(2006:2024)){#unique(debourr_Alps$year)){
  print(annee)
  # reconstruc_Tday_sites = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
  # reconstruc_Tday_sites$julian_day = yday(as.Date(reconstruc_Tday_sites$date))
  # Tcalcs = apply(debourr_Alps[debourr_Alps$year == annee, c("julian_day", "id_base_site")],1,
  #                FUN=function(x){Tperiod(jourcible = x[1], site=x[2], dureej = 30, tabT=reconstruc_Tday_sites)})
  # # ---------------------------------------------------------------------------------------------------------------------------*
  # # /!\ pour certaines espèces précoces, pour 30 jours avant la date de débourrement on passe à l'année précédente !!!
  # #     => corriger le calcul pour le cas de ces espèces précoces !!
  # # ---------------------------------------------------------------------------------------------------------------------------*
  if (annee != 2006){ # si on est en 2006, c'est la première année où on a des infos météo...
    reconstruc_Tday_sites_an0 = reconstruc_Tday_sites_an1
    # on renumérote les jours de l'année précédente en négatif
    reconstruc_Tday_sites_an0$julian_day = reconstruc_Tday_sites_an0$julian_day - yday(as.Date(paste0(annee-1,"/12/31"))) # pour prendre en compte les années bisextiles
    }
  reconstruc_Tday_sites_an1 = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
  reconstruc_Tday_sites_an1$julian_day = yday(as.Date(reconstruc_Tday_sites_an1$date))
  if (annee != 2006){
    reconstruc_Tday_sites = rbind(reconstruc_Tday_sites_an0, reconstruc_Tday_sites_an1)}
  else {
    reconstruc_Tday_sites = reconstruc_Tday_sites_an1
  }
  tabT = rename(reconstruc_Tday_sites, c(jul_day=vartabT[1], sitePheno=vartabT[2], T=vartabT[3], , Tmin=vartabT[4]))
  
  # Tcalcs_per = apply(debourr_Alps[debourr_Alps$year == annee, c("id_base_site", "species")],1,
  #                    FUN=function(x){Tperiod_med(esp = x[2], site=x[1], dates_medianes = dates_deb_med, dureej = 30,
  #                                            tabT=reconstruc_Tday_sites,
  #                                            tempcalc = varTselec)})
  for (esp in unique(debourr_Alps$species[debourr_Alps$year == annee])){
    jourscibles = dates_deb_med[dates_deb_med$species==esp,-1]
    for (periode in colnames(jourscibles)){
      jourcible = as.numeric(jourscibles[periode])
      Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-30) %>% group_by(sitePheno) %>%
        summarise(Tmoy30j = mean(T, na.rm=T),
                  gel30j = length(Tmin[Tmin<0]))
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-40) %>% group_by(sitePheno) %>%
                      summarise(Tmoy40j = mean(T, na.rm=T),
                                gel40j = length(Tmin[Tmin<0])))
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-50) %>% group_by(sitePheno) %>%
                      summarise(Tmoy50j = mean(T, na.rm=T),
                                gel50j = length(Tmin[Tmin<0])))
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= jourcible & jul_day >= 0) %>% group_by(sitePheno) %>%
                      summarise(GDD0 = sum(T[T > 0], na.rm=T),
                                GDD5 = sum(T[T >=5]-5, na.rm=T)), 
                    by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= 31 & jul_day >= -60) %>% group_by(sitePheno) %>%
                      summarise(Tmoyhiv = mean(T, na.rm=T)), 
                    by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= jourcible & jul_day >= -60) %>% group_by(sitePheno) %>%
                      summarise(dChill = length(T[T <=8 & T>=0]) / length(T)), 
                    by="sitePheno",all.x=T, all.y=T)
      
      Tcalc$periode = periode
      Tcalc$esp = esp
      Tcalc$annee = annee
      
      output = rbind(output, Tcalc)
      
    }
    
    
  }

}

output = output[!is.na(output$sitePheno),]
Tperiodes = data.table::dcast(data=setDT(output), formula=sitePheno+esp+annee ~ periode, value.var=varTselec)

# ggplot(Tperiodes, aes(y=Tmoy30j_med_0624, x=annee, col=sitePheno)) + geom_point() + geom_smooth(method='lm')

write.csv(Tperiodes, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/varTprintemps_allphenosites.csv")


debourr_Alps = merge(debourr_Alps, Tperiodes, by.x=c("id_base_site","species","year"), by.y=c("sitePheno","esp","annee"), all.x=T, all.y=F)

write.csv(debourr_Alps, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Debourr_Alps_T.csv")



############################################################################################-
# MODÈLES DE DÉBOURREMENT                                                                ----
############################################################################################-

debourr_Alps = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Debourr_Alps_T.csv", row.names=1)
debourr_Alps$yearQ = factor(debourr_Alps$year)
debourr_Alps$cl_alt = factor(debourr_Alps$cl_alt, levels=c("150-450","450-750" ,  "750-1050"   ,"1050-1350" , "1350-1650" ,"1650-1950", "1950-2250" ), 
                             ordered = T)

# On renomme les variables de température calculées sur la date médiane de débourrement de la période 2006-2024 ('par défaut')
debourr_Alps = debourr_Alps %>% rename(Tmoy30j = Tmoy30j_med_0624,
                                       Tmoy40j = Tmoy40j_med_0624,
                                       Tmoy50j = Tmoy50j_med_0624,
                                       gel30j = gel30j_med_0624,
                                       gel40j = gel40j_med_0624,
                                       gel50j = gel50j_med_0624,
                                       GDD0 = GDD0_med_0624,
                                       GDD5 = GDD5_med_0624,
                                       Tmoyhiv = Tmoyhiv_med_0624,
                                       dChill = dChill_med_0624)

# Remarque : au vu des corrélations entre ces variables de température, il est peu logique d'en intégrer plusieurs si elles sont trop corrélées
varTselec = c("Tmoy30j","Tmoy40j","Tmoy50j","gel30j","gel40j","gel50j","GDD0","GDD5","Tmoyhiv","dChill")
corrplot::corrplot(cor(na.omit(debourr_Alps[,varTselec])))
# => on choisit une parmi Tmoy30, 40, 50j, GDD0, GDD5, qu'on combine avec Tmoyhiv et/ou dChill

# Aperçu des données de débourrement
ggplot(debourr_Alps[!is.na(debourr_Alps$Tmoy30j),], 
       aes(x=altitude, y=julian_day, col=nom_cite)) + geom_point() + theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
ggplot(debourr_Alps[!is.na(debourr_Alps$Tmoy30j),], 
       aes(x=Tmoy30j, y=julian_day, col=nom_cite)) + geom_point() + theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
# À première vue, la relation entre date de débourrement et altitude est plus nette que la relation entre date et température (notamment pour 
# l'épicéa) --> tester dans les modèles si l'ajout de cette variable température améliore les modèles par rapport à des modèles avec seulement
# l'altitude.


# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce, pour les 2 périodes)
resultats = data.frame(species = NA, periode = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(debourr_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()

##############################################-
#*---- Bouleau verruqueux ----

deb_Bpen = debourr_Alps[debourr_Alps$species == "Bouleau_verruqueux",]
# ggplot(deb_Bpen, aes(x=year)) + geom_histogram(binwidth=1) 
# # on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Bpen = deb_Bpen[!is.na(deb_Bpen$Tmoy30j),]
# ggplot(deb_Bpen, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Bpen_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Bpen, REML=F)
summary(mod_deb_Bpen_Alt)
altyearmods = c(list("Bouleau_verruqueux"=mod_deb_Bpen_Alt), altyearmods)
deb_Bpen_10ans = deb_Bpen[deb_Bpen$year <= 2016,]
mod_deb_Bpen_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Bpen_10ans, REML=F)
summary(mod_deb_Bpen_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Bpen_Alt)
# + Validation croisée
n = nrow(deb_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Bpen[trainIndex ,]
test <- deb_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Bpen$ID_zone[drop=T])[!unique(deb_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Bpen_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6689
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Bpen_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)

AIC(mod_deb_Bpen_Tmoy30j) ; AIC(mod_deb_Bpen_Tmoy40j) ; AIC(mod_deb_Bpen_Tmoy50j) ; AIC(mod_deb_Bpen_gel30j) ; AIC(mod_deb_Bpen_gel40j) ; AIC(mod_deb_Bpen_gel50j) ; AIC(mod_deb_Bpen_GDD0) ; AIC(mod_deb_Bpen_GDD5); AIC(mod_deb_Bpen_Tmoyhiv) ; AIC(mod_deb_Bpen_dChill)
# La variable Tmoy40j donne les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
# /!\ ce qui serait le plus logique serait de ne garder que la meilleure des 5 variables de Tmoy et GDD (cf corrélations) et éventuellement Tmoyhiv et/ou dChill
mod_deb_Bpen_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
# On enlève les variables une par une : dChill puis Tmoy30j puis Tmoyhiv --> c'est le modèle avec le meilleur AIC... MAIS GDD5 n'a pas un effet 
# significatif au seuil 0.05... et les variables de température Tmoy+GDD sont super corrélées ! 
# En combinant ces 3 critères (AIC, significativité des effets, variables peu corrélées), on retient donc :
mod_deb_Bpen_multiT = lmer(julian_day ~ Tmoy40j + gel50j + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F) # AIC = 12583.2
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Bpen, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Bpen, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Bpen_Tmoy40j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Bpen_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j + year + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Bpen, REML=F)
mod_deb_Bpen_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Bpen, REML=F)
AIC(mod_deb_Bpen_Tmoy40j) ; AIC(mod_deb_Bpen_Tmoy40j_Alt) ; AIC(mod_deb_Bpen_Tmoy40j_year) ; AIC(mod_deb_Bpen_Tmoy40j_clAlt) ; AIC(mod_deb_Bpen_Tmoy30j_Alt_) ; AIC(mod_deb_Bpen_Tmoy40j_Alt_) ; AIC(mod_deb_Bpen_Tmoy40j_clAlt_)
# Mieux avec l'altitude en plus, et en ayant l'effet aléatoire 'année' sur la température (AIC = 12453.3)

#========= MEILLEUR MODELE :
bestmod_deb_Bpen = lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Bpen, REML=F) # AIC = 12453.3
bestmods = c(list("Bouleau_verruqueux"=bestmod_deb_Bpen), bestmods)


summary(bestmod_deb_Bpen)
# gradient altitudinal : 1.8 jour de retard quand on monte de 100m
# effet de la température moyenne dans les 30 jours précédant le débourrement : 1.4 jour d'avance quand on gagne 1°C

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Bpen)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Bpen)
# - Validation croisée
n = nrow(deb_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Bpen[trainIndex ,]
test <- deb_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Bpen$ID_zone[drop=T])[!unique(deb_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Bpen), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Bpen))),
                                        coef = coef(summary(bestmod_deb_Bpen))[,1],
                                        std = coef(summary(bestmod_deb_Bpen))[,2],
                                        pval = coef(summary(bestmod_deb_Bpen))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Bpen)[rownames(coef(summary(bestmod_deb_Bpen)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Bpen_10ans = deb_Bpen[deb_Bpen$year <= 2016,]
bestmod_deb_Bpen_10ans <- lmer(julian_day ~ Tmoy30j_med_0616 + altitude + (1|ID_zone) + (Tmoy30j_med_0616|yearQ), deb_Bpen_10ans, REML=F)
summary(bestmod_deb_Bpen_10ans)
# On a les mêmes tendances sur les deux périodes (même en prenant la température calculée avec la médiane des débourrements 2006-2016)


resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Bpen_10ans))),
                                        coef = coef(summary(bestmod_deb_Bpen_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Bpen_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Bpen_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Bpen_10ans)[rownames(coef(summary(bestmod_deb_Bpen_10ans)))], error = function(e) return(NA))))


##############################################-
#*---- Noisetier ----

deb_Cave = debourr_Alps[debourr_Alps$species == "Noisetier",]
# ggplot(deb_Cave, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Cave = deb_Cave[!is.na(deb_Cave$Tmoy30j),]
# ggplot(deb_Cave, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Noisetier", x="année",y="jour julien")


#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Cave_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Cave, REML=F)
summary(mod_deb_Cave_Alt)
altyearmods = c(list("Noisetier"=mod_deb_Cave_Alt), altyearmods)
deb_Cave_10ans = deb_Cave[deb_Cave$year <= 2016,]
mod_deb_Cave_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Cave_10ans, REML=F)
summary(mod_deb_Cave_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019 (mais effet année plus faible)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Noisetier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Cave_Alt)
# + Validation croisée
n = nrow(deb_Cave)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Cave[trainIndex ,]
test <- deb_Cave[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Cave$ID_zone[drop=T])[!unique(deb_Cave$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Cave_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.5817
R2_models[R2_models$species == "Noisetier", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Cave_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)

AIC(mod_deb_Cave_Tmoy30j) ; AIC(mod_deb_Cave_Tmoy40j) ; AIC(mod_deb_Cave_Tmoy50j) ;AIC(mod_deb_Cave_gel30j) ; AIC(mod_deb_Cave_gel40j) ; AIC(mod_deb_Cave_gel50j) ; AIC(mod_deb_Cave_GDD0) ; AIC(mod_deb_Cave_GDD5); AIC(mod_deb_Cave_Tmoyhiv) ; AIC(mod_deb_Cave_dChill)
# La variable Tmoy50j donne les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Cave_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
# On enlève les variables une par une : Tmoy30j puis dChill puis Tmoyhiv puis Tmoy50j--> c'est le modèle avec le meilleur AIC... MAIS les variables
# restantes sont quand même fortement corrélées...!
mod_deb_Cave_multiT = lmer(julian_day ~ Tmoy40j + GDD5 + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F) # AIC = 12845.5
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Cave, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Cave, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Cave_Tmoy40j, "Tmoy40j")
# ggplot(deb_Bpen, aes(x=Tmoy40j, y=julian_day)) + geom_point() + stat_ellipse(aes(col=yearQ), type="norm", level=0.15, lwd=1.5)

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Cave_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j + year + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40jGDD5_Alt_ <- lmer(julian_day ~ Tmoy40j + GDD5 + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Cave, REML=F)
mod_deb_Cave_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Cave, REML=F)
AIC(mod_deb_Cave_Tmoy40j) ; AIC(mod_deb_Cave_Tmoy40j_Alt) ; AIC(mod_deb_Cave_Tmoy40j_year) ; AIC(mod_deb_Cave_Tmoy40j_clAlt) ; AIC(mod_deb_Cave_Tmoy40j_Alt_) ; AIC(mod_deb_Cave_Tmoy40jGDD5_Alt_) ; AIC(mod_deb_Cave_Tmoy40j_clAlt_)
# Mieux avec l'altitude en plus, en ayant une combinaison GDD5 + Tmoy40j, et en ayant l'effet aléatoire 'année' sur la température (AIC = 12798.2)

#========= MEILLEUR MODELE :
bestmod_deb_Cave = lmer(julian_day ~ Tmoy40j + GDD5 + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Cave, REML=F) # AIC = 12798.2
bestmods = c(list("Noisetier"=bestmod_deb_Cave), bestmods)

summary(bestmod_deb_Cave)
# gradient altitudinal : 1.9 jour de retard quand on monte de 100m
# effet de la température moyenne dans les 40 jours précédant le débourrement : 3.5 jours d'avance quand on gagne 1°C

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Noisetier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Cave)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Cave)
# - Validation croisée
n = nrow(deb_Cave)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Cave[trainIndex ,]
test <- deb_Cave[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Cave$ID_zone[drop=T])[!unique(deb_Cave$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Cave), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6328
R2_models[R2_models$species == "Noisetier", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Noisetier",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Cave))),
                                        coef = coef(summary(bestmod_deb_Cave))[,1],
                                        std = coef(summary(bestmod_deb_Cave))[,2],
                                        pval = coef(summary(bestmod_deb_Cave))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Cave)[rownames(coef(summary(bestmod_deb_Cave)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Cave_10ans = deb_Cave[deb_Cave$year <= 2016,]
bestmod_deb_Cave_10ans <- lmer(julian_day ~ Tmoy40j_med_0616 + GDD5_med_0616 + altitude + (1|ID_zone) + (Tmoy40j_med_0616|yearQ), deb_Cave_10ans, REML=F)
summary(bestmod_deb_Cave_10ans)
# On a les mêmes tendances sur les deux périodes (même en prenant la température calculée avec la médiane des débourrements 2006-2016)


resultats = rbind(resultats, data.frame(species = "Noisetier",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Cave_10ans))),
                                        coef = coef(summary(bestmod_deb_Cave_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Cave_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Cave_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Cave_10ans)[rownames(coef(summary(bestmod_deb_Cave_10ans)))], error = function(e) return(NA))))



##############################################-
#*---- Frêne ----
deb_Fexc = debourr_Alps[debourr_Alps$species == "Frene",]
# ggplot(deb_Fexc, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Fexc = deb_Fexc[!is.na(deb_Fexc$Tmoy30j),]
# ggplot(deb_Fexc, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Frêne", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Fexc_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Fexc, REML=F)
summary(mod_deb_Fexc_Alt)
altyearmods = c(list("Frene"=mod_deb_Fexc_Alt), altyearmods)
deb_Fexc_10ans = deb_Fexc[deb_Fexc$year <= 2016,]
mod_deb_Fexc_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Fexc_10ans, REML=F)
summary(mod_deb_Fexc_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019 (mais effet année plus faible)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Frene", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Fexc_Alt)
# + Validation croisée
n = nrow(deb_Fexc)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Fexc[trainIndex ,]
test <- deb_Fexc[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Fexc$ID_zone[drop=T])[!unique(deb_Fexc$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Fexc_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6114
R2_models[R2_models$species == "Frene", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Fexc_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)

AIC(mod_deb_Fexc_Tmoy30j) ; AIC(mod_deb_Fexc_Tmoy40j) ; AIC(mod_deb_Fexc_Tmoy50j) ; AIC(mod_deb_Fexc_gel30j) ; AIC(mod_deb_Fexc_gel40j) ; AIC(mod_deb_Fexc_gel50j) ; AIC(mod_deb_Fexc_GDD0) ; AIC(mod_deb_Fexc_GDD5); AIC(mod_deb_Fexc_Tmoyhiv) ; AIC(mod_deb_Fexc_dChill)
# La variable Tmoy30j donne les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Fexc_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
# On enlève les variables une par une : Tmoy30j puis dChill puis Tmoyhiv puis Tmoy50j--> c'est le modèle avec le meilleur AIC... MAIS les variables
# restantes sont quand même fortement corrélées...!
mod_deb_Fexc_multiT = lmer(julian_day ~ Tmoy30j + GDD0 + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F) # AIC = 13924.1
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Fexc, aes(x=year, y=Tmoy30j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Fexc, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Fexc_Tmoy30j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Fexc_Tmoy30j_Alt <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy30j_year <- lmer(julian_day ~ Tmoy30j + year + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy30j_clAlt <- lmer(julian_day ~ Tmoy30j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy30jGDD0_Alt_ <- lmer(julian_day ~ Tmoy30j + GDD0 + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_GDD0_Alt_ <- lmer(julian_day ~ GDD0 + altitude + (1|ID_zone) + (GDD0|yearQ), deb_Fexc, REML=F)
mod_deb_Fexc_Tmoy30j_clAlt_ <- lmer(julian_day ~ Tmoy30j * cl_alt + (1|ID_zone) + (Tmoy30j|yearQ), deb_Fexc, REML=F)
AIC(mod_deb_Fexc_Tmoy30j) ; AIC(mod_deb_Fexc_Tmoy30j_Alt) ; AIC(mod_deb_Fexc_Tmoy30j_year) ; AIC(mod_deb_Fexc_Tmoy30j_clAlt) ; AIC(mod_deb_Fexc_Tmoy30j_Alt_); AIC(mod_deb_Fexc_Tmoy30jGDD0_Alt_); AIC(mod_deb_Fexc_GDD0_Alt_) ; AIC(mod_deb_Fexc_Tmoy30j_clAlt_)


#========= MEILLEUR MODELE :
bestmod_deb_Fexc = lmer(julian_day ~ GDD0 + Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Fexc, REML=F) # AIC = 13761.28
bestmods = c(list("Frene"=bestmod_deb_Fexc), bestmods)

summary(bestmod_deb_Fexc)
# gradient altitudinal : 3.3 jours de retard quand on monte de 100m
# effet de la température moyenne dans les 30 jours précédant le débourrement : 0.5 jours d'avance quand on gagne 1°C MAIS non significatif (pval = 0.49)

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Frene", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Fexc)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Fexc)
# - Validation croisée
n = nrow(deb_Fexc)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Fexc[trainIndex ,]
test <- deb_Fexc[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Fexc$ID_zone[drop=T])[!unique(deb_Fexc$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Fexc), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7564
R2_models[R2_models$species == "Frene", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Frene",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Fexc))),
                                        coef = coef(summary(bestmod_deb_Fexc))[,1],
                                        std = coef(summary(bestmod_deb_Fexc))[,2],
                                        pval = coef(summary(bestmod_deb_Fexc))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Fexc)[rownames(coef(summary(bestmod_deb_Fexc)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Fexc_10ans = deb_Fexc[deb_Fexc$year <= 2016,]
bestmod_deb_Fexc_10ans <- lmer(julian_day ~ Tmoy30j_med_0616 + GDD0_med_0616 + altitude + (1|ID_zone) + (Tmoy30j_med_0616|yearQ), deb_Fexc_10ans, REML=F)
summary(bestmod_deb_Fexc_10ans)
# On a les mêmes tendances sur les deux périodes (même en prenant la température calculée avec la médiane des débourrements 2006-2016)


resultats = rbind(resultats, data.frame(species = "Frene",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Fexc_10ans))),
                                        coef = coef(summary(bestmod_deb_Fexc_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Fexc_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Fexc_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Fexc_10ans)[rownames(coef(summary(bestmod_deb_Fexc_10ans)))], error = function(e) return(NA))))


##############################################-
#*---- Mélèze ----

deb_Ldec = debourr_Alps[debourr_Alps$species == "Meleze",]
# ggplot(deb_Ldec, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Ldec = deb_Ldec[!is.na(deb_Ldec$Tmoy30j),]
# ggplot(deb_Ldec, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Mélèze", x="année",y="jour julien")
# On retire les points de basse altitude (il n'y en a que 2 et ça pose problème pour la suite des modèles)
deb_Ldec = deb_Ldec[deb_Ldec$cl_alt != "150-450",]

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Ldec_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Ldec, REML=F)
summary(mod_deb_Ldec_Alt)
altyearmods = c(list("Meleze"=mod_deb_Ldec_Alt), altyearmods)
deb_Ldec_10ans = deb_Ldec[deb_Ldec$year <= 2016,]
mod_deb_Ldec_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Ldec_10ans, REML=F)
summary(mod_deb_Ldec_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019 (mais effet année plus fort)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Ldec_Alt)
# + Validation croisée
n = nrow(deb_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Ldec[trainIndex ,]
test <- deb_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Ldec$ID_zone[drop=T])[!unique(deb_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Ldec_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6565
R2_models[R2_models$species == "Meleze", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Ldec_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)

AIC(mod_deb_Ldec_Tmoy30j) ; AIC(mod_deb_Ldec_Tmoy40j) ; AIC(mod_deb_Ldec_Tmoy50j) ; AIC(mod_deb_Ldec_gel30j) ; AIC(mod_deb_Ldec_gel40j) ; AIC(mod_deb_Ldec_gel50j) ; AIC(mod_deb_Ldec_GDD0) ; AIC(mod_deb_Ldec_GDD5); AIC(mod_deb_Ldec_Tmoyhiv) ; AIC(mod_deb_Ldec_dChill)
# La variable Tmoy40j donne les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Ldec_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
# On enlève les variables une par une : Tmoy40j reste le meilleur modèle
mod_deb_Ldec_multiT = lmer(julian_day ~ Tmoy40j + gel30j + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F) # AIC = 13370.29

# ggplot(deb_Ldec, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Ldec, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Ldec_Tmoy40j, "Tmoy40j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Ldec_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j  + year + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Ldec, REML=F)
mod_deb_Ldec_Tmoy40j_clAlt_x <- lmer(julian_day ~ Tmoy40j * altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Ldec, REML=F)
AIC(mod_deb_Ldec_Tmoy40j) ; AIC(mod_deb_Ldec_Tmoy40j_Alt) ; AIC(mod_deb_Ldec_Tmoy40j_year) ; AIC(mod_deb_Ldec_Tmoy40j_clAlt) ; AIC(mod_deb_Ldec_Tmoy40j_Alt_) ; AIC(mod_deb_Ldec_Tmoy40j_clAlt_) ; AIC(mod_deb_Ldec_Tmoy40j_clAlt_x)


#========= MEILLEUR MODELE :
bestmod_deb_Ldec = lmer(julian_day ~ Tmoy40j * gel30j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Ldec, REML=F) # AIC = 13271.05
bestmods = c(list("Meleze"=bestmod_deb_Ldec), bestmods)

summary(bestmod_deb_Ldec)
visreg(bestmod_deb_Ldec, "Tmoy40j", by="gel30j")
# gradient altitudinal : 2.3 jours de retard par 100m
# effet de la température moyenne dans les 40 jours précédant le débourrement : 1.6 jour d'avance quand on gagne 1°C MAIS effet moins fort s'il y a des jours de gel

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Ldec)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Ldec)
# - Validation croisée
n = nrow(deb_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Ldec[trainIndex ,]
test <- deb_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Ldec$ID_zone[drop=T])[!unique(deb_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Ldec), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.820
R2_models[R2_models$species == "Meleze", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Ldec))),
                                        coef = coef(summary(bestmod_deb_Ldec))[,1],
                                        std = coef(summary(bestmod_deb_Ldec))[,2],
                                        pval = coef(summary(bestmod_deb_Ldec))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Ldec)[rownames(coef(summary(bestmod_deb_Ldec)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Ldec_10ans = deb_Ldec[deb_Ldec$year <= 2016,]
bestmod_deb_Ldec_10ans <- lmer(julian_day ~ Tmoy40j_med_0616 * gel30j_med_0616 + altitude + (1|ID_zone) + (Tmoy40j_med_0616|yearQ),  deb_Ldec_10ans, REML=F)
summary(bestmod_deb_Ldec_10ans)
# On a les mêmes tendances sur les deux périodes 


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Ldec_10ans))),
                                        coef = coef(summary(bestmod_deb_Ldec_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Ldec_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Ldec_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Ldec_10ans)[rownames(coef(summary(bestmod_deb_Ldec_10ans)))], error = function(e) return(NA))))

##############################################-
#*---- Epicéa ----

deb_Pabi = debourr_Alps[debourr_Alps$species == "Epicea",]
# ggplot(deb_Pabi, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Pabi = deb_Pabi[!is.na(deb_Pabi$Tmoy30j),]
# ggplot(deb_Pabi, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Epicéa", x="année",y="jour julien")
deb_Pabi = deb_Pabi[deb_Pabi$cl_alt != "150-450",] # on retire la classe d'altitude la plus basse parce qu'il n'y a pas assez de données pour la suite

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Pabi_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Pabi, REML=F)
summary(mod_deb_Pabi_Alt)
altyearmods = c(list("Epicea"=mod_deb_Pabi_Alt), altyearmods)
deb_Pabi_10ans = deb_Pabi[deb_Pabi$year <= 2016,]
mod_deb_Pabi_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Pabi_10ans, REML=F)
summary(mod_deb_Pabi_all_10ans)
# => L'effet altitude est beaucoup plus faible que Bison et al 2019 ! =============================       /!\ /!\ /!\

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Epicea", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Pabi_Alt)
# + Validation croisée
n = nrow(deb_Pabi)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Pabi[trainIndex ,]
test <- deb_Pabi[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Pabi$ID_zone[drop=T])[!unique(deb_Pabi$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Pabi_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 
R2_models[R2_models$species == "Epicea", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Pabi_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)

AIC(mod_deb_Pabi_Tmoy30j) ; AIC(mod_deb_Pabi_Tmoy40j) ; AIC(mod_deb_Pabi_Tmoy50j) ; AIC(mod_deb_Pabi_gel30j) ; AIC(mod_deb_Pabi_gel40j) ; AIC(mod_deb_Pabi_gel50j) ; AIC(mod_deb_Pabi_GDD0) ; AIC(mod_deb_Pabi_GDD5); AIC(mod_deb_Pabi_Tmoyhiv) ; AIC(mod_deb_Pabi_dChill)
# La variable Tmoy50j donne les meilleurs résultats : c'est celle qu'on garde pour la suite
# /!\ Quand on ajoute d'autres variables, notamment l'altitude, Tmoyhiv devient meilleure !! ======================================== /!\ /!\ /!\

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Pabi_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
# On enlève les variables une par une.. jusqu'à avoir le meilleur AIC
mod_deb_Pabi_multiT = lmer(julian_day ~ Tmoy50j + gel30j + GDD5 + Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F) # AIC = 10217.52
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Pabi, aes(x=year, y=GDD0, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Pabi, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Pabi_GDD0, "GDD0")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Pabi_gel30j_Alt <- lmer(julian_day ~ gel30j + altitude + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel30j_clAlt <- lmer(julian_day ~ gel30j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel30j_Alt_ <- lmer(julian_day ~ gel30j + altitude + (1|ID_zone) + (gel30j|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_gel30j_clAlt_ <- lmer(julian_day ~ gel30j * cl_alt + (1|ID_zone) + (gel30j|yearQ), deb_Pabi, REML=F)
AIC(mod_deb_Pabi_gel30j) ; AIC(mod_deb_Pabi_gel30j_Alt) ; AIC(mod_deb_Pabi_gel30j_clAlt) ; AIC(mod_deb_Pabi_gel30j_Alt_) ; AIC(mod_deb_Pabi_gel30j_clAlt_)
mod_deb_Pabi_Tmoyhiv_Alt <- lmer(julian_day ~ Tmoyhiv + altitude + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv_year <- lmer(julian_day ~ Tmoyhiv + year + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv_clAlt <- lmer(julian_day ~ Tmoyhiv * cl_alt + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv_Alt_ <- lmer(julian_day ~ Tmoyhiv + altitude + (1|ID_zone) + (Tmoyhiv|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv_clAlt_ <- lmer(julian_day ~ Tmoyhiv * cl_alt + (1|ID_zone) + (Tmoyhiv|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoyhiv_clAlt_year <- lmer(julian_day ~ Tmoyhiv * cl_alt + year + (1|ID_zone) + (Tmoyhiv|yearQ), deb_Pabi, REML=F)
AIC(mod_deb_Pabi_Tmoyhiv) ; AIC(mod_deb_Pabi_Tmoyhiv_Alt) ; AIC(mod_deb_Pabi_Tmoyhiv_year) ;AIC(mod_deb_Pabi_Tmoyhiv_clAlt) ; AIC(mod_deb_Pabi_Tmoyhiv_Alt_) ; AIC(mod_deb_Pabi_Tmoyhiv_clAlt_) ; AIC(mod_deb_Pabi_Tmoyhiv_clAlt_year)
mod_deb_Pabi_Tmoy50j_Alt <- lmer(julian_day ~ Tmoy50j + altitude + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j_year <- lmer(julian_day ~ Tmoy50j + year + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j_clAlt <- lmer(julian_day ~ Tmoy50j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j_Alt_ <- lmer(julian_day ~ Tmoy50j + altitude + (1|ID_zone) + (Tmoy50j|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j_clAlt_ <- lmer(julian_day ~ Tmoy50j * cl_alt + (1|ID_zone) + (Tmoy50j|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50j_clAlt_year <- lmer(julian_day ~ Tmoy50j * cl_alt + year + (1|ID_zone) + (Tmoy50j|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50jTmoyhiv_clAlt_year <- lmer(julian_day ~ (Tmoy50j + Tmoyhiv) * cl_alt + year + (1|ID_zone) + (Tmoy50j|yearQ), deb_Pabi, REML=F)
mod_deb_Pabi_Tmoy50jTmoyhivgel30j_clAlt_year <- lmer(julian_day ~ (Tmoy50j + gel30j + Tmoyhiv) * cl_alt + year + (1|ID_zone) + (Tmoy50j|yearQ), deb_Pabi, REML=F)
AIC(mod_deb_Pabi_Tmoy50j) ; AIC(mod_deb_Pabi_Tmoy50j_Alt) ; AIC(mod_deb_Pabi_Tmoy50j_year) ;AIC(mod_deb_Pabi_Tmoy50j_clAlt) ; AIC(mod_deb_Pabi_Tmoy50j_Alt_) ; AIC(mod_deb_Pabi_Tmoy50j_clAlt_) ; AIC(mod_deb_Pabi_Tmoy50j_clAlt_year) ; AIC(mod_deb_Pabi_Tmoy50jTmoyhiv_clAlt_year) ; AIC(mod_deb_Pabi_Tmoy50jTmoyhivgel30j_clAlt_year)
# Mieux avec l'altitude en plus, en prenant Tmoyhiv + Tmoy50j + gel30j, et en ayant l'effet aléatoire 'année' sur la température et une interaction (AIC = 10194.72)

#========= MEILLEUR MODELE :
bestmod_deb_Pabi = lmer(julian_day ~ Tmoy30j * altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Pabi, REML=F) # AIC = 10189.21
bestmods = c(list("Epicea"=bestmod_deb_Pabi), bestmods)

summary(bestmod_deb_Pabi)
visreg(bestmod_deb_Pabi, "Tmoy30j", by="altitude")


# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Epicea", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Pabi)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Pabi)
# - Validation croisée
n = nrow(deb_Pabi)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Pabi[trainIndex ,]
test <- deb_Pabi[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Pabi$ID_zone[drop=T])[!unique(deb_Pabi$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Pabi), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7706
R2_models[R2_models$species == "Epicea", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Epicea",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Pabi))),
                                        coef = coef(summary(bestmod_deb_Pabi))[,1],
                                        std = coef(summary(bestmod_deb_Pabi))[,2],
                                        pval = coef(summary(bestmod_deb_Pabi))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Pabi)[rownames(coef(summary(bestmod_deb_Pabi)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Pabi_10ans = deb_Pabi[deb_Pabi$year <= 2016,]
bestmod_deb_Pabi_10ans <- lmer(julian_day ~ Tmoy30j_med_0616 * altitude + (1|ID_zone) + (Tmoy30j_med_0616|yearQ), deb_Pabi_10ans, REML=F)
summary(bestmod_deb_Pabi_10ans)
# On a les mêmes tendances sur les deux périodes (même en prenant la température calculée avec la médiane des débourrements 2006-2016)


resultats = rbind(resultats, data.frame(species = "Epicea",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Pabi_10ans))),
                                        coef = coef(summary(bestmod_deb_Pabi_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Pabi_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Pabi_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Pabi_10ans)[rownames(coef(summary(bestmod_deb_Pabi_10ans)))], error = function(e) return(NA))))



##############################################-
#*---- Sorbier ----

deb_Sacu = debourr_Alps[debourr_Alps$species == "Sorbier",]
# ggplot(deb_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Sacu = deb_Sacu[!is.na(deb_Sacu$Tmoy30j),]
deb_Sacu = deb_Sacu[!is.na(deb_Sacu$altitude),]
# ggplot(deb_Sacu, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Sorbier", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Sacu_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Sacu, REML=F)
summary(mod_deb_Sacu_Alt)
altyearmods = c(list("Sorbier"=mod_deb_Sacu_Alt), altyearmods)
deb_Sacu_10ans = deb_Sacu[deb_Sacu$year <= 2016,]
mod_deb_Sacu_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Sacu_10ans, REML=F)
summary(mod_deb_Sacu_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019 (mais effet année plus faible)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Sacu_Alt)
# + Validation croisée
n = nrow(deb_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Sacu[trainIndex ,]
test <- deb_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Sacu$ID_zone[drop=T])[!unique(deb_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Sacu_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.613
R2_models[R2_models$species == "Sorbier", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Sacu_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)

AIC(mod_deb_Sacu_Tmoy30j) ; AIC(mod_deb_Sacu_Tmoy40j) ; AIC(mod_deb_Sacu_Tmoy50j) ; AIC(mod_deb_Sacu_gel30j) ; AIC(mod_deb_Sacu_gel40j) ; AIC(mod_deb_Sacu_gel50j) ; AIC(mod_deb_Sacu_GDD0) ; AIC(mod_deb_Sacu_GDD5); AIC(mod_deb_Sacu_Tmoyhiv) ; AIC(mod_deb_Sacu_dChill)
# Les variables Tmoy40j et GDD0 donnent les meilleurs résultats 

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Sacu_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
# On enlève les variables une par une : Il reste 4 variables dans le modèle avec le meilleur AIC... MAIS les variables restantes sont quand même 
# fortement corrélées, et pas toutes significatives, donc on en retire encore 2
mod_deb_Sacu_multiT = lmer(julian_day ~ Tmoy40j + gel50j + dChill + GDD0 + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F) # AIC = 5728.5
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Sacu, aes(x=year, y=Tmoy30j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Sacu, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Sacu_Tmoy30j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Sacu_GDD0_Alt <- lmer(julian_day ~ GDD0 + altitude + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0dChill_Alt <- lmer(julian_day ~ GDD0 + dChill + altitude + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0_year <- lmer(julian_day ~ GDD0 + altitude + year + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0_clAlt <- lmer(julian_day ~ GDD0 * cl_alt + (1|ID_zone) + (1|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0_Alt_ <- lmer(julian_day ~ GDD0 + altitude + (1|ID_zone) + (GDD0|yearQ), deb_Sacu, REML=F)
mod_deb_Sacu_GDD0_clAlt_ <- lmer(julian_day ~ GDD0 * cl_alt + (1|ID_zone) + (GDD0|yearQ), deb_Sacu, REML=F)
AIC(mod_deb_Sacu_GDD0) ; AIC(mod_deb_Sacu_GDD0_Alt) ; AIC(mod_deb_Sacu_GDD0dChill_Alt) ; AIC(mod_deb_Sacu_GDD0_year) ; AIC(mod_deb_Sacu_GDD0_clAlt) ; AIC(mod_deb_Sacu_GDD0_Alt_) ; AIC(mod_deb_Sacu_GDD0_clAlt_)
# Mieux avec l'altitude en plus (AIC = 5701.37)

#========= MEILLEUR MODELE :
bestmod_deb_Sacu = lmer(julian_day ~ GDD0 + dChill + altitude + (1|ID_zone) + (dChill|yearQ), deb_Sacu, REML=F) # AIC = 5712.1
bestmods = c(list("Sorbier"=bestmod_deb_Sacu), bestmods)

summary(bestmod_deb_Sacu)
# gradient altitudinal : 1.9 jours de retard quand on monte de 100m
# effet des GDD0 cumulés avant la date médiane de débourrement : 0.04 jours d'avance quand on gagne 1°C.J

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Sacu)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Sacu)
# - Validation croisée
n = nrow(deb_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Sacu[trainIndex ,]
test <- deb_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Sacu$ID_zone[drop=T])[!unique(deb_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Sacu), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6873
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Sacu))),
                                        coef = coef(summary(bestmod_deb_Sacu))[,1],
                                        std = coef(summary(bestmod_deb_Sacu))[,2],
                                        pval = coef(summary(bestmod_deb_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Sacu)[rownames(coef(summary(bestmod_deb_Sacu)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Sacu_10ans = deb_Sacu[deb_Sacu$year <= 2016,]
bestmod_deb_Sacu_10ans <- lmer(julian_day ~ GDD0_med_0616 + dChill_med_0616 + altitude + (1|ID_zone) + (dChill_med_0616|yearQ), deb_Sacu_10ans, REML=F)
summary(bestmod_deb_Sacu_10ans)
# On a les mêmes tendances sur les deux périodes mais le gradient altitudinale est plus faible


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Sacu_10ans))),
                                        coef = coef(summary(bestmod_deb_Sacu_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Sacu_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Sacu_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Sacu_10ans)[rownames(coef(summary(bestmod_deb_Sacu_10ans)))], error = function(e) return(NA))))

##############################################-
#*---- Lilas ----

deb_Svul = debourr_Alps[debourr_Alps$species == "Lilas",]
# ggplot(deb_Svul, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Svul = deb_Svul[!is.na(deb_Svul$Tmoy30j),]
# ggplot(deb_Svul, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Lilas", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Svul_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Svul, REML=F)
summary(mod_deb_Svul_Alt)
altyearmods = c(list("Lilas"=mod_deb_Svul_Alt), altyearmods)
# deb_Svul_10ans = deb_Svul[deb_Svul$year <= 2016,]
# mod_deb_Svul_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Svul_10ans, REML=F)
# summary(mod_deb_Svul_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019 (mais résultat pas indiqué pour la première période !)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Lilas", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Svul_Alt)
# + Validation croisée
n = nrow(deb_Svul)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Svul[trainIndex ,]
test <- deb_Svul[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Svul$ID_zone[drop=T])[!unique(deb_Svul$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Svul_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.5828
R2_models[R2_models$species == "Lilas", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Svul_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)

AIC(mod_deb_Svul_Tmoy30j) ; AIC(mod_deb_Svul_Tmoy40j) ; AIC(mod_deb_Svul_Tmoy50j) ; AIC(mod_deb_Svul_gel30j) ; AIC(mod_deb_Svul_gel40j) ; AIC(mod_deb_Svul_gel50j) ; AIC(mod_deb_Svul_GDD0) ; AIC(mod_deb_Svul_GDD5); AIC(mod_deb_Svul_Tmoyhiv) ; AIC(mod_deb_Svul_dChill)
# Les variables Tmoy30j et gel30j donnent les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
mod_deb_Svul_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
# On enlève les variables une par une : Tmoy30j puis dChill puis Tmoyhiv puis Tmoy50j--> c'est le modèle avec le meilleur AIC... MAIS les variables
# restantes sont quand même fortement corrélées...!
mod_deb_Svul_multiT = lmer(julian_day ~ Tmoy30j + gel30j + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F) # AIC = 4273.4
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Svul, aes(x=year, y=Tmoy30j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Svul, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Svul_Tmoy30j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
mod_deb_Svul_Tmoy30j_Alt <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy30j_year <- lmer(julian_day ~ Tmoy30j + year + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy30j_clAlt <- lmer(julian_day ~ Tmoy30j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy30jgel30j_Alt_ <- lmer(julian_day ~ Tmoy30j + gel30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Svul, REML=F)
mod_deb_Svul_gel30j_Alt_ <- lmer(julian_day ~ gel30j + altitude + (1|ID_zone) + (gel30j|yearQ), deb_Svul, REML=F)
mod_deb_Svul_Tmoy30j_clAlt_ <- lmer(julian_day ~ Tmoy30j * cl_alt + (1|ID_zone) + (Tmoy30j|yearQ), deb_Svul, REML=F)
AIC(mod_deb_Svul_Tmoy30j) ; AIC(mod_deb_Svul_Tmoy30j_Alt) ; AIC(mod_deb_Svul_Tmoy30j_year) ; AIC(mod_deb_Svul_Tmoy30j_clAlt) ; AIC(mod_deb_Svul_Tmoy30j_Alt_) ; AIC(mod_deb_Svul_Tmoy30jgel30j_Alt_) ; AIC(mod_deb_Svul_gel30j_Alt_) ; AIC(mod_deb_Svul_Tmoy30j_clAlt_)
# Mieux avec l'ajout altitude (AIC = 4241.98)

#========= MEILLEUR MODELE :
bestmod_deb_Svul = lmer(julian_day ~ altitude + gel30j + (1|ID_zone) + (gel30j|yearQ), deb_Svul, REML=F) # AIC = 4259.7
bestmods = c(list("Lilas"=bestmod_deb_Svul), bestmods)

summary(bestmod_deb_Svul)
visreg(bestmod_deb_Svul, "gel30j", by="altitude")
# gradient altitudinal : 2.0 jours de retard quand on monte de 100m
# effet du nombre de jours de gel 30 j avant la date médiane de débourrement : 0.5 jours de retard par jour de gel

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Lilas", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Svul)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Svul)
# - Validation croisée
n = nrow(deb_Svul)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Svul[trainIndex ,]
test <- deb_Svul[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Svul$ID_zone[drop=T])[!unique(deb_Svul$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Svul), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.715
R2_models[R2_models$species == "Lilas", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Lilas",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Svul))),
                                        coef = coef(summary(bestmod_deb_Svul))[,1],
                                        std = coef(summary(bestmod_deb_Svul))[,2],
                                        pval = coef(summary(bestmod_deb_Svul))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Svul)[rownames(coef(summary(bestmod_deb_Svul)))], error = function(e) return(NA))))


# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Svul_10ans = deb_Svul[deb_Svul$year <= 2016,]
bestmod_deb_Svul_10ans <- lmer(julian_day ~ altitude + gel30j_med_0616+ (1|ID_zone) + (gel30j_med_0616|yearQ), deb_Svul_10ans, REML=F)
summary(bestmod_deb_Svul_10ans)
# On a les mêmes tendances sur les deux périodes mais le gradient de température est un peu plus faible


resultats = rbind(resultats, data.frame(species = "Lilas",
                                        periode = "2006-2016",
                                        variable = rownames(coef(summary(bestmod_deb_Svul_10ans))),
                                        coef = coef(summary(bestmod_deb_Svul_10ans))[,1],
                                        std = coef(summary(bestmod_deb_Svul_10ans))[,2],
                                        pval = coef(summary(bestmod_deb_Svul_10ans))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Svul_10ans)[rownames(coef(summary(bestmod_deb_Svul_10ans)))], error = function(e) return(NA))))






##############################################-
#*---- Hetre ----

# /!\ Espèce ajoutée en 2018 !! 
#     => Pas incluse dans l'étude de Bison et al. 2019
#     => Ne pas l'inclure dans le calcul des indices

deb_Fsyl = debourr_Alps[debourr_Alps$species == "Hetre",]
# ggplot(deb_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Fsyl = deb_Fsyl[!is.na(deb_Fsyl$Tmoy30j),]
# ggplot(deb_Fsyl, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Fsyl_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Fsyl, REML=F)
summary(mod_deb_Fsyl_Alt)
altyearmods = c(list("Hetre"=mod_deb_Fsyl_Alt), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Fsyl_Alt)
# + Validation croisée
n = nrow(deb_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Fsyl[trainIndex ,]
test <- deb_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Fsyl$ID_zone[drop=T])[!unique(deb_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Fsyl_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.4343799
R2_models[R2_models$species == "Hetre", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Fsyl_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
mod_deb_Fsyl_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)

AIC(mod_deb_Fsyl_Tmoy30j) ; AIC(mod_deb_Fsyl_Tmoy40j) ; AIC(mod_deb_Fsyl_Tmoy50j) ; AIC(mod_deb_Fsyl_gel30j) ; AIC(mod_deb_Fsyl_gel40j) ; AIC(mod_deb_Fsyl_gel50j) ; AIC(mod_deb_Fsyl_GDD0) ; AIC(mod_deb_Fsyl_GDD5); AIC(mod_deb_Fsyl_Tmoyhiv) ; AIC(mod_deb_Fsyl_dChill)
# La variable Tmoyhiv donne les meilleurs résultats 

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
# /!\ ce qui serait le plus logique serait de ne garder que la meilleure des 5 variables de Tmoy et GDD (cf corrélations) et éventuellement Tmoyhiv et/ou dChill
mod_deb_Fsyl_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
# On enlève les variables une par une : dChill puis Tmoy30j puis Tmoyhiv --> c'est le modèle avec le meilleur AIC... MAIS GDD5 n'a pas un effet 
# significatif au seuil 0.05... et les variables de température Tmoy+GDD sont super corrélées ! 
# En combinant ces 3 critères (AIC, significativité des effets, variables peu corrélées), on retient donc :
mod_deb_Fsyl_multiT = lmer(julian_day ~ Tmoy50j + gel40j +GDD5 + Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F) # AIC = 543.5
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Fsyl, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Fsyl, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Fsyl_Tmoy40j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)

# mod_deb_Fsyl_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
# mod_deb_Fsyl_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j + year + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
# mod_deb_Fsyl_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F)
# mod_deb_Fsyl_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Fsyl, REML=F)
# mod_deb_Fsyl_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Fsyl, REML=F)
# mod_deb_Fsyl_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Fsyl, REML=F)
# AIC(mod_deb_Fsyl_Tmoy40j) ; AIC(mod_deb_Fsyl_Tmoy40j_Alt) ; AIC(mod_deb_Fsyl_Tmoy40j_year) ; AIC(mod_deb_Fsyl_Tmoy40j_clAlt) ; AIC(mod_deb_Fsyl_Tmoy30j_Alt_) ; AIC(mod_deb_Fsyl_Tmoy40j_Alt_) ; AIC(mod_deb_Fsyl_Tmoy40j_clAlt_)
# # Mieux avec l'altitude en plus, et en ayant l'effet aléatoire 'année' sur la température 

#========= MEILLEUR MODELE (en limitant le nombre de variables explicatives) :
bestmod_deb_Fsyl = lmer(julian_day ~ Tmoyhiv + Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Fsyl, REML=F) # AIC = 545.3
bestmods = c(list("Hetre"=bestmod_deb_Fsyl), bestmods)


summary(bestmod_deb_Fsyl)
# effet de la température moyenne hivernale : 4.9 jour d'avance quand il a fait plus chaud l'hiver
# effet de la température moyenne dans les 30 jours précédant le débourrement : 2.9 jour d'avance quand on gagne 1°C

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Fsyl)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Fsyl)
# - Validation croisée
n = nrow(deb_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Fsyl[trainIndex ,]
test <- deb_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Fsyl$ID_zone[drop=T])[!unique(deb_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Fsyl), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Fsyl))),
                                        coef = coef(summary(bestmod_deb_Fsyl))[,1],
                                        std = coef(summary(bestmod_deb_Fsyl))[,2],
                                        pval = coef(summary(bestmod_deb_Fsyl))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Fsyl)[rownames(coef(summary(bestmod_deb_Fsyl)))], error = function(e) return(NA))))


##############################################-
#*---- Pin ----

# /!\ Espèce ajoutée en 2018 !! 
#     => Pas incluse dans l'étude de Bison et al. 2019
#     => Ne pas l'inclure dans le calcul des indices

deb_Psyl = debourr_Alps[debourr_Alps$species == "Pin",]
# ggplot(deb_Psyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Psyl = deb_Psyl[!is.na(deb_Psyl$Tmoy30j),]
# ggplot(deb_Psyl, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Bouleau", x="année",y="jour julien")
deb_Psyl = deb_Psyl[deb_Psyl$julian_day>50,] # On supprime une valeur bizarre de débourrement à 47 jours (plus d'un mois avant toutes les autres)

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Psyl_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Psyl, REML=F)
summary(mod_deb_Psyl_Alt)
altyearmods = c(list("Pin"=mod_deb_Psyl_Alt), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Pin", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Psyl_Alt)
# + Validation croisée
n = nrow(deb_Psyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Psyl[trainIndex ,]
test <- deb_Psyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Psyl$ID_zone[drop=T])[!unique(deb_Psyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Psyl_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6469443
R2_models[R2_models$species == "Pin", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Psyl_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
mod_deb_Psyl_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)

AIC(mod_deb_Psyl_Tmoy30j) ; AIC(mod_deb_Psyl_Tmoy40j) ; AIC(mod_deb_Psyl_Tmoy50j) ; AIC(mod_deb_Psyl_gel30j) ; AIC(mod_deb_Psyl_gel40j) ; AIC(mod_deb_Psyl_gel50j) ; AIC(mod_deb_Psyl_GDD0) ; AIC(mod_deb_Psyl_GDD5); AIC(mod_deb_Psyl_Tmoyhiv) ; AIC(mod_deb_Psyl_dChill)
# La variable dChill donne les meilleurs résultats 

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
# /!\ ce qui serait le plus logique serait de ne garder que la meilleure des 5 variables de Tmoy et GDD (cf corrélations) et éventuellement Tmoyhiv et/ou dChill
mod_deb_Psyl_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
# On enlève les variables une par une, en regardant leur significativité, leurs corrélations + l'AIC
# En combinant ces 3 critères (AIC, significativité des effets, variables peu corrélées), on retient donc :
mod_deb_Psyl_multiT = lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F) # AIC = 331.3
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Psyl, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Psyl, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Psyl_Tmoy40j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)

# mod_deb_Psyl_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
# mod_deb_Psyl_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j + year + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
# mod_deb_Psyl_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F)
# mod_deb_Psyl_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Psyl, REML=F)
# mod_deb_Psyl_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Psyl, REML=F)
# mod_deb_Psyl_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Psyl, REML=F)
# AIC(mod_deb_Psyl_Tmoy40j) ; AIC(mod_deb_Psyl_Tmoy40j_Alt) ; AIC(mod_deb_Psyl_Tmoy40j_year) ; AIC(mod_deb_Psyl_Tmoy40j_clAlt) ; AIC(mod_deb_Psyl_Tmoy30j_Alt_) ; AIC(mod_deb_Psyl_Tmoy40j_Alt_) ; AIC(mod_deb_Psyl_Tmoy40j_clAlt_)
# # Mieux avec l'altitude en plus, et en ayant l'effet aléatoire 'année' sur la température 

#========= MEILLEUR MODELE (en limitant le nombre de variables explicatives) :
bestmod_deb_Psyl = lmer(julian_day ~ dChill + GDD5 + (1|ID_zone) + (1|yearQ), deb_Psyl, REML=F) # AIC = 327.8
bestmods = c(list("Pin"=bestmod_deb_Psyl), bestmods)


summary(bestmod_deb_Psyl)
# effet du nombre de la proportion de jours de chilling (entre 0 et 8°C) : 8.6 jours d'avance pour 10% de jours de chilling en plus
# effet de l'accumulation de chaleur (GDD5) : 0.1 jour d'avance quand on gagne 1°C

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Pin", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Psyl)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Psyl)
# - Validation croisée
n = nrow(deb_Psyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Psyl[trainIndex ,]
test <- deb_Psyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Psyl$ID_zone[drop=T])[!unique(deb_Psyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Psyl), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7922845
R2_models[R2_models$species == "Pin", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Pin",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Psyl))),
                                        coef = coef(summary(bestmod_deb_Psyl))[,1],
                                        std = coef(summary(bestmod_deb_Psyl))[,2],
                                        pval = coef(summary(bestmod_deb_Psyl))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Psyl)[rownames(coef(summary(bestmod_deb_Psyl)))], error = function(e) return(NA))))


##############################################-
#*---- Sapin ----

# /!\ Espèce ajoutée en 2018 !! 
#     => Pas incluse dans l'étude de Bison et al. 2019
#     => Ne pas l'inclure dans le calcul des indices

# /!\ Il n'y a que 12 obs au total pour cette espèce !!! => pas de sens de caler des modèles avec 12 points seulement (on peut visualiser les tendances qd même)

deb_Aalb = debourr_Alps[debourr_Alps$species == "Sapin",]
deb_Aalb = deb_Aalb[!is.na(deb_Aalb$Tmoy30j),]
ggplot(deb_Aalb, aes(x=altitude, y=julian_day)) + geom_point() + labs(title = "Débourrement 10% - Sapin", x="altitude",y="jour julien")
ggplot(deb_Aalb, aes(x=Tmoy30j, y=julian_day)) + geom_point() + labs(title = "Débourrement 10% - Sapin", x="Tmoy30j",y="jour julien")
ggplot(deb_Aalb, aes(x=Tmoyhiv, y=julian_day)) + geom_point() + labs(title = "Débourrement 10% - Sapin", x="Tmoyhiv",y="jour julien")


##############################################-
#*---- Bouleau pubescent ----

deb_Bpub = debourr_Alps[debourr_Alps$species == "Bouleau_pubescent",]
# ggplot(deb_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Bpub = deb_Bpub[!is.na(deb_Bpub$Tmoy30j),]
# ggplot(deb_Bpub, aes(x=julian_day)) + geom_histogram()
# ggplot(deb_Bpub, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019 ---- 
mod_deb_Bpub_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|ID_zone), deb_Bpub, REML=F)
summary(mod_deb_Bpub_Alt)
altyearmods = c(list("Bouleau_pubescent"=mod_deb_Bpub_Alt), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_deb_Bpub_Alt)
# + Validation croisée
n = nrow(deb_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Bpub[trainIndex ,]
test <- deb_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Bpub$ID_zone[drop=T])[!unique(deb_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_deb_Bpub_Alt), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6222258
R2_models[R2_models$species == "Bouleau_pubescent", "R2_altyear_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Bpub_Tmoy30j <- lmer(julian_day ~ Tmoy30j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_Tmoy40j <- lmer(julian_day ~ Tmoy40j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_Tmoy50j <- lmer(julian_day ~ Tmoy50j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_gel30j <- lmer(julian_day ~ gel30j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_gel40j <- lmer(julian_day ~ gel40j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_gel50j <- lmer(julian_day ~ gel50j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_GDD0 <- lmer(julian_day ~ GDD0 + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_GDD5 <- lmer(julian_day ~ GDD5 + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_Tmoyhiv <- lmer(julian_day ~ Tmoyhiv + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
mod_deb_Bpub_dChill <- lmer(julian_day ~ dChill + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)

AIC(mod_deb_Bpub_Tmoy30j) ; AIC(mod_deb_Bpub_Tmoy40j) ; AIC(mod_deb_Bpub_Tmoy50j) ; AIC(mod_deb_Bpub_gel30j) ; AIC(mod_deb_Bpub_gel40j) ; AIC(mod_deb_Bpub_gel50j) ; AIC(mod_deb_Bpub_GDD0) ; AIC(mod_deb_Bpub_GDD5); AIC(mod_deb_Bpub_Tmoyhiv) ; AIC(mod_deb_Bpub_dChill)
# La variable Tmoy40j donne les meilleurs résultats : c'est celle qu'on garde pour la suite

# NB : on peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
# /!\ ce qui serait le plus logique serait de ne garder que la meilleure des 5 variables de Tmoy et GDD (cf corrélations) et éventuellement Tmoyhiv et/ou dChill
mod_deb_Bpub_multiT = lmer(julian_day ~ Tmoy30j + Tmoy40j + Tmoy50j + gel30j + gel40j + gel50j + GDD0 + GDD5 + Tmoyhiv + dChill + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
# On enlève les variables une par une
# En combinant ces 3 critères (AIC, significativité des effets, variables peu corrélées), on retient donc :
mod_deb_Bpub_multiT = lmer(julian_day ~ Tmoy30j + gel30j + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F) # AIC = 637.8
# On garde en tête ces multivariables de température pour la suite

# ggplot(deb_Bpub, aes(x=year, y=Tmoy40j, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Bpub, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm")
# visreg(mod_deb_Bpub_Tmoy40j, "Tmoy30j")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
# mod_deb_Bpub_Tmoy40j_Alt <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
# mod_deb_Bpub_Tmoy40j_year <- lmer(julian_day ~ Tmoy40j + year + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
# mod_deb_Bpub_Tmoy40j_clAlt <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F)
# mod_deb_Bpub_Tmoy30j_Alt_ <- lmer(julian_day ~ Tmoy30j + altitude + (1|ID_zone) + (Tmoy30j|yearQ), deb_Bpub, REML=F)
# mod_deb_Bpub_Tmoy40j_Alt_ <- lmer(julian_day ~ Tmoy40j + altitude + (1|ID_zone) + (Tmoy40j|yearQ), deb_Bpub, REML=F)
# mod_deb_Bpub_Tmoy40j_clAlt_ <- lmer(julian_day ~ Tmoy40j * cl_alt + (1|ID_zone) + (Tmoy40j|yearQ), deb_Bpub, REML=F)
# AIC(mod_deb_Bpub_Tmoy40j) ; AIC(mod_deb_Bpub_Tmoy40j_Alt) ; AIC(mod_deb_Bpub_Tmoy40j_year) ; AIC(mod_deb_Bpub_Tmoy40j_clAlt) ; AIC(mod_deb_Bpub_Tmoy30j_Alt_) ; AIC(mod_deb_Bpub_Tmoy40j_Alt_) ; AIC(mod_deb_Bpub_Tmoy40j_clAlt_)
# # Mieux avec l'altitude en plus, et en ayant l'effet aléatoire 'année' sur la température (AIC = 12453.3)

#========= MEILLEUR MODELE :
bestmod_deb_Bpub = lmer(julian_day ~ Tmoy30j + gel30j + altitude + (1|ID_zone) + (1|yearQ), deb_Bpub, REML=F) # AIC = 631.0
bestmods = c(list("Bouleau_pubescent"=bestmod_deb_Bpub), bestmods)


summary(bestmod_deb_Bpub)
# gradient altitudinal : 2.3 jour de retard quand on monte de 100m
# effet de la température moyenne dans les 30 jours précédant le débourrement : 4.4 jours d'avance quand on gagne 1°C
# effet du gel dans les 30j précédant le débourrement : 1 jour d'avance quand il y a 1 jour de gel en plus

# EVALUATION DU MODELE
# - Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_deb_Bpub)
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_deb_Bpub)
# - Validation croisée
n = nrow(deb_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- deb_Bpub[trainIndex ,]
test <- deb_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(deb_Bpub$ID_zone[drop=T])[!unique(deb_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
test1 <- test1[!test1$yearQ%in%(unique(deb_Bpub$yearQ[drop=T])[!unique(deb_Bpub$yearQ[drop=T])%in%unique(train$yearQ[drop=T])]),]
mod_train <- lmer(formula(bestmod_deb_Bpub), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.6386649
R2_models[R2_models$species == "Bouleau_pubescent", "R2_bestmod_calibval"] = summary(lm(predictions~julian_day, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        periode = "2006-2024",
                                        variable = rownames(coef(summary(bestmod_deb_Bpub))),
                                        coef = coef(summary(bestmod_deb_Bpub))[,1],
                                        std = coef(summary(bestmod_deb_Bpub))[,2],
                                        pval = coef(summary(bestmod_deb_Bpub))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_deb_Bpub)[rownames(coef(summary(bestmod_deb_Bpub)))], error = function(e) return(NA))))








#*---- Enregistrement des résultats ----

write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_deb__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_deb__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_deb.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_deb_v1altyear.Rdata")



#*---- Diagnostic / aperçu des résidus et diagnostic des effets aléatoires ----
DIAG = data.frame(species = unique(debourr_Alps$species), 
                  sd_res_altyear=NA, sd_pred_altyear=NA, sd_alea_ID_zone_altyear=NA, 
                  sd_res_bestmod=NA, sd_pred_bestmod=NA, sd_alea_ID_zone_bestmod=NA, sd_alea_yearQ_bestmod=NA)
for (esp in names(bestmods)){
  DIAG[DIAG$species == esp,c("sd_res_bestmod","sd_pred_bestmod","sd_alea_ID_zone_bestmod","sd_alea_yearQ_bestmod")] = 
    c(sd(bestmods[[esp]]@resp$wtres), sd(bestmods[[esp]]@resp$mu), attr(summary(bestmods[[esp]])$varcor$ID_zone,"stddev")[1], attr(summary(bestmods[[esp]])$varcor$yearQ,"stddev")[1])
  DIAG[DIAG$species == esp,c("sd_res_altyear","sd_pred_altyear","sd_alea_ID_zone_altyear")] = 
    c(sd(altyearmods[[esp]]@resp$wtres), sd(altyearmods[[esp]]@resp$mu), attr(summary(altyearmods[[esp]])$varcor$ID_zone,"stddev")[1], attr(summary(altyearmods[[esp]])$varcor$yearQ,"stddev")[1])

}

############################################################################################-
# INDICE PHÉNOCLIM BASÉ SUR CES MODÈLES                                                  ----
############################################################################################-

##############################################-
#*----- Par département ----

# L'idée serait de calculer pour chaque année x espèce x dept x altitude (ou classe d'altitude) la date de débourrement, en fonction de la 
# température moyenne mesurée une année donnée dans la classe d'altitude x dept.
# Pour calculer ces variables de température, on se base sur les reconstructions de température.

Tperiodes = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/varTprintemps_allphenosites.csv")
# phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/data/_CLEANED_data_pheno.csv")

Tperiodes = merge(Tperiodes, 
                  phenoclim[!duplicated(phenoclim$id_base_site),c("id_base_site","altitude","cl_alt","cl_alt2","cl_alt3",
                                                                  "coord_x_2154", "coord_y_2154",
                                                                  "dept1","region1","pays1","nom_massif","nom_massif_v2019")],
                  by.x="sitePheno", by.y="id_base_site", all.x=T)

Tperiodes$dept2 = Tperiodes$dept1
Tperiodes[is.element(Tperiodes$dept1, "District d'Aigle"),]$dept2 <- "Vaud"
Tperiodes[is.element(Tperiodes$dept1, "District de la Riviera-Pays-d’Enhaut"),]$dept2 <- "Vaud"
Tperiodes[is.element(Tperiodes$dept1, "District de la Veveyse"),]$dept2 <- "Fribourg"
Tperiodes[is.element(Tperiodes$dept1, "District de Lausanne"),]$dept2 <- "Vaud"
Tperiodes[is.element(Tperiodes$dept1, "District de Lavaux-Oron"),]$dept2 <- "Vaud"
Tperiodes[is.element(Tperiodes$dept1, "Hérens"),]$dept2 <- "Valais"
Tperiodes[is.element(Tperiodes$dept1, "Entremont"),]$dept2 <- "Valais"
Tperiodes[is.element(Tperiodes$dept1, "Wahlkreis Sarganserland"),]$dept2 <- "Saint Gall"
Tperiodes[is.element(Tperiodes$region1, "Valle d'Aosta / Vallée d'Aoste"),]$dept2 <- "Aoste"
#         --------------------------------------------------------------------------------------
#           /!\ À voir si on groupe aussi Vercelli avec un département italien, et Wahlkreis See-Gaster avec Saint Gall ou autre
#         --------------------------------------------------------------------------------------

# Aperçu des zones (département x classe d'altitude) où on a des données
table(Tperiodes$dept2, Tperiodes$cl_alt)



T_pourpred = Tperiodes %>% group_by(dept2, annee, esp, cl_alt,
                                          cl_alt2, cl_alt3, region1, pays1, nom_massif) %>% summarise(across(ends_with("0624"), mean))
T_pourpred$altitude = as.numeric(as.character(factor(T_pourpred$cl_alt,
                                                  levels = c("150-450", "450-750", "750-1050", "1050-1350", "1350-1650", "1650-1950", "1950-2250"),
                                                  labels = c(300, 600, 900, 1200, 1500, 1800, 2100))))


# On peut donc utiliser les modèles de chaque espèce pour faire les prédictions de débourrement
#         --------------------------------------------------------------------------------------
#           /!\ je n'ai pas fait de modèle pour le tussilage, la primevère et le pin sylvestre !! À voir s'ils sont dans l'indice calculé par Marjo & co
#         --------------------------------------------------------------------------------------

indice_deb = T_pourpred %>% rename(Tmoy30j = Tmoy30j_med_0624,
                                   Tmoy40j = Tmoy40j_med_0624,
                                   Tmoy50j = Tmoy50j_med_0624,
                                   gel30j = gel30j_med_0624,
                                   gel40j = gel40j_med_0624,
                                   gel50j = gel50j_med_0624,
                                   GDD0 = GDD0_med_0624,
                                   GDD5 = GDD5_med_0624,
                                   Tmoyhiv = Tmoyhiv_med_0624,
                                   dChill = dChill_med_0624,
                                   year = annee)
indice_deb$cl_alt = factor(indice_deb$cl_alt, levels=c("150-450","450-750" ,  "750-1050"   ,"1050-1350" , "1350-1650" ,"1650-1950", "1950-2250" ), 
                            ordered = T)
indice_deb$yearQ = factor(indice_deb$year)

# Chargement des best models
load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_deb.Rdata") 
# Prédictions
# RQ : on peut choisir de prédire en prenant en compte les effets aléatoires ou sans. Ici les effets aléatoires sont ID_zone (effet local, qu'on
#      cherche à lisser en agrégeant par département) et yearQ (effet annuel lié aux spécificités climatiques non prises en compte dans le modèle,
#      ce qu'on peut vouloir garder). Dans les prédictions V1, l'effet ID_zone n'est pas pris en compte, et l'effet year est pris en compte comme 
#      effet fixe dans le modèle. Dans la V2, on choisit de faire avec et sans l'effet aléatoire yearQ (en factor).
for (esp in c("Bouleau_verruqueux","Noisetier","Frene","Meleze","Epicea","Sorbier","Lilas")){#names(bestmods)){
    print(esp)
  if ("cl_alt" %in% colnames(bestmods[[esp]]@frame)){
    # Prédiction
    indice_deb$pred_deb_fixef[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp,], re.form=NA)
    indice_deb$pred_deb_allef[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp,], re.form=~(1|yearQ))
  } else {
    # Prédiction
    indice_deb$pred_deb_fixef[indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$esp == esp,], re.form=NA)
    indice_deb$pred_deb_allef[indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$esp == esp,], re.form=~(1|yearQ))
  }
}




indice_deb_dept = indice_deb %>% group_by(year, dept2, cl_alt, cl_alt2) %>% summarise(pred_deb_fixef = mean(pred_deb_fixef, na.rm=T),
                                                                                      pred_deb_allef = mean(pred_deb_allef, na.rm=T))

write.csv(indice_deb_dept, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_dept.csv", row.names = F)


# ggplot(indice_deb_dept, aes(y=year + 0.1*as.numeric(as.character(factor(cl_alt, levels=levels(cl_alt), labels=1:7))), 
#                     x=pred_deb_fixef,
#                     col=cl_alt)) + 
#   # geom_rect(aes(ymax = year + 0.5, 
#   #               ymin = year - 0.5, 
#   #               xmin = -Inf, 
#   #               xmax = Inf, color=NULL, 
#   #               fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
#   geom_point(shape=15, size=3) + 
#   geom_segment(aes(x=pred_deb_fixef, xend=+Inf, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_alt, levels=levels(cl_alt), labels=1:7))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_alt, levels=levels(cl_alt), labels=1:7))),
#                    col=cl_alt ), lty=2, lwd=0.7)+
#   # scale_color_manual(values=c("darkorange2","gold2")) +
#   scale_color_discrete(type=terrain.colors(7))+
#   # geom_vline(xintercept = 0, col="black") +
#   scale_y_continuous(breaks=seq(2005,2024, by=1), labels=seq(2005,2024, by=1))+
#   # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
#   labs(x="", y="")   + facet_wrap(~dept2) +
#   theme(legend.position = "none", 
#         panel.border = element_blank(),  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank())
# 
# 
# 
# 
# 
# 


# ##############################################-
# #*----- Par maille dans la zone alpine ----
# 
# # Pour chaque maille on veut une donnée de température et l'altitude moyenne (ou classe d'altitude)
# 
# # Pour choisir la taille de la maille, on regarde la distribution spatiale des points Phénoclim, et le nombre qu'on peut en avoir en fonction du
# # maillage sélectionné
# rast_pred_deb_fixef = vect(phenoclim[phenoclim$nom_massif_v2019 == "Alpes",], geom=c("coord_x_2154", "coord_y_2154"), crs="epsg:2154")
# rast_pred_deb_fixef = rasterize(x=rast_pred_deb_fixef, y=rast(ext(rast_pred_deb_fixef), resolution=10000), field="julian_day",fun="count")
# plot(rast_pred_deb_fixef)
# 
# #                    => il y a pas mal de trous, et après discussion avec Marjo l'approche par département convient


##############################################-
#*----- Par espace protégé ----

# L'idée serait de calculer pour chaque année x espèce x espprot x altitude (ou classe d'altitude) la date de débourrement, en fonction de la 
# température moyenne mesurée une année donnée dans la classe d'altitude x espprot.
# Pour calculer ces variables de température, on se base sur les reconstructions de température.

Tperiodes = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/varTprintemps_allphenosites.csv")
# phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/data/_CLEANED_data_pheno.csv")


# Espaces protégés considérés : PN Vanoise, PN Écrins, RNs + Parcs italiens du Gran Paradisio (ID = 555580376) et Mont Avic (ID = 555528148) + Espace Mt Blanc
# => data France : https://data.naturefrance.fr/geonetwork/srv/fre/catalog.search#/metadata/cc1bfe07-625e-4df6-86a2-49fd07d193d7
RNNs = vect("data/_lim_admin/espaces_proteges/N_ENP_RNN_S_000.shp")
RNRs = vect("data/_lim_admin/espaces_proteges/N_ENP_RNR_S_000.shp")
PNs = vect("data/_lim_admin/espaces_proteges/N_ENP_PN_S_000.shp")
PNRs = vect("data/_lim_admin/espaces_proteges/pnr_polygonPolygon.shp", opts="ENCODING=LATIN1")
# => data Italie : UNEP-WCMC and IUCN (2025), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], May 2025, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net. 
ITprot = vect("data/_lim_admin/espaces_proteges/WDPA_WDOECM_May2025_Public_EU_shp_1/WDPA_WDOECM_May2025_Public_EU_shp-polygons.shp")
ITprot = ITprot[ITprot$WDPAID == 555580376 | ITprot$WDPAID == 555528148,] # On sélectionne uniquement les 2 parcs qui nous intéressent
# => Espace Mont Blanc : https://www.espace-mont-blanc.com/asset/img/carte_emb_dl.jpg
EMB = vect("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/05. DONNEES/Shapefiles/Perim_administratif/EMB/Perimetre_EMB_2017_poly.shp")

espprot = phenoclim[,c("id_base_site","coord_x_2154","coord_y_2154")]
espprot$espprot_RNN = extract(RNNs, vect(espprot, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$NOM_SITE
espprot$espprot_RNR = extract(RNRs, vect(espprot, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$NOM_SITE
espprot$espprot_PN = extract(PNs, vect(espprot, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$NOM_SITE
espprot$espprot_PNR = extract(project(PNRs, "epsg:2154"), vect(espprot, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$name
espprot$espprot_ITprot = extract(project(ITprot, "epsg:2154"), vect(espprot, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$NAME

# On groupe tout en une seule variable, puisqu'en théorie il n'y a pas de recoupement RNR / RN / PN /ITprot
# /!\ l'espace Mont-Blanc doit être considéré séparément !!
phenoclim$espprot = apply(espprot[,grep("espprot", colnames(espprot))], 1, 
                          function(x){ifelse(length(x[!is.na(x)])==0,NA,paste(c("RNN","RNR","PN","","")[!is.na(x)],x[!is.na(x)],sep=" "))})
phenoclim$EMB = extract(project(EMB, "epsg:2154"), vect(phenoclim, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))$EMB
phenoclim$EMB[!is.na(phenoclim$EMB)] = "EMB"


Tperiodes = merge(Tperiodes, 
                  phenoclim[!duplicated(phenoclim$id_base_site),c("id_base_site","altitude","cl_alt","cl_alt2","cl_alt3",
                                                                  "coord_x_2154", "coord_y_2154",
                                                                  "dept1","region1","pays1","nom_massif","nom_massif_v2019","espprot","EMB")],
                  by.x="sitePheno", by.y="id_base_site", all.x=T)



# # Aperçu des zones (espace protégé x classe d'altitude) où on a des données
# table(Tperiodes$espprot, Tperiodes$cl_alt)
# table(phenoclim$espprot, phenoclim$year)
# table(Tperiodes$espprot, Tperiodes$annee)

##########################################################################################################################################-
# /!\ il n'y a pas des données pour tous les espaces et toutes les années --> ajuster les prédictions !!!
##########################################################################################################################################-



T_pourpred = Tperiodes %>% group_by(espprot, annee, esp, cl_alt,
                                    cl_alt2, cl_alt3, region1, pays1, nom_massif) %>% summarise(across(ends_with("0624"), mean))
T_pourpred$altitude = as.numeric(as.character(factor(T_pourpred$cl_alt,
                                                     levels = c("150-450", "450-750", "750-1050", "1050-1350", "1350-1650", "1650-1950", "1950-2250"),
                                                     labels = c(300, 600, 900, 1200, 1500, 1800, 2100))))


# On peut donc utiliser les modèles de chaque espèce pour faire les prédictions de débourrement
#         --------------------------------------------------------------------------------------
#           /!\ je n'ai pas fait de modèle pour le tussilage, la primevère et le pin sylvestre !! À voir s'ils sont dans l'indice calculé par Marjo & co
#         --------------------------------------------------------------------------------------

indice_deb = T_pourpred %>% rename(Tmoy30j = Tmoy30j_med_0624,
                                   Tmoy40j = Tmoy40j_med_0624,
                                   Tmoy50j = Tmoy50j_med_0624,
                                   gel30j = gel30j_med_0624,
                                   gel40j = gel40j_med_0624,
                                   gel50j = gel50j_med_0624,
                                   GDD0 = GDD0_med_0624,
                                   GDD5 = GDD5_med_0624,
                                   Tmoyhiv = Tmoyhiv_med_0624,
                                   dChill = dChill_med_0624,
                                   year = annee)
indice_deb$cl_alt = factor(indice_deb$cl_alt, levels=c("150-450","450-750" ,  "750-1050"   ,"1050-1350" , "1350-1650" ,"1650-1950", "1950-2250" ), 
                           ordered = T)
indice_deb$yearQ = factor(indice_deb$year)

# Chargement des best models
load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_deb.Rdata") 
# Prédictions
# RQ : on peut choisir de prédire en prenant en compte les effets aléatoires ou sans. Ici les effets aléatoires sont ID_zone (effet local, qu'on
#      cherche à lisser en agrégeant par département) et yearQ (effet annuel lié aux spécificités climatiques non prises en compte dans le modèle,
#      ce qu'on peut vouloir garder). Dans les prédictions V1, l'effet ID_zone n'est pas pris en compte, et l'effet year est pris en compte comme 
#      effet fixe dans le modèle. Dans la V2, on choisit de faire avec et sans l'effet aléatoire yearQ (en factor).
for (esp in c("Bouleau_verruqueux","Noisetier","Frene","Meleze","Epicea","Sorbier","Lilas")){#names(bestmods)){
  print(esp)
  if ("cl_alt" %in% colnames(bestmods[[esp]]@frame)){ # On sépare le cas des classes d'altitude comme variable explicative (variable factorielle dont tous les niveaux ne sont pas nécessairement présents dans le tableau pour les prédictions)
    # Prédiction
    indice_deb$pred_deb_fixef[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp,], re.form=NA)
    indice_deb$pred_deb_allef[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$cl_alt %in% unique(bestmods[[esp]]@frame$cl_alt)  & indice_deb$esp == esp,], re.form=~(1|yearQ))
  } else {
    # Prédiction
    indice_deb$pred_deb_fixef[indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$esp == esp,], re.form=NA)
    indice_deb$pred_deb_allef[indice_deb$esp == esp] = predict(bestmods[[esp]], indice_deb[indice_deb$esp == esp,], re.form=~(1|yearQ))
  }
}




indice_deb_espprot__all = indice_deb %>% group_by(year, espprot, cl_alt, cl_alt2) %>% summarise(pred_deb_fixef = mean(pred_deb_fixef, na.rm=T),
                                                                                      pred_deb_allef = mean(pred_deb_allef, na.rm=T),
                                                                                      nb_esp = length(unique(esp)))
# # Si on veut sélectionner uniquement les prédictions où il y a des données d'observation (pour espèce x année x classe d'altitude) :
# indice_deb_espprot__selecobs = indice_deb %>% filter(esp %in% phenoclim$species[phenoclim$espprot[phenoclim$cl_alt == cl_alt] == espprot]) %>% 
#   group_by(year, espprot, cl_alt, cl_alt2) %>% 
#   summarise(pred_deb_fixef = mean(pred_deb_fixef, na.rm=T),
#             pred_deb_allef = mean(pred_deb_allef, na.rm=T),
#             nb_esp = length(unique(esp[esp %in% phenoclim$species[phenoclim$espprot[phenoclim$cl_alt == cl_alt] == espprot]])))
# ============================= /!\ À CORRIGER POUR QUE ÇA MARCHE !! ============================================================================================-

write.csv(indice_deb_espprot__all, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_espprot__extra.csv", row.names = F)
# write.csv(indice_deb_espprot__selecobs, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_espprot__ouobs.csv", row.names = F)


#----- Visualisation des indices par espace protégé

indice_deb_espprot$cl_2alt = ifelse(indice_deb_espprot$cl_alt2 == "<1050", "Inf1050", "Sup1050")

indice_deb_espprot_V2 = indice_deb_espprot[!is.na(indice_deb_espprot$espprot),] %>% group_by(espprot,cl_2alt,year) %>% summarise(pred_deb_fixef = mean(pred_deb_fixef),
                                                                                    pred_deb_allef = mean(pred_deb_allef))
indice_deb_espprot_V2 = merge(indice_deb_espprot_V2,
                           indice_deb_espprot %>% group_by(espprot,cl_2alt) %>% summarise(pred_deb_fixef_ref = mean(pred_deb_fixef),
                                                                                     pred_deb_allef_ref = mean(pred_deb_allef)),
                           by=c("espprot", "cl_2alt"))
indice_deb_espprot_V2$diff_allef = indice_deb_espprot_V2$pred_deb_allef - indice_deb_espprot_V2$pred_deb_allef_ref
indice_deb_espprot_V2$diff_fixef = indice_deb_espprot_V2$pred_deb_fixef - indice_deb_espprot_V2$pred_deb_fixef_ref

plot_indice_espprot_fixef = ggplot(indice_deb_espprot_V2, 
                              aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                  x=diff_fixef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diff_fixef, xend=+Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  scale_color_manual(values=c("darkgreen","yellowgreen"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2024, by=1), labels=seq(2006,2024, by=1))+
  # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="", title = "Indice de débourrement - V2 fixed effects")   + facet_wrap(~espprot) + #xlim(50,170)+
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
plot_indice_espprot_allef = ggplot(indice_deb_espprot_V2, 
                              aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                  x=diff_allef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diff_allef, xend=+Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  scale_color_manual(values=c("darkgreen","yellowgreen"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2024, by=1), labels=seq(2006,2024, by=1))+
  # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="", title = "Indice de débourrement - V2 all effects")   + facet_wrap(~espprot) + #xlim(50,170)+
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())



pdf("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_espprot.pdf", width=12, height = 10)
plot_indice_espprot_fixef ; plot_indice_espprot_allef
dev.off()





##############################################-
#*----- Comparaison de l'indice obtenu avec les indices précédemment calculés (selon le script de Marjorie et Colin) ----

indice_deb_V1_dept = read.csv2("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/Analyse/Indice_phenoclim/pheno_year_deb6.csv", sep=",", dec=".")
indice_deb_V1 = read.csv2("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/Analyse/Indice_phenoclim/pheno_year_global.csv", sep=",", dec=".", row.names=1)

indice_deb$cl_2alt = ifelse(indice_deb$cl_alt2 == "<1050", "Inf1050", "Sup1050")
indice_deb_dept$cl_2alt = ifelse(indice_deb_dept$cl_alt2 == "<1050", "Inf1050", "Sup1050")

indice_deb_V2_dept = indice_deb_dept %>% group_by(dept2,cl_2alt,year) %>% summarise(pred_deb_fixef = mean(pred_deb_fixef),
                                                                                    pred_deb_allef = mean(pred_deb_allef))
indice_deb_V2_dept = merge(indice_deb_V2_dept,
                           indice_deb_dept %>% group_by(dept2,cl_2alt) %>% summarise(pred_deb_fixef_ref = mean(pred_deb_fixef),
                                                                                     pred_deb_allef_ref = mean(pred_deb_allef)),
                           by=c("dept2", "cl_2alt"))
indice_deb_V2_dept$diff_allef = indice_deb_V2_dept$pred_deb_allef - indice_deb_V2_dept$pred_deb_allef_ref
indice_deb_V2_dept$diff_fixef = indice_deb_V2_dept$pred_deb_fixef - indice_deb_V2_dept$pred_deb_fixef_ref

  
#----- Visualisation de l'écart entre les différentes méthodes de calcul des indices de débourrement

test_indice = merge(indice_deb_V2_dept, indice_deb_V1_dept[,c("year","dept1.x","cl_2alt.x","pred.x")],  
                    by.x=c("year","dept2","cl_2alt"), by.y=c("year","dept1.x","cl_2alt.x"))

ggplot(test_indice, aes(y=pred.x, x=pred_deb_fixef)) + geom_point(size=2, aes(col=dept2)) + #, shape=cl_2alt
  geom_abline(intercept = 0, slope=1, lty=2) + scale_colour_brewer(palette = "Paired") + 
  facet_wrap(~cl_2alt) +
  theme(legend.position="bottom")
ggplot(test_indice, aes(y=pred.x, x=pred_deb_allef)) + geom_point(size=2, aes(col=dept2)) + #, shape=cl_2alt
  geom_abline(intercept = 0, slope=1, lty=2) + scale_colour_brewer(palette = "Paired") + 
  facet_wrap(~cl_2alt) +
  theme(legend.position="bottom")

# cor.test(test_indice$pred_deb_fixef, test_indice$pred.x)
# cor.test(test_indice$pred_deb_allef, test_indice$pred.x)


#----- Visualisation des indices

plot_indice_v1 = ggplot(indice_deb_V1_dept, 
                        aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt.x, levels=unique(cl_2alt.x), labels=c(-1,1)))), 
                            x=diff, col=cl_2alt.x)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diff, xend=+Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt.x, levels=unique(cl_2alt.x), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt.x, levels=unique(cl_2alt.x), labels=c(-1,1)))),
                   col=cl_2alt.x ), lty=3, lwd=0.7)+
  scale_color_manual(values=c("darkgreen","yellowgreen"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2024, by=1), labels=seq(2006,2024, by=1))+
  # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="", title = "Indice de débourrement - V1")   + facet_wrap(~dept1.x) + #xlim(50,170) +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
plot_indice_v2_fixef = ggplot(indice_deb_V2_dept[indice_deb_V2_dept$dept2 %in% indice_deb_V1_dept$dept1.x,], 
                              aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                  x=diff_fixef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diff_fixef, xend=+Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  scale_color_manual(values=c("darkgreen","yellowgreen"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2024, by=1), labels=seq(2006,2024, by=1))+
  # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="", title = "Indice de débourrement - V2 fixed effects")   + facet_wrap(~dept2) + #xlim(50,170)+
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
plot_indice_v2_allef = ggplot(indice_deb_V2_dept[indice_deb_V2_dept$dept2 %in% indice_deb_V1_dept$dept1.x,], 
                              aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                  x=diff_allef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diff_allef, xend=+Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  scale_color_manual(values=c("darkgreen","yellowgreen"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2024, by=1), labels=seq(2006,2024, by=1))+
  # scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="", title = "Indice de débourrement - V2 all effects")   + facet_wrap(~dept2) + #xlim(50,170)+
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())

pdf("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_dept.pdf", width=12, height = 10)
plot_indice_v1 ; plot_indice_v2_fixef ; plot_indice_v2_allef
dev.off()






############################################################################################-
# DIFFÉRENCES SELON LES CATÉGORIES DE PARTICIPANT·ES                                     ----
############################################################################################-

load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_deb.Rdata")

coeff_cat = data.frame(species = NA, variable = NA, coef_all = NA, coef_PROF = NA, coef_SCOL = NA, coef_PART = NA)

##############################################-
#*---- Bouleau verruqueux ----

deb_Bpen = debourr_Alps[debourr_Alps$species == "Bouleau_verruqueux",]
deb_Bpen = deb_Bpen[!is.na(deb_Bpen$Tmoy30j),]

coef_all = coef(summary(lmer(formula(bestmods$Bouleau_verruqueux), deb_Bpen, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Bouleau_verruqueux), deb_Bpen[deb_Bpen$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Bouleau_verruqueux), deb_Bpen[deb_Bpen$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Bouleau_verruqueux), deb_Bpen[deb_Bpen$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Bouleau_verruqueux"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])

##############################################-
#*---- Noisetier ----

deb_Cave = debourr_Alps[debourr_Alps$species == "Noisetier",]
deb_Cave = deb_Cave[!is.na(deb_Cave$Tmoy30j),]

coef_all = coef(summary(lmer(formula(bestmods$Noisetier), deb_Cave, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Noisetier), deb_Cave[deb_Cave$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Noisetier), deb_Cave[deb_Cave$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Noisetier), deb_Cave[deb_Cave$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Noisetier"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])

##############################################-
#*---- Frêne ----

deb_Fexc = debourr_Alps[debourr_Alps$species == "Frene",]
deb_Fexc = deb_Fexc[!is.na(deb_Fexc$Tmoy30j),]

coef_all = coef(summary(lmer(formula(bestmods$Frene), deb_Fexc, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Frene), deb_Fexc[deb_Fexc$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Frene), deb_Fexc[deb_Fexc$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Frene), deb_Fexc[deb_Fexc$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Frene"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])
##############################################-
#*---- Mélèze ----    

deb_Ldec = debourr_Alps[debourr_Alps$species == "Meleze",]
deb_Ldec = deb_Ldec[!is.na(deb_Ldec$Tmoy30j),]
deb_Ldec = deb_Ldec[deb_Ldec$cl_alt != "150-450",]

# table(deb_Ldec$cl_alt, deb_Ldec$new_cat)
# # /!\ la classe d'altitude la plus élevée n'est pas présente pour les scolaires et les particuliers !!
coef_all = coef(summary(lmer(formula(bestmods$Meleze), deb_Ldec, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Meleze), deb_Ldec[deb_Ldec$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Meleze), deb_Ldec[deb_Ldec$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Meleze), deb_Ldec[deb_Ldec$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
              merge(coef_PROF[,1], 
                    merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                    by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Meleze"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])
  

##############################################-
#*---- Épicéa ----

deb_Pabi = debourr_Alps[debourr_Alps$species == "Epicea",]
deb_Pabi = deb_Pabi[!is.na(deb_Pabi$Tmoy30j),]
deb_Pabi = deb_Pabi[deb_Pabi$cl_alt != "150-450",] # on retire la classe d'altitude la plus basse parce qu'il n'y a pas assez de données pour la suite

# table(deb_Pabi$cl_alt, deb_Pabi$new_cat)
# # /!\ la classe d'altitude la plus élevée n'est pas présente pour les scolaires et les professionnels !!
coef_all = coef(summary(lmer(formula(bestmods$Epicea), deb_Pabi, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Epicea), deb_Pabi[deb_Pabi$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Epicea), deb_Pabi[deb_Pabi$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Epicea), deb_Pabi[deb_Pabi$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Epicea"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])


##############################################-
#*---- Sorbier ----

deb_Sacu = debourr_Alps[debourr_Alps$species == "Sorbier",]
deb_Sacu = deb_Sacu[!is.na(deb_Sacu$Tmoy30j),]

coef_all = coef(summary(lmer(formula(bestmods$Sorbier), deb_Sacu, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Sorbier), deb_Sacu[deb_Sacu$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Sorbier), deb_Sacu[deb_Sacu$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Sorbier), deb_Sacu[deb_Sacu$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Sorbier"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])

##############################################-
#*---- Lilas ----

deb_Svul = debourr_Alps[debourr_Alps$species == "Lilas",]
deb_Svul = deb_Svul[!is.na(deb_Svul$Tmoy30j),]

coef_all = coef(summary(lmer(formula(bestmods$Lilas), deb_Svul, REML=F)))
coef_PROF = coef(summary(lmer(formula(bestmods$Lilas), deb_Svul[deb_Svul$new_cat=="Professionnels",], REML=F)))
coef_SCOL = coef(summary(lmer(formula(bestmods$Lilas), deb_Svul[deb_Svul$new_cat=="Multiscolaires",], REML=F)))
coef_PART = coef(summary(lmer(formula(bestmods$Lilas), deb_Svul[deb_Svul$new_cat=="Particuliers",], REML=F)))
coef = merge(coef_all[,1],
             merge(coef_PROF[,1], 
                   merge(coef_SCOL[,1], coef_PART[,1], by="row.names", all.x=T, all.y=T), 
                   by.x="row.names", by.y="Row.names",all.x=T, all.y=T), 
             by.x="row.names", by.y="Row.names",all.x=T, all.y=T)
colnames(coef) = c("variable","coef_all","coef_PROF","coef_SCOL","coef_PART")
coef$species = "Lilas"

coeff_cat = rbind(coeff_cat, coef[,colnames(coeff_cat)])



write.csv(coeff_cat[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_deb__coeff_entrecateg.csv", row.names = F)

