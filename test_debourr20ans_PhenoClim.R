##############################################-
# ANALYSE DES 20 ANS DE DONNÉES PHÉNOCLIM  ----
#  Ninon Fontaine - février 2025              -
##############################################-


library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
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
# 


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
phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
# Données obtenues avec le script 1_mise_en_forme_BDD.R



##############################################-
# INTÉGRATION D'UNE VARIABLE EXPLICATIVE 'TEMPÉRATURE' ----
##############################################-

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

#*------------ Récupération des données météo reconstruites et calcul de variable pertinente pour la phénologie ----

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


Tperiod <- function(jourcible, site,  tabT, 
                    vartabT = c("julian_day","station_name","Tair_moy"),
                    dureej, tempcalc="Tmoy_Xj"){
  # NB : le tableau de températures en entrée est un tableau avec une seule année, donc on peut fonctionner en "julian day"
  tabT = rename(tabT, c(jul_day=vartabT[1], sitePheno=vartabT[2], T=vartabT[3]))
  output = matrix(ncol=length(tempcalc))
  colnames(output)=tempcalc
  # tabT$date = as.Date(tabT$date)
  # jourcible = as.Date(jourcible)
  if ("Tmoy_Xj" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
    output[,"Tmoy_Xj"] = mean(Tcalc$T, na.rm=T)
  }
  if ("Tmin" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
    output[,"Tmin"]=min(Tcalc$T, na.rm=T)
  }
  if ("Tmax" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= jourcible-dureej & sitePheno ==site)
    output[,"Tmax"]=max(Tcalc$T, na.rm=T)
  }
  if ("GDD0" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
    output[,"GDD0"]=sum(Tcalc$T[Tcalc$T > 0], na.rm=T)
  }
  if ("GDD5" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= 0 & sitePheno ==site) #on calcule les GDD sur la période du 1er janvier au débourrement
    output[,"GDD5"]=sum(Tcalc$T[Tcalc$T >=5]-5, na.rm=T)
  }
  if ("Tmoyhiv" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= 31 & jul_day >= -60 & sitePheno ==site) #on définit l'hiver comme la période novembre-janvier (Vitasse)
    output[,"Tmoyhiv"]=mean(Tcalc$T, na.rm=T)
  }
  if ("dChill" %in% tempcalc){
    Tcalc = tabT %>% filter(jul_day <= jourcible & jul_day >= -60 & sitePheno ==site) #on calcule les jours de chilling du 1/11 au débourrement (Vitasse) - en proportion
    output[,"dChill"]=nrow(Tcalc[Tcalc$T <=8 & Tcalc$T>=0,]) / nrow(Tcalc)
  }
  return(output)
}

# test=debourrBpen_Alps[debourrBpen_Alps$year ==2023,c("visit_date", "julian_day", "id_base_site", "altitude", "nom_zone")]
# reconstruc_sites_2023$date = as.Date(reconstruc_sites_2023$date)
# 
# test$Tmoy_30j = apply(test[,c("julian_day", "id_base_site")],1, 
#                     FUN=function(x){Tperiod(jourcible = x[1], site=x[2], dureej = 30, tabT=reconstruc_sites_2023)})
# # plot(as.Date(test$visit_date[!is.na(test$visit_date)]), test$Tcalc[!is.na(test$visit_date)])
# ggplot(test, aes(y=as.Date(visit_date), x=Tmoy_30j, col=altitude)) + geom_point()
# 
# 
# global_deb1 <- lmer(julian_day~I(altitude-1100)+I(Tmoy_30j)+(I(year-2014)|nom_zone), test, REML=FALSE)
# summary(global_deb1)
# global_debT <- lmer(julian_day~I(Tmoy_30j)+(I(year-2014)|nom_zone), test, REML=FALSE)
# summary(global_debT)
# global_debAlt <- lmer(julian_day~I(altitude-1100)+(I(year-2014)|nom_zone), test, REML=FALSE)
# summary(global_debAlt)
# 
# visreg(global_deb1, "Tmoy_30j")
# visreg(global_deb1, "altitude")
# visreg(global_debAlt, "altitude")


# Construction d'un tableau sur le débourrement dans les Alpes, où chaque donnée (site x année) est associée à une valeur de température supposée
# pertinente pour la phénologie : la température moyenne pendant les 30 jours précédant la date de débourrement.

debourr_Alps = phenoclim[phenoclim$pheno_etape_value == "Debourrement - Ok 10%" & 
                           phenoclim$nom_massif == "Alpes",]
# /!\ il y a beaucoup de lignes remplies de NA -> à supprimer :
debourr_Alps = debourr_Alps[!is.na(debourr_Alps$ids_observers),]

varTselec = c("Tmoy_Xj","GDD0","GDD5","Tmoyhiv","dChill")
debourr_Alps[,varTselec] = NA

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
  
  Tcalcs = apply(debourr_Alps[debourr_Alps$year == annee, c("julian_day", "id_base_site")],1, 
                 FUN=function(x){Tperiod(jourcible = x[1], site=x[2], dureej = 30, 
                                         tabT=reconstruc_Tday_sites,
                                         tempcalc = varTselec)})
  
  debourr_Alps[debourr_Alps$year == annee, varTselec] = matrix(unlist(Tcalcs),ncol=length(varTselec),byrow=T)
}

debourr_Alps = debourr_Alps %>% rename("Tmoy_30j"="Tmoy_Xj")

write.csv(debourr_Alps, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Debourr_Alps_T.csv")



##############################################-
# MODÈLES DE DÉBOURREMENT ----
##############################################-

debourr_Alps = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Debourr_Alps.csv", row.names=1)
debourr_Alps$yearQ = factor(debourr_Alps$year)

# Aperçu des données de débourrement
ggplot(debourr_Alps[!is.na(debourr_Alps$Tmoy_30j),], 
       aes(x=altitude, y=julian_day, col=nom_cite)) + geom_point() + theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
ggplot(debourr_Alps[!is.na(debourr_Alps$Tmoy_30j),], 
       aes(x=Tmoy_30j, y=julian_day, col=nom_cite)) + geom_point() + theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
# À première vue, la relation entre date de débourrement et altitude est plus nette que la relation entre date et température (notamment pour 
# l'épicéa) --> tester dans les modèles si l'ajout de cette variable température améliore les modèles par rapport à des modèles avec seulement
# l'altitude.




##############################################-
#*---- Bouleau verruqueux ----

deb_Bpen = debourr_Alps[debourr_Alps$species == "Bouleau_verruqueux",]
# ggplot(deb_Bpen, aes(x=year)) + geom_histogram(binwidth=1) 
# # on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Bpen = deb_Bpen[!is.na(deb_Bpen$Tmoy_30j),]
# ggplot(deb_Bpen, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On reprend le modèle de Bison et al. 2019, pour voir si on a à peu près les mêmes résultats ---- 
mod_deb_Bpen_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Bpen, REML=F)
summary(mod_deb_Bpen_Alt)
deb_Bpen_10ans = deb_Bpen[deb_Bpen$year <= 2016,]
mod_deb_Bpen_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Bpen_10ans, REML=F)
summary(mod_deb_Bpen_all_10ans)
# => OK : on retrouve à peu près les mêmes résultats que Bison et al 2019

#*-------- 2) On teste l'effet de différentes variables de température à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)
mod_deb_Bpen_Tmoy <- lmer(julian_day ~ I(Tmoy_30j) + (1|nom_zone) + (1|year), deb_Bpen, REML=F)
mod_deb_Bpen_GDD0 <- lmer(julian_day ~ I(GDD0) + (1|nom_zone) + (1|year), deb_Bpen, REML=F)
mod_deb_Bpen_GDD5 <- lmer(julian_day ~ I(GDD5) + (1|nom_zone) + (1|year), deb_Bpen, REML=F)

AIC(mod_deb_Bpen_Tmoy) ; AIC(mod_deb_Bpen_GDD0) ; AIC(mod_deb_Bpen_GDD5)
# La variable GDD0 donne largement les meilleurs résultats : c'est celle qu'on garde pour la suite

# On notera que les GDD0 augmentent au cours du temps : la végétation peut accumuler plus de degrés-jour
ggplot(deb_Bpen, aes(x=year, y=GDD0, col=cl_alt)) + geom_point() + geom_smooth(method="lm")

# # Visualisation
# ggplot(deb_Bpen, aes(x=GDD0, y=julian_day)) + geom_point() + geom_smooth(method = "lm") 
# visreg(mod_deb_Bpen_GDD0, "GDD0")

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle, en ajoutant un effet altitude
mod_deb_Bpen_GDD0_Alt <- lmer(julian_day ~ I(GDD0) +I(altitude) + (1|nom_zone) + (1|yearQ), deb_Bpen, REML=F)
AIC(mod_deb_Bpen_GDD0) ; AIC(mod_deb_Bpen_GDD0_Alt)
# Mieux avec l'altitude en plus
summary(mod_deb_Bpen_GDD0_Alt)
# gradient altitudinal : 3.3 jours de retard quand on monte de 100m (colle avec la règle 'classique')
# visreg(mod_deb_Bpen_GDD0_Alt, "altitude")

# Effets fixes VS tous les effets inclus
r.squaredGLMM(mod_deb_Bpen_GDD0_Alt)
# Part de variance expliquée par les différentes variables
calcVarPart(mod_deb_Bpen_GDD0_Alt)



mod_deb_Bpen_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Bpen, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Bpen_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Bpen, REML=F) 
# mod_deb_Bpen_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Bpen, REML=F) 

AIC(mod_deb_Bpen_Tmoy) ; AIC(mod_deb_Bpen_GDD0) ; AIC(mod_deb_Bpen_GDD5) ; AIC(mod_deb_Bpen_Alt) ;AIC(mod_deb_Bpen_all) #; AIC(mod_deb_Bpen_interT) ; AIC(mod_deb_Bpen_interAlt)

summary(mod_deb_Bpen_all)
summary(mod_deb_Bpen_Tmoy)
# summary(mod_deb_Bpen_interT)
# summary(mod_deb_Bpen_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a une AVANCÉE dans le débourrement au fil des années, de l'ordre de -0.6j par décennie.
# La température a au contraire un effet RETARD, de l'ordre de +1.8j par degré supplémentaire... (??)

visreg(mod_deb_Bpen_all, "altitude")
visreg(mod_deb_Bpen_all, "Tmoy_30j")
visreg(mod_deb_Bpen_all, "year")
# visreg(mod_deb_Bpen_interAlt, "year")
# visreg(mod_deb_Bpen_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
deb_Bpen_10ans = deb_Bpen[deb_Bpen$year <= 2016,]
mod_deb_Bpen_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Bpen_10ans, REML=F)
summary(mod_deb_Bpen_all_10ans)
# Sur cette période, on a aussi une AVANCÉE dans le débourrement au fil des années, mais de l'ordre de -4j par décennie.



##############################################-
#*---- Noisetier ----

deb_Cave = debourr_Alps[debourr_Alps$species == "Noisetier",]
ggplot(deb_Cave, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Cave = deb_Cave[!is.na(deb_Cave$Tmoy_30j),]
ggplot(deb_Cave, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Noisetier", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Cave_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Cave, REML=F)
mod_deb_Cave_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Cave, REML=F)
mod_deb_Cave_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Cave, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Cave_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Cave, REML=F) 
# mod_deb_Cave_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Cave, REML=F) 

AIC(mod_deb_Cave_all) ; AIC(mod_deb_Cave_T) ; AIC(mod_deb_Cave_Alt) #; AIC(mod_deb_Cave_interT) ; AIC(mod_deb_Cave_interAlt)

summary(mod_deb_Cave_all)
summary(mod_deb_Cave_T)
summary(mod_deb_Cave_Alt)
# summary(mod_deb_Cave_interT)
# summary(mod_deb_Cave_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a une AVANCÉE dans le débourrement au fil des années, de l'ordre de -3.9j par décennie.
# La température a au contraire un effet RETARD, de l'ordre de +3.0j par degré supplémentaire... (??)

visreg(mod_deb_Cave_all, "altitude")
visreg(mod_deb_Cave_all, "Tmoy_30j")
visreg(mod_deb_Cave_all, "year")
# visreg(mod_deb_Cave_interAlt, "year")
# visreg(mod_deb_Cave_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Cave_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Cave[deb_Cave$year <= 2016,], REML=F)
summary(mod_deb_Cave_all_10ans)
# Sur cette période, on a aussi une AVANCÉE dans le débourrement au fil des années, mais de l'ordre de -2.3j par décennie.



##############################################-
#*---- Frêne ----
deb_Fexc = debourr_Alps[debourr_Alps$species == "Frene",]
ggplot(deb_Fexc, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Fexc = deb_Fexc[!is.na(deb_Fexc$Tmoy_30j),]
ggplot(deb_Fexc, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Frêne", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Fexc_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Fexc, REML=F)
mod_deb_Fexc_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Fexc, REML=F)
mod_deb_Fexc_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Fexc, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Fexc_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Fexc, REML=F) 
# mod_deb_Fexc_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Fexc, REML=F) 

AIC(mod_deb_Fexc_all) ; AIC(mod_deb_Fexc_T) ; AIC(mod_deb_Fexc_Alt) #; AIC(mod_deb_Fexc_interT) ; AIC(mod_deb_Fexc_interAlt)

summary(mod_deb_Fexc_all)
summary(mod_deb_Fexc_T)
summary(mod_deb_Fexc_Alt)
# summary(mod_deb_Fexc_interT)
# summary(mod_deb_Fexc_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a un RETARD dans le débourrement au fil des années, de l'ordre de +3j par décennie.
# La température a aussi un effet retard, de l'ordre de +3j par degré supplémentaire... (??)

visreg(mod_deb_Fexc_all, "altitude")
visreg(mod_deb_Fexc_all, "Tmoy_30j")
visreg(mod_deb_Fexc_all, "year")
# visreg(mod_deb_Fexc_interAlt, "year")
# visreg(mod_deb_Fexc_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Fexc_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Fexc[deb_Fexc$year <= 2016,], REML=F)
summary(mod_deb_Fexc_all_10ans)
# Sur cette période, on a aussi un RETARD dans le débourrement au fil des années, mais de l'ordre de +2j par décennie.

# REMARQUE : 
r.squaredGLMM(mod_deb_Fexc_all)
# l'effet aléatoire représente quand même une bonne partie de la variance expliquée par le modèle... (~ 20 %)



##############################################-
#*---- Mélèze ----

deb_Ldec = debourr_Alps[debourr_Alps$species == "Meleze",]
ggplot(deb_Ldec, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Ldec = deb_Ldec[!is.na(deb_Ldec$Tmoy_30j),]
ggplot(deb_Ldec, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Mélèze", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Ldec_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Ldec, REML=F)
mod_deb_Ldec_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Ldec, REML=F)
mod_deb_Ldec_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Ldec, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Ldec_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Ldec, REML=F) 
# mod_deb_Ldec_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Ldec, REML=F) 

AIC(mod_deb_Ldec_all) ; AIC(mod_deb_Ldec_T) ; AIC(mod_deb_Ldec_Alt) #; AIC(mod_deb_Ldec_interT) ; AIC(mod_deb_Ldec_interAlt)

summary(mod_deb_Ldec_all)
summary(mod_deb_Ldec_T)
summary(mod_deb_Ldec_Alt)
# summary(mod_deb_Ldec_interT)
# summary(mod_deb_Ldec_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a une AVANCÉE dans le débourrement au fil des années, de l'ordre de -1.5j par décennie.
# La température a au contraire un effet retard, de l'ordre de +1.4j par degré supplémentaire... (??)

visreg(mod_deb_Ldec_all, "altitude")
visreg(mod_deb_Ldec_all, "Tmoy_30j")
visreg(mod_deb_Ldec_all, "year")
# visreg(mod_deb_Ldec_interAlt, "year")
# visreg(mod_deb_Ldec_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Ldec_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Ldec[deb_Ldec$year <= 2016,], REML=F)
summary(mod_deb_Ldec_all_10ans)
# Sur cette période, on a aussi une AVANCÉE dans le débourrement au fil des années, mais de l'ordre de -2.5j par décennie.


##############################################-
#*---- Epicéa ----

deb_Pabi = debourr_Alps[debourr_Alps$species == "Epicea",]
ggplot(deb_Pabi, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Pabi = deb_Pabi[!is.na(deb_Pabi$Tmoy_30j),]
ggplot(deb_Pabi, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Epicéa", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Pabi_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Pabi, REML=F)
mod_deb_Pabi_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Pabi, REML=F)
mod_deb_Pabi_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Pabi, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Pabi_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Pabi, REML=F) 
# mod_deb_Pabi_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Pabi, REML=F) 

AIC(mod_deb_Pabi_all) ; AIC(mod_deb_Pabi_T) ; AIC(mod_deb_Pabi_Alt) #; AIC(mod_deb_Pabi_interT) ; AIC(mod_deb_Pabi_interAlt)

summary(mod_deb_Pabi_all)
summary(mod_deb_Pabi_T)
summary(mod_deb_Pabi_Alt)
# summary(mod_deb_Pabi_interT)
# summary(mod_deb_Pabi_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a un RETARD dans le débourrement au fil des années, de l'ordre de +4.0j par décennie.
# La température a aussi un effet retard, de l'ordre de +1j par degré supplémentaire... (??)

visreg(mod_deb_Pabi_all, "altitude")
visreg(mod_deb_Pabi_all, "Tmoy_30j")
visreg(mod_deb_Pabi_all, "year")
# visreg(mod_deb_Pabi_interAlt, "year")
# visreg(mod_deb_Pabi_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Pabi_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Pabi[deb_Pabi$year <= 2016,], REML=F)
summary(mod_deb_Pabi_all_10ans)
# Sur cette période, on a aussi un RETARD dans le débourrement au fil des années, mais de l'ordre de +1j par décennie.



##############################################-
#*---- Sorbier ----

deb_Sacu = debourr_Alps[debourr_Alps$species == "Sorbier",]
ggplot(deb_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Sacu = deb_Sacu[!is.na(deb_Sacu$Tmoy_30j),]
ggplot(deb_Sacu, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Sorbier", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Sacu_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Sacu, REML=F)
mod_deb_Sacu_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Sacu, REML=F)
mod_deb_Sacu_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Sacu, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Sacu_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Sacu, REML=F) 
# mod_deb_Sacu_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Sacu, REML=F) 

AIC(mod_deb_Sacu_all) ; AIC(mod_deb_Sacu_T) ; AIC(mod_deb_Sacu_Alt) #; AIC(mod_deb_Sacu_interT) ; AIC(mod_deb_Sacu_interAlt)

summary(mod_deb_Sacu_all)
summary(mod_deb_Sacu_T)
summary(mod_deb_Sacu_Alt)
# summary(mod_deb_Sacu_interT)
# summary(mod_deb_Sacu_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a une AVANCÉE dans le débourrement au fil des années, de l'ordre de -0.8j par décennie.
# La température a au contraire un effet retard, de l'ordre de +3j par degré supplémentaire... (??)

visreg(mod_deb_Sacu_all, "altitude")
visreg(mod_deb_Sacu_all, "Tmoy_30j")
visreg(mod_deb_Sacu_all, "year")
# visreg(mod_deb_Sacu_interAlt, "year")
# visreg(mod_deb_Sacu_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Sacu_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Sacu[deb_Sacu$year <= 2016,], REML=F)
summary(mod_deb_Sacu_all_10ans)
# Sur cette période, on a plutôt une AVANCÉE dans le débourrement au fil des années, de l'ordre de -1j par décennie.


##############################################-
#*---- Lilas ----

deb_Svul = debourr_Alps[debourr_Alps$species == "Lilas",]
ggplot(deb_Svul, aes(x=year)) + geom_histogram(binwidth=1) 
# on voit qu'il y a peu de données en 2005, les retirer (parce qu'on n'a pas d'information de température) ne devrait pas changer grand-chose
deb_Svul = deb_Svul[!is.na(deb_Svul$Tmoy_30j),]
ggplot(deb_Svul, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Débourrement 10% - Lilas", x="année",y="jour julien")

# On reprend le modèle de Bison et al. 2019, en le complétant par les données température : 
# /!\ dans cet article, le terme aléatoire porte sur l'année, mais cela pose problème parce qu'il n'y a parfois qu'une seule mesure par année x zone
#     ça me semble plus logique de mettre le terme aléatoire sur l'intercept uniquement, en considérant qu'un·e observateur·rice a un entraînement
#     qui fait qu'il ou elle observe de la même manière d'une année à l'autre. Par contre il peut y avoir un effet général du site et/ou de l'obs.
mod_deb_Svul_all <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Svul, REML=F)
mod_deb_Svul_T <- lmer(julian_day ~ I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Svul, REML=F)
mod_deb_Svul_Alt <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) + (I(year-2014)|nom_zone), deb_Svul, REML=F)
# # pour tester un effet d'avancement inégal le long du gradient altitudinal, on ajoute un terme d'interaction :
# # /!\ difficile à interpréter les interactions, à voir si je les garde !!
# mod_deb_Svul_interT <- lmer(julian_day ~ I(altitude-1100) + I(year-2014) * I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Svul, REML=F) 
# mod_deb_Svul_interAlt <- lmer(julian_day ~ I(altitude-1100) * I(year-2014) + I(Tmoy_30j) + (I(year-2014)|nom_zone), deb_Svul, REML=F) 

AIC(mod_deb_Svul_all) ; AIC(mod_deb_Svul_T) ; AIC(mod_deb_Svul_Alt) #; AIC(mod_deb_Svul_interT) ; AIC(mod_deb_Svul_interAlt)

summary(mod_deb_Svul_all)
summary(mod_deb_Svul_T)
summary(mod_deb_Svul_Alt)
# summary(mod_deb_Svul_interT)
# summary(mod_deb_Svul_interAlt)

# Selon le modèle 'complet' (meilleur score AIC), on a une AVANCÉE dans le débourrement au fil des années, de l'ordre de -6.6j par décennie (beaucoup !).
# La température a au contraire un effet retard, de l'ordre de +3j par degré supplémentaire... (??)

visreg(mod_deb_Svul_all, "altitude")
visreg(mod_deb_Svul_all, "Tmoy_30j")
visreg(mod_deb_Svul_all, "year")
# visreg(mod_deb_Svul_interAlt, "year")
# visreg(mod_deb_Svul_interAlt, "year", by="altitude")

# En utilisant ce même modèle mais pour la période 2006-2016 (cf article Bison et al), on obtient :
mod_deb_Svul_all_10ans <- lmer(julian_day ~ I(altitude-1100) + I(Tmoy_30j) + I(year-2014) + (I(year-2014)|nom_zone), deb_Svul[deb_Svul$year <= 2016,], REML=F)
summary(mod_deb_Svul_all_10ans)
# Sur cette période, on a au contraire un RETARD dans le débourrement au fil des années, de l'ordre de +4.6j par décennie.









