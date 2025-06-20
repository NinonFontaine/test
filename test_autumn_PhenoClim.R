############################################################################################-
# ANALYSE DES 20 ANS DE DONNÉES PHÉNOCLIM - Phénologie automnale                        ----
#  Ninon Fontaine - printemps 2025                                                         -
############################################################################################-


library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(data.table)
library(reshape2)
library(dplyr)
library(terra)
library(raster)
library(lme4)
library(lmerTest) # permet d'avoir les pvalues qui s'affichent pour les lmer
library(visreg)
library(MuMIn) # fonction r.squaredGLMM notamment
library(variancePartition) # fonction calcVarPart
library(brms)
library(sjPlot)
library(meteor) # pour le calcul de la photopériode
library(leaflet)


##############################################-
# OBJECTIFS ----
##############################################-

# Concernant les données PhenoClim et leur analyse, certains éléments ont déjà été explorés (notamment par Marjorie Bison).
# Les analyses se sont principalement concentrées sur les Alpes, la phase de débourrement, et sur une première série de données, de 
# 2004 à 2016. Les phases de floraison, de feuillaison et de *sénescence* restent à explorer, tout comme le cas des autres massifs 
# (Pyrénées, Massif Central, Vosges...).
# 
# La phénologie automnale (sénescence) est moins étudiée que la phénologie printanière (débourrement, floraison, feuillaison) parce que
# ses déterminants sont moins bien identifiés. PhénoClim peut offrir un jeu de données conséquent pour explorer ces questions.
#
# Il peut être intéressant de considérer à la fois le début de la sénescence (sen10%) et sa durée (sen50% - sen10%), les déterminants
# pouvant être différents (Zohner et al 2023).
# Les facteurs explicatifs étudiés jusqu'à présents sont :
# - la photopériode, 
# - les températures minimales (~Tnuit), maximales (~Tjour),
# - une accumulation de froid (sorte de GDD inversé), 
# - les températures moyennes avant ou après le solstice (maijuin vs aoûtseptoct), 
# - la température moyenne estivale (en lien avec la photosynthèse et la productivité estivale)
# - la date de début de saison de végétation, 
# - le déficit hydrique estival (en lien avec la sécheresse, qui limite les transferts de sucres entre organes).
# (Zohner et al 2023, Zani et al 2020, Wu et al 2018, Delpierre et al 2009, Jolly et al 2005)


# Avec l'exploration des effets climatiques, d'autres questions peuvent être posées : 
#   (1) l'ajout de variable climatique améliore-t-il les modèles phénologiques, par rapport aux modèles altitude + année ?
#   (2) y a-t-il une tendance au retard de la sénescence des arbres (et donc à l'allongement de la saison de végétation) ?
#   (3) cette tendance est-elle différente selon les périodes étudiées (2006-2016 vs 2006-2024 vs 2016-2024) ?
#   (4) les tendances dans les décalages phénologiques sont-elles différentes selon les altitudes considérées ?
#   (5) les tendances sont-elles différentes selon les catégories d'observateur·rices (scolaires, professionnels, particuliers) ?



##############################################-
# DATA & APERÇU                           ----
##############################################-

# phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 
phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/data/_CLEANED_data_pheno.csv")
# Données obtenues avec le script 1_mise_en_forme_BDD.R



# pheno_aut10 = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 10%",]
# ggplot(pheno_aut10, aes(x=yearQ, y=julian_day, col=altitude)) +
#   geom_point(size=0.5)+
#   geom_smooth(method="lm", formula=y~x) +
#   facet_wrap(.~ nom_cite)
# 
# 
# pheno_aut50 = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 50%",]
# ggplot(pheno_aut50, aes(x=yearQ, y=julian_day)) +
#   geom_point(size=0.5)+
#   geom_smooth(method="lm", formula=y~x) +
#   facet_wrap(.~ nom_cite)
# 
# 
# # TESTS SUR LE MÉLÈZE
# pheno_aut10_meleze = pheno_aut10[pheno_aut10$species == "Meleze",]
# model_meleze <- brm(. # /!\ bayésien fait bien chauffer l'ordi ! Privilégier d'autres méthodes comme pour le débourrement ?
#   bf(julian_day ~ altitude+yearQ+(yearQ|ID_zone)),
#   data = pheno_aut10_meleze,
#   init = "0",
#   chains = 4, iter = 5000, warmup = 1000,
#   cores = 4
# ) # /!\ PB convergence !!! --> faire plus d'itérations ? Modifier l'initialisation ? Modifier le '0' en faisant 'altitude-1100' comme MB ?
# 
# # plot(model_meleze) # to look at chain mixing
# # summary(model_meleze)
# summary(model_meleze)$fixed
# ## + 2.8 j / décennie
# 
# # Visualisation des effets aléatoires
# # plot_model(model_meleze, type="re")
# # plot_model(model_meleze, type="re", axis.lim=c(-2,2))
# plot_model(model_meleze, type="emm", terms="yearQ")





######################################################################################################-
#   VISUALISATION DE L'INDICE D'AUTOMNE PROPOSÉ PAR COLIN VR ET MARJORIE B                        ----
######################################################################################################-

# DATA 
resume_automne = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/Analyse/Indice_phenoclim/pheno_year_global.csv")

# Visual
ggplot(resume_automne, aes(y=year + 0.2*ifelse(cl_2alt=="Inf1050",-1,1), 
                           x=diff,
                           col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, 
                ymin = year - 0.5, 
                xmin = -Inf, 
                xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_point(shape=15, size=3) + 
  geom_segment(aes(x=-Inf, xend=diff, y=year+ 0.2*ifelse(cl_2alt=="Inf1050",-1,1) , 
                   yend=year+ 0.2*ifelse(cl_2alt=="Inf1050",-1,1), col=cl_2alt ), lty=2, lwd=0.7)+
  scale_color_manual(values=c("darkorange2","gold2")) +
  geom_vline(xintercept = 0, col="black") +
  scale_y_continuous(breaks=seq(2005,2024, by=1), labels=seq(2005,2024, by=1))+
  scale_x_continuous(breaks=seq(-12,12, by=2), labels=seq(-12,12, by=2)) + 
  labs(x="", y="")   +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())


############################################################################################-
# INTÉGRATION DE VARIABLES EXPLICATIVES (principalement climatiques)                     ----
############################################################################################-

# On essaie d'expliquer la sénescence par la date de débourrement, les précipitations, les températures

sites_pheno_Alps = phenoclim %>% filter(nom_massif_v2019 == "Alpes") %>% 
  distinct(id_base_site, species, coord_x_4326, coord_y_4326, coord_x_2154, coord_y_2154,altitude, cl_alt, cl_alt2, cl_alt3, 
           ID_zone, ID_zone, dept1, region1, pays1, circumference, year #new_cat, 
           )

#=============================================================================================================================*
# On associe les données de sénescence aux données de débourrement des mêmes arbres

sites_pheno_Alps = merge(sites_pheno_Alps, 
                         merge(phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Debourrement - Ok 10%",c("id_base_site","julian_day","year")],
                               merge(phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Changement de couleur - Ok 10%",c("id_base_site","julian_day","year")] %>% distinct(id_base_site, year, .keep_all = T),
                                     phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Changement de couleur - Ok 50%",c("id_base_site","julian_day","year")] %>% distinct(id_base_site, year, .keep_all = T), 
                                     by=c("id_base_site","year"), all.x=T, all.y=T),
                               by=c("id_base_site","year"), all.x=T, all.y=T),
                         by=c("id_base_site","year"), all.x=T, all.y=T)
colnames(sites_pheno_Alps)[(ncol(sites_pheno_Alps)-2):ncol(sites_pheno_Alps)] = c("debou10","senes10","senes50")

sites_pheno_Alps = sites_pheno_Alps[!is.na(sites_pheno_Alps$year),]
# on ne garde que les lignes où il y a des infos de sénescence
sites_pheno_Alps = sites_pheno_Alps[!(is.na(sites_pheno_Alps$senes10) & is.na(sites_pheno_Alps$senes50)),]

# dim(sites_pheno_Alps[!is.na(sites_pheno_Alps$debou10) & !is.na(sites_pheno_Alps$senes10) & !is.na(sites_pheno_Alps$senes50),])

#=============================================================================================================================*
# On récupère les précipitations de ERA5-land (Monthly averaged reanalysis)
# --> données mensuelles (précipitations totales en m) de 2004 à 2024, à une résolution de 9km (https://cds.climate.copernicus.eu/requests?tab=all)
# --> données de précipitation --> https://codes.ecmwf.int/grib/param-db/228
GRIB<-brick("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/ERA5land_precip.grib") 
PRECIP = GRIB[[2*(1:252)]]
names(PRECIP) = paste("P",rep(2004:2024, each=12), rep(1:12,21),sep="_")
# # Visualisation des précipitations du mois de juin 2023
# plot(PRECIP[["P_2023_6"]])

# On extrait les précipitations des mois d'été (juin, juillet, août), puis on somme par année
PRECIP_estiv = extract(PRECIP[[grep("_6|_7|_8", names(PRECIP))]], 
                       crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP))))
PRECIP_estiv = t(as.data.frame(t(PRECIP_estiv)) %>% group_by(rep(2004:2024, each=3)) %>% summarise(across(1:nrow(PRECIP_estiv), sum)))
colnames(PRECIP_estiv) = paste0("Psummer",PRECIP_estiv[1,])
sites_pheno_Alps = cbind(sites_pheno_Alps, PRECIP_estiv[-1,])

# # Visualisation
# par(mfrow=c(2,2))
# for (annee in c(2004, 2010, 2023, 2024)){hist(sites_pheno_Alps[,paste0("Psummer",annee)], main=annee, xlim=c(0,0.03), breaks=10)}
# # --> année 2023 particulièrement sèche !

sites_pheno_Alps$P_summer = unlist(apply(sites_pheno_Alps[,c("year",grep("Psummer",colnames(sites_pheno_Alps), value=T))],1,
                                  function(x){x[x[1]-2002]}))
sites_pheno_Alps = sites_pheno_Alps[,-grep("Psummer", colnames(sites_pheno_Alps))]


#=============================================================================================================================*
# On récupère les températures via les reconstructions basées sur les stations météo Phénoclim (cf scripts reconstruction G. Klein)
# On calcule des températures basées sur la date *médiane* de sénescence et la période précédente, et/ou en lien avec la photopériode

# senes10_Alps = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 10%" &
#                            phenoclim$nom_massif == "Alpes",]
# senes10_Alps = senes10_Alps[!is.na(senes10_Alps$ids_observers),] # suppression des lignes remplies de NA
# 
# senes50_Alps = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 50%" &
#                            phenoclim$nom_massif == "Alpes",]
# senes50_Alps = senes50_Alps[!is.na(senes50_Alps$ids_observers),]


# 1) Calcul de la date de sénescence médiane par espèce, sur différentes périodes
dates_med = sites_pheno_Alps %>% group_by(species) %>% 
  summarise(debou10_med_0616 = median(debou10[year %in% 2006:2016], na.rm=T),
            debou10_med_1624 = median(debou10[year %in% 2016:2024], na.rm=T),
            debou10_med_0624 = median(debou10[year %in% 2006:2024], na.rm=T),
            senes10_med_0616 = median(senes10[year %in% 2006:2016], na.rm=T),
            senes10_med_1624 = median(senes10[year %in% 2016:2024], na.rm=T),
            senes10_med_0624 = median(senes10[year %in% 2006:2024], na.rm=T),
            senes50_med_0616 = median(senes50[year %in% 2006:2016], na.rm=T),
            senes50_med_1624 = median(senes50[year %in% 2016:2024], na.rm=T),
            senes50_med_0624 = median(senes50[year %in% 2006:2024], na.rm=T))


# 2) Calcul de la date correspondant à une photopériode seuil, en chaque point où il y a une observation
dates_photoper = phenoclim %>% distinct(id_base_site, coord_y_4326)
dates_photoper$P12.5 = sapply(dates_photoper$coord_y_4326, function(x){max(c(1:366)[photoperiod(1:366,x) > 12.5])})
# /!\ PB : ce seuil de 12.5 h est atteint autour du 16/09 (doy=260), alors que la sénescence a déjà commencé à cette date pour les 5 espèces !
#         => privilégier une date fixe, par exemple le 1/07 (doy=183), comme dans Wu et al 2018 ?


# 3) Calcul des différentes variables de température
# Construction d'un tableau sur la sénescence dans les Alpes, où chaque donnée (site x année) est associée à une valeur de température supposée
# pertinente pour la phénologie (cf liste des variables potentielles ci-dessus)

varTselec = c("Tmoy_GS",                                                          #Zani et al 2020 (growing season = mediane debourr - mediane senes10)
              "Tmoy_MJ","Tmoy_ASO","Tmoy_AS",                                     #Zohner et al 2023
              "Tnight21j_jsenes10", "Tnight30j_jsenes10","Tnight40j_jsenes10",
              "Tnight21j_jsenes50", "Tnight30j_jsenes50","Tnight40j_jsenes50",    #Jolly et al 2005, Wu et al 2018, Zohner et al 2023 (Tnight = moyenne des minimums)
              "GDDinv25_jsenes10","GDDinv20_jsenes10","GDDinv15_jsenes10",
              "GDDinv25_jsenes50","GDDinv20_jsenes50","GDDinv15_jsenes50")        #Delpierre et al 2009 (accumulation de froid, depuis le solstice = j173)
              

output = data.frame(matrix(ncol=4+length(varTselec)))
colnames(output) = c("sitePheno",varTselec, "periode","esp","annee")

vartabT = c("julian_day","station_name","Tair_moy","Tair_min","Tair_max")



for(annee in c(2006:2024)){
  
  print(annee)
  
  reconstruc_Tday_sites = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
  reconstruc_Tday_sites = reconstruc_Tday_sites[!is.na(reconstruc_Tday_sites$station_name),]
  reconstruc_Tday_sites$julian_day = yday(as.Date(reconstruc_Tday_sites$date))
  tabT = rename(reconstruc_Tday_sites, c(jul_day=vartabT[1], sitePheno=vartabT[2], Tmoy=vartabT[3], Tmin=vartabT[4], Tmax=vartabT[5]))
  
  for (esp in unique(sites_pheno_Alps$species[sites_pheno_Alps$year == annee & !(is.na(sites_pheno_Alps$senes10) & is.na(sites_pheno_Alps$senes50))])){

    for (periode in c("med_0616","med_1624","med_0624")){
      med_debou10 = as.numeric(dates_med[dates_med$species==esp,paste0("debou10_", periode)])
      med_senes10 = as.numeric(dates_med[dates_med$species==esp,paste0("senes10_", periode)])
      med_senes50 = as.numeric(dates_med[dates_med$species==esp,paste0("senes50_", periode)])
      
      Tcalc = tabT %>% filter(jul_day <= med_senes10 & jul_day >= med_debou10) %>% group_by(sitePheno) %>% summarise(Tmoy_GS = mean(Tmoy, na.rm=T))
      Tcalc = merge(Tcalc, tabT %>% filter(month==5 | month==6) %>% group_by(sitePheno) %>% summarise(Tmoy_MJ = mean(Tmoy, na.rm=T)), by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(month==8 | month==9) %>% group_by(sitePheno) %>% summarise(Tmoy_AS = mean(Tmoy, na.rm=T)), by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(month==8 | month==9 | month==10) %>% group_by(sitePheno) %>% summarise(Tmoy_ASO = mean(Tmoy, na.rm=T)), by="sitePheno",all.x=T, all.y=T)

      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes10 & jul_day >= med_senes10-21) %>% group_by(sitePheno) %>%
                      summarise(Tnight21j_jsenes10 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes10 & jul_day >= med_senes10-30) %>% group_by(sitePheno) %>%
                      summarise(Tnight30j_jsenes10 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes10 & jul_day >= med_senes10-40) %>% group_by(sitePheno) %>%
                      summarise(Tnight40j_jsenes10 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes50 & jul_day >= med_senes50-21) %>% group_by(sitePheno) %>%
                      summarise(Tnight21j_jsenes50 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes50 & jul_day >= med_senes50-30) %>% group_by(sitePheno) %>%
                      summarise(Tnight30j_jsenes50 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes50 & jul_day >= med_senes50-40) %>% group_by(sitePheno) %>%
                      summarise(Tnight40j_jsenes50 = mean(Tmin, na.rm=T)),  by="sitePheno",all.x=T, all.y=T)
      
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes10 & jul_day >= 173) %>% group_by(sitePheno) %>%
                      summarise(GDDinv25_jsenes10 = sum(max(25-Tmoy, 0), na.rm=T),
                                GDDinv20_jsenes10 = sum(max(20-Tmoy, 0), na.rm=T),
                                GDDinv15_jsenes10 = sum(max(15-Tmoy, 0), na.rm=T)), by="sitePheno",all.x=T, all.y=T)
      Tcalc = merge(Tcalc, tabT %>% filter(jul_day <= med_senes50 & jul_day >= 173) %>% group_by(sitePheno) %>%
                      summarise(GDDinv25_jsenes50 = sum(max(25-Tmoy, 0), na.rm=T),
                                GDDinv20_jsenes50 = sum(max(20-Tmoy, 0), na.rm=T),
                                GDDinv15_jsenes50 = sum(max(15-Tmoy, 0), na.rm=T)), by="sitePheno",all.x=T, all.y=T)
      
      Tcalc$periode = periode
      Tcalc$esp = esp
      Tcalc$annee = annee
      
      output = rbind(output, Tcalc)
      
    }
    
    
  }
  
}

output = output[!is.na(output$sitePheno),]
Tperiodes = data.table::dcast(data=data.table::setDT(output), formula=sitePheno+esp+annee ~ periode, value.var=varTselec)

# ggplot(Tperiodes, aes(y=Tmoy30j_med_0624, x=annee, col=sitePheno)) + geom_point() + geom_smooth(method='lm')

write.csv(Tperiodes, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/varTautomne_allphenosites.csv")



data_Senes_Alps = merge(sites_pheno_Alps, 
                         Tperiodes, by.x=c("id_base_site","species","year"), by.y=c("sitePheno","esp","annee"), all.x=T, all.y=F)

write.csv(data_Senes_Alps, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Senes_Alps_T.csv")



##############################################-
# MODÈLES DE SÉNESCENCE                    ----
##############################################-

senes_Alps_all = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/data_Senes_Alps_T.csv", row.names=1)
senes_Alps_all$yearQ = factor(senes_Alps_all$year)
senes_Alps_all$cl_alt = factor(senes_Alps_all$cl_alt, levels=c("150-450","450-750" ,  "750-1050"   ,"1050-1350" , "1350-1650" ,"1650-1950", "1950-2250" ), 
                             ordered = T)

# On crée une variable associée à la durée de la sénescence 
# (/!\ ce n'est pas vraiment une durée puisqu'on regarde le nombre de jours entre 10 et 50% de changement de couleur)
senes_Alps_all$senesduree = senes_Alps_all$senes50 - senes_Alps_all$senes10

# On renomme les variables de température calculées sur la date médiane de débourrement de la période 2006-2024 ('par défaut')
senes_Alps_all = senes_Alps_all %>% rename(Tmoy_GS = Tmoy_GS_med_0624,
                                   Tmoy_MJ = Tmoy_MJ_med_0624,
                                   Tmoy_AS = Tmoy_AS_med_0624,
                                   Tmoy_ASO = Tmoy_ASO_med_0624,
                                   Tnight21j_jsenes10 = Tnight21j_jsenes10_med_0624,
                                   Tnight30j_jsenes10 = Tnight30j_jsenes10_med_0624,
                                   Tnight40j_jsenes10 = Tnight40j_jsenes10_med_0624,
                                   Tnight21j_jsenes50 = Tnight21j_jsenes50_med_0624,
                                   Tnight30j_jsenes50 = Tnight30j_jsenes50_med_0624,
                                   Tnight40j_jsenes50 = Tnight40j_jsenes50_med_0624,
                                   GDDinv25_jsenes10 = GDDinv25_jsenes10_med_0624,
                                   GDDinv20_jsenes10 = GDDinv20_jsenes10_med_0624,
                                   GDDinv15_jsenes10 = GDDinv15_jsenes10_med_0624,
                                   GDDinv25_jsenes50 = GDDinv25_jsenes50_med_0624,
                                   GDDinv20_jsenes50 = GDDinv20_jsenes50_med_0624,
                                   GDDinv15_jsenes50 = GDDinv15_jsenes50_med_0624)

# Remarque : au vu des corrélations entre ces variables de température, il est peu logique d'en intégrer plusieurs si elles sont trop corrélées
corrplot::corrplot(cor(na.omit(senes_Alps_all[,varTselec])))


#-----------------------------------------------------------------------
#*-- Aperçu de la localisation des données phénologiques automnales ----
senes_Alps_teledec = senes_Alps_all %>% group_by(id_base_site) %>% mutate(nb_years = length(unique(year))) %>% ungroup()
senes_Alps_leaflet = vect(senes_Alps_teledec, geom=c("coord_x_4326","coord_y_4326"), "epsg:4326")
  
colors = data.frame(ID = c("Bouleau_pubescent","Bouleau_verruqueux","Hetre","Meleze","Sorbier"), 
                    col = c("darkgreen","lightgreen","brown","yellow","orange"))
# carto = leaflet() %>% #senes_Alps_leaflet[senes_Alps_leaflet$year >= 2017,]) %>%
#   addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>% 
#   addCircleMarkers(lng = jitter(senes_Alps_teledec$coord_x_4326[senes_Alps_teledec$year >= 2017], factor = 2), lat = jitter(senes_Alps_teledec$coord_y_4326[senes_Alps_teledec$year >= 2017], factor = 2), # si on veut limiter la superposition : jitter
#                    label = paste(senes_Alps_teledec$id_base_site[senes_Alps_teledec$year >= 2017], " : ",senes_Alps_teledec$nb_years[senes_Alps_teledec$year >= 2017], " obs.", sep=""),
#                    fillColor =as.character(factor(senes_Alps_teledec$species,
#                                                   levels=colors$ID,
#                                                   labels=colors$col)), 
#                    color="black",opacity=1,fillOpacity = 1, weight=0.1, radius=5) %>%
#   addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Espece", opacity = 1)
# carto = leaflet(senes_Alps_leaflet) %>%
#   addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>% 
#   addCircleMarkers(label = paste(senes_Alps_leaflet$id_base_site, " : ",senes_Alps_leaflet$nb_years, " obs.", sep=""),
#                    fillColor =as.character(factor(senes_Alps_leaflet$species,
#                                                   levels=colors$ID,
#                                                   labels=colors$col)), 
#                    color="black",opacity=1,fillOpacity = 1, weight=0.1, radius=5) %>%
#   addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Espece", opacity = 1)
#   
# carto
# 
# Pour le croisement avec les données de télédétection (cf projet IUT Digne), on peut croiser les observations Phénoclim et les mailles / la couverture des données 
# télédétection
# Les données de télédétection Pléiades couvrent des tuiles de 20km --> on identifie ainsi les zones d'une vingtaine de kilomètres où on a suffisamment de données
# Phénoclim
# rast20km = rast(xmin=860000, xmax=1080000, ymin=6300000, ymax=6600000, resolution = 20000, crs="epsg:2154")
# values(rast20km) = 1:length(values(rast20km))
# names(rast20km) = "IDcell"
# senes_Alps_teledec$IDcell = extract(rast20km, vect(senes_Alps_teledec, geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"))[,"IDcell"]
senes_Alps_teledec$dept1[is.na(senes_Alps_teledec$dept1)] = senes_Alps_teledec$region1[is.na(senes_Alps_teledec$dept1)]

# recap = senes_Alps_teledec %>% group_by(IDcell, species) %>% summarise(nb_arbres = length(unique(id_base_site)),
#                                                                        nb_obssenes = length(id_base_site),
#                                                                        nb_annees = length(unique(year)),
#                                                                        anneemin_obssenes = min(year),
#                                                                        anneemax_obssenes = max(year),
#                                                                        dept = paste(unique(dept1), collapse = ", "))
# 
# recap = senes_Alps_teledec %>% group_by(IDcell, species) %>% 
#                                 summarise(nb_arbres = length(unique(id_base_site)),
#                                           nb_obssenes = length(id_base_site),
#                                           nb_annees = length(unique(year)),
#                                           nb_annees_par_arbre_moy = mean(nb_years),
#                                           anneemin_obssenes = min(year),
#                                           anneemax_obssenes = max(year),
#                                           dept = paste(unique(dept1), collapse = ", "))
# 
# senes_Alps_teledec$selection = apply(senes_Alps_teledec[,c("IDcell","species")], 1,
#                                      function(X, rec=recap){rec = rec[rec$species == X[2],]
#                                      return(ifelse(X[1] %in% rec$IDcell[rec$nb_arbres >=5 & rec$nb_obssenes >=10], X[2],"-"))})
# senes_Alps_teledec$selection = ifelse(apply(senes_Alps_teledec[,c("IDcell","species")], 1,
#                                      function(X, rec=recap){rec = rec[rec$species == X[2],]
#                                      return(ifelse(X[1] %in% rec$IDcell[rec$nb_arbres >=5 & rec$nb_obssenes >=10], X[2],"-"))}) == senes_Alps_teledec$species, "oui","non")
# 
# 
# senes_Alps_leaflet = vect(senes_Alps_teledec, geom=c("coord_x_4326","coord_y_4326"), "epsg:4326")
# 
# # carto_selec = leaflet(senes_Alps_leaflet[grep("oui",senes_Alps_leaflet$selection),]) %>%
# carto_selec = leaflet(senes_Alps_leaflet) %>%
#   addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>%
#   addCircleMarkers(label = paste("nbyr: ",senes_Alps_leaflet$nb_years,"cell: ",senes_Alps_leaflet$IDcell),#paste(senes_Alps_leaflet$id_base_site[grep("oui",senes_Alps_leaflet$selection)], " : ",senes_Alps_leaflet$nb_years[grep("oui",senes_Alps_leaflet$selection)], " obs.", sep=""),
#                    fillColor =as.character(factor(senes_Alps_leaflet$species,#[grep("oui",senes_Alps_leaflet$selection)],
#                                                   levels=colors$ID,
#                                                   labels=colors$col)),
#                    radius = as.character(factor(senes_Alps_leaflet$nb_years >5,
#                                                 levels=c(T,F),
#                                                 labels=c(5,2))),
#                    fillOpacity = as.character(factor(senes_Alps_leaflet$IDcell %in% c(29,30,53,54, 127,128,138,139, 94,115,116),
#                                                 levels=c(T,F),
#                                                 labels=c(1,0.7))),
#                    color="black",opacity=1, weight=0.1) %>%
#   addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Espece", opacity = 1)
# 
# carto_selec
  


recap2 = senes_Alps_teledec %>% group_by(ID_zone, species) %>% 
  summarise(x = mean(coord_x_4326), y = mean(coord_y_4326),
            nb_arbres = length(unique(id_base_site)),
            nb_obssenes = length(id_base_site),
            nb_annees = length(unique(year)),
            nb_annees_par_arbre_moy = mean(nb_years),
            anneemin_obssenes = min(year),
            anneemax_obssenes = max(year),
            dept = paste(unique(dept1), collapse = ", "))
leaflet(recap2) %>%
  addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>%
  addCircleMarkers(lng=~x, lat=~y,
                   label = paste("IDzone: ", recap2$ID_zone, ", nbobs: ",recap2$nb_obssenes,", nbyears: ",recap2$nb_annees_par_arbre_moy),
                   fillColor =as.character(factor(recap2$species,#[grep("oui",senes_Alps_leaflet$selection)],
                                                  levels=colors$ID,
                                                  labels=colors$col)),
                   radius = ~nb_obssenes*0.1,
                   fillOpacity = 1,
                   color="black",opacity=1, weight=0.1) %>%
  addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Espece", opacity = 1)

# Les ID_zone intéressants sont :
# - Chamonix : zone692
# - Gran Paradisio E : zone782, zone783, zone784
# - Gran Paradisio W : zone750 (nombreuses zones alentours, mais avec moins d'observations)
# - Vanoise : zone521
# - Mercantour : zone687
# - Vercors : zone279
# - Ecrins : zone554, zone552, zone615
# - Digne-Auzet : zone512
# On crée un buffer de 10km autour de ces zones, pour capter les points supplémentaires, où il y a moins de données mais où ça peut apporter des points complémentaires
senes_Alps_teledec$selec = as.character(factor(senes_Alps_teledec$ID_zone,
                                  levels = c("zone692", "zone782", "zone783", "zone784", "zone750", "zone521","zone687","zone279","zone554","zone552","zone615", "zone512"),
                                  labels = c("Chamonix",rep("GranParadisioE",3), "GranParadisioW", "Vanoise","Mercantour","Vercors",rep("Ecrins",3),"Digne-Auzet")))
vectselec = aggregate(buffer(vect(senes_Alps_teledec[!is.na(senes_Alps_teledec$selec),], geom=c("coord_x_2154","coord_y_2154"), "epsg:2154"), 10000))
senes_Alps_teledec$selec[is.na(senes_Alps_teledec$selec)] = ifelse(as.vector(extract(vectselec, senes_Alps_teledec[is.na(senes_Alps_teledec$selec),c("coord_x_2154","coord_y_2154")])[,2])==1,"extra",NA)

# Si c'est une autre espèce que mélèze, bouleau verruqueux ou sorbier, on les met en "extra"
senes_Alps_teledec$selec[!is.na(senes_Alps_teledec$selec) & !(senes_Alps_teledec$species %in% c("Bouleau_verruqueux","Sorbier","Meleze"))] = "extra"


write.csv(senes_Alps_teledec[!is.na(senes_Alps_teledec$selec),c(1:19, 71,72)], "/Users/ninonfontaine/Desktop/projetsR/TEST/data/PhenoClim/data_projetIUT_senes.csv")
writeVector(vectselec, "/Users/ninonfontaine/Desktop/projetsR/TEST/data/PhenoClim/data_projetIUT_zones_senes.shp")

leaflet(vect(senes_Alps_teledec[!is.na(senes_Alps_teledec$selec),], geom=c("coord_x_4326","coord_y_4326"), "epsg:4326")) %>%
  addTiles() %>% #'http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}')%>%
  addCircleMarkers(label = paste("IDzone: ", senes_Alps_teledec$ID_zone),#, ", nbobs: ",recap2$nb_obssenes,", nbyears: ",recap2$nb_annees_par_arbre_moy),
                   fillColor =as.character(factor(senes_Alps_teledec$species[!is.na(senes_Alps_teledec$selec)],#[grep("oui",senes_Alps_leaflet$selection)],
                                                  levels=colors$ID,
                                                  labels=colors$col)),
                   radius = ifelse(senes_Alps_teledec$selec[!is.na(senes_Alps_teledec$selec)]=="extra",3,7),
                   fillOpacity = 1,
                   color="black",opacity=1, weight=0.1) %>%
  addLegend(position="bottomright",colors=colors$col, labels = colors$ID, title = "Espece", opacity = 1)



##############################################-
# *-- ONSET DE SÉNESCENCE                  ----
##############################################-

senes10_Alps = senes_Alps_all[!is.na(senes_Alps_all$senes10),]

# # Aperçu des données de sénescence
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=altitude, y=senes10, col=species)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes10, col=species)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes10, col=yearQ)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=GDDinv20_jsenes10, y=senes10, col=yearQ)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=P_summer, y=senes10)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes10_Alps[!is.na(senes10_Alps$Tmoy_GS),], aes(x=debou10, y=senes10)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)



# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2024, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes10_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()


##############################################-
#*---- Bouleau verruqueux ----

senes10_Bpen = senes10_Alps[senes10_Alps$species == "Bouleau_verruqueux",]
# # on retire les lignes où on n'a pas de données de température (2005)
senes10_Bpen = senes10_Bpen[!is.na(senes10_Bpen$Tmoy_GS),]
ggplot(senes10_Bpen, aes(x=senes10)) + geom_histogram(binwidth=5)
# On retire une valeur extrême (sénescence notée à 78jours = 18/03 !)
senes10_Bpen = senes10_Bpen[senes10_Bpen$senes10 > 100,]
# ggplot(senes10_Bpen, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Changement de couleur 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Bpen_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Bpen, REML=F)
summary(mod_senes10_Bpen_Altyear)
altyearmods = c(list("Bouleau_verruqueux"=mod_senes10_Bpen_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes10_Bpen_Altyear)
# + Validation croisée
n = nrow(senes10_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Bpen[trainIndex ,]
test <- senes10_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Bpen$ID_zone[drop=T])[!unique(senes10_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes10_Bpen_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_altyear_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes10_Bpen_all = senes10_Bpen
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes10_Bpen = senes10_Bpen_all# senes10_Bpen[!is.na(senes10_Bpen$debou10),] #senes10_Bpen_all

mod_senes10_Bpen_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) 
mod_senes10_Bpen_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)

AIC(mod_senes10_Bpen_debou10, 
    mod_senes10_Bpen_P_summer,mod_senes10_Bpen_Tmoy_GS,mod_senes10_Bpen_Tmoy_MJ,mod_senes10_Bpen_Tmoy_AS,mod_senes10_Bpen_Tmoy_ASO,
    mod_senes10_Bpen_Tnight21j_jsenes10,mod_senes10_Bpen_Tnight30j_jsenes10,mod_senes10_Bpen_Tnight40j_jsenes10,
    mod_senes10_Bpen_GDDinv25_jsenes10,mod_senes10_Bpen_GDDinv20_jsenes10,mod_senes10_Bpen_GDDinv15_jsenes10)
# Les variables de GDDinverse à 25 et 20°C donnent les meilleurs résultats (AIC = 13502)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                               (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
mod_senes10_Bpen_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
                               (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AIC = 13422.9
# visreg(mod_senes10_Bpen_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes10_Bpen_multivar = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + 
                                 (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AIC = 13417.5
# visreg(mod_senes10_Bpen_multivar, "GDDinv25_jsenes10", by="P_summer")
visreg(mod_senes10_Bpen_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes10_Bpen = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + altitude +
                                 (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AIC = 13411.8
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes10_Bpen), bestmods)


summary(bestmod_senes10_Bpen)
# gradient altitudinal : sénescence 0.9 jour plus tôt quand on monte de 100m
# précipitations : sénescence 2.3 jours plus tôt quand on gagne 1mm de pluie 
# accumulation de "froid" : sénescence 2.2 jours plus tôt quand on gagne 1°C de froid // effet d'avancée moins marqué s'il y a plus de pluie l'été
qqnorm(resid(bestmod_senes10_Bpen))
qqline(resid(bestmod_senes10_Bpen))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Bpen)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Bpen)
# - Validation croisée
n = nrow(senes10_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Bpen[trainIndex ,]
test <- senes10_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Bpen$ID_zone[drop=T])[!unique(senes10_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Bpen), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        variable = rownames(coef(summary(bestmod_senes10_Bpen))),
                                        coef = coef(summary(bestmod_senes10_Bpen))[,1],
                                        std = coef(summary(bestmod_senes10_Bpen))[,2],
                                        pval = coef(summary(bestmod_senes10_Bpen))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Bpen)[rownames(coef(summary(bestmod_senes10_Bpen)))], error = function(e) return(NA))))



##############################################-
#*---- Meleze ----

senes10_Ldec = senes10_Alps[senes10_Alps$species == "Meleze",]
# # on retire les lignes où on n'a pas de données de température (2005)
senes10_Ldec = senes10_Ldec[!is.na(senes10_Ldec$Tmoy_GS),]
ggplot(senes10_Ldec, aes(x=senes10)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Ldec_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Ldec, REML=F)
summary(mod_senes10_Ldec_Altyear)
altyearmods = c(list("Meleze"=mod_senes10_Ldec_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes10_Ldec_Altyear)
# + Validation croisée
n = nrow(senes10_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Ldec[trainIndex ,]
test <- senes10_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Ldec$ID_zone[drop=T])[!unique(senes10_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes10_Ldec_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_altyear_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes10_Ldec_all = senes10_Ldec
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes10_Ldec = senes10_Ldec_all#senes10_Ldec[!is.na(senes10_Ldec$debou10),] #senes10_Ldec_all

mod_senes10_Ldec_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) 
mod_senes10_Ldec_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)

AIC(mod_senes10_Ldec_debou10, 
    mod_senes10_Ldec_P_summer,mod_senes10_Ldec_Tmoy_GS,mod_senes10_Ldec_Tmoy_MJ,mod_senes10_Ldec_Tmoy_AS,mod_senes10_Ldec_Tmoy_ASO,
    mod_senes10_Ldec_Tnight21j_jsenes10,mod_senes10_Ldec_Tnight30j_jsenes10,mod_senes10_Ldec_Tnight40j_jsenes10,
    mod_senes10_Ldec_GDDinv25_jsenes10,mod_senes10_Ldec_GDDinv20_jsenes10,mod_senes10_Ldec_GDDinv15_jsenes10)
# Les températures moyennes estivales (août-septembre) donnent les meilleurs résultats (AIC = 13870.48)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on reste sur le modèle à une variable :
mod_senes10_Ldec_multivar = lmer(senes10 ~ Tmoy_GS + 
                                 (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) # AIC = 13835.2
# visreg(mod_senes10_Ldec_multivar, "Tmoy_GS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes10_Ldec = lmer(senes10 ~ altitude  +
                            (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) # AIC = 13840.1
# Finalement l'altitude donne de meilleurs scores que les variables climatiques...!
bestmods = c(list("Meleze"=bestmod_senes10_Ldec), bestmods)


summary(bestmod_senes10_Ldec)
# Sénescence qui commence 2.4 jours plus tôt par 100m d'altitude
qqnorm(resid(bestmod_senes10_Ldec))
qqline(resid(bestmod_senes10_Ldec))
# /!\ résidus !!!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Ldec)
# R2 des effets fixes plutôt nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Ldec)
# - Validation croisée
n = nrow(senes10_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Ldec[trainIndex ,]
test <- senes10_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Ldec$ID_zone[drop=T])[!unique(senes10_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Ldec), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        variable = rownames(coef(summary(bestmod_senes10_Ldec))),
                                        coef = coef(summary(bestmod_senes10_Ldec))[,1],
                                        std = coef(summary(bestmod_senes10_Ldec))[,2],
                                        pval = coef(summary(bestmod_senes10_Ldec))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Ldec)[rownames(coef(summary(bestmod_senes10_Ldec)))], error = function(e) return(NA))))



##############################################-
#*---- Bouleau pubescent ----

senes10_Bpub = senes10_Alps[senes10_Alps$species == "Bouleau_pubescent",]
# ggplot(senes10_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes10_Bpub = senes10_Bpub[!is.na(senes10_Bpub$Tmoy_GS),]
ggplot(senes10_Bpub, aes(x=senes10)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Bpub_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Bpub, REML=F)
summary(mod_senes10_Bpub_Altyear)
altyearmods = c(list("Bouleau_pubescent"=mod_senes10_Bpub_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes10_Bpub_Altyear)
# + Validation croisée
n = nrow(senes10_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Bpub[trainIndex ,]
test <- senes10_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Bpub$ID_zone[drop=T])[!unique(senes10_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes10_Bpub_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_altyear_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes10_Bpub_all = senes10_Bpub
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes10_Bpub = senes10_Bpub_all# senes10_Bpub[!is.na(senes10_Bpub$debou10),] #senes10_Bpub_all

mod_senes10_Bpub_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) 
mod_senes10_Bpub_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)

AIC(mod_senes10_Bpub_debou10, 
    mod_senes10_Bpub_P_summer,mod_senes10_Bpub_Tmoy_GS,mod_senes10_Bpub_Tmoy_MJ,mod_senes10_Bpub_Tmoy_AS,mod_senes10_Bpub_Tmoy_ASO,
    mod_senes10_Bpub_Tnight21j_jsenes10,mod_senes10_Bpub_Tnight30j_jsenes10,mod_senes10_Bpub_Tnight40j_jsenes10,
    mod_senes10_Bpub_GDDinv25_jsenes10,mod_senes10_Bpub_GDDinv20_jsenes10,mod_senes10_Bpub_GDDinv15_jsenes10)
# Les variables de GDDinverse à 25 et 20°C donnent les meilleurs résultats (AIC = 646.6667)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F), corr=F) # /!\ en intégrant debou10 on divise le jeu de données par 2 ! (de 80 à 40 points)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes10_Bpub_multivar = lmer(senes10 ~ GDDinv25_jsenes10 + 
                                 (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AIC = 615.6
# visreg(mod_senes10_Bpub_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes10_Bpub_multivar = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + 
                                 (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AIC = 611.2
# visreg(mod_senes10_Bpub_multivar, "GDDinv25_jsenes10", by="P_summer")
visreg(mod_senes10_Bpub_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes10_Bpub = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + altitude +
                            (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AIC = 606.6
bestmods = c(list("Bouleau_pubescent"=bestmod_senes10_Bpub), bestmods)
visreg(bestmod_senes10_Bpub, "P_summer", by="GDDinv25_jsenes10")
visreg(bestmod_senes10_Bpub, "GDDinv25_jsenes10", by="P_summer")


summary(bestmod_senes10_Bpub)
# gradient altitudinal : sénescence 2.0 jours plus tôt quand on monte de 100m
# précipitations : sénescence 5.8 jours plus tôt quand on gagne 1mm de pluie
# accumulation de "froid" : sénescence 6.2 jours plus tôt quand on gagne 1°C de froid (précocité d'autant plus forte qu'il fait sec)

qqnorm(resid(bestmod_senes10_Bpub))
qqline(resid(bestmod_senes10_Bpub))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Bpub)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Bpub)
# - Validation croisée
n = nrow(senes10_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Bpub[trainIndex ,]
test <- senes10_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Bpub$ID_zone[drop=T])[!unique(senes10_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Bpub), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        variable = rownames(coef(summary(bestmod_senes10_Bpub))),
                                        coef = coef(summary(bestmod_senes10_Bpub))[,1],
                                        std = coef(summary(bestmod_senes10_Bpub))[,2],
                                        pval = coef(summary(bestmod_senes10_Bpub))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Bpub)[rownames(coef(summary(bestmod_senes10_Bpub)))], error = function(e) return(NA))))



##############################################-
#*---- Hetre ----

senes10_Fsyl = senes10_Alps[senes10_Alps$species == "Hetre",]
# ggplot(senes10_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes10_Fsyl = senes10_Fsyl[!is.na(senes10_Fsyl$Tmoy_GS),]
ggplot(senes10_Fsyl, aes(x=senes10)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Fsyl_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Fsyl, REML=F)
summary(mod_senes10_Fsyl_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes10_Fsyl_Altyear)
# + Validation croisée
n = nrow(senes10_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Fsyl[trainIndex ,]
test <- senes10_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Fsyl$ID_zone[drop=T])[!unique(senes10_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes10_Fsyl_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Hetre", "R2_altyear_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes10_Fsyl_all = senes10_Fsyl
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes10_Fsyl = senes10_Fsyl_all# senes10_Fsyl[!is.na(senes10_Fsyl$debou10),] #senes10_Fsyl_all

mod_senes10_Fsyl_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) 
mod_senes10_Fsyl_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)

AIC(mod_senes10_Fsyl_debou10, 
    mod_senes10_Fsyl_P_summer,mod_senes10_Fsyl_Tmoy_GS,mod_senes10_Fsyl_Tmoy_MJ,mod_senes10_Fsyl_Tmoy_AS,mod_senes10_Fsyl_Tmoy_ASO,
    mod_senes10_Fsyl_Tnight21j_jsenes10,mod_senes10_Fsyl_Tnight30j_jsenes10,mod_senes10_Fsyl_Tnight40j_jsenes10,
    mod_senes10_Fsyl_GDDinv25_jsenes10,mod_senes10_Fsyl_GDDinv20_jsenes10,mod_senes10_Fsyl_GDDinv15_jsenes10)
# Les variables de température estivales donnent les meilleurs résultats (AIC = 528...)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes10_Fsyl_multivar = lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_ASO + 
                                 Tnight21j_jsenes10 + 
                                 GDDinv25_jsenes10 +
                                 (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AIC = 486.3
# visreg(mod_senes10_Fsyl_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes10_Fsyl_multivar = lmer(senes10 ~ P_summer * Tmoy_AS + 
                                 (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AIC = 481.9
# visreg(mod_senes10_Fsyl_multivar, "Tmoy_AS", by="P_summer")
# visreg(mod_senes10_Fsyl_multivar, "P_summer", by="Tmoy_AS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes10_Fsyl = lmer(senes10 ~ P_summer * Tmoy_AS + 
                            (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AIC = 481.9
bestmods = c(list("Hetre"=bestmod_senes10_Fsyl), bestmods)
# visreg(bestmod_senes10_Fsyl, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes10_Fsyl, "Tmoy_AS", by="P_summer")


summary(bestmod_senes10_Fsyl)
# précipitations : sénescence 28 jours plus tôt quand on gagne 1mm de pluie
# chaleur de fin de saison (août-septembre) : sénescence 18 jours plus tard quand on gagne 1°C (retard d'autant plus fort qu'il fait sec)
# /!\ ça me semble beaucoup ces chiffres !!!

qqnorm(resid(bestmod_senes10_Fsyl))
qqline(resid(bestmod_senes10_Fsyl))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Fsyl)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Fsyl)
# - Validation croisée
n = nrow(senes10_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Fsyl[trainIndex ,]
test <- senes10_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Fsyl$ID_zone[drop=T])[!unique(senes10_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Fsyl), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        variable = rownames(coef(summary(bestmod_senes10_Fsyl))),
                                        coef = coef(summary(bestmod_senes10_Fsyl))[,1],
                                        std = coef(summary(bestmod_senes10_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senes10_Fsyl))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Fsyl)[rownames(coef(summary(bestmod_senes10_Fsyl)))], error = function(e) return(NA))))



##############################################-
#*---- Sorbier ----

senes10_Sacu = senes10_Alps[senes10_Alps$species == "Sorbier",]
# ggplot(senes10_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes10_Sacu = senes10_Sacu[!is.na(senes10_Sacu$Tmoy_GS),]
ggplot(senes10_Sacu, aes(x=senes10)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Sacu_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Sacu, REML=F)
summary(mod_senes10_Sacu_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes10_Sacu_Altyear)
# + Validation croisée
n = nrow(senes10_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Sacu[trainIndex ,]
test <- senes10_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Sacu$ID_zone[drop=T])[!unique(senes10_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes10_Sacu_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_altyear_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes10_Sacu_all = senes10_Sacu
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes10_Sacu = senes10_Sacu_all# senes10_Sacu[!is.na(senes10_Sacu$debou10),] #senes10_Sacu_all

mod_senes10_Sacu_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) 
mod_senes10_Sacu_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)

AIC(mod_senes10_Sacu_debou10, 
    mod_senes10_Sacu_P_summer,mod_senes10_Sacu_Tmoy_GS,mod_senes10_Sacu_Tmoy_MJ,mod_senes10_Sacu_Tmoy_AS,mod_senes10_Sacu_Tmoy_ASO,
    mod_senes10_Sacu_Tnight21j_jsenes10,mod_senes10_Sacu_Tnight30j_jsenes10,mod_senes10_Sacu_Tnight40j_jsenes10,
    mod_senes10_Sacu_GDDinv25_jsenes10,mod_senes10_Sacu_GDDinv20_jsenes10,mod_senes10_Sacu_GDDinv15_jsenes10)
# Les variables de température nocturne donnent les meilleurs résultats (AIC = 4940.168)...
# ... ainsi que la variable de débourrement !
# on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes10_Sacu_multivar = lmer(senes10 ~ debou10 + 
                                 Tnight21j_jsenes10 + 
                                 GDDinv25_jsenes10 +
                                 (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) # AIC = 4035.9 # /!\ en ayant ajouté debou10, on tronque la base de données de 20% !
# visreg(mod_senes10_Sacu_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions (ici ça n'améliore rien) 


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes10_Sacu = lmer(senes10 ~ debou10 + 
                            Tnight21j_jsenes10 + 
                            GDDinv25_jsenes10 +
                            (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) # AIC = 4035.9
# Ici l'altitude n'améliore rien, MAIS en mettant la date de débourrement en variable explicative on tronque la base de données
bestmods = c(list("Sorbier"=bestmod_senes10_Sacu), bestmods)
# visreg(bestmod_senes10_Sacu, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes10_Sacu, "Tmoy_AS", by="P_summer")


summary(bestmod_senes10_Sacu)
# températures nocturnes un peu avant la sénescence : sénescence 1.5 jours plus tard quand on gagne 1°C en moyenne la nuit
# accumulation de froid (pas vraiment significatif !) : sénescence 1.0 jour plus tard quand on gagne 1°C de froid
# quand le débourrement est plus tardif d'une semaine, la sénescence est retardée d'1 jour seulement

qqnorm(resid(bestmod_senes10_Sacu))
qqline(resid(bestmod_senes10_Sacu))
# résidus quasi ok !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Sacu)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Sacu)
# - Validation croisée
n = nrow(senes10_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Sacu[trainIndex ,]
test <- senes10_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Sacu$ID_zone[drop=T])[!unique(senes10_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Sacu), train)
predictions <- mod_train %>% predict(test1) # BUG !!!
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        variable = rownames(coef(summary(bestmod_senes10_Sacu))),
                                        coef = coef(summary(bestmod_senes10_Sacu))[,1],
                                        std = coef(summary(bestmod_senes10_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes10_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Sacu)[rownames(coef(summary(bestmod_senes10_Sacu)))], error = function(e) return(NA))))



write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes10.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senes10_v1altyear.Rdata")




##############################################-
# *-- 50% DE SÉNESCENCE                  ----
##############################################-

senes50_Alps = senes_Alps_all[!is.na(senes_Alps_all$senes50),]

# # Aperçu des données de sénescence
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=altitude, y=senes50, col=species)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes50, col=species)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes50, col=yearQ)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=GDDinv20_jsenes50, y=senes50, col=yearQ)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=P_summer, y=senes50)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes50_Alps[!is.na(senes50_Alps$Tmoy_GS),], aes(x=debou10, y=senes50)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)



# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2024, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes50_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()


##############################################-
#*---- Bouleau verruqueux ----

senes50_Bpen = senes50_Alps[senes50_Alps$species == "Bouleau_verruqueux",]
# # on retire les lignes où on n'a pas de données de température (2005)
senes50_Bpen = senes50_Bpen[!is.na(senes50_Bpen$Tmoy_GS),]
ggplot(senes50_Bpen, aes(x=senes50)) + geom_histogram(binwidth=5)
# On retire une valeur extrême (sénescence notée à 75jours = 15/03 !)
senes50_Bpen = senes50_Bpen[senes50_Bpen$senes50 > 100,]
# ggplot(senes50_Bpen, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Changement de couleur 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Bpen_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Bpen, REML=F)
summary(mod_senes50_Bpen_Altyear)
altyearmods = c(list("Bouleau_verruqueux"=mod_senes50_Bpen_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes50_Bpen_Altyear)
# + Validation croisée
n = nrow(senes50_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Bpen[trainIndex ,]
test <- senes50_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Bpen$ID_zone[drop=T])[!unique(senes50_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes50_Bpen_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_altyear_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes50_Bpen_all = senes50_Bpen
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes50_Bpen = senes50_Bpen_all# senes50_Bpen[!is.na(senes50_Bpen$debou10),] #senes50_Bpen_all

mod_senes50_Bpen_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) 
mod_senes50_Bpen_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tmoy_GS <- lmer(senes50 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tmoy_MJ <- lmer(senes50 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tmoy_AS <- lmer(senes50 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tmoy_ASO <- lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tnight21j_jsenes50 <- lmer(senes50 ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tnight30j_jsenes50 <- lmer(senes50 ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_Tnight40j_jsenes50 <- lmer(senes50 ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_GDDinv25_jsenes50 <- lmer(senes50 ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_GDDinv20_jsenes50 <- lmer(senes50 ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_GDDinv15_jsenes50 <- lmer(senes50 ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)

AIC(mod_senes50_Bpen_debou10, 
    mod_senes50_Bpen_P_summer,mod_senes50_Bpen_Tmoy_GS,mod_senes50_Bpen_Tmoy_MJ,mod_senes50_Bpen_Tmoy_AS,mod_senes50_Bpen_Tmoy_ASO,
    mod_senes50_Bpen_Tnight21j_jsenes50,mod_senes50_Bpen_Tnight30j_jsenes50,mod_senes50_Bpen_Tnight40j_jsenes50,
    mod_senes50_Bpen_GDDinv25_jsenes50,mod_senes50_Bpen_GDDinv20_jsenes50,mod_senes50_Bpen_GDDinv15_jsenes50)
# Les variables de GDDinverse à 25 et 20°C donnent les meilleurs résultats (AIC = 13047.74)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
mod_senes50_Bpen_multivar = lmer(senes50 ~ GDDinv25_jsenes50 + 
                                   (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AIC = 13005.7
# visreg(mod_senes50_Bpen_multivar, "GDDinv25_jsenes50")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senes50_Bpen_multivar = lmer(senes50 ~ P_summer * GDDinv25_jsenes50 + 
#                                    (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AIC = 13491.7
# # visreg(mod_senes50_Bpen_multivar, "GDDinv25_jsenes50", by="P_summer")
# visreg(mod_senes50_Bpen_multivar, "P_summer", by="GDDinv25_jsenes50")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes50_Bpen = lmer(senes50 ~ altitude + GDDinv25_jsenes50 + 
                              (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AIC = 12996.6
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes50_Bpen), bestmods)


summary(bestmod_senes50_Bpen)
# gradient altitudinal : sénescence 1.0 jour plus tôt quand on monte de 100m
# accumulation de "froid" : sénescence 0.9 jours plus tôt quand on gagne 1°C de froid 
qqnorm(resid(bestmod_senes50_Bpen))
qqline(resid(bestmod_senes50_Bpen))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Bpen)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Bpen)
# - Validation croisée
n = nrow(senes50_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Bpen[trainIndex ,]
test <- senes50_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Bpen$ID_zone[drop=T])[!unique(senes50_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Bpen), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        variable = rownames(coef(summary(bestmod_senes50_Bpen))),
                                        coef = coef(summary(bestmod_senes50_Bpen))[,1],
                                        std = coef(summary(bestmod_senes50_Bpen))[,2],
                                        pval = coef(summary(bestmod_senes50_Bpen))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Bpen)[rownames(coef(summary(bestmod_senes50_Bpen)))], error = function(e) return(NA))))



##############################################-
#*---- Meleze ----

senes50_Ldec = senes50_Alps[senes50_Alps$species == "Meleze",]
# # on retire les lignes où on n'a pas de données de température (2005)
senes50_Ldec = senes50_Ldec[!is.na(senes50_Ldec$Tmoy_GS),]
ggplot(senes50_Ldec, aes(x=senes50)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Ldec_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Ldec, REML=F)
summary(mod_senes50_Ldec_Altyear)
altyearmods = c(list("Meleze"=mod_senes50_Ldec_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes50_Ldec_Altyear)
# + Validation croisée
n = nrow(senes50_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Ldec[trainIndex ,]
test <- senes50_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Ldec$ID_zone[drop=T])[!unique(senes50_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes50_Ldec_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_altyear_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes50_Ldec_all = senes50_Ldec
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes50_Ldec = senes50_Ldec_all#senes50_Ldec[!is.na(senes50_Ldec$debou10),] #senes50_Ldec_all

mod_senes50_Ldec_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) 
mod_senes50_Ldec_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tmoy_GS <- lmer(senes50 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tmoy_MJ <- lmer(senes50 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tmoy_AS <- lmer(senes50 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tmoy_ASO <- lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tnight21j_jsenes50 <- lmer(senes50 ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tnight30j_jsenes50 <- lmer(senes50 ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_Tnight40j_jsenes50 <- lmer(senes50 ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_GDDinv25_jsenes50 <- lmer(senes50 ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_GDDinv20_jsenes50 <- lmer(senes50 ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_GDDinv15_jsenes50 <- lmer(senes50 ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)

AIC(mod_senes50_Ldec_debou10, 
    mod_senes50_Ldec_P_summer,mod_senes50_Ldec_Tmoy_GS,mod_senes50_Ldec_Tmoy_MJ,mod_senes50_Ldec_Tmoy_AS,mod_senes50_Ldec_Tmoy_ASO,
    mod_senes50_Ldec_Tnight21j_jsenes50,mod_senes50_Ldec_Tnight30j_jsenes50,mod_senes50_Ldec_Tnight40j_jsenes50,
    mod_senes50_Ldec_GDDinv25_jsenes50,mod_senes50_Ldec_GDDinv20_jsenes50,mod_senes50_Ldec_GDDinv15_jsenes50)
# Les températures moyennes estivales (août-septembre) donnent les meilleurs résultats (AIC = 12604.73)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on reste sur le modèle à une variable :
mod_senes50_Ldec_multivar = lmer(senes50 ~ Tmoy_GS + 
                                   (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) # AIC = 12570.8
# visreg(mod_senes50_Ldec_multivar, "Tmoy_GS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes50_Ldec = lmer(senes50 ~ altitude * Tmoy_GS  +
                              (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) # AIC = 12513.7
# Finalement l'altitude donne de meilleurs scores que les variables climatiques...!
bestmods = c(list("Meleze"=bestmod_senes50_Ldec), bestmods)


summary(bestmod_senes50_Ldec)
visreg(bestmod_senes50_Ldec, "altitude", by="Tmoy_GS")
# MiSénescence 5.8 jours plus tôt par 100m d'altitude 
# MiSénescence 3.8 jours plus tôt par 1°C de température moyenne pendant la saison de végétation (gradient plus marqué à plus basse altitude, qui s'inverse quand on est plus haut)
qqnorm(resid(bestmod_senes50_Ldec))
qqline(resid(bestmod_senes50_Ldec))
# /!\ résidus !!!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Ldec)
# R2 des effets fixes plutôt nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Ldec)
# - Validation croisée
n = nrow(senes50_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Ldec[trainIndex ,]
test <- senes50_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Ldec$ID_zone[drop=T])[!unique(senes50_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Ldec), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        variable = rownames(coef(summary(bestmod_senes50_Ldec))),
                                        coef = coef(summary(bestmod_senes50_Ldec))[,1],
                                        std = coef(summary(bestmod_senes50_Ldec))[,2],
                                        pval = coef(summary(bestmod_senes50_Ldec))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Ldec)[rownames(coef(summary(bestmod_senes50_Ldec)))], error = function(e) return(NA))))



##############################################-
#*---- Bouleau pubescent ----

senes50_Bpub = senes50_Alps[senes50_Alps$species == "Bouleau_pubescent",]
# ggplot(senes50_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes50_Bpub = senes50_Bpub[!is.na(senes50_Bpub$Tmoy_GS),]
ggplot(senes50_Bpub, aes(x=senes50)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Bpub_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Bpub, REML=F)
summary(mod_senes50_Bpub_Altyear)
altyearmods = c(list("Bouleau_pubescent"=mod_senes50_Bpub_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes50_Bpub_Altyear)
# + Validation croisée
n = nrow(senes50_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Bpub[trainIndex ,]
test <- senes50_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Bpub$ID_zone[drop=T])[!unique(senes50_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes50_Bpub_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_altyear_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes50_Bpub_all = senes50_Bpub
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes50_Bpub = senes50_Bpub_all# senes50_Bpub[!is.na(senes50_Bpub$debou10),] #senes50_Bpub_all

mod_senes50_Bpub_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) 
mod_senes50_Bpub_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tmoy_GS <- lmer(senes50 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tmoy_MJ <- lmer(senes50 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tmoy_AS <- lmer(senes50 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tmoy_ASO <- lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tnight21j_jsenes50 <- lmer(senes50 ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tnight30j_jsenes50 <- lmer(senes50 ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_Tnight40j_jsenes50 <- lmer(senes50 ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_GDDinv25_jsenes50 <- lmer(senes50 ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_GDDinv20_jsenes50 <- lmer(senes50 ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_GDDinv15_jsenes50 <- lmer(senes50 ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)

AIC(mod_senes50_Bpub_debou10, 
    mod_senes50_Bpub_P_summer,mod_senes50_Bpub_Tmoy_GS,mod_senes50_Bpub_Tmoy_MJ,mod_senes50_Bpub_Tmoy_AS,mod_senes50_Bpub_Tmoy_ASO,
    mod_senes50_Bpub_Tnight21j_jsenes50,mod_senes50_Bpub_Tnight30j_jsenes50,mod_senes50_Bpub_Tnight40j_jsenes50,
    mod_senes50_Bpub_GDDinv25_jsenes50,mod_senes50_Bpub_GDDinv20_jsenes50,mod_senes50_Bpub_GDDinv15_jsenes50)
# Les températures estivales donnent les meilleurs résultats... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F), corr=F) # 
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes50_Bpub_multivar = lmer(senes50 ~ Tmoy_GS + 
                                   (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) # AIC = 482.4
# visreg(mod_senes50_Bpub_multivar, "GDDinv25_jsenes50")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senes50_Bpub_multivar = lmer(senes50 ~ Tmoy_GS + 
#                                    (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) # AIC = 641.3
# # visreg(mod_senes50_Bpub_multivar, "GDDinv25_jsenes50", by="P_summer")
# visreg(mod_senes50_Bpub_multivar, "P_summer", by="GDDinv25_jsenes50")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes50_Bpub = lmer(senes50 ~ Tmoy_GS +
                              (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) # AIC = 482.4
bestmods = c(list("Bouleau_pubescent"=bestmod_senes50_Bpub), bestmods)
visreg(bestmod_senes50_Bpub, "Tmoy_GS")


summary(bestmod_senes50_Bpub)
# températures pendant la saison de végétation : sénescence 4 jours plus tard quand on gagne 1°C

qqnorm(resid(bestmod_senes50_Bpub))
qqline(resid(bestmod_senes50_Bpub))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Bpub)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Bpub)
# - Validation croisée
n = nrow(senes50_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Bpub[trainIndex ,]
test <- senes50_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Bpub$ID_zone[drop=T])[!unique(senes50_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Bpub), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        variable = rownames(coef(summary(bestmod_senes50_Bpub))),
                                        coef = coef(summary(bestmod_senes50_Bpub))[,1],
                                        std = coef(summary(bestmod_senes50_Bpub))[,2],
                                        pval = coef(summary(bestmod_senes50_Bpub))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Bpub)[rownames(coef(summary(bestmod_senes50_Bpub)))], error = function(e) return(NA))))



##############################################-
#*---- Hetre ----

senes50_Fsyl = senes50_Alps[senes50_Alps$species == "Hetre",]
# ggplot(senes50_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes50_Fsyl = senes50_Fsyl[!is.na(senes50_Fsyl$Tmoy_GS),]
ggplot(senes50_Fsyl, aes(x=senes50)) + geom_histogram(binwidth=5)
ggplot(senes50_Fsyl, aes(x=year, y=senes50)) + geom_point()
# Une valeur sort du lot... on peut la retirer ? C'est une obs du CREA !! (zone Fouine, arbre 14011 en 2023... MAIS c'est cohérent entre senes10 et senes50,
# et on ne l'a pas retirée pour senes10 même si plus précoce de 20 jours environ... (30 jours pour senes50))
senes50_Fsyl = senes50_Fsyl[senes50_Fsyl$senes50 > 250,]

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Fsyl_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Fsyl, REML=F)
summary(mod_senes50_Fsyl_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes50_Fsyl_Altyear)
# + Validation croisée
n = nrow(senes50_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Fsyl[trainIndex ,]
test <- senes50_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Fsyl$ID_zone[drop=T])[!unique(senes50_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes50_Fsyl_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Hetre", "R2_altyear_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes50_Fsyl_all = senes50_Fsyl
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes50_Fsyl = senes50_Fsyl_all# senes50_Fsyl[!is.na(senes50_Fsyl$debou10),] #senes50_Fsyl_all

mod_senes50_Fsyl_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) 
mod_senes50_Fsyl_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tmoy_GS <- lmer(senes50 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tmoy_MJ <- lmer(senes50 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tmoy_AS <- lmer(senes50 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tmoy_ASO <- lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tnight21j_jsenes50 <- lmer(senes50 ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tnight30j_jsenes50 <- lmer(senes50 ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_Tnight40j_jsenes50 <- lmer(senes50 ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_GDDinv25_jsenes50 <- lmer(senes50 ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_GDDinv20_jsenes50 <- lmer(senes50 ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_GDDinv15_jsenes50 <- lmer(senes50 ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)

AIC(mod_senes50_Fsyl_debou10, 
    mod_senes50_Fsyl_P_summer,mod_senes50_Fsyl_Tmoy_GS,mod_senes50_Fsyl_Tmoy_MJ,mod_senes50_Fsyl_Tmoy_AS,mod_senes50_Fsyl_Tmoy_ASO,
    mod_senes50_Fsyl_Tnight21j_jsenes50,mod_senes50_Fsyl_Tnight30j_jsenes50,mod_senes50_Fsyl_Tnight40j_jsenes50,
    mod_senes50_Fsyl_GDDinv25_jsenes50,mod_senes50_Fsyl_GDDinv20_jsenes50,mod_senes50_Fsyl_GDDinv15_jsenes50)
# Tout donne à peu près la même chose... bizarre...


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 

# NUL NUL NUL !!!!


# # En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
# mod_senes50_Fsyl_multivar = lmer(senes50 ~ P_summer + Tmoy_AS + Tmoy_ASO + 
#                                    GDDinv25_jsenes50 +
#                                    (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) # AIC = 526.1
# # visreg(mod_senes50_Fsyl_multivar, "GDDinv25_jsenes50")
#
# # On peut aussi complexifier en ajoutant des interactions :
# mod_senes50_Fsyl_multivar = lmer(senes50 ~ P_summer * Tmoy_AS + 
#                                    (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) # AIC = 521.4
# # visreg(mod_senes50_Fsyl_multivar, "Tmoy_AS", by="P_summer")
# # visreg(mod_senes50_Fsyl_multivar, "P_summer", by="Tmoy_AS")
# # TOUJOURS NUL !!!


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes50_Fsyl = lmer(senes50 ~ altitude + 
                              (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) # AIC = 392.6
bestmods = c(list("Hetre"=bestmod_senes50_Fsyl), bestmods)
# visreg(bestmod_senes50_Fsyl, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes50_Fsyl, "Tmoy_AS", by="P_summer")


summary(bestmod_senes50_Fsyl)
# mi-sénescence 3.4 jours plus tôt paour 100m d'altitude

qqnorm(resid(bestmod_senes50_Fsyl))
qqline(resid(bestmod_senes50_Fsyl))
# résidus à peu près ok 


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Fsyl)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Fsyl)
# - Validation croisée
n = nrow(senes50_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Fsyl[trainIndex ,]
test <- senes50_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Fsyl$ID_zone[drop=T])[!unique(senes50_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Fsyl), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        variable = rownames(coef(summary(bestmod_senes50_Fsyl))),
                                        coef = coef(summary(bestmod_senes50_Fsyl))[,1],
                                        std = coef(summary(bestmod_senes50_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senes50_Fsyl))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Fsyl)[rownames(coef(summary(bestmod_senes50_Fsyl)))], error = function(e) return(NA))))



##############################################-
#*---- Sorbier ----

senes50_Sacu = senes50_Alps[senes50_Alps$species == "Sorbier",]
# ggplot(senes50_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes50_Sacu = senes50_Sacu[!is.na(senes50_Sacu$Tmoy_GS),]
ggplot(senes50_Sacu, aes(x=senes50)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Sacu_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Sacu, REML=F)
summary(mod_senes50_Sacu_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes50_Sacu_Altyear)
# + Validation croisée
n = nrow(senes50_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Sacu[trainIndex ,]
test <- senes50_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Sacu$ID_zone[drop=T])[!unique(senes50_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senes50_Sacu_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_altyear_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senes50_Sacu_all = senes50_Sacu
# ICI le timing de débourrement ne semble pas affecter la sénescence
senes50_Sacu = senes50_Sacu_all# senes50_Sacu[!is.na(senes50_Sacu$debou10),] #senes50_Sacu_all

mod_senes50_Sacu_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) 
mod_senes50_Sacu_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tmoy_GS <- lmer(senes50 ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tmoy_MJ <- lmer(senes50 ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tmoy_AS <- lmer(senes50 ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tmoy_ASO <- lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tnight21j_jsenes50 <- lmer(senes50 ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tnight30j_jsenes50 <- lmer(senes50 ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_Tnight40j_jsenes50 <- lmer(senes50 ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_GDDinv25_jsenes50 <- lmer(senes50 ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_GDDinv20_jsenes50 <- lmer(senes50 ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_GDDinv15_jsenes50 <- lmer(senes50 ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)

AIC(mod_senes50_Sacu_debou10, 
    mod_senes50_Sacu_P_summer,mod_senes50_Sacu_Tmoy_GS,mod_senes50_Sacu_Tmoy_MJ,mod_senes50_Sacu_Tmoy_AS,mod_senes50_Sacu_Tmoy_ASO,
    mod_senes50_Sacu_Tnight21j_jsenes50,mod_senes50_Sacu_Tnight30j_jsenes50,mod_senes50_Sacu_Tnight40j_jsenes50,
    mod_senes50_Sacu_GDDinv25_jsenes50,mod_senes50_Sacu_GDDinv20_jsenes50,mod_senes50_Sacu_GDDinv15_jsenes50)
# Les variables de température AS - ASO donnent les meilleurs résultats (AIC = 4940.168)...
# on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F), corr=F)  # /!\ en ayant ajouté debou10, on tronque la base de données de 20% !
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes50_Sacu_multivar = lmer(senes50 ~ Tnight30j_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) # AIC = 5095.0
# visreg(mod_senes50_Sacu_multivar, "GDDinv25_jsenes50")

# On peut aussi complexifier en ajoutant des interactions 
mod_senes50_Sacu_multivar = lmer(senes50 ~ Tnight30j_jsenes50 * P_summer +
                                   (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) # AIC = 5077.0


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes50_Sacu = lmer(senes50 ~ Tnight30j_jsenes50 * P_summer +
                              (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) # AIC = 5077.0
# Ici l'altitude n'améliore rien, MAIS en mettant la date de débourrement en variable explicative on tronque la base de données
bestmods = c(list("Sorbier"=bestmod_senes50_Sacu), bestmods)
# visreg(bestmod_senes50_Sacu, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes50_Sacu, "Tmoy_AS", by="P_summer")


summary(bestmod_senes50_Sacu)
# températures nocturnes un peu avant la sénescence : sénescence 1.5 jours plus tard quand on gagne 1°C en moyenne la nuit
# accumulation de froid (pas vraiment significatif !) : sénescence 1.1 jour plus tard quand on gagne 1°C de froid
# quand le débourrement est plus tardif d'une semaine, la sénescence est retardée d'1 jour seulement

qqnorm(resid(bestmod_senes50_Sacu))
qqline(resid(bestmod_senes50_Sacu))
# résidus quasi ok !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Sacu)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Sacu)
# - Validation croisée
n = nrow(senes50_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Sacu[trainIndex ,]
test <- senes50_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Sacu$ID_zone[drop=T])[!unique(senes50_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Sacu), train)
predictions <- mod_train %>% predict(test1) # BUG !!!
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        variable = rownames(coef(summary(bestmod_senes50_Sacu))),
                                        coef = coef(summary(bestmod_senes50_Sacu))[,1],
                                        std = coef(summary(bestmod_senes50_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes50_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Sacu)[rownames(coef(summary(bestmod_senes50_Sacu)))], error = function(e) return(NA))))



write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes50__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes50__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes50.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senes50_v1altyear.Rdata")



##############################################-
# *-- DURÉE DE SÉNESCENCE                 ----
##############################################-


senesduree_Alps = senes_Alps_all[!is.na(senes_Alps_all$senesduree),]
# on considère qu'il n'est pas possible de passer de 10 à 50% de feuilles qui ont changé de couleur en moins de 2 jours (erreur de saisie ?)
senesduree_Alps = senesduree_Alps[senesduree_Alps$senesduree >= 2,] 

# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2024, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senesduree_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()


##############################################-
#*---- Bouleau verruqueux ----

senesduree_Bpen = senesduree_Alps[senesduree_Alps$species == "Bouleau_verruqueux",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Bpen = senesduree_Bpen[!is.na(senesduree_Bpen$Tmoy_GS),]
ggplot(senesduree_Bpen, aes(x=senesduree)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Bpen_Altyear <- lmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Bpen, REML=F)
summary(mod_senesduree_Bpen_Altyear)
altyearmods = c(list("Bouleau_verruqueux"=mod_senesduree_Bpen_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Bpen_Altyear)
# + Validation croisée
n = nrow(senesduree_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Bpen[trainIndex ,]
test <- senesduree_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Bpen$ID_zone[drop=T])[!unique(senesduree_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senesduree_Bpen_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Bpen_all = senesduree_Bpen
# ICI le timing de débourrement ne semble pas affecter la sénescence
senesduree_Bpen = senesduree_Bpen_all#senesduree_Bpen[!is.na(senesduree_Bpen$debou10),] #senesduree_Bpen_all

mod_senesduree_Bpen_debou10 <- lmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F) 
mod_senesduree_Bpen_P_summer <- lmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)
mod_senesduree_Bpen_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F)

AIC(mod_senesduree_Bpen_debou10, 
    mod_senesduree_Bpen_P_summer,mod_senesduree_Bpen_Tmoy_GS,mod_senesduree_Bpen_Tmoy_MJ,mod_senesduree_Bpen_Tmoy_AS,mod_senesduree_Bpen_Tmoy_ASO,
    mod_senesduree_Bpen_Tnight21j_jsenes10,mod_senesduree_Bpen_Tnight30j_jsenes10,mod_senesduree_Bpen_Tnight40j_jsenes10,
    mod_senesduree_Bpen_GDDinv25_jsenes10,mod_senesduree_Bpen_GDDinv20_jsenes10,mod_senesduree_Bpen_GDDinv15_jsenes10)
# Pas une grosse différence entre les variables... et la date de débourrement n'apporte pas grandchose (mais tronque 20% des données)

# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
mod_senesduree_Bpen_multivar = lmer(senesduree ~ P_summer + Tnight40j_jsenes50 + 
                                   (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F) # AIC = 11841.5
# visreg(mod_senesduree_Bpen_multivar, "Tnight40j_jsenes50")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senesduree_Bpen_multivar = lmer(senesduree ~ P_summer * GDDinv25_jsenes10 + 
#                                    (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F) 
# # visreg(mod_senesduree_Bpen_multivar, "GDDinv25_jsenes10", by="P_summer")
# visreg(mod_senesduree_Bpen_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senesduree_Bpen = lmer(senesduree ~ P_summer + Tnight40j_jsenes50 +  
                              (1|ID_zone) + (1|yearQ), senesduree_Bpen, REML=F) # AIC = 11841.5
bestmods = c(list("Bouleau_verruqueux"=bestmod_senesduree_Bpen), bestmods)


summary(bestmod_senesduree_Bpen)
# sénescence rallongée de 1 jour par 5mm de pluie supplémentaire pendant l'été
# sénescence rallongée de 0.7 jour quand il fait plus froid la nuit de 1°C
qqnorm(resid(bestmod_senesduree_Bpen))
qqline(resid(bestmod_senesduree_Bpen))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Bpen)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Bpen)
# - Validation croisée
n = nrow(senesduree_Bpen)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Bpen[trainIndex ,]
test <- senesduree_Bpen[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Bpen$ID_zone[drop=T])[!unique(senesduree_Bpen$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senesduree_Bpen), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        variable = rownames(coef(summary(bestmod_senesduree_Bpen))),
                                        coef = coef(summary(bestmod_senesduree_Bpen))[,1],
                                        std = coef(summary(bestmod_senesduree_Bpen))[,2],
                                        pval = coef(summary(bestmod_senesduree_Bpen))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Bpen)[rownames(coef(summary(bestmod_senesduree_Bpen)))], error = function(e) return(NA))))



##############################################-
#*---- Meleze ----

senesduree_Ldec = senesduree_Alps[senesduree_Alps$species == "Meleze",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Ldec = senesduree_Ldec[!is.na(senesduree_Ldec$Tmoy_GS),]
ggplot(senesduree_Ldec, aes(x=senesduree)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Ldec_Altyear <- lmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Ldec, REML=F)
summary(mod_senesduree_Ldec_Altyear)
altyearmods = c(list("Meleze"=mod_senesduree_Ldec_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Ldec_Altyear)
# + Validation croisée
n = nrow(senesduree_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Ldec[trainIndex ,]
test <- senesduree_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Ldec$ID_zone[drop=T])[!unique(senesduree_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senesduree_Ldec_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Ldec_all = senesduree_Ldec
# ICI le timing de débourrement ne semble pas affecter la sénescence
senesduree_Ldec = senesduree_Ldec_all#senesduree_Ldec[!is.na(senesduree_Ldec$debou10),] #senesduree_Ldec_all

mod_senesduree_Ldec_debou10 <- lmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F) 
mod_senesduree_Ldec_P_summer <- lmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)
mod_senesduree_Ldec_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F)

AIC(mod_senesduree_Ldec_debou10, 
    mod_senesduree_Ldec_P_summer,mod_senesduree_Ldec_Tmoy_GS,mod_senesduree_Ldec_Tmoy_MJ,mod_senesduree_Ldec_Tmoy_AS,mod_senesduree_Ldec_Tmoy_ASO,
    mod_senesduree_Ldec_Tnight21j_jsenes10,mod_senesduree_Ldec_Tnight30j_jsenes10,mod_senesduree_Ldec_Tnight40j_jsenes10,
    mod_senesduree_Ldec_GDDinv25_jsenes10,mod_senesduree_Ldec_GDDinv20_jsenes10,mod_senesduree_Ldec_GDDinv15_jsenes10)
# Pas de grande différence entre les variables testées... et le débourrement n'apporte pas grand-chose (et ça tronque 15% des données)

# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on reste sur le modèle à une variable :
mod_senesduree_Ldec_multivar = lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_AS + Tnight40j_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F) # AIC = 12125.3
# visreg(mod_senesduree_Ldec_multivar, "Tmoy_GS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senesduree_Ldec = lmer(senesduree ~ P_summer + Tnight21j_jsenes50  +
                              (1|ID_zone) + (1|yearQ), senesduree_Ldec, REML=F) # AIC = 12126.8
# Finalement retirer 2 variables supplémentaires n'améliore pas le modèle, mais c'est plus facile à interpréter et on limite la corrélation entre variables
bestmods = c(list("Meleze"=bestmod_senesduree_Ldec), bestmods)


summary(bestmod_senesduree_Ldec)
# Sénescence allongée de 1 jour pour 4 mm de pluie supplémentaire sur la période estivale
# Sénescence allongée de 1 jour pour 2°C de plus dans les températures nocturnes
qqnorm(resid(bestmod_senesduree_Ldec))
qqline(resid(bestmod_senesduree_Ldec))
# /!\ résidus !!!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Ldec)
# R2 des effets fixes plutôt nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Ldec)
# - Validation croisée
n = nrow(senesduree_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Ldec[trainIndex ,]
test <- senesduree_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Ldec$ID_zone[drop=T])[!unique(senesduree_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senesduree_Ldec), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        variable = rownames(coef(summary(bestmod_senesduree_Ldec))),
                                        coef = coef(summary(bestmod_senesduree_Ldec))[,1],
                                        std = coef(summary(bestmod_senesduree_Ldec))[,2],
                                        pval = coef(summary(bestmod_senesduree_Ldec))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Ldec)[rownames(coef(summary(bestmod_senesduree_Ldec)))], error = function(e) return(NA))))



##############################################-
#*---- Bouleau pubescent ----

senesduree_Bpub = senesduree_Alps[senesduree_Alps$species == "Bouleau_pubescent",]
# ggplot(senesduree_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Bpub = senesduree_Bpub[!is.na(senesduree_Bpub$Tmoy_GS),]
ggplot(senesduree_Bpub, aes(x=senesduree)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Bpub_Altyear <- lmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Bpub, REML=F)
summary(mod_senesduree_Bpub_Altyear)
altyearmods = c(list("Bouleau_pubescent"=mod_senesduree_Bpub_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Bpub_Altyear)
# + Validation croisée
n = nrow(senesduree_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Bpub[trainIndex ,]
test <- senesduree_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Bpub$ID_zone[drop=T])[!unique(senesduree_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senesduree_Bpub_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Bpub_all = senesduree_Bpub
# ICI le timing de débourrement ne semble pas affecter la sénescence
senesduree_Bpub = senesduree_Bpub_all# senesduree_Bpub[!is.na(senesduree_Bpub$debou10),] #senesduree_Bpub_all

mod_senesduree_Bpub_debou10 <- lmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F) 
mod_senesduree_Bpub_P_summer <- lmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)
mod_senesduree_Bpub_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F)

AIC(mod_senesduree_Bpub_debou10, 
    mod_senesduree_Bpub_P_summer,mod_senesduree_Bpub_Tmoy_GS,mod_senesduree_Bpub_Tmoy_MJ,mod_senesduree_Bpub_Tmoy_AS,mod_senesduree_Bpub_Tmoy_ASO,
    mod_senesduree_Bpub_Tnight21j_jsenes10,mod_senesduree_Bpub_Tnight30j_jsenes10,mod_senesduree_Bpub_Tnight40j_jsenes10,
    mod_senesduree_Bpub_GDDinv25_jsenes10,mod_senesduree_Bpub_GDDinv20_jsenes10,mod_senesduree_Bpub_GDDinv15_jsenes10)
# précipitations expliquent mieux


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F), corr=F) # 
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senesduree_Bpub_multivar = lmer(senesduree ~ P_summer + 
                                   (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F) # AIC = 344.0
# visreg(mod_senesduree_Bpub_multivar, "P_summer")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senesduree_Bpub_multivar = lmer(senesduree ~ Tmoy_GS + 
#                                    (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F) # AIC = 641.3
# # visreg(mod_senesduree_Bpub_multivar, "GDDinv25_jsenes10", by="P_summer")
# visreg(mod_senesduree_Bpub_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senesduree_Bpub = lmer(senesduree ~ P_summer +
                              (1|ID_zone) + (1|yearQ), senesduree_Bpub, REML=F) # AIC = 344.0
bestmods = c(list("Bouleau_pubescent"=bestmod_senesduree_Bpub), bestmods)


summary(bestmod_senesduree_Bpub)
# sénescence raccourcie de 0.8 jour quand il y a 1mm de précipitations en plus dans l'été

qqnorm(resid(bestmod_senesduree_Bpub))
qqline(resid(bestmod_senesduree_Bpub))
# résidus à peu près ok

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Bpub)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Bpub)
# - Validation croisée
n = nrow(senesduree_Bpub)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Bpub[trainIndex ,]
test <- senesduree_Bpub[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Bpub$ID_zone[drop=T])[!unique(senesduree_Bpub$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senesduree_Bpub), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        variable = rownames(coef(summary(bestmod_senesduree_Bpub))),
                                        coef = coef(summary(bestmod_senesduree_Bpub))[,1],
                                        std = coef(summary(bestmod_senesduree_Bpub))[,2],
                                        pval = coef(summary(bestmod_senesduree_Bpub))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Bpub)[rownames(coef(summary(bestmod_senesduree_Bpub)))], error = function(e) return(NA))))



##############################################-
#*---- Hetre ----

senesduree_Fsyl = senesduree_Alps[senesduree_Alps$species == "Hetre",]
# ggplot(senesduree_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Fsyl = senesduree_Fsyl[!is.na(senesduree_Fsyl$Tmoy_GS),]
ggplot(senesduree_Fsyl, aes(x=senesduree)) + geom_histogram(binwidth=5)
ggplot(senesduree_Fsyl, aes(x=year, y=senesduree)) + geom_point()
# Une valeur sort du lot... on peut la retirer ? C'est une obs du CREA !! (zone Fouine - arbre 14006 en 2023)
senesduree_Fsyl = senesduree_Fsyl[senesduree_Fsyl$senesduree < 65,]

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Fsyl_Altyear <- lmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Fsyl, REML=F)
summary(mod_senesduree_Fsyl_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Fsyl_Altyear)
# + Validation croisée
n = nrow(senesduree_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Fsyl[trainIndex ,]
test <- senesduree_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Fsyl$ID_zone[drop=T])[!unique(senesduree_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senesduree_Fsyl_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Hetre", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Fsyl_all = senesduree_Fsyl
# ICI le timing de débourrement ne semble pas affecter la sénescence
senesduree_Fsyl = senesduree_Fsyl_all#senesduree_Fsyl[!is.na(senesduree_Fsyl$debou10),] #senesduree_Fsyl_all

mod_senesduree_Fsyl_debou10 <- lmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F) 
mod_senesduree_Fsyl_P_summer <- lmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)
mod_senesduree_Fsyl_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F)

AIC(mod_senesduree_Fsyl_debou10, 
    mod_senesduree_Fsyl_P_summer,mod_senesduree_Fsyl_Tmoy_GS,mod_senesduree_Fsyl_Tmoy_MJ,mod_senesduree_Fsyl_Tmoy_AS,mod_senesduree_Fsyl_Tmoy_ASO,
    mod_senesduree_Fsyl_Tnight21j_jsenes10,mod_senesduree_Fsyl_Tnight30j_jsenes10,mod_senesduree_Fsyl_Tnight40j_jsenes10,
    mod_senesduree_Fsyl_GDDinv25_jsenes10,mod_senesduree_Fsyl_GDDinv20_jsenes10,mod_senesduree_Fsyl_GDDinv15_jsenes10)
# Précip meilleur + debou intéressant aussi


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En enlevant les points où on n'a pas d'info de débourrement, on n'enlève que 4 données (8%) = acceptable

# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senesduree_Fsyl_multivar = lmer(senesduree ~ debou10 + P_summer +
                                   (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F) # AIC = 342.7
# visreg(mod_senesduree_Fsyl_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senesduree_Fsyl_multivar = lmer(senesduree ~ P_summer * debou10 +
                                   (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F) # AIC = 339.7
# visreg(mod_senesduree_Fsyl_multivar, "debou10", by="P_summer")
# visreg(mod_senesduree_Fsyl_multivar, "P_summer", by="Tmoy_AS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senesduree_Fsyl = lmer(senesduree ~ P_summer * debou10 +
                              (1|ID_zone) + (1|yearQ), senesduree_Fsyl, REML=F) # AIC = 339.7
bestmods = c(list("Hetre"=bestmod_senesduree_Fsyl), bestmods)
# visreg(bestmod_senesduree_Fsyl, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senesduree_Fsyl, "Tmoy_AS", by="P_summer")


summary(bestmod_senesduree_Fsyl)
# précipitations : sénescence 28 jours plus tôt quand on gagne 1mm de pluie
# chaleur de fin de saison (août-septembre) : sénescence 20 jours plus tard quand on gagne 1°C (retard d'autant plus fort qu'il fait sec)
# /!\ ça me semble beaucoup ces chiffres !!!

qqnorm(resid(bestmod_senesduree_Fsyl))
qqline(resid(bestmod_senesduree_Fsyl))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Fsyl)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Fsyl)
# - Validation croisée
n = nrow(senesduree_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Fsyl[trainIndex ,]
test <- senesduree_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Fsyl$ID_zone[drop=T])[!unique(senesduree_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senesduree_Fsyl), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        variable = rownames(coef(summary(bestmod_senesduree_Fsyl))),
                                        coef = coef(summary(bestmod_senesduree_Fsyl))[,1],
                                        std = coef(summary(bestmod_senesduree_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senesduree_Fsyl))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Fsyl)[rownames(coef(summary(bestmod_senesduree_Fsyl)))], error = function(e) return(NA))))



##############################################-
#*---- Sorbier ----

senesduree_Sacu = senesduree_Alps[senesduree_Alps$species == "Sorbier",]
# ggplot(senesduree_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Sacu = senesduree_Sacu[!is.na(senesduree_Sacu$Tmoy_GS),]
ggplot(senesduree_Sacu, aes(x=senesduree)) + geom_histogram(binwidth=5)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Sacu_Altyear <- lmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Sacu, REML=F)
summary(mod_senesduree_Sacu_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Sacu_Altyear)
# + Validation croisée
n = nrow(senesduree_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Sacu[trainIndex ,]
test <- senesduree_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Sacu$ID_zone[drop=T])[!unique(senesduree_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(mod_senesduree_Sacu_Altyear), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Sacu_all = senesduree_Sacu
# ICI le timing de débourrement ne semble pas affecter la sénescence
senesduree_Sacu = senesduree_Sacu_all# senesduree_Sacu[!is.na(senesduree_Sacu$debou10),] #senesduree_Sacu_all

mod_senesduree_Sacu_debou10 <- lmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F) 
mod_senesduree_Sacu_P_summer <- lmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)
mod_senesduree_Sacu_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F)

AIC(mod_senesduree_Sacu_debou10, 
    mod_senesduree_Sacu_P_summer,mod_senesduree_Sacu_Tmoy_GS,mod_senesduree_Sacu_Tmoy_MJ,mod_senesduree_Sacu_Tmoy_AS,mod_senesduree_Sacu_Tmoy_ASO,
    mod_senesduree_Sacu_Tnight21j_jsenes10,mod_senesduree_Sacu_Tnight30j_jsenes10,mod_senesduree_Sacu_Tnight40j_jsenes10,
    mod_senesduree_Sacu_GDDinv25_jsenes10,mod_senesduree_Sacu_GDDinv20_jsenes10,mod_senesduree_Sacu_GDDinv15_jsenes10)
# debou semble intéressant, ainsi que T


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ debou10 + P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senesduree_Sacu_multivar = lmer(senesduree ~ debou10 +
                                   (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F) # AIC = 3640.4 # /!\ en ayant ajouté debou10, on tronque la base de données de 14% !
# visreg(mod_senesduree_Sacu_multivar, "GDDinv25_jsenes10")

# # On peut aussi complexifier en ajoutant des interactions 
# mod_senesduree_Sacu_multivar = lmer(senesduree ~ debou10 +
#                                    (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F) # AIC = 3640.4


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senesduree_Sacu = lmer(senesduree ~ debou10 +
                              (1|ID_zone) + (1|yearQ), senesduree_Sacu, REML=F) # AIC = 3640.4
# Ici l'altitude n'améliore rien, MAIS en mettant la date de débourrement en variable explicative on tronque la base de données
bestmods = c(list("Sorbier"=bestmod_senesduree_Sacu), bestmods)
# visreg(bestmod_senesduree_Sacu, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senesduree_Sacu, "Tmoy_AS", by="P_summer")


summary(bestmod_senesduree_Sacu)
# températures nocturnes un peu avant la sénescence : sénescence 1.5 jours plus tard quand on gagne 1°C en moyenne la nuit
# accumulation de froid (pas vraiment significatif !) : sénescence 1.1 jour plus tard quand on gagne 1°C de froid
# quand le débourrement est plus tardif d'une semaine, la sénescence est retardée d'1 jour seulement

qqnorm(resid(bestmod_senesduree_Sacu))
qqline(resid(bestmod_senesduree_Sacu))
# résidus quasi ok !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Sacu)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Sacu)
# - Validation croisée
n = nrow(senesduree_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Sacu[trainIndex ,]
test <- senesduree_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Sacu$ID_zone[drop=T])[!unique(senesduree_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senesduree_Sacu), train)
predictions <- mod_train %>% predict(test1) # BUG !!!
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        variable = rownames(coef(summary(bestmod_senesduree_Sacu))),
                                        coef = coef(summary(bestmod_senesduree_Sacu))[,1],
                                        std = coef(summary(bestmod_senesduree_Sacu))[,2],
                                        pval = coef(summary(bestmod_senesduree_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Sacu)[rownames(coef(summary(bestmod_senesduree_Sacu)))], error = function(e) return(NA))))



write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senesduree.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senesduree_v1altyear.Rdata")








