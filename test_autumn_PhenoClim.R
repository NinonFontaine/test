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
library(lme4)
library(lmerTest) # permet d'avoir les pvalues qui s'affichent pour les lmer
library(visreg)
library(MuMIn) # fonction r.squaredGLMM notamment
library(variancePartition) # fonction calcVarPart
library(brms)
library(sjPlot)
library(meteor) # pour le calcul de la photopériode


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
phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Phenoclim/data/_CLEANED_data_pheno.csv")
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
#   bf(julian_day ~ altitude+yearQ+(yearQ|nom_zone)),
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
resume_automne = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/prototool/Phenoclim/Analyse/Indice_phenoclim/pheno_year_global.csv")

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
           nom_zone, dept1, region1, pays1, circumference, new_cat, year
           )

#=============================================================================================================================*
# On associe les données de sénescence aux données de débourrement des mêmes arbres

sites_pheno_Alps = merge(sites_pheno_Alps, 
                         merge(phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Debourrement - Ok 10%",c("id_base_site","julian_day","year")],
                               merge(phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Changement de couleur - Ok 10%",c("id_base_site","julian_day","year")],
                                     phenoclim[phenoclim$nom_massif_v2019=="Alpes" & phenoclim$pheno_etape_value == "Changement de couleur - Ok 50%",c("id_base_site","julian_day","year")], 
                                     by=c("id_base_site","year"), all.x=T, all.y=T),
                               by=c("id_base_site","year"), all.x=T, all.y=T),
                         by=c("id_base_site","year"), all.x=T, all.y=T)
colnames(sites_pheno_Alps)[(ncol(sites_pheno_Alps)-2):ncol(sites_pheno_Alps)] = c("debou10","senes10","senes50")

sites_pheno_Alps = sites_pheno_Alps[!is.na(sites_pheno_Alps$year),]

# dim(sites_pheno_Alps[!is.na(sites_pheno_Alps$debou10) & !is.na(sites_pheno_Alps$senes10) & !is.na(sites_pheno_Alps$senes50),])

#=============================================================================================================================*
# On récupère les précipitations de ERA5-land (Monthly averaged reanalysis)
# --> données mensuelles de 2004 à 2024, à une résolution de 9km (https://cds.climate.copernicus.eu/requests?tab=all)
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

sites_pheno_Alps$P_summer = unlist(apply(sites_pheno_Alps[,c("year",grep("Psummer",colnames(sites_pheno_Alps), value=T))],1,
                                  function(x){x[x[1]-2002]}))
sites_pheno_Alps = sites_pheno_Alps[,-grep("Psummer", colnames(sites_pheno_Alps))]

# # Visualisation
# par(mfrow=c(2,2))
# for (annee in c(2004, 2010, 2023, 2024)){hist(sites_pheno_Alps[,paste0("P_summer",annee)], main=annee, xlim=c(0,0.03), breaks=10)}
# # --> année 2023 particulièrement sèche !

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
  
  for (esp in unique(sites_pheno_Alps$species[sites_pheno_Alps$year == annee & (!is.na(sites_pheno_Alps$senes10) | !is.na(sites_pheno_Alps$senes50))])){

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
# (/!\ ce n'est pas vraiment une durée puisqu'on regarde le nombre de jours entre 1à et 50% de changement de couleur)
senes_Alps_all$duree_senes = senes_Alps_all$senes50 - senes_Alps_all$senes10

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
# => on choisit une parmi Tmoy30, 40, 50j, GDD0, GDD5, qu'on combine avec Tmoyhiv et/ou dChill


##############################################-
# *-- ONSET DE SÉNESCENCE                  ----
##############################################-

senes_Alps = senes_Alps_all[!is.na(senes_Alps_all$senes10),]

# # Aperçu des données de sénescence
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=altitude, y=senes10, col=species)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes10, col=species)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senes10, col=yearQ)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=GDDinv20_jsenes10, y=senes10, col=yearQ)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=P_summer, y=senes10)) + geom_point() + 
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)
# ggplot(senes_Alps[!is.na(senes_Alps$Tmoy_GS),], aes(x=debou10, y=senes10)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) + ylim(200,325)



# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2024, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes_Alps$species), R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_bestmod_fixef = NA, R2_bestmod_allef = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()


##############################################-
#*---- Bouleau verruqueux ----

senes_Bpen = senes_Alps[senes_Alps$species == "Bouleau_verruqueux",]
# ggplot(senes_Bpen, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Bpen = senes_Bpen[!is.na(senes_Bpen$Tmoy_GS),]
# On retire une valeur extrême (sénescence notée à 78jours = 18/03 !)
senes_Bpen = senes_Bpen[senes_Bpen$senes10 > 100,]
# ggplot(senes_Bpen, aes(x=year, y=julian_day)) + geom_point() + geom_smooth(method=lm) + labs(title = "Changement de couleur 10% - Bouleau", x="année",y="jour julien")

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Bpen_Altyear <- lmer(senes10 ~ altitude + year + (year|nom_zone), senes_Bpen, REML=F)
summary(mod_senes_Bpen_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Bpen_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# mod_senes_Bpen_debou10 <- lmer(senes10 ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) 
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Bpen_P_summer <- lmer(senes10 ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)

AIC(#mod_senes_Bpen_debou10, 
    mod_senes_Bpen_P_summer,mod_senes_Bpen_Tmoy_GS,mod_senes_Bpen_Tmoy_MJ,mod_senes_Bpen_Tmoy_AS,mod_senes_Bpen_Tmoy_ASO,
    mod_senes_Bpen_Tnight21j_jsenes10,mod_senes_Bpen_Tnight30j_jsenes10,mod_senes_Bpen_Tnight40j_jsenes10,
    mod_senes_Bpen_GDDinv25_jsenes10,mod_senes_Bpen_GDDinv20_jsenes10,mod_senes_Bpen_GDDinv15_jsenes10)
# Les variables de GDDinverse à 25 et 20°C donnent les meilleurs résultats (AIC = 13346.73)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                               (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Bpen_multivar = lmer(senes10 ~ P_summer + Tmoy_ASO + Tnight30j_jsenes10 + GDDinv25_jsenes10 + 
                               (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) # AIC = 13338.5
# visreg(mod_senes_Bpen_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes_Bpen_multivar = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + 
                                 (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) # AIC = 13334.6
# visreg(mod_senes_Bpen_multivar, "GDDinv25_jsenes10", by="P_summer")
visreg(mod_senes_Bpen_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Bpen = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + altitude +
                                 (1|nom_zone) + (P_summer|yearQ), senes_Bpen, REML=F) # AIC = 13328.7
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes_Bpen), bestmods)


summary(bestmod_senes_Bpen)
# gradient altitudinal : sénescence 0.9 jour plus tôt quand on monte de 100m
# précipitations : sénescence 1.8 jour plus tôt quand on gagne 1mm de pluie
# accumulation de "froid" : sénescence 2.0 jours plus tôt quand on gagne 1°C de froid
qqnorm(resid(bestmod_senes_Bpen))
qqline(resid(bestmod_senes_Bpen))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Bpen)
# R2 des effets fixes très très nul !!

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Bpen)


resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        variable = rownames(coef(summary(bestmod_senes_Bpen))),
                                        coef = coef(summary(bestmod_senes_Bpen))[,1],
                                        std = coef(summary(bestmod_senes_Bpen))[,2],
                                        pval = coef(summary(bestmod_senes_Bpen))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Bpen)[rownames(coef(summary(bestmod_senes_Bpen)))]))



##############################################-
#*---- Meleze ----

senes_Ldec = senes_Alps[senes_Alps$species == "Meleze",]
# ggplot(senes_Ldec, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Ldec = senes_Ldec[!is.na(senes_Ldec$Tmoy_GS),]


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Ldec_Altyear <- lmer(senes10 ~ altitude + year + (year|nom_zone), senes_Ldec, REML=F)
summary(mod_senes_Ldec_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Ldec_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# mod_senes_Ldec_debou10 <- lmer(senes10 ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
# /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
#     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Ldec_P_summer <- lmer(senes10 ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)

AIC(#mod_senes_Ldec_debou10, 
  mod_senes_Ldec_P_summer,mod_senes_Ldec_Tmoy_GS,mod_senes_Ldec_Tmoy_MJ,mod_senes_Ldec_Tmoy_AS,mod_senes_Ldec_Tmoy_ASO,
  mod_senes_Ldec_Tnight21j_jsenes10,mod_senes_Ldec_Tnight30j_jsenes10,mod_senes_Ldec_Tnight40j_jsenes10,
  mod_senes_Ldec_GDDinv25_jsenes10,mod_senes_Ldec_GDDinv20_jsenes10,mod_senes_Ldec_GDDinv15_jsenes10)
# Les températures moyennes estivales (août-septembre) donnent les meilleurs résultats (AIC = 13716.53)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on reste sur le modèle à une variable :
mod_senes_Ldec_multivar = lmer(senes10 ~ Tmoy_AS + 
                                 (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) # AIC = 13716.5
# visreg(mod_senes_Ldec_multivar, "Tmoy_AS")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senes_Ldec_multivar = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + 
#                                  (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) 
# # -> ça n'améliore pas les modèles

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Ldec = lmer(senes10 ~ Tmoy_AS  +
                            (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) # AIC = 13716.5
bestmods = c(list("Meleze"=bestmod_senes_Ldec), bestmods)


summary(bestmod_senes_Ldec)
# Quand il fait plus chaud en été : sénescence 1.2 jours plus tard quand on gagne 1°C l'été
qqnorm(resid(bestmod_senes_Ldec))
qqline(resid(bestmod_senes_Ldec))
# /!\ résidus !!!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Ldec)
# R2 des effets fixes très très nul !! Voire moins bien que le modèle altitude + année !

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Ldec)


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        variable = rownames(coef(summary(bestmod_senes_Ldec))),
                                        coef = coef(summary(bestmod_senes_Ldec))[,1],
                                        std = coef(summary(bestmod_senes_Ldec))[,2],
                                        pval = coef(summary(bestmod_senes_Ldec))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Ldec)[rownames(coef(summary(bestmod_senes_Ldec)))]))



##############################################-
#*---- Bouleau pubescent ----

senes_Bpub = senes_Alps[senes_Alps$species == "Bouleau_pubescent",]
# ggplot(senes_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Bpub = senes_Bpub[!is.na(senes_Bpub$Tmoy_GS),]


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Bpub_Altyear <- lmer(senes10 ~ altitude + year + (year|nom_zone), senes_Bpub, REML=F)
summary(mod_senes_Bpub_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Bpub_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Bpub = senes_Bpub[!is.na(senes_Bpub$debou10),]
# mod_senes_Bpub_debou10 <- lmer(senes10 ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) 
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Bpub_P_summer <- lmer(senes10 ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)

AIC(#mod_senes_Bpub_debou10, 
  mod_senes_Bpub_P_summer,mod_senes_Bpub_Tmoy_GS,mod_senes_Bpub_Tmoy_MJ,mod_senes_Bpub_Tmoy_AS,mod_senes_Bpub_Tmoy_ASO,
  mod_senes_Bpub_Tnight21j_jsenes10,mod_senes_Bpub_Tnight30j_jsenes10,mod_senes_Bpub_Tnight40j_jsenes10,
  mod_senes_Bpub_GDDinv25_jsenes10,mod_senes_Bpub_GDDinv20_jsenes10,mod_senes_Bpub_GDDinv15_jsenes10)
# Les variables de GDDinverse à 25 et 20°C donnent les meilleurs résultats (AIC = 636.2224)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Bpub_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
                                 (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) # AIC = 634.9
# visreg(mod_senes_Bpub_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes_Bpub_multivar = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + 
                                 (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) # AIC = 631.4
# visreg(mod_senes_Bpub_multivar, "GDDinv25_jsenes10", by="P_summer")
visreg(mod_senes_Bpub_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Bpub = lmer(senes10 ~ P_summer * GDDinv25_jsenes10 + altitude +
                            (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) # AIC = 627.1
bestmods = c(list("Bouleau_pubescent"=bestmod_senes_Bpub), bestmods)
visreg(bestmod_senes_Bpub, "P_summer", by="GDDinv25_jsenes10")
visreg(bestmod_senes_Bpub, "GDDinv25_jsenes10", by="P_summer")


summary(bestmod_senes_Bpub)
# gradient altitudinal : sénescence 2.0 jours plus tôt quand on monte de 100m
# précipitations : sénescence 5.8 jour plus tôt quand on gagne 1mm de pluie
# accumulation de "froid" : sénescence 6.3 jours plus tôt quand on gagne 1°C de froid (précocité d'autant plus forte qu'il fait sec)

qqnorm(resid(bestmod_senes_Bpub))
qqline(resid(bestmod_senes_Bpub))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Bpub)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Bpub)


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        variable = rownames(coef(summary(bestmod_senes_Bpub))),
                                        coef = coef(summary(bestmod_senes_Bpub))[,1],
                                        std = coef(summary(bestmod_senes_Bpub))[,2],
                                        pval = coef(summary(bestmod_senes_Bpub))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Bpub)[rownames(coef(summary(bestmod_senes_Bpub)))]))



##############################################-
#*---- Hetre ----

senes_Fsyl = senes_Alps[senes_Alps$species == "Hetre",]
# ggplot(senes_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Fsyl = senes_Fsyl[!is.na(senes_Fsyl$Tmoy_GS),]


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Fsyl_Altyear <- lmer(senes10 ~ altitude + year + (year|nom_zone), senes_Fsyl, REML=F)
summary(mod_senes_Fsyl_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Fsyl_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Fsyl = senes_Fsyl[!is.na(senes_Fsyl$debou10),]
# mod_senes_Fsyl_debou10 <- lmer(senes10 ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Fsyl_P_summer <- lmer(senes10 ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)

AIC(#mod_senes_Fsyl_debou10, 
  mod_senes_Fsyl_P_summer,mod_senes_Fsyl_Tmoy_GS,mod_senes_Fsyl_Tmoy_MJ,mod_senes_Fsyl_Tmoy_AS,mod_senes_Fsyl_Tmoy_ASO,
  mod_senes_Fsyl_Tnight21j_jsenes10,mod_senes_Fsyl_Tnight30j_jsenes10,mod_senes_Fsyl_Tnight40j_jsenes10,
  mod_senes_Fsyl_GDDinv25_jsenes10,mod_senes_Fsyl_GDDinv20_jsenes10,mod_senes_Fsyl_GDDinv15_jsenes10)
# Les variables de température estivales donnent les meilleurs résultats (AIC = 528.2573)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Fsyl_multivar = lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_ASO + 
                                 Tnight21j_jsenes10 + 
                                 GDDinv25_jsenes10 +
                                 (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F) # AIC = 526.6
# visreg(mod_senes_Fsyl_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions :
mod_senes_Fsyl_multivar = lmer(senes10 ~ P_summer * Tmoy_AS + 
                                 (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F) # AIC = 522.1
# visreg(mod_senes_Fsyl_multivar, "Tmoy_AS", by="P_summer")
# visreg(mod_senes_Fsyl_multivar, "P_summer", by="Tmoy_AS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Fsyl = lmer(senes10 ~ P_summer * Tmoy_AS + 
                            (1|nom_zone) + (P_summer|yearQ), senes_Fsyl, REML=F) # AIC = 520.1
bestmods = c(list("Hetre"=bestmod_senes_Fsyl), bestmods)
# visreg(bestmod_senes_Fsyl, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes_Fsyl, "Tmoy_AS", by="P_summer")


summary(bestmod_senes_Fsyl)
# précipitations : sénescence 28 jours plus tôt quand on gagne 1mm de pluie
# chaleur de fin de saison (août-septembre) : sénescence 20 jours plus tard quand on gagne 1°C (retard d'autant plus fort qu'il fait sec)

qqnorm(resid(bestmod_senes_Fsyl))
qqline(resid(bestmod_senes_Fsyl))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Fsyl)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Fsyl)


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        variable = rownames(coef(summary(bestmod_senes_Fsyl))),
                                        coef = coef(summary(bestmod_senes_Fsyl))[,1],
                                        std = coef(summary(bestmod_senes_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senes_Fsyl))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Fsyl)[rownames(coef(summary(bestmod_senes_Fsyl)))]))



##############################################-
#*---- Sorbier ----

senes_Sacu = senes_Alps[senes_Alps$species == "Sorbier",]
# ggplot(senes_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Sacu = senes_Sacu[!is.na(senes_Sacu$Tmoy_GS),]
# hist(senes_Sacu$senes10)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Sacu_Altyear <- lmer(senes10 ~ altitude + year + (year|nom_zone), senes_Sacu, REML=F)
summary(mod_senes_Sacu_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Sacu_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Sacu = senes_Sacu[!is.na(senes_Sacu$debou10),]
# mod_senes_Sacu_debou10 <- lmer(senes10 ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Sacu_P_summer <- lmer(senes10 ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_GS <- lmer(senes10 ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_MJ <- lmer(senes10 ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_AS <- lmer(senes10 ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_ASO <- lmer(senes10 ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight21j_jsenes10 <- lmer(senes10 ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight30j_jsenes10 <- lmer(senes10 ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight40j_jsenes10 <- lmer(senes10 ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv25_jsenes10 <- lmer(senes10 ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv20_jsenes10 <- lmer(senes10 ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv15_jsenes10 <- lmer(senes10 ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)

AIC(#mod_senes_Sacu_debou10, 
  mod_senes_Sacu_P_summer,mod_senes_Sacu_Tmoy_GS,mod_senes_Sacu_Tmoy_MJ,mod_senes_Sacu_Tmoy_AS,mod_senes_Sacu_Tmoy_ASO,
  mod_senes_Sacu_Tnight21j_jsenes10,mod_senes_Sacu_Tnight30j_jsenes10,mod_senes_Sacu_Tnight40j_jsenes10,
  mod_senes_Sacu_GDDinv25_jsenes10,mod_senes_Sacu_GDDinv20_jsenes10,mod_senes_Sacu_GDDinv15_jsenes10)
# Les variables de température nocturne donnent les meilleurs résultats (AIC = 4940.168)...
# ... ainsi que la variable de débourrement !
# on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Sacu_multivar = lmer(senes10 ~ debou10 + 
                                 Tnight21j_jsenes10 + 
                                 GDDinv25_jsenes10 +
                                 (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F) # AIC = 4057.1 # /!\ en ayant ajouté debou10, on tronque la base de données de 20% !
# visreg(mod_senes_Sacu_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions (ici ça n'améliore rien) 


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Sacu = lmer(senes10 ~ debou10 + 
                            Tnight21j_jsenes10 + 
                            GDDinv25_jsenes10 +
                            (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F) # AIC = 4057.1
# Ici l'altitude n'améliore rien, et l'interaction température/année via l'effet aléatoire ne permet plus au modèle de converger
bestmods = c(list("Sorbier"=bestmod_senes_Sacu), bestmods)
# visreg(bestmod_senes_Sacu, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes_Sacu, "Tmoy_AS", by="P_summer")


summary(bestmod_senes_Sacu)
# températures nocturnes un peu avant la sénescence : sénescence 1.6 jours plus tard quand on gagne 1°C en moyenne la nuit
# accumulation de froid (pas vraiment significatif !) : sénescence 1.0 jour plus tard quand on gagne 1°C de froid

qqnorm(resid(bestmod_senes_Sacu))
qqline(resid(bestmod_senes_Sacu))
# résidus ok !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Sacu)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Sacu)


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        variable = rownames(coef(summary(bestmod_senes_Sacu))),
                                        coef = coef(summary(bestmod_senes_Sacu))[,1],
                                        std = coef(summary(bestmod_senes_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes_Sacu))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Sacu)[rownames(coef(summary(bestmod_senes_Sacu)))]))



write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes10.Rdata")



##############################################-
# *-- DURÉE DE SÉNESCENCE                 ----
##############################################-


senes_Alps = senes_Alps_all[!is.na(senes_Alps_all$senesduree),]
# on considère qu'il n'est pas possible de passer de 10 à 50% de feuilles qui ont changé de couleur en moins de 2 jours (erreur de saisie ?)
senes_Alps = senes_Alps_all[senes_Alps_all$senesduree >= 2,] 

# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2024, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes_Alps$species), R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_bestmod_fixef = NA, R2_bestmod_allef = NA)

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()


##############################################-
#*---- Bouleau verruqueux ----

senes_Bpen = senes_Alps[senes_Alps$species == "Bouleau_verruqueux",]
# ggplot(senes_Bpen, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Bpen = senes_Bpen[!is.na(senes_Bpen$Tmoy_GS),]
# hist(senes_Bpen$senesduree)
# ggplot(senes_Bpen, aes(x=year, y=senesduree)) + geom_point() + geom_smooth(method=lm) + labs(title = "Durée de la sénescence - Bouleau", x="Année",y="Durée (jours)")

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Bpen_Altyear <- lmer(senesduree ~ altitude + year + (year|nom_zone), senes_Bpen, REML=F)
summary(mod_senes_Bpen_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Bpen_Altyear)
# Pas ouf du tout !! 


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# mod_senes_Bpen_debou10 <- lmer(senesduree ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) 
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Bpen_P_summer <- lmer(senesduree ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight21j_jsenes50 <- lmer(senesduree ~ Tnight21j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight30j_jsenes50 <- lmer(senesduree ~ Tnight30j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_Tnight40j_jsenes50 <- lmer(senesduree ~ Tnight40j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)
mod_senes_Bpen_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F)

AIC(#mod_senes_Bpen_debou10, 
  mod_senes_Bpen_P_summer,mod_senes_Bpen_Tmoy_GS,mod_senes_Bpen_Tmoy_MJ,mod_senes_Bpen_Tmoy_AS,mod_senes_Bpen_Tmoy_ASO,
  mod_senes_Bpen_Tnight21j_jsenes10,mod_senes_Bpen_Tnight30j_jsenes10,mod_senes_Bpen_Tnight40j_jsenes10,
  mod_senes_Bpen_Tnight21j_jsenes50,mod_senes_Bpen_Tnight30j_jsenes50,mod_senes_Bpen_Tnight40j_jsenes50,
  mod_senes_Bpen_GDDinv25_jsenes10,mod_senes_Bpen_GDDinv20_jsenes10,mod_senes_Bpen_GDDinv15_jsenes10)
# Les températures nocturnes pendant et un peu avant sénescence donnent les meilleurs résultats (AIC = 11755.76)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Bpen_multivar = lmer(senesduree ~ Tmoy_ASO + Tnight40j_jsenes50 + 
                                 (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) # AIC = 11753.7
# visreg(mod_senes_Bpen_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions (ici ça n'améliore pas le modèle) :
# mod_senes_Bpen_multivar = lmer(senesduree ~ P_summer * GDDinv25_jsenes10 + 
#                                  (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) 
# # visreg(mod_senes_Bpen_multivar, "GDDinv25_jsenes10", by="P_summer")
# visreg(mod_senes_Bpen_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
# Ici ça n'améliore pas le modèle
bestmod_senes_Bpen = lmer(senesduree ~ Tmoy_ASO + Tnight40j_jsenes50 + 
                            (1|nom_zone) + (1|yearQ), senes_Bpen, REML=F) # AIC = 11753.7
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes_Bpen), bestmods)


summary(bestmod_senes_Bpen)
# températures nocturnes : sénescence 1.5 jours plus longue / ralentie quand il fait 1°C plus chaud
qqnorm(resid(bestmod_senes_Bpen))
qqline(resid(bestmod_senes_Bpen))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_verruqueux", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Bpen)
# R2 mieux que sans variable climatique, mais ça reste assez nul !!

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Bpen)


resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        variable = rownames(coef(summary(bestmod_senes_Bpen))),
                                        coef = coef(summary(bestmod_senes_Bpen))[,1],
                                        std = coef(summary(bestmod_senes_Bpen))[,2],
                                        pval = coef(summary(bestmod_senes_Bpen))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Bpen)[rownames(coef(summary(bestmod_senes_Bpen)))]))



##############################################-
#*---- Meleze ----

senes_Ldec = senes_Alps[senes_Alps$species == "Meleze",]
# ggplot(senes_Ldec, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Ldec = senes_Ldec[!is.na(senes_Ldec$Tmoy_GS),]
# hist(senes_Ldec$senesduree)
# ggplot(senes_Bpen, aes(x=year, y=senesduree)) + geom_point() + geom_smooth(method=lm) + labs(title = "Durée de la sénescence - Meleze", x="Année",y="Durée (jours)")

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Ldec_Altyear <- lmer(senesduree ~ altitude + year + (year|nom_zone), senes_Ldec, REML=F)
summary(mod_senes_Ldec_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Meleze", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Ldec_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# mod_senes_Ldec_debou10 <- lmer(senesduree ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
# /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
#     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Ldec_P_summer <- lmer(senesduree ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight21j_jsenes50 <- lmer(senesduree ~ Tnight21j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight30j_jsenes50 <- lmer(senesduree ~ Tnight30j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_Tnight40j_jsenes50 <- lmer(senesduree ~ Tnight40j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)
mod_senes_Ldec_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F)

AIC(#mod_senes_Ldec_debou10, 
  mod_senes_Ldec_P_summer,mod_senes_Ldec_Tmoy_GS,mod_senes_Ldec_Tmoy_MJ,mod_senes_Ldec_Tmoy_AS,mod_senes_Ldec_Tmoy_ASO,
  mod_senes_Ldec_Tnight21j_jsenes10,mod_senes_Ldec_Tnight30j_jsenes10,mod_senes_Ldec_Tnight40j_jsenes10,
  mod_senes_Ldec_Tnight21j_jsenes50,mod_senes_Ldec_Tnight30j_jsenes50,mod_senes_Ldec_Tnight40j_jsenes50,
  mod_senes_Ldec_GDDinv25_jsenes10,mod_senes_Ldec_GDDinv20_jsenes10,mod_senes_Ldec_GDDinv15_jsenes10)
# Les températures nocturnes pendant et un peu avant sénescence donnent les meilleurs résultats (AIC = 12158.60)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets) :
mod_senes_Ldec_multivar = lmer(senesduree ~ P_summer + Tnight21j_jsenes50 + 
                                 (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) # AIC = 12155.8
# visreg(mod_senes_Ldec_multivar, "Tmoy_AS")

# # On peut aussi complexifier en ajoutant des interactions :
# mod_senes_Ldec_multivar = lmer(senesduree ~ P_summer * Tnight21j_jsenes50 + 
#                                  (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) 
# # -> ça n'améliore pas les modèles

#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Ldec = lmer(senesduree ~ P_summer + Tnight21j_jsenes50  +
                            (1|nom_zone) + (1|yearQ), senes_Ldec, REML=F) # AIC = 12152.7
bestmods = c(list("Meleze"=bestmod_senes_Ldec), bestmods)


summary(bestmod_senes_Ldec)
# températures nocturnes : sénescence 0.6 jours plus courte / accélérée quand il fait 1°C plus froid
# précipitation estivales : sénescence 0.3 jours plus longue / ralentie quand il y a 1mm de plus de précipitations.
qqnorm(resid(bestmod_senes_Ldec))
qqline(resid(bestmod_senes_Ldec))
# /!\ résidus !!!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Ldec)
# R2 des effets fixes très très nul !! Voire moins bien que le modèle altitude + année !

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Ldec)


resultats = rbind(resultats, data.frame(species = "Meleze",
                                        variable = rownames(coef(summary(bestmod_senes_Ldec))),
                                        coef = coef(summary(bestmod_senes_Ldec))[,1],
                                        std = coef(summary(bestmod_senes_Ldec))[,2],
                                        pval = coef(summary(bestmod_senes_Ldec))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Ldec)[rownames(coef(summary(bestmod_senes_Ldec)))]))



##############################################-
#*---- Bouleau pubescent ----

senes_Bpub = senes_Alps[senes_Alps$species == "Bouleau_pubescent",]
# ggplot(senes_Bpub, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Bpub = senes_Bpub[!is.na(senes_Bpub$Tmoy_GS),]
# hist(senes_Bpub$senesduree)
# ggplot(senes_Bpub, aes(x=year, y=senesduree)) + geom_point() + geom_smooth(method=lm) + labs(title = "Durée de la sénescence - Bouleau pubescent", x="Année",y="Durée (jours)")


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Bpub_Altyear <- lmer(senesduree ~ altitude + year + (year|nom_zone), senes_Bpub, REML=F)
summary(mod_senes_Bpub_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Bpub_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Bpub = senes_Bpub[!is.na(senes_Bpub$debou10),]
# mod_senes_Bpub_debou10 <- lmer(senesduree ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) 
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Bpub_P_summer <- lmer(senesduree ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight21j_jsenes50 <- lmer(senesduree ~ Tnight21j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight30j_jsenes50 <- lmer(senesduree ~ Tnight30j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_Tnight40j_jsenes50 <- lmer(senesduree ~ Tnight40j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)
mod_senes_Bpub_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F)

AIC(#mod_senes_Bpub_debou10, 
  mod_senes_Bpub_P_summer,mod_senes_Bpub_Tmoy_GS,mod_senes_Bpub_Tmoy_MJ,mod_senes_Bpub_Tmoy_AS,mod_senes_Bpub_Tmoy_ASO,
  mod_senes_Bpub_Tnight21j_jsenes10,mod_senes_Bpub_Tnight30j_jsenes10,mod_senes_Bpub_Tnight40j_jsenes10,
  mod_senes_Bpub_Tnight21j_jsenes50,mod_senes_Bpub_Tnight30j_jsenes50,mod_senes_Bpub_Tnight40j_jsenes50,
  mod_senes_Bpub_GDDinv25_jsenes10,mod_senes_Bpub_GDDinv20_jsenes10,mod_senes_Bpub_GDDinv15_jsenes10)
# La variable précip summer donne les meilleurs résultats (AIC = 344.3094)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Bpub_multivar = lmer(senesduree ~ P_summer  + 
                                 (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) # AIC = 344.3
# visreg(mod_senes_Bpub_multivar, "GDDinv25_jsenes10")

# # On peut aussi complexifier en ajoutant des interactions :
# # Ici ça n'améliore rien
# mod_senes_Bpub_multivar = lmer(senesduree ~ P_summer * GDDinv25_jsenes10 + 
#                                  (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) 
# # visreg(mod_senes_Bpub_multivar, "GDDinv25_jsenes10", by="P_summer")
# visreg(mod_senes_Bpub_multivar, "P_summer", by="GDDinv25_jsenes10")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Bpub = lmer(senesduree ~ P_summer +
                            (1|nom_zone) + (1|yearQ), senes_Bpub, REML=F) # AIC = 344.3
bestmods = c(list("Bouleau_pubescent"=bestmod_senes_Bpub), bestmods)
# visreg(bestmod_senes_Bpub, "P_summer")


summary(bestmod_senes_Bpub)
# précipitations : sénescence 0.8 jour plus courte quand on gagne 1mm de pluie (les précipitations estivales accélèrent la sénescence)
# ----------- NOTE : quid des précipitations des mois suivants (août-septembre, ou que septembre, ou que octobre ?) ---------------*

qqnorm(resid(bestmod_senes_Bpub))
qqline(resid(bestmod_senes_Bpub))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Bpub)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Bpub)


resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        variable = rownames(coef(summary(bestmod_senes_Bpub))),
                                        coef = coef(summary(bestmod_senes_Bpub))[,1],
                                        std = coef(summary(bestmod_senes_Bpub))[,2],
                                        pval = coef(summary(bestmod_senes_Bpub))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Bpub)[rownames(coef(summary(bestmod_senes_Bpub)))]))



##############################################-
#*---- Hetre ----

senes_Fsyl = senes_Alps[senes_Alps$species == "Hetre",]
# ggplot(senes_Fsyl, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Fsyl = senes_Fsyl[!is.na(senes_Fsyl$Tmoy_GS),]
# hist(senes_Fsyl$senesduree)
# ggplot(senes_Fsyl, aes(x=year, y=senesduree)) + geom_point() + geom_smooth(method=lm) + labs(title = "Durée de la sénescence - Hetre", x="Année",y="Durée (jours)")


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Fsyl_Altyear <- lmer(senesduree ~ altitude + year + (year|nom_zone), senes_Fsyl, REML=F)
summary(mod_senes_Fsyl_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Fsyl_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Fsyl = senes_Fsyl[!is.na(senes_Fsyl$debou10),]
# mod_senes_Fsyl_debou10 <- lmer(senesduree ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Fsyl_P_summer <- lmer(senesduree ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight21j_jsenes50 <- lmer(senesduree ~ Tnight21j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight30j_jsenes50 <- lmer(senesduree ~ Tnight30j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_Tnight40j_jsenes50 <- lmer(senesduree ~ Tnight40j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)
mod_senes_Fsyl_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F)

AIC(#mod_senes_Fsyl_debou10, 
  mod_senes_Fsyl_P_summer,mod_senes_Fsyl_Tmoy_GS,mod_senes_Fsyl_Tmoy_MJ,mod_senes_Fsyl_Tmoy_AS,mod_senes_Fsyl_Tmoy_ASO,
  mod_senes_Fsyl_Tnight21j_jsenes10,mod_senes_Fsyl_Tnight30j_jsenes10,mod_senes_Fsyl_Tnight40j_jsenes10,
  mod_senes_Fsyl_Tnight21j_jsenes50,mod_senes_Fsyl_Tnight30j_jsenes50,mod_senes_Fsyl_Tnight40j_jsenes50,
  mod_senes_Fsyl_GDDinv25_jsenes10,mod_senes_Fsyl_GDDinv20_jsenes10,mod_senes_Fsyl_GDDinv15_jsenes10)
# La variable précipitations estivales donne les meilleurs résultats (AIC = 434.0872)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Fsyl_multivar = lmer(senesduree ~ Tmoy_ASO + 
                                 Tnight21j_jsenes50 + 
                                 (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F) # AIC = 426.6

# # On peut aussi complexifier en ajoutant des interactions :
# # Ici ça n'améliore rien
# mod_senes_Fsyl_multivar = lmer(senesduree ~ Tnight21j_jsenes50 * Tmoy_ASO + 
#                                  (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F) 
# # visreg(mod_senes_Fsyl_multivar, "Tmoy_AS", by="P_summer")
# # visreg(mod_senes_Fsyl_multivar, "P_summer", by="Tmoy_AS")


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
# Ici ça n'améliore pas grand chose, on reste sur le modèle plus simple
bestmod_senes_Fsyl = lmer(senesduree ~ Tnight21j_jsenes50 + Tmoy_ASO + 
                            (1|nom_zone) + (1|yearQ), senes_Fsyl, REML=F) # AIC = 426.6
bestmods = c(list("Hetre"=bestmod_senes_Fsyl), bestmods)
# visreg(bestmod_senes_Fsyl, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes_Fsyl, "Tmoy_AS", by="P_summer")


summary(bestmod_senes_Fsyl)
# températures nocturnes : sénescence 9.2 jours plus longue / ralentie quand on gagne 1°C
# chaleur de fin de saison (août-septembre-octobre) : sénescence raccourcie / accélérée de 11.1 jours quand on gagne 1°C

qqnorm(resid(bestmod_senes_Fsyl))
qqline(resid(bestmod_senes_Fsyl))
# résidus pas oufs !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Fsyl)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Fsyl)


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        variable = rownames(coef(summary(bestmod_senes_Fsyl))),
                                        coef = coef(summary(bestmod_senes_Fsyl))[,1],
                                        std = coef(summary(bestmod_senes_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senes_Fsyl))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Fsyl)[rownames(coef(summary(bestmod_senes_Fsyl)))]))



##############################################-
#*---- Sorbier ----

senes_Sacu = senes_Alps[senes_Alps$species == "Sorbier",]
# ggplot(senes_Sacu, aes(x=year)) + geom_histogram(binwidth=1) 
# # on retire les lignes où on n'a pas de données de température (2005)
senes_Sacu = senes_Sacu[!is.na(senes_Sacu$Tmoy_GS),]
# hist(senes_Sacu$senesduree)

#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes_Sacu_Altyear <- lmer(senesduree ~ altitude + year + (year|nom_zone), senes_Sacu, REML=F)
summary(mod_senes_Sacu_Altyear)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senes_Sacu_Altyear)


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# senes_Sacu = senes_Sacu[!is.na(senes_Sacu$debou10),]
# mod_senes_Sacu_debou10 <- lmer(senesduree ~ debou10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
# # /!\ on ne peut pas comparer ce modèle avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
mod_senes_Sacu_P_summer <- lmer(senesduree ~ P_summer + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_GS <- lmer(senesduree ~ Tmoy_GS + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_MJ <- lmer(senesduree ~ Tmoy_MJ + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_AS <- lmer(senesduree ~ Tmoy_AS + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tmoy_ASO <- lmer(senesduree ~ Tmoy_ASO + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight21j_jsenes10 <- lmer(senesduree ~ Tnight21j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight30j_jsenes10 <- lmer(senesduree ~ Tnight30j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight40j_jsenes10 <- lmer(senesduree ~ Tnight40j_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight21j_jsenes50 <- lmer(senesduree ~ Tnight21j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight30j_jsenes50 <- lmer(senesduree ~ Tnight30j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_Tnight40j_jsenes50 <- lmer(senesduree ~ Tnight40j_jsenes50 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv25_jsenes10 <- lmer(senesduree ~ GDDinv25_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv20_jsenes10 <- lmer(senesduree ~ GDDinv20_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)
mod_senes_Sacu_GDDinv15_jsenes10 <- lmer(senesduree ~ GDDinv15_jsenes10 + (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F)

AIC(#mod_senes_Sacu_debou10, 
  mod_senes_Sacu_P_summer,mod_senes_Sacu_Tmoy_GS,mod_senes_Sacu_Tmoy_MJ,mod_senes_Sacu_Tmoy_AS,mod_senes_Sacu_Tmoy_ASO,
  mod_senes_Sacu_Tnight21j_jsenes10,mod_senes_Sacu_Tnight30j_jsenes10,mod_senes_Sacu_Tnight40j_jsenes10,
  mod_senes_Sacu_Tnight21j_jsenes50,mod_senes_Sacu_Tnight30j_jsenes50,mod_senes_Sacu_Tnight40j_jsenes50,
  mod_senes_Sacu_GDDinv25_jsenes10,mod_senes_Sacu_GDDinv20_jsenes10,mod_senes_Sacu_GDDinv15_jsenes10)
# Les variables de température nocturne donnent les meilleurs résultats (AIC = 4229.620)...
# on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senesduree ~ P_summer + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F), corr=F)
# On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC 
# En combinant ces 2 critères (AIC et significativité des effets), on retient donc :
mod_senes_Sacu_multivar = lmer(senesduree ~  P_summer +
                                 Tnight21j_jsenes10 +
                                 (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F) # AIC = 4225.5 
# visreg(mod_senes_Sacu_multivar, "GDDinv25_jsenes10")

# On peut aussi complexifier en ajoutant des interactions (ici ça n'améliore rien) 


#*-------- 3) On complexifie le modèle ---- 
# On teste en complexifiant le modèle [différents test à faire] :
# - en ajoutant un effet altitude, en interaction avec la température (ex. adaptation au fait qu'il y a plus de risque de gel tardif ?) 
# - en regardant l'interaction température / année pour l'effet aléatoire (en considérant que l'effet température ne sera pas le même tous les ans
#   s'il y a aussi des effets précipitations par exemple)
bestmod_senes_Sacu = lmer(senesduree ~ P_summer +
                            Tnight21j_jsenes10 +
                            (1|nom_zone) + (1|yearQ), senes_Sacu, REML=F) # AIC = 4225.5
# Ici l'altitude n'améliore rien, et l'interaction température/année via l'effet aléatoire ne permet plus au modèle de converger
bestmods = c(list("Sorbier"=bestmod_senes_Sacu), bestmods)
# visreg(bestmod_senes_Sacu, "P_summer", by="Tmoy_AS")
# visreg(bestmod_senes_Sacu, "Tmoy_AS", by="P_summer")


summary(bestmod_senes_Sacu)
# températures nocturnes un peu avant la sénescence : sénescence rallongée / ralentie de 0.7 jours quand on gagne 1°C en moyenne la nuit
# précipitations estivales : sénescence rallongée / ralentie de 0.3 jour quand on gagne 1mm de pluie en été

qqnorm(resid(bestmod_senes_Sacu))
qqline(resid(bestmod_senes_Sacu))
# résidus bifbof !


# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes_Sacu)

# Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes_Sacu)


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        variable = rownames(coef(summary(bestmod_senes_Sacu))),
                                        coef = coef(summary(bestmod_senes_Sacu))[,1],
                                        std = coef(summary(bestmod_senes_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes_Sacu))[,5],
                                        varexpli = calcVarPart(bestmod_senes_Sacu)[rownames(coef(summary(bestmod_senes_Sacu)))]))



write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senesduree.Rdata")

