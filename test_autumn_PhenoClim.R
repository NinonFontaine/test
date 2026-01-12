############################################################################################-
# ANALYSE DES 20 ANS DE DONNÉES PHÉNOCLIM - Phénologie automnale                        ----
#  Ninon Fontaine - printemps 2025                                                         -
############################################################################################-


library(ggmap)
library(ggplot2)
library(tidyr)
library(stringr)
library(lubridate)
library(data.table)
library(reshape2)
library(dplyr)
library(terra)
library(raster)
library(lme4)
library(lmerTest) # permet d'avoir les pvalues qui s'affichent pour les lmer
library(visreg)
library(MuMIn) # fonction r.squaredGLMM notamment + dredge pour la sélection de variables
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
#   (3) cette tendance est-elle différente selon les périodes étudiées (2006-2016 vs 2006-2025 vs 2016-2025) ?
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
  scale_y_continuous(breaks=seq(2005,2025, by=1), labels=seq(2005,2025, by=1))+
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
           ID_zone, nom_zone, dept1, region1, pays1, year #new_cat, 
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
# --> données mensuelles (précipitations totales en m) de 2004 à 2025, à une résolution de 9km (https://cds.climate.copernicus.eu/requests?tab=all)
# --> données de précipitation --> https://codes.ecmwf.int/grib/param-db/228
PRECIP<-brick("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/ERA5land_precip.grib") 
# /!\ téléchargement au 27 novembre --> on n'a pas novembre et décembre 2025 ! Donc on retire ces deux derniers mois
names(PRECIP) = paste("P",rep(2004:2025, each=12), rep(1:12,length(2004:2025)),sep="_")[1:length(names(PRECIP))]
# # Visualisation des précipitations du mois de juin 2023
# plot(PRECIP[["P_2023_6"]])

# PRECIPsummer = stackApply(PRECIP[[grep("_6|_7|_8", names(PRECIP))]], rep(2004:2025, each=3), fun = sum)
# par(mfrow = c(4,6), mar = rep(0, 4))
# for(i in 1:22){plot(PRECIPsummer[[i]], col=colorRampPalette(c("white","lightblue","blue"))(20), zlim=c(0.0001,0.05), xlim=c(6,10),ylim=c(44,48))}

# # Visualisation
# par(mfrow=c(2,2))
# for (annee in c(2004, 2010, 2023, 2024)){hist(sites_pheno_Alps[,paste0("Psummer",annee)], main=annee, xlim=c(0,0.03), breaks=10)}
# # --> année 2023 particulièrement sèche !

# sites_pheno_Alps$P_summer = apply(cbind(sites_pheno_Alps$year, 
#                                         extract(PRECIPsummer, crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIPsummer))))),
#                                   1, function(x){x[x[1]-2002]})

sites_pheno_Alps[,c("P_juin","P_juillet","P_aout", "P_septembre","P_octobre")] = 
      t(apply(cbind(sites_pheno_Alps$year, 
                    extract(PRECIP[[grep("_6", names(PRECIP))]], crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP)))),
                    extract(PRECIP[[grep("_7", names(PRECIP))]], crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP)))),
                    extract(PRECIP[[grep("_8", names(PRECIP))]], crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP)))),
                    extract(PRECIP[[grep("_9", names(PRECIP))]], crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP)))),
                    extract(PRECIP[[grep("_10", names(PRECIP))]], crds(project(vect(sites_pheno_Alps,geom=c("coord_x_4326", "coord_y_4326"), "epsg:4326"), crs(PRECIP))))),
              1, function(x){x[x[1]-c(2002,1980,1958,1936,1914)]}))

sites_pheno_Alps$P_summer = sites_pheno_Alps$P_juin + sites_pheno_Alps$P_juillet + sites_pheno_Alps$P_aout
sites_pheno_Alps$P_autumn = sites_pheno_Alps$P_septembre + sites_pheno_Alps$P_octobre



#=============================================================================================================================*
# On récupère les températures via les reconstructions basées sur les stations météo Phénoclim (cf scripts reconstruction G. Klein)
# On calcule des températures basées sur la date *médiane* de sénescence de l'année en cours


# 1) Calcul de la date de sénescence médiane par espèce, sur différentes périodes
dates_med = sites_pheno_Alps %>% group_by(species) %>%
  summarise(med_debou10_0616 = median(debou10[year %in% 2006:2016], na.rm=T),
            med_debou10_1625 = median(debou10[year %in% 2016:2025], na.rm=T),
            med_debou10_0625 = median(debou10[year %in% 2006:2025], na.rm=T),
            med_senes10_0616 = median(senes10[year %in% 2006:2016], na.rm=T),
            med_senes10_1625 = median(senes10[year %in% 2016:2025], na.rm=T),
            med_senes10_0625 = median(senes10[year %in% 2006:2025], na.rm=T),
            med_senes50_0616 = median(senes50[year %in% 2006:2016], na.rm=T),
            med_senes50_1625 = median(senes50[year %in% 2016:2025], na.rm=T),
            med_senes50_0625 = median(senes50[year %in% 2006:2025], na.rm=T))
dates_med_ann = sites_pheno_Alps %>% group_by(species, year) %>% 
  summarise(med_debou10 = median(debou10, na.rm=T),
            med_senes10 = median(senes10, na.rm=T),
            med_senes50 = median(senes50, na.rm=T))


# # 2) Calcul de la date correspondant à une photopériode seuil, en chaque point où il y a une observation
# dates_photoper = phenoclim %>% distinct(id_base_site, coord_y_4326)
# dates_photoper$P12.5 = sapply(dates_photoper$coord_y_4326, function(x){max(c(1:366)[photoperiod(1:366,x) > 12.5])})
# # /!\ PB : ce seuil de 12.5 h est atteint autour du 16/09 (doy=260), alors que la sénescence a déjà commencé à cette date pour les 5 espèces !
# #         => privilégier une date fixe, par exemple le 1/07 (doy=183), comme dans Wu et al 2018 ?


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



for(annee in c(2006:2025)){
  
  print(annee)
  
  reconstruc_Tday_sites = read.csv(paste0("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/",annee,"_Tairday_SitesPheno_reconstruit.csv"))
  reconstruc_Tday_sites = reconstruc_Tday_sites[!is.na(reconstruc_Tday_sites$station_name),]
  reconstruc_Tday_sites$julian_day = yday(as.Date(reconstruc_Tday_sites$date))
  tabT = rename(reconstruc_Tday_sites, c(jul_day=vartabT[1], sitePheno=vartabT[2], Tmoy=vartabT[3], Tmin=vartabT[4], Tmax=vartabT[5]))
  
  for (esp in unique(sites_pheno_Alps$species[sites_pheno_Alps$year == annee & !(is.na(sites_pheno_Alps$senes10) & is.na(sites_pheno_Alps$senes50))])){

    jourscibles = dates_med[dates_med$species==esp,-1]
    jourscibles_ann = dates_med_ann[dates_med_ann$species==esp & dates_med_ann$year==annee,-1]
    
    for (periode in c("0616","1625","0625","ann")){
      med_debou10 = ifelse(periode=="ann",as.numeric(jourscibles_ann$med_debou10),as.numeric(jourscibles[paste("med_debou10",periode,sep="_")]))
      med_senes10 = ifelse(periode=="ann",as.numeric(jourscibles_ann$med_senes10),as.numeric(jourscibles[paste("med_senes10",periode,sep="_")]))
      med_senes50 = ifelse(periode=="ann",as.numeric(jourscibles_ann$med_senes50),as.numeric(jourscibles[paste("med_senes50",periode,sep="_")]))
      

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
# senes_Alps_all$P_autumn = senes_Alps_all$P_septembre + senes_Alps_all$P_octobre

# /!\ VIGILANCE : les précipitations d'octobre (et d'automne) sont partiellement sur la période de sénescence, parfois après !!       ----
ggplot(senes_Alps_all) + geom_histogram(aes(x=senes10),fill="orange",alpha=0.7,binwidth=5) + geom_histogram(aes(senes50),fill="brown",alpha=0.7,binwidth=5) + facet_wrap(~species) + 
  geom_vline(xintercept = c(yday(as.Date("01/09/2025",format="%d/%m/%Y")),yday(as.Date("01/10/2025",format="%d/%m/%Y")),yday(as.Date("01/11/2025",format="%d/%m/%Y"))),lty=2)

# On renomme les variables de température calculées sur la date médiane de débourrement calculée année par année ('par défaut')
senes_Alps_all = senes_Alps_all %>% rename_with(~str_replace(.,"_ann","")) 
# senes_Alps_all = senes_Alps_all %>% rename_with(~str_replace(.,"_0625","")) # ... ou la médiane de débourrement de la période 2006-2025


# Remarque : au vu des corrélations entre ces variables de température, il est peu logique d'en intégrer plusieurs si elles sont trop corrélées // idem pour les variables de précipitations
corrplot::corrplot(cor(na.omit(senes_Alps_all[,varTselec])))
corrplot::corrplot(cor(na.omit(senes_Alps_all[,grep("P_",colnames(senes_Alps_all))])))


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
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2025, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, periode = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes10_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)
# (RQ : On pourra aussi regarder la RMSE comme indice de qualité des modèles (dans un second temps))

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()

options(na.action = "na.fail")



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
senes10_Bpen = senes10_Bpen[!is.na(senes10_Bpen$debou10),] #senes10_Bpen_all

mod_senes10_Bpen_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) 
mod_senes10_Bpen_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_P_juin <- lmer(senes10 ~ P_juin + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_P_juillet <- lmer(senes10 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_P_aout <- lmer(senes10 ~ P_aout + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
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
    mod_senes10_Bpen_P_summer,mod_senes10_Bpen_P_juin,mod_senes10_Bpen_P_juillet,mod_senes10_Bpen_P_aout,
    mod_senes10_Bpen_Tmoy_GS,mod_senes10_Bpen_Tmoy_MJ,mod_senes10_Bpen_Tmoy_AS,mod_senes10_Bpen_Tmoy_ASO,
    mod_senes10_Bpen_Tnight21j_jsenes10,mod_senes10_Bpen_Tnight30j_jsenes10,mod_senes10_Bpen_Tnight40j_jsenes10,
    mod_senes10_Bpen_GDDinv25_jsenes10,mod_senes10_Bpen_GDDinv20_jsenes10,mod_senes10_Bpen_GDDinv15_jsenes10)
# Les températures moyennes en août-septembre ou août-septembre-octobre donnent les meilleurs résultats sur le jeu de données tronqué (78% des données)...
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                               (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes10_Bpen_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
#                                (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AIC = 13422.9
# # visreg(mod_senes10_Bpen_multivar, "GDDinv25_jsenes10")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes10_Bpen_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
mod_senes10_Bpen_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F)
dredge(mod_senes10_Bpen_multivar)
mod_senes10_Bpen_multivar = lmer(senes10 ~ debou10 + Tmoy_ASO + Tnight40j_jsenes10 +
                               (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AICc = 10969.9
# La date de débourrement est retenue comme variable importante donc on considère le jeu de données tronqué


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes10 ~ debou10 * Tmoy_ASO + Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F))
summary(lmer(senes10 ~ debou10 + Tmoy_ASO * Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F))
summary(lmer(senes10 ~ debou10  * Tnight40j_jsenes10 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes10_Bpen = lmer(senes10 ~ debou10 + Tmoy_ASO + Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Bpen, REML=F) # AIC = 10974.22
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes10_Bpen), bestmods)


summary(bestmod_senes10_Bpen)
# la sénescence est un peu plus précoce lorque le débourrement l'a été, mais pas 1j = 1j --> plutôt 10j = 1j
# les températures de la fin de saison estivale / début d'automne ont un effet retard, SAUF s'il ne fait pas assez froid la nuit

# qqnorm(resid(bestmod_senes10_Bpen))
# qqline(resid(bestmod_senes10_Bpen))

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
                                        periode = "2006-2025",
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
senes10_Ldec = senes10_Ldec[!is.na(senes10_Ldec$debou10),] #senes10_Ldec_all

mod_senes10_Ldec_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) 
mod_senes10_Ldec_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_P_juin <- lmer(senes10 ~ P_juin + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_P_juillet <- lmer(senes10 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_P_aout <- lmer(senes10 ~ P_aout + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
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

AIC(mod_senes10_Ldec_debou10, mod_senes10_Ldec_P_juin,mod_senes10_Ldec_P_juillet,mod_senes10_Ldec_P_aout,
    mod_senes10_Ldec_Tmoy_GS,mod_senes10_Ldec_Tmoy_MJ,mod_senes10_Ldec_Tmoy_AS,mod_senes10_Ldec_Tmoy_ASO,
    mod_senes10_Ldec_Tnight21j_jsenes10,mod_senes10_Ldec_Tnight30j_jsenes10,mod_senes10_Ldec_Tnight40j_jsenes10,
    mod_senes10_Ldec_GDDinv25_jsenes10,mod_senes10_Ldec_GDDinv20_jsenes10,mod_senes10_Ldec_GDDinv15_jsenes10)
# Les températures moyennes de la growing season donnent les meilleurs résultats (AIC = 12227.82)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes10_Ldec_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
#                                (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) # AIC = 13422.9
# # visreg(mod_senes10_Ldec_multivar, "GDDinv25_jsenes10")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes10_Ldec_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
mod_senes10_Ldec_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
dredge(mod_senes10_Ldec_multivar)
mod_senes10_Ldec_multivar = lmer(senes10 ~ P_juillet + Tmoy_GS +
                                   (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) # AICc = 10969.9
# La date de débourrement n'est pas retenue comme variable importante donc on relance la sélection sur le jeu de données complet
senes10_Ldec = senes10_Ldec_all
mod_senes10_Ldec_multivar = lmer(senes10 ~ P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F)
dredge(mod_senes10_Ldec_multivar)
# Ce même modèle (P_juillet + Tmoy_GS) n'est pas sélectionné en premier, mais le delta-AICc est très faible (seulement 2.1 d'écart avec le meilleur modèle) --> on reste sur ce modèle-là

#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes10 ~ P_juillet * Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes10_Ldec = lmer(senes10 ~ P_juillet + Tmoy_GS + (1|ID_zone) + (1|yearQ), senes10_Ldec, REML=F) # AIC = 12224.53 / 14608.26 avec le jeu de données complet
bestmods = c(list("Meleze"=bestmod_senes10_Ldec), bestmods)


summary(bestmod_senes10_Ldec)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes10_Ldec))
# qqline(resid(bestmod_senes10_Ldec))
# # /!\ résidus !!

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
                                        periode = "2006-2025",
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
senes10_Bpub = senes10_Bpub[!is.na(senes10_Bpub$debou10),] #senes10_Bpub_all

mod_senes10_Bpub_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) 
mod_senes10_Bpub_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_P_juin <- lmer(senes10 ~ P_juin + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_P_juillet <- lmer(senes10 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_P_aout <- lmer(senes10 ~ P_aout + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
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

AIC(mod_senes10_Bpub_debou10, mod_senes10_Bpub_P_summer,mod_senes10_Bpub_P_juin,mod_senes10_Bpub_P_juillet,mod_senes10_Bpub_P_aout,
    mod_senes10_Bpub_Tmoy_GS,mod_senes10_Bpub_Tmoy_MJ,mod_senes10_Bpub_Tmoy_AS,mod_senes10_Bpub_Tmoy_ASO,
    mod_senes10_Bpub_Tnight21j_jsenes10,mod_senes10_Bpub_Tnight30j_jsenes10,mod_senes10_Bpub_Tnight40j_jsenes10,
    mod_senes10_Bpub_GDDinv25_jsenes10,mod_senes10_Bpub_GDDinv20_jsenes10,mod_senes10_Bpub_GDDinv15_jsenes10)
# Les températures moyennes de mai-juin, aout-septembre ou les GDDinverse donnent les meilleurs résultats (AIC = 304.6)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes10_Bpub_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
#                                (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AIC = 13422.9
# # visreg(mod_senes10_Bpub_multivar, "GDDinv25_jsenes10")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes10_Bpub_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
mod_senes10_Bpub_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS +
                                   Tnight30j_jsenes10  +
                                   GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
dredge(mod_senes10_Bpub_multivar)
mod_senes10_Bpub_multivar = lmer(senes10 ~ GDDinv20_jsenes10 + P_juin + Tmoy_AS +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AICc = 295.2915
# La date de débourrement n'est pas retenue comme variable importante donc on relance la sélection sur le jeu de données complet
senes10_Bpub = senes10_Bpub_all
mod_senes10_Bpub_multivar = lmer(senes10 ~ P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS +
                                   Tnight30j_jsenes10  +
                                   GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F)
dredge(mod_senes10_Bpub_multivar)
# Ce même modèle (GDDinv20_jsenes10 + P_juin + Tmoy_AS) n'est pas sélectionné en premier, mais GDDinv20_jsenes10 + P_summer + Tmoy_AS oui, et c'était le second dans la version avec jeu de donnée tronqué --> on le sélectionne


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes10 ~ GDDinv20_jsenes10 * P_summer + Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F))
summary(lmer(senes10 ~ GDDinv20_jsenes10 + P_summer * Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F))
summary(lmer(senes10 ~ GDDinv20_jsenes10 * Tmoy_AS + P_summer + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes10_Bpub = lmer(senes10 ~ GDDinv20_jsenes10 + P_summer + Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Bpub, REML=F) # AIC = 291.9859 / 595.8356 avec le jeu de données complet
bestmods = c(list("Bouleau_pubescent"=bestmod_senes10_Bpub), bestmods)


summary(bestmod_senes10_Bpub)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes10_Bpub))
# qqline(resid(bestmod_senes10_Bpub))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Bpub)
# R2 des effets fixes plutôt nul !!
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
                                        periode = "2006-2025",
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
senes10_Fsyl = senes10_Fsyl[!is.na(senes10_Fsyl$Tmoy_GS),] %>% filter(!is.na(altitude)) #on élimine le point sans altitude
ggplot(senes10_Fsyl, aes(x=senes10)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes10_Fsyl_Altyear <- lmer(senes10 ~ altitude + year + (year|ID_zone), senes10_Fsyl, REML=F)
summary(mod_senes10_Fsyl_Altyear)
altyearmods = c(list("Hetre"=mod_senes10_Fsyl_Altyear), altyearmods)

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
senes10_Fsyl = senes10_Fsyl[!is.na(senes10_Fsyl$debou10),] #senes10_Fsyl_all

mod_senes10_Fsyl_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) 
mod_senes10_Fsyl_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_P_juin <- lmer(senes10 ~ P_juin + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_P_juillet <- lmer(senes10 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_P_aout <- lmer(senes10 ~ P_aout + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
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

AIC(mod_senes10_Fsyl_debou10, mod_senes10_Fsyl_P_summer,mod_senes10_Fsyl_P_juin,mod_senes10_Fsyl_P_juillet,mod_senes10_Fsyl_P_aout,
    mod_senes10_Fsyl_Tmoy_GS,mod_senes10_Fsyl_Tmoy_MJ,mod_senes10_Fsyl_Tmoy_AS,mod_senes10_Fsyl_Tmoy_ASO,
    mod_senes10_Fsyl_Tnight21j_jsenes10,mod_senes10_Fsyl_Tnight30j_jsenes10,mod_senes10_Fsyl_Tnight40j_jsenes10,
    mod_senes10_Fsyl_GDDinv25_jsenes10,mod_senes10_Fsyl_GDDinv20_jsenes10,mod_senes10_Fsyl_GDDinv15_jsenes10)
# Les températures minimales donnent les meilleurs résultats (AIC = 304.6)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes10_Fsyl_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
#                                (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AIC = 13422.9
# # visreg(mod_senes10_Fsyl_multivar, "GDDinv25_jsenes10")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes10_Fsyl_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
mod_senes10_Fsyl_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS +
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10  +
                                    GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
dredge(mod_senes10_Fsyl_multivar)
# La date de débourrement n'est pas retenue comme variable importante donc on relance la sélection sur le jeu de données complet
senes10_Fsyl = senes10_Fsyl_all
mod_senes10_Fsyl_multivar = lmer(senes10 ~ P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS +
                                   Tnight30j_jsenes10  +
                                   (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F)
dredge(mod_senes10_Fsyl_multivar)
mod_senes10_Fsyl_multivar = lmer(senes10 ~ Tmoy_MJ + Tmoy_AS + Tnight30j_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AICc = 569.4133
# /!\ difficile de faire une sélection convenable...!


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes10 ~ Tmoy_MJ * Tmoy_AS + Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F))
summary(lmer(senes10 ~ Tmoy_MJ + Tmoy_AS * Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F))
summary(lmer(senes10 ~ Tmoy_MJ * Tnight30j_jsenes10 + Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes10_Fsyl = lmer(senes10 ~ Tmoy_MJ * Tnight30j_jsenes10 + Tmoy_AS + (1|ID_zone) + (1|yearQ), senes10_Fsyl, REML=F) # AIC = 555.4 avec le jeu de données complet
bestmods = c(list("Hetre"=bestmod_senes10_Fsyl), bestmods)


summary(bestmod_senes10_Fsyl)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes10_Fsyl))
# qqline(resid(bestmod_senes10_Fsyl))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Fsyl)
# R2 des effets fixes plutôt nul !!
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
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        periode = "2006-2025",
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
altyearmods = c(list("Sorbier"=mod_senes10_Sacu_Altyear), altyearmods)

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
senes10_Sacu = senes10_Sacu[!is.na(senes10_Sacu$debou10),] #senes10_Sacu_all

mod_senes10_Sacu_debou10 <- lmer(senes10 ~ debou10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) 
mod_senes10_Sacu_P_summer <- lmer(senes10 ~ P_summer + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_P_juin <- lmer(senes10 ~ P_juin + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_P_juillet <- lmer(senes10 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_P_aout <- lmer(senes10 ~ P_aout + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
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

AIC(mod_senes10_Sacu_debou10, mod_senes10_Sacu_P_summer,mod_senes10_Sacu_P_juin,mod_senes10_Sacu_P_juillet,mod_senes10_Sacu_P_aout,
    mod_senes10_Sacu_Tmoy_GS,mod_senes10_Sacu_Tmoy_MJ,mod_senes10_Sacu_Tmoy_AS,mod_senes10_Sacu_Tmoy_ASO,
    mod_senes10_Sacu_Tnight21j_jsenes10,mod_senes10_Sacu_Tnight30j_jsenes10,mod_senes10_Sacu_Tnight40j_jsenes10,
    mod_senes10_Sacu_GDDinv25_jsenes10,mod_senes10_Sacu_GDDinv20_jsenes10,mod_senes10_Sacu_GDDinv15_jsenes10)
# La date de débourrement et les précipitations du mois d'août donnent les meilleurs résultats (AIC = 4302)... on regarde en version multivariables


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
               (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes10_Sacu_multivar = lmer(senes10 ~ P_summer + GDDinv25_jsenes10 + 
#                                (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) # AIC = 13422.9
# # visreg(mod_senes10_Sacu_multivar, "GDDinv25_jsenes10")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes10_Sacu_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + P_aout + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                                   GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
mod_senes10_Sacu_multivar = lmer(senes10 ~ debou10 + P_summer + P_juin + P_juillet + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO +
                                   Tnight30j_jsenes10  + Tnight40j_jsenes10  +
                                   GDDinv20_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F)
dredge(mod_senes10_Sacu_multivar)
mod_senes10_Sacu_multivar = lmer(senes10 ~ debou10 + GDDinv20_jsenes10 + Tnight30j_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) # AICc = 4300.2
# La date de débourrement est retenue comme importante donc on fonctionne sur les données tronquées

#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes10 ~ debou10 * GDDinv20_jsenes10 + Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F))
summary(lmer(senes10 ~ debou10 + GDDinv20_jsenes10 * Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F))
summary(lmer(senes10 ~ debou10 * Tnight30j_jsenes10 + GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes10_Sacu = lmer(senes10 ~ debou10 + GDDinv20_jsenes10 + Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senes10_Sacu, REML=F) # AIC = 291.9859 / 595.8356 avec le jeu de données complet
bestmods = c(list("Sorbier"=bestmod_senes10_Sacu), bestmods)


summary(bestmod_senes10_Sacu)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes10_Sacu))
# qqline(resid(bestmod_senes10_Sacu))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes10_Sacu)
# R2 des effets fixes plutôt nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes10_Sacu)
# - Validation croisée
n = nrow(senes10_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes10_Sacu[trainIndex ,]
test <- senes10_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes10_Sacu$ID_zone[drop=T])[!unique(senes10_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes10_Sacu), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senes10, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senes10_Sacu))),
                                        coef = coef(summary(bestmod_senes10_Sacu))[,1],
                                        std = coef(summary(bestmod_senes10_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes10_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes10_Sacu)[rownames(coef(summary(bestmod_senes10_Sacu)))], error = function(e) return(NA))))




##############################################-
#*---- Enregistrement des résultats ----


write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10_20062025__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10_20062025__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes10_20062025.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senes10_20062025_v1altyear.Rdata")




# ##############################################-
# # *-- 50% DE SÉNESCENCE                  ----
# ##############################################-
# 
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
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2025, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, periode = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senes50_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)
# (RQ : On pourra aussi regarder la RMSE comme indice de qualité des modèles (dans un second temps))

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()

options(na.action = "na.fail")



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
senes50_Bpen = senes50_Bpen[!is.na(senes50_Bpen$debou10),] #senes50_Bpen_all

mod_senes50_Bpen_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) 
mod_senes50_Bpen_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_juin <- lmer(senes50 ~ P_juin + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_juillet <- lmer(senes50 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_aout <- lmer(senes50 ~ P_aout + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_septembre <- lmer(senes50 ~ P_septembre + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_octobre <- lmer(senes50 ~ P_octobre + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_P_autumn <- lmer(senes50 ~ P_autumn + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
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
    mod_senes50_Bpen_P_summer, mod_senes50_Bpen_P_autumn,mod_senes50_Bpen_P_juin,mod_senes50_Bpen_P_juillet,mod_senes50_Bpen_P_aout,mod_senes50_Bpen_P_septembre,mod_senes50_Bpen_P_octobre,
    mod_senes50_Bpen_Tmoy_GS,mod_senes50_Bpen_Tmoy_MJ,mod_senes50_Bpen_Tmoy_AS,mod_senes50_Bpen_Tmoy_ASO,
    mod_senes50_Bpen_Tnight21j_jsenes50,mod_senes50_Bpen_Tnight30j_jsenes50,mod_senes50_Bpen_Tnight40j_jsenes50,
    mod_senes50_Bpen_GDDinv25_jsenes50,mod_senes50_Bpen_GDDinv20_jsenes50,mod_senes50_Bpen_GDDinv15_jsenes50)
# Les précipitations d'octobre donnent les meilleurs résultats sur le jeu de données tronqué... MAIS 22% des senes50 sont avant le mois d'octobre !
# sinon, ce osnt les GDD inversés qui donnent les meilleurs résultats
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes50_Bpen_multivar = lmer(senes50 ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AIC = 13422.9
# # visreg(mod_senes50_Bpen_multivar, "GDDinv25_jsenes50")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes50_Bpen_multivar = lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
mod_senes50_Bpen_multivar = lmer(senes50 ~ debou10 + P_summer + P_autumn + P_octobre + Tmoy_MJ + Tmoy_AS + #Tmoy_ASO + #P_juin + P_aout + P_septembre +
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 +
                                   GDDinv25_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F)
dredge(mod_senes50_Bpen_multivar)
mod_senes50_Bpen_multivar = lmer(senes50 ~ debou10 + Tmoy_AS + P_octobre + Tnight21j_jsenes50 +GDDinv25_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AICc = 10969.9
# La date de débourrement est retenue comme variable importante donc on considère le jeu de données tronqué


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
AICc(lmer(senes50 ~ debou10 * Tmoy_AS + P_octobre + Tnight21j_jsenes50 + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 * P_octobre + Tmoy_AS + Tnight21j_jsenes50 + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 * Tnight21j_jsenes50 + Tmoy_AS + P_octobre + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 * GDDinv25_jsenes50 + Tmoy_AS + P_octobre + Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS * P_octobre + Tnight21j_jsenes50 + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS * Tnight21j_jsenes50 + P_octobre + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS * GDDinv25_jsenes50 + P_octobre + Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS + P_octobre * Tnight21j_jsenes50 + GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS + P_octobre * GDDinv25_jsenes50 + Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
AICc(lmer(senes50 ~ debou10 + Tmoy_AS + P_octobre + Tnight21j_jsenes50 * GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F))
# ajouter des interactions améliore le modèle (AICc)

bestmod_senes50_Bpen = lmer(senes50 ~ debou10 + Tmoy_AS + P_octobre * GDDinv25_jsenes50 + Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpen, REML=F) # AIC = 10728.3
visreg(bestmod_senes50_Bpen, "P_octobre", by="GDDinv25_jsenes50")
bestmods = c(list("Bouleau_verruqueux"=bestmod_senes50_Bpen), bestmods)


summary(bestmod_senes50_Bpen)
# la sénescence est un peu plus précoce lorque le débourrement l'a été, mais pas 1j = 1j --> plutôt 10j = 0.8j
# les températures de la fin de saison estivale / début d'automne ont un effet retard
# les précipitations d'octobre ont un effet avancement, d'autant plus marqué s'il fait plus chaud (le froid ralenti le processus)
# les GDD inversés ont un effet avancement --> plus il fait froid = plus il y a accumulation de froid = plus la sénescence est précoce (mais effet inverse pour
#     les températures minimales = nocturnes : plus il fait froid la nuit, plus la sénescence est tardive !)

# qqnorm(resid(bestmod_senes50_Bpen))
# qqline(resid(bestmod_senes50_Bpen))

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
                                        periode = "2006-2025",
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
# senes50_Ldec = senes50_Ldec[!is.na(senes50_Ldec$debou10),] #senes50_Ldec_all
senes50_Ldec = senes50_Ldec_all

mod_senes50_Ldec_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) 
mod_senes50_Ldec_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_juin <- lmer(senes50 ~ P_juin + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_juillet <- lmer(senes50 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_aout <- lmer(senes50 ~ P_aout + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_septembre <- lmer(senes50 ~ P_septembre + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_octobre <- lmer(senes50 ~ P_octobre + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_P_autumn <- lmer(senes50 ~ P_autumn + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
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
    mod_senes50_Ldec_P_summer, mod_senes50_Ldec_P_autumn,mod_senes50_Ldec_P_juin,mod_senes50_Ldec_P_juillet,mod_senes50_Ldec_P_aout,mod_senes50_Ldec_P_septembre,mod_senes50_Ldec_P_octobre,
    mod_senes50_Ldec_Tmoy_GS,mod_senes50_Ldec_Tmoy_MJ,mod_senes50_Ldec_Tmoy_AS,mod_senes50_Ldec_Tmoy_ASO,
    mod_senes50_Ldec_Tnight21j_jsenes50,mod_senes50_Ldec_Tnight30j_jsenes50,mod_senes50_Ldec_Tnight40j_jsenes50,
    mod_senes50_Ldec_GDDinv25_jsenes50,mod_senes50_Ldec_GDDinv20_jsenes50,mod_senes50_Ldec_GDDinv15_jsenes50)
# Les températures estivales donnent les meilleurs résultats sur le jeu de données tronqué...
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet (parce que debou10 ne semble pas ressortir)
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senes50_Ldec_multivar = lmer(senes50 ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
# # visreg(mod_senes50_Ldec_multivar, "GDDinv25_jsenes50")

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# La date de débourrement n'est pas retenue comme variable importante donc on relance la sélection sur le jeu de données complet
mod_senes50_Ldec_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
mod_senes50_Ldec_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_septembre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv20_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F)
dredge(mod_senes50_Ldec_multivar)
mod_senes50_Ldec_multivar = lmer(senes50 ~ P_autumn + P_juin + Tmoy_MJ +
                                   (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) # AICc = 13176.8



#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes50 ~ P_autumn * P_juin + Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F))
summary(lmer(senes50 ~ P_autumn + P_juin * Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F))
summary(lmer(senes50 ~ P_autumn * Tmoy_MJ + P_juin + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes50_Ldec = lmer(senes50 ~ P_autumn + P_juin + Tmoy_MJ + (1|ID_zone) + (1|yearQ), senes50_Ldec, REML=F) # AIC = 12224.53 / 14608.26 avec le jeu de données complet
bestmods = c(list("Meleze"=bestmod_senes50_Ldec), bestmods)


summary(bestmod_senes50_Ldec)
# la sénescence est un plus précoce quand il fait plus sec en juin (1.1j d'avance par mm de pluie en moins), et plus précoce quand c'est plus humide à l'automne
# les températures de début de saison (mai-juin) ont un effet retard --> 2.4j de plus s'il fait 1°C de plus

# qqnorm(resid(bestmod_senes50_Ldec))
# qqline(resid(bestmod_senes50_Ldec))
# # /!\ résidus !!

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
                                        periode = "2006-2025",
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
senes50_Bpub = senes50_Bpub[!is.na(senes50_Bpub$debou10),] #senes50_Bpub_all

mod_senes50_Bpub_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) 
mod_senes50_Bpub_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_juin <- lmer(senes50 ~ P_juin + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_juillet <- lmer(senes50 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_aout <- lmer(senes50 ~ P_aout + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_septembre <- lmer(senes50 ~ P_septembre + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_octobre <- lmer(senes50 ~ P_octobre + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_P_autumn <- lmer(senes50 ~ P_autumn + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
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
    mod_senes50_Bpub_P_summer, mod_senes50_Bpub_P_autumn,mod_senes50_Bpub_P_juin,mod_senes50_Bpub_P_juillet,mod_senes50_Bpub_P_aout,mod_senes50_Bpub_P_septembre,mod_senes50_Bpub_P_octobre,
    mod_senes50_Bpub_Tmoy_GS,mod_senes50_Bpub_Tmoy_MJ,mod_senes50_Bpub_Tmoy_AS,mod_senes50_Bpub_Tmoy_ASO,
    mod_senes50_Bpub_Tnight21j_jsenes50,mod_senes50_Bpub_Tnight30j_jsenes50,mod_senes50_Bpub_Tnight40j_jsenes50,
    mod_senes50_Bpub_GDDinv25_jsenes50,mod_senes50_Bpub_GDDinv20_jsenes50,mod_senes50_Bpub_GDDinv15_jsenes50)
# Les températures moyennes pendant la saison de végétation donnent les meilleurs résultats sur le jeu de données tronqué...
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F), corr=F)
# debou0 ne sort pas, on travaille avec le jeu de données complet
senes50_Bpub = senes50_Bpub_all

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes50_Bpub_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
mod_senes50_Bpub_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_septembre + Tmoy_GS + Tmoy_AS +
                                  Tnight40j_jsenes50  +
                                   GDDinv20_jsenes50  +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F)
dredge(mod_senes50_Bpub_multivar)
mod_senes50_Bpub_multivar = lmer(senes50 ~ P_autumn + Tnight40j_jsenes50 + Tmoy_ASO +
                                   (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) # AICc = 457.5


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes50 ~ P_autumn * Tnight40j_jsenes50 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F))
summary(lmer(senes50 ~ P_autumn + Tnight40j_jsenes50 * Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F))
summary(lmer(senes50 ~ P_autumn * Tmoy_ASO + Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes50_Bpub = lmer(senes50 ~ P_autumn + Tnight40j_jsenes50 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Bpub, REML=F) # AIC = 291.9859 / 595.8356 avec le jeu de données complet
bestmods = c(list("Bouleau_pubescent"=bestmod_senes50_Bpub), bestmods)


summary(bestmod_senes50_Bpub)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes50_Bpub))
# qqline(resid(bestmod_senes50_Bpub))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Bouleau_pubescent", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Bpub)
# R2 des effets fixes plutôt nul !!
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
                                        periode = "2006-2025",
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
senes50_Fsyl = senes50_Fsyl[!is.na(senes50_Fsyl$Tmoy_GS),] %>% filter(!is.na(altitude)) #on élimine le point sans altitude
ggplot(senes50_Fsyl, aes(x=senes50)) + geom_histogram(binwidth=5)


#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senes50_Fsyl_Altyear <- lmer(senes50 ~ altitude + year + (year|ID_zone), senes50_Fsyl, REML=F)
summary(mod_senes50_Fsyl_Altyear)
altyearmods = c(list("Hetre"=mod_senes50_Fsyl_Altyear), altyearmods)

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
senes50_Fsyl = senes50_Fsyl[!is.na(senes50_Fsyl$debou10),] #senes50_Fsyl_all

mod_senes50_Fsyl_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) 
mod_senes50_Fsyl_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_juin <- lmer(senes50 ~ P_juin + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_juillet <- lmer(senes50 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_aout <- lmer(senes50 ~ P_aout + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_septembre <- lmer(senes50 ~ P_septembre + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_octobre <- lmer(senes50 ~ P_octobre + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_P_autumn <- lmer(senes50 ~ P_autumn + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
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
    mod_senes50_Fsyl_P_summer, mod_senes50_Fsyl_P_autumn,mod_senes50_Fsyl_P_juin,mod_senes50_Fsyl_P_juillet,mod_senes50_Fsyl_P_aout,mod_senes50_Fsyl_P_septembre,mod_senes50_Fsyl_P_octobre,
    mod_senes50_Fsyl_Tmoy_GS,mod_senes50_Fsyl_Tmoy_MJ,mod_senes50_Fsyl_Tmoy_AS,mod_senes50_Fsyl_Tmoy_ASO,
    mod_senes50_Fsyl_Tnight21j_jsenes50,mod_senes50_Fsyl_Tnight30j_jsenes50,mod_senes50_Fsyl_Tnight40j_jsenes50,
    mod_senes50_Fsyl_GDDinv25_jsenes50,mod_senes50_Fsyl_GDDinv20_jsenes50,mod_senes50_Fsyl_GDDinv15_jsenes50)
# Les précipitations du mois de septembre donnent les meilleurs résultats sur le jeu de données tronqué...
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F), corr=F)

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes50_Fsyl_multivar = lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
mod_senes50_Fsyl_multivar = lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_septembre + Tmoy_GS + Tmoy_MJ + Tmoy_AS +
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50  +
                                   GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
dredge(mod_senes50_Fsyl_multivar)
# La date de débourrement n'est pas retenue comme variable importante donc on relance la sélection sur le jeu de données complet
senes50_Fsyl = senes50_Fsyl_all
mod_senes50_Fsyl_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_septembre + Tmoy_MJ + Tmoy_AS +
                                   Tnight30j_jsenes50  +
                                   GDDinv15_jsenes50  +
                                   (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F)
dredge(mod_senes50_Fsyl_multivar)
mod_senes50_Fsyl_multivar = lmer(senes50 ~ P_septembre + Tmoy_AS  + # avec Tmoy_MJ en plus l'AICc est mieux, mais la corrélation est trop importante entre MJ et AS
                                   (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) # AICc = 492.6
# /!\ difficile de faire une sélection convenable...!


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(lmer(senes50 ~ P_septembre * Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F))
# ajouter des interactions n'améliore pas le modèle

bestmod_senes50_Fsyl = lmer(senes50 ~ P_septembre + Tmoy_AS + (1|ID_zone) + (1|yearQ), senes50_Fsyl, REML=F) # AIC = 555.4 avec le jeu de données complet
bestmods = c(list("Hetre"=bestmod_senes50_Fsyl), bestmods)


summary(bestmod_senes50_Fsyl)
# la sénescence est un plus précoce quand il fait plus humide en septembre (2.1j d'avance par mm de pluie en plus)
# les températures de fin d'été (août-septembre) ont un effet retard --> 2.5j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes50_Fsyl))
# qqline(resid(bestmod_senes50_Fsyl))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Hetre", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Fsyl)
# R2 des effets fixes plutôt nul !!
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
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Hetre",
                                        periode = "2006-2025",
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
altyearmods = c(list("Sorbier"=mod_senes50_Sacu_Altyear), altyearmods)

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
senes50_Sacu = senes50_Sacu[!is.na(senes50_Sacu$debou10),] #senes50_Sacu_all

mod_senes50_Sacu_debou10 <- lmer(senes50 ~ debou10 + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) 
mod_senes50_Sacu_P_summer <- lmer(senes50 ~ P_summer + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_juin <- lmer(senes50 ~ P_juin + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_juillet <- lmer(senes50 ~ P_juillet + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_aout <- lmer(senes50 ~ P_aout + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_septembre <- lmer(senes50 ~ P_septembre + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_octobre <- lmer(senes50 ~ P_octobre + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
mod_senes50_Sacu_P_autumn <- lmer(senes50 ~ P_autumn + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
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
    mod_senes50_Sacu_P_summer, mod_senes50_Sacu_P_autumn,mod_senes50_Sacu_P_juin,mod_senes50_Sacu_P_juillet,mod_senes50_Sacu_P_aout,mod_senes50_Sacu_P_septembre,mod_senes50_Sacu_P_octobre,
    mod_senes50_Sacu_Tmoy_GS,mod_senes50_Sacu_Tmoy_MJ,mod_senes50_Sacu_Tmoy_AS,mod_senes50_Sacu_Tmoy_ASO,
    mod_senes50_Sacu_Tnight21j_jsenes50,mod_senes50_Sacu_Tnight30j_jsenes50,mod_senes50_Sacu_Tnight40j_jsenes50,
    mod_senes50_Sacu_GDDinv25_jsenes50,mod_senes50_Sacu_GDDinv20_jsenes50,mod_senes50_Sacu_GDDinv15_jsenes50)
# Les Tmoy_ASO donnent les meilleurs résultats sur le jeu de données tronqué...
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F), corr=F)

# /!\ on peut aussi utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
mod_senes50_Sacu_multivar = lmer(senes50 ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                   Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                                   GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
# La date de débourrement n'est pas retenue comme importante donc on fonctionne sur les données complètes
senes50_Sacu = senes50_Sacu_all
mod_senes50_Sacu_multivar = lmer(senes50 ~ P_summer + P_autumn + P_juin + P_juillet + P_septembre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO +
                                   Tnight30j_jsenes50 +
                                   GDDinv20_jsenes50 +
                                   (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F)
dredge(mod_senes50_Sacu_multivar)
mod_senes50_Sacu_multivar = lmer(senes50 ~ Tmoy_ASO +
                                   (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) # AICc = 4300.2

#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
# pas de sens avec une seule variable

bestmod_senes50_Sacu = lmer(senes50 ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senes50_Sacu, REML=F) # AIC = 291.9859 / 595.8356 avec le jeu de données complet
bestmods = c(list("Sorbier"=bestmod_senes50_Sacu), bestmods)


summary(bestmod_senes50_Sacu)
# la sénescence est un plus précoce quand il fait plus sec en juillet (0.8j d'avance par mm de pluie en moins)
# les températures de la saison de croissance (entre débourrement et sénescence) ont un effet retard --> 2j de plus s'il fait en moyenne 1°C de plus

# qqnorm(resid(bestmod_senes50_Sacu))
# qqline(resid(bestmod_senes50_Sacu))
# # /!\ résidus !!

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Sorbier", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senes50_Sacu)
# R2 des effets fixes plutôt nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senes50_Sacu)
# - Validation croisée
n = nrow(senes50_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senes50_Sacu[trainIndex ,]
test <- senes50_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senes50_Sacu$ID_zone[drop=T])[!unique(senes50_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- lmer(formula(bestmod_senes50_Sacu), train)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senes50, test2))[["adj.r.squared"]]


resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senes50_Sacu))),
                                        coef = coef(summary(bestmod_senes50_Sacu))[,1],
                                        std = coef(summary(bestmod_senes50_Sacu))[,2],
                                        pval = coef(summary(bestmod_senes50_Sacu))[,5],
                                        varexpli = tryCatch(calcVarPart(bestmod_senes50_Sacu)[rownames(coef(summary(bestmod_senes50_Sacu)))], error = function(e) return(NA))))





##############################################-
#*---- Enregistrement des résultats ----


write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes50_20062025__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes50_20062025__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes50_20062025.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senes50_20062025_v1altyear.Rdata")




##############################################-
# *-- DURÉE DE SÉNESCENCE                 ----
##############################################-


senesduree_Alps = senes_Alps_all[!is.na(senes_Alps_all$senesduree),]
# on considère qu'il n'est pas possible de passer de 10 à 50% de feuilles qui ont changé de couleur en moins de 2 jours (erreur de saisie ?)
senesduree_Alps = senesduree_Alps[senesduree_Alps$senesduree >= 2,] 

# On met les précipitations en mm plutôt qu'en m, sinon ça fait trop de différence d'échelle entre les variables
senesduree_Alps = senesduree_Alps %>% mutate(across(starts_with("P_"), ~.x*1000))

ggplot(senesduree_Alps, aes(x=senesduree)) + geom_histogram() + facet_wrap(~species)
# C'est quand même étonnant d'avoir de la sénescence qui dure plus de 2 mois !! (hêtre, bouleau verruqueux, mélèze, sorbier)
# => OK d'après Colin (ex d'un été sec, où la sénescence commence très tôt, puis se poursuit plus tard)
#
# Pour modéliser cette "durée avant un événement", on utilise une loi de Poisson (cf distribution des données non-normale)



# # Aperçu des données de sénescence
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=altitude, y=senesduree, col=species)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senesduree, col=species)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~year) + geom_smooth(method=lm)
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=Tmoy_GS, y=senesduree, col=yearQ)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm) 
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=GDDinv20_jsenes10, y=senesduree, col=yearQ)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm)
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=P_summer, y=senesduree)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm)
# ggplot(senesduree_Alps[!is.na(senesduree_Alps$Tmoy_GS),], aes(x=debou10, y=senesduree)) + geom_point() +
#   theme(legend.position = "none") + facet_wrap(.~species) + geom_smooth(method=lm)

corrplot::corrplot(cor(na.omit(senes_Alps_all[,varTselec])))

# Initialisation d'un tableau de résultat (coefficients du meilleur modèle pour chaque espèce)
# (RQ : pour le débourrement, on avait comparé les résultats sur les périodes 2006-2016 et 2006-2025, en lien avec le papier de Bison et al 2019. Pour 
#       le changement de couleur des feuilles, on ne regarde pas cet aspect (dans un premier temps))
resultats = data.frame(species = NA, periode = NA, variable = NA, coef = NA, std=NA, pval=NA, varexpli=NA)[0,]

# Initialisation d'un tableau pour comparer l'intérêt d'avoir des modèles intégrant la température, par rapport aux modèles avec année et altitude
R2_models = data.frame(species = unique(senesduree_Alps$species), 
                       R2_altyear_fixef = NA, R2_altyear_allef = NA, R2_altyear_calibval = NA, 
                       R2_bestmod_fixef = NA, R2_bestmod_allef = NA, R2_bestmod_calibval = NA)
# (RQ : On pourra aussi regarder la RMSE comme indice de qualité des modèles (dans un second temps))

# Initialisation d'une liste avec tous les meilleurs modèles
bestmods = list()
altyearmods = list()

options(na.action = "na.fail")



##############################################-
#*---- Bouleau verruqueux ----

senesduree_Bpen = senesduree_Alps[senesduree_Alps$species == "Bouleau_verruqueux",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Bpen = senesduree_Bpen[!is.na(senesduree_Bpen$Tmoy_GS),]
ggplot(senesduree_Bpen, aes(x=senesduree)) + geom_histogram(binwidth=5)



#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Bpen_Altyear <- glmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Bpen, family=poisson)
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
mod_train <- glmer(formula(mod_senesduree_Bpen_Altyear), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Bpen_all = senesduree_Bpen
senesduree_Bpen = senesduree_Bpen[!is.na(senesduree_Bpen$debou10),] #senesduree_Bpen_all

mod_senesduree_Bpen_debou10 <- glmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson) 
mod_senesduree_Bpen_P_summer <- glmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_juin <- glmer(senesduree ~ P_juin + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_juillet <- glmer(senesduree ~ P_juillet + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_aout <- glmer(senesduree ~ P_aout + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_septembre <- glmer(senesduree ~ P_septembre+ (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_octobre <- glmer(senesduree ~ P_octobre + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_P_autumn <- glmer(senesduree ~ P_autumn + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tmoy_GS <- glmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tmoy_MJ <- glmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tmoy_AS <- glmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tmoy_ASO <- glmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight21j_jsenes10 <- glmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight30j_jsenes10 <- glmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight40j_jsenes10 <- glmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv25_jsenes10 <- glmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv20_jsenes10 <- glmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv15_jsenes10 <- glmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight21j_jsenes50 <- glmer(senesduree ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight30j_jsenes50 <- glmer(senesduree ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_Tnight40j_jsenes50 <- glmer(senesduree ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv25_jsenes50 <- glmer(senesduree ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv20_jsenes50 <- glmer(senesduree ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
mod_senesduree_Bpen_GDDinv15_jsenes50 <- glmer(senesduree ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)

AIC(mod_senesduree_Bpen_debou10, 
    mod_senesduree_Bpen_P_summer,mod_senesduree_Bpen_P_juin,mod_senesduree_Bpen_P_juillet,mod_senesduree_Bpen_P_aout,mod_senesduree_Bpen_P_septembre,mod_senesduree_Bpen_P_octobre,mod_senesduree_Bpen_P_autumn,
    mod_senesduree_Bpen_Tmoy_GS,mod_senesduree_Bpen_Tmoy_MJ,mod_senesduree_Bpen_Tmoy_AS,mod_senesduree_Bpen_Tmoy_ASO,
    mod_senesduree_Bpen_Tnight21j_jsenes10,mod_senesduree_Bpen_Tnight30j_jsenes10,mod_senesduree_Bpen_Tnight40j_jsenes10,
    mod_senesduree_Bpen_GDDinv25_jsenes10,mod_senesduree_Bpen_GDDinv20_jsenes10,mod_senesduree_Bpen_GDDinv15_jsenes10,
    mod_senesduree_Bpen_Tnight21j_jsenes50,mod_senesduree_Bpen_Tnight30j_jsenes50,mod_senesduree_Bpen_Tnight40j_jsenes50,
    mod_senesduree_Bpen_GDDinv25_jsenes50,mod_senesduree_Bpen_GDDinv20_jsenes50,mod_senesduree_Bpen_GDDinv15_jsenes50)
# les précipitations d'octobre donnent les meilleurs résultats sur le jeu de données tronqué, mais il y a très peu de différences entre les AIC des modèles simples...
# Ce sont aussi les meilleures variables quand on regarde le jeu de données complet
# -> on regarde en version multivariables, avec ou sans debou10 (i.e. jeu de données complet ou non)


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(glmer(senesduree ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
               Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
               GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 + 
               Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
               GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
               (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senesduree_Bpen_multivar = glmer(senesduree ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson) # AIC = 13422.9
# # visreg(mod_senesduree_Bpen_multivar, "GDDinv25_jsenes50")

# Plutôt que de faire "a la mano", on peut utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# Une première sélection basée sur le résumé du modèle global (summary) permet de limiter la lourdeur/durée de la sélection dredge
# La date de débourrement n'est pas retenue comme variable importante donc on considère le jeu de données complet
senesduree_Bpen = senesduree_Bpen_all
mod_senesduree_Bpen_multivar = glmer(senesduree ~ P_autumn + P_juillet + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                    Tnight30j_jsenes10 +
                                    GDDinv25_jsenes10 + 
                                    Tnight30j_jsenes50 +
                                    GDDinv25_jsenes50 + 
                                    (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson)
dredge(mod_senesduree_Bpen_multivar)
mod_senesduree_Bpen_multivar = glmer(senesduree ~ P_autumn + Tmoy_ASO + Tnight30j_jsenes10 +
                                   (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson) # AICc = 153331


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
summary(glmer(senesduree ~ P_autumn * Tmoy_ASO + Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson))
summary(glmer(senesduree ~ P_autumn + Tmoy_ASO * Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson))
summary(glmer(senesduree ~ P_autumn * Tnight30j_jsenes10 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson))
# les interactions qui améliorent l'AIC font que le modèle ne converge pas

bestmod_senesduree_Bpen = glmer(senesduree ~ P_autumn + Tmoy_ASO + Tnight30j_jsenes10 +
                                  (1|ID_zone) + (1|yearQ), senesduree_Bpen, family=poisson) # AICc = 15333.5
bestmods = c(list("Bouleau_verruqueux"=bestmod_senesduree_Bpen), bestmods)


summary(bestmod_senesduree_Bpen)
# la sénescence est dure plus longtemps si l'automne est plus sec (i.e. les précipitations d'automne ont un effet accélérateur)
# les températures de la fin de saison estivale / début d'automne ont un effet accélérateur, SAUF s'il ne fait pas assez froid la nuit

# qqnorm(resid(bestmod_senesduree_Bpen))
# qqline(resid(bestmod_senesduree_Bpen))

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
mod_train <- glmer(formula(bestmod_senesduree_Bpen), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_verruqueux", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Bouleau_verruqueux",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senesduree_Bpen))),
                                        coef = coef(summary(bestmod_senesduree_Bpen))[,1],
                                        std = coef(summary(bestmod_senesduree_Bpen))[,2],
                                        pval = coef(summary(bestmod_senesduree_Bpen))[,4],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Bpen)[rownames(coef(summary(bestmod_senesduree_Bpen)))], error = function(e) return(NA))))



##############################################-
#*---- Meleze ----

senesduree_Ldec = senesduree_Alps[senesduree_Alps$species == "Meleze",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Ldec = senesduree_Ldec[!is.na(senesduree_Ldec$Tmoy_GS),]
ggplot(senesduree_Ldec, aes(x=senesduree)) + geom_histogram(binwidth=5)



#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Ldec_Altyear <- glmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Ldec, family=poisson)
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
mod_train <- glmer(formula(mod_senesduree_Ldec_Altyear), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Meleze", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Ldec_all = senesduree_Ldec
senesduree_Ldec = senesduree_Ldec[!is.na(senesduree_Ldec$debou10),] #senesduree_Ldec_all

mod_senesduree_Ldec_debou10 <- glmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson) 
mod_senesduree_Ldec_P_summer <- glmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_juin <- glmer(senesduree ~ P_juin + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_juillet <- glmer(senesduree ~ P_juillet + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_aout <- glmer(senesduree ~ P_aout + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_septembre <- glmer(senesduree ~ P_septembre+ (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_octobre <- glmer(senesduree ~ P_octobre + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_P_autumn <- glmer(senesduree ~ P_autumn + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tmoy_GS <- glmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tmoy_MJ <- glmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tmoy_AS <- glmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tmoy_ASO <- glmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight21j_jsenes10 <- glmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight30j_jsenes10 <- glmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight40j_jsenes10 <- glmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv25_jsenes10 <- glmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv20_jsenes10 <- glmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv15_jsenes10 <- glmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight21j_jsenes50 <- glmer(senesduree ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight30j_jsenes50 <- glmer(senesduree ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_Tnight40j_jsenes50 <- glmer(senesduree ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv25_jsenes50 <- glmer(senesduree ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv20_jsenes50 <- glmer(senesduree ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
mod_senesduree_Ldec_GDDinv15_jsenes50 <- glmer(senesduree ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)

AIC(mod_senesduree_Ldec_debou10, 
    mod_senesduree_Ldec_P_summer,mod_senesduree_Ldec_P_juin,mod_senesduree_Ldec_P_juillet,mod_senesduree_Ldec_P_aout,mod_senesduree_Ldec_P_septembre,mod_senesduree_Ldec_P_octobre,mod_senesduree_Ldec_P_autumn,
    mod_senesduree_Ldec_Tmoy_GS,mod_senesduree_Ldec_Tmoy_MJ,mod_senesduree_Ldec_Tmoy_AS,mod_senesduree_Ldec_Tmoy_ASO,
    mod_senesduree_Ldec_Tnight21j_jsenes10,mod_senesduree_Ldec_Tnight30j_jsenes10,mod_senesduree_Ldec_Tnight40j_jsenes10,
    mod_senesduree_Ldec_GDDinv25_jsenes10,mod_senesduree_Ldec_GDDinv20_jsenes10,mod_senesduree_Ldec_GDDinv15_jsenes10,
    mod_senesduree_Ldec_Tnight21j_jsenes50,mod_senesduree_Ldec_Tnight30j_jsenes50,mod_senesduree_Ldec_Tnight40j_jsenes50,
    mod_senesduree_Ldec_GDDinv25_jsenes50,mod_senesduree_Ldec_GDDinv20_jsenes50,mod_senesduree_Ldec_GDDinv15_jsenes50)
# la date de débourrement donne qqch d'intéressant à cette étape !
# => fonctionner avec le jeu de données tronqué



# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(glmer(senesduree ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 + 
                Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senesduree_Ldec_multivar = glmer(senesduree ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson) # AIC = 13422.9
# # visreg(mod_senesduree_Ldec_multivar, "GDDinv25_jsenes50")

# Plutôt que de faire "a la mano", on peut utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# Une première sélection basée sur le résumé du modèle global (summary) permet de limiter la lourdeur/durée de la sélection dredge
# /!\ trop de variables incluses et trop de corrélations dans leurs effets
mod_senesduree_Ldec_multivar = glmer(senesduree ~ debou10 + P_summer + P_autumn + Tmoy_GS + Tmoy_MJ + 
                                       Tnight30j_jsenes10 +
                                       GDDinv25_jsenes10 + 
                                       Tnight30j_jsenes50 +
                                       GDDinv25_jsenes50 + 
                                       (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson)
dredge(mod_senesduree_Ldec_multivar)
mod_senesduree_Ldec_multivar = glmer(senesduree ~ debou10 + P_juin + P_autumn + Tmoy_GS + Tnight30j_jsenes50 +
                                       (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson) # AICc = 13008.8


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
AICc(glmer(senesduree ~ debou10 * P_juin + P_autumn + Tmoy_GS + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 * P_autumn + P_juin + Tmoy_GS + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 * Tmoy_GS + P_juin + P_autumn + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 * Tnight30j_jsenes50 + P_juin + P_autumn + Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin * P_autumn + Tmoy_GS + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin * Tmoy_GS + P_autumn + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin * Tnight30j_jsenes50 + P_autumn + Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin + P_autumn * Tmoy_GS + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin + P_autumn * Tnight30j_jsenes50 + Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
AICc(glmer(senesduree ~ debou10 + P_juin + P_autumn + Tmoy_GS * Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Ldec, family=poisson))
# les interactions qui améliorent l'AIC font que le modèle ne converge pas

bestmod_senesduree_Ldec = glmer(senesduree ~ debou10 + P_juin + P_autumn + Tmoy_GS + Tnight30j_jsenes50 +
                                  (1|ID_zone) + (1|yearQ), senesduree_Ldec %>% mutate(debou10=debou10/100), family=poisson) # AICc = 13008.8 #on rescale la variable de débourrement pour pouvoir calculer les indices
bestmods = c(list("Meleze"=bestmod_senesduree_Ldec), bestmods)


summary(bestmod_senesduree_Ldec)
# la sénescence est accélérée lorque le débourrement a été plus précoce
# les précipitations de début de saison plus importantes ralentissent la sénescence, mais celles d'automne l'accélèrent
# les températures plus chaudes pendant la saison de croissance accélèrent la sénescence
# des températures nocturnes juste avant et pendant la sénescence plus froides accélèrent la sénescence

# qqnorm(resid(bestmod_senesduree_Ldec))
# qqline(resid(bestmod_senesduree_Ldec))

# Effets fixes VS tous les effets inclus
R2_models[R2_models$species == "Meleze", c("R2_bestmod_fixef","R2_bestmod_allef")] = r.squaredGLMM(bestmod_senesduree_Ldec)
# R2 des effets fixes très très nul !!
# - Part de variance expliquée par les différentes variables
calcVarPart(bestmod_senesduree_Ldec)
# - Validation croisée
n = nrow(senesduree_Ldec)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Ldec[trainIndex ,]
test <- senesduree_Ldec[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Ldec$ID_zone[drop=T])[!unique(senesduree_Ldec$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- glmer(formula(bestmod_senesduree_Ldec), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Meleze", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Meleze",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senesduree_Ldec))),
                                        coef = coef(summary(bestmod_senesduree_Ldec))[,1],
                                        std = coef(summary(bestmod_senesduree_Ldec))[,2],
                                        pval = coef(summary(bestmod_senesduree_Ldec))[,4],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Ldec)[rownames(coef(summary(bestmod_senesduree_Ldec)))], error = function(e) return(NA))))



##############################################-
#*---- Bouleau pubescent ----

senesduree_Bpub = senesduree_Alps[senesduree_Alps$species == "Bouleau_pubescent",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Bpub = senesduree_Bpub[!is.na(senesduree_Bpub$Tmoy_GS),]
ggplot(senesduree_Bpub, aes(x=senesduree)) + geom_histogram(binwidth=5)



#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Bpub_Altyear <- glmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Bpub, family=poisson)
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
mod_train <- glmer(formula(mod_senesduree_Bpub_Altyear), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Bouleau_pubescent", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Bpub_all = senesduree_Bpub
senesduree_Bpub = senesduree_Bpub[!is.na(senesduree_Bpub$debou10),] #senesduree_Bpub_all

mod_senesduree_Bpub_debou10 <- glmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson) 
mod_senesduree_Bpub_P_summer <- glmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_juin <- glmer(senesduree ~ P_juin + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_juillet <- glmer(senesduree ~ P_juillet + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_aout <- glmer(senesduree ~ P_aout + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_septembre <- glmer(senesduree ~ P_septembre+ (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_octobre <- glmer(senesduree ~ P_octobre + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_P_autumn <- glmer(senesduree ~ P_autumn + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tmoy_GS <- glmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tmoy_MJ <- glmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tmoy_AS <- glmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tmoy_ASO <- glmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight21j_jsenes10 <- glmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight30j_jsenes10 <- glmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight40j_jsenes10 <- glmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv25_jsenes10 <- glmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv20_jsenes10 <- glmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv15_jsenes10 <- glmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight21j_jsenes50 <- glmer(senesduree ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight30j_jsenes50 <- glmer(senesduree ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_Tnight40j_jsenes50 <- glmer(senesduree ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv25_jsenes50 <- glmer(senesduree ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv20_jsenes50 <- glmer(senesduree ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
mod_senesduree_Bpub_GDDinv15_jsenes50 <- glmer(senesduree ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)

AIC(mod_senesduree_Bpub_debou10, 
    mod_senesduree_Bpub_P_summer,mod_senesduree_Bpub_P_juin,mod_senesduree_Bpub_P_juillet,mod_senesduree_Bpub_P_aout,mod_senesduree_Bpub_P_septembre,mod_senesduree_Bpub_P_octobre,mod_senesduree_Bpub_P_autumn,
    mod_senesduree_Bpub_Tmoy_GS,mod_senesduree_Bpub_Tmoy_MJ,mod_senesduree_Bpub_Tmoy_AS,mod_senesduree_Bpub_Tmoy_ASO,
    mod_senesduree_Bpub_Tnight21j_jsenes10,mod_senesduree_Bpub_Tnight30j_jsenes10,mod_senesduree_Bpub_Tnight40j_jsenes10,
    mod_senesduree_Bpub_GDDinv25_jsenes10,mod_senesduree_Bpub_GDDinv20_jsenes10,mod_senesduree_Bpub_GDDinv15_jsenes10,
    mod_senesduree_Bpub_Tnight21j_jsenes50,mod_senesduree_Bpub_Tnight30j_jsenes50,mod_senesduree_Bpub_Tnight40j_jsenes50,
    mod_senesduree_Bpub_GDDinv25_jsenes50,mod_senesduree_Bpub_GDDinv20_jsenes50,mod_senesduree_Bpub_GDDinv15_jsenes50)
# la date de débourrement ne ressort pas vraiment à cette étape, on peut tester avec le jeu de données complet
# En fonction du jeu de données utilisé, soit la précipiations de septembre et les températures moyennes d'août-septembre donnent les meilleurs AIC


# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(glmer(senesduree ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 + 
                Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senesduree_Bpub_multivar = glmer(senesduree ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson) # AIC = 13422.9
# # visreg(mod_senesduree_Bpub_multivar, "GDDinv25_jsenes50")

# Plutôt que de faire "a la mano", on peut utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# Une première sélection basée sur le résumé du modèle global (summary) permet de limiter la lourdeur/durée de la sélection dredge
# /!\ trop de variables incluses et trop de corrélations dans leurs effets
# Le débourrement ne ressort pas donc on compare les modèles sur le jeu de données complet
senesduree_Bpub = senesduree_Bpub_all
mod_senesduree_Bpub_multivar = glmer(senesduree ~ P_summer + P_autumn + P_juillet + P_septembre + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                                       (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson)
dredge(mod_senesduree_Bpub_multivar)
mod_senesduree_Bpub_multivar = glmer(senesduree ~ Tmoy_AS +
                                       (1|ID_zone) + (1|yearQ), senesduree_Bpub, family=poisson) # AICc = 13008.8


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions (pas de sens avec une seule variable explicative)

bestmod_senesduree_Bpub = glmer(senesduree ~ Tmoy_AS +
                                  (1|ID_zone) + (1|yearQ), senesduree_Bpub %>% mutate(debou10=debou10/100), family=poisson) # AICc = 13008.8 #on rescale la variable de débourrement pour pouvoir calculer les indices
bestmods = c(list("Bouleau_pubescent"=bestmod_senesduree_Bpub), bestmods)


summary(bestmod_senesduree_Bpub)
# la sénescence est accélérée quand il fait plus froid en août-septembre

# qqnorm(resid(bestmod_senesduree_Bpub))
# qqline(resid(bestmod_senesduree_Bpub))

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
mod_train <- glmer(formula(bestmod_senesduree_Bpub), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Bouleau_pubescent", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Bouleau_pubescent",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senesduree_Bpub))),
                                        coef = coef(summary(bestmod_senesduree_Bpub))[,1],
                                        std = coef(summary(bestmod_senesduree_Bpub))[,2],
                                        pval = coef(summary(bestmod_senesduree_Bpub))[,4],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Bpub)[rownames(coef(summary(bestmod_senesduree_Bpub)))], error = function(e) return(NA))))





##############################################-
#*---- Hetre ----

senesduree_Fsyl = senesduree_Alps[senesduree_Alps$species == "Hetre",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Fsyl = senesduree_Fsyl[!is.na(senesduree_Fsyl$Tmoy_GS),]
ggplot(senesduree_Fsyl, aes(x=senesduree)) + geom_histogram(binwidth=5)



#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Fsyl_Altyear <- glmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Fsyl, family=poisson)
summary(mod_senesduree_Fsyl_Altyear)
altyearmods = c(list("Hetre"=mod_senesduree_Fsyl_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Hetre", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Fsyl_Altyear)
# + Validation croisée
n = nrow(senesduree_Fsyl)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Fsyl[trainIndex ,]
test <- senesduree_Fsyl[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Fsyl$ID_zone[drop=T])[!unique(senesduree_Fsyl$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- glmer(formula(mod_senesduree_Fsyl_Altyear), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Hetre", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Fsyl_all = senesduree_Fsyl
senesduree_Fsyl = senesduree_Fsyl[!is.na(senesduree_Fsyl$debou10),] #senesduree_Fsyl_all

mod_senesduree_Fsyl_debou10 <- glmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson) 
mod_senesduree_Fsyl_P_summer <- glmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_juin <- glmer(senesduree ~ P_juin + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_juillet <- glmer(senesduree ~ P_juillet + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_aout <- glmer(senesduree ~ P_aout + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_septembre <- glmer(senesduree ~ P_septembre+ (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_octobre <- glmer(senesduree ~ P_octobre + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_P_autumn <- glmer(senesduree ~ P_autumn + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tmoy_GS <- glmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tmoy_MJ <- glmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tmoy_AS <- glmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tmoy_ASO <- glmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight21j_jsenes10 <- glmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight30j_jsenes10 <- glmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight40j_jsenes10 <- glmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv25_jsenes10 <- glmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv20_jsenes10 <- glmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv15_jsenes10 <- glmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight21j_jsenes50 <- glmer(senesduree ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight30j_jsenes50 <- glmer(senesduree ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_Tnight40j_jsenes50 <- glmer(senesduree ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv25_jsenes50 <- glmer(senesduree ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv20_jsenes50 <- glmer(senesduree ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
mod_senesduree_Fsyl_GDDinv15_jsenes50 <- glmer(senesduree ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)

AIC(mod_senesduree_Fsyl_debou10, 
    mod_senesduree_Fsyl_P_summer,mod_senesduree_Fsyl_P_juin,mod_senesduree_Fsyl_P_juillet,mod_senesduree_Fsyl_P_aout,mod_senesduree_Fsyl_P_septembre,mod_senesduree_Fsyl_P_octobre,mod_senesduree_Fsyl_P_autumn,
    mod_senesduree_Fsyl_Tmoy_GS,mod_senesduree_Fsyl_Tmoy_MJ,mod_senesduree_Fsyl_Tmoy_AS,mod_senesduree_Fsyl_Tmoy_ASO,
    mod_senesduree_Fsyl_Tnight21j_jsenes10,mod_senesduree_Fsyl_Tnight30j_jsenes10,mod_senesduree_Fsyl_Tnight40j_jsenes10,
    mod_senesduree_Fsyl_GDDinv25_jsenes10,mod_senesduree_Fsyl_GDDinv20_jsenes10,mod_senesduree_Fsyl_GDDinv15_jsenes10,
    mod_senesduree_Fsyl_Tnight21j_jsenes50,mod_senesduree_Fsyl_Tnight30j_jsenes50,mod_senesduree_Fsyl_Tnight40j_jsenes50,
    mod_senesduree_Fsyl_GDDinv25_jsenes50,mod_senesduree_Fsyl_GDDinv20_jsenes50,mod_senesduree_Fsyl_GDDinv15_jsenes50)
# les températures moyennes d'août-septembre donnent les meilleurs AIC



# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(glmer(senesduree ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 + 
                Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senesduree_Fsyl_multivar = glmer(senesduree ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson) # AIC = 13422.9
# # visreg(mod_senesduree_Fsyl_multivar, "GDDinv25_jsenes50")

# Plutôt que de faire "a la mano", on peut utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# Une première sélection basée sur le résumé du modèle global (summary) permet de limiter la lourdeur/durée de la sélection dredge
# /!\ trop de variables incluses et trop de corrélations dans leurs effets
# Le débourrement ne ressort pas donc on compare les modèles sur le jeu de données complet
senesduree_Fsyl = senesduree_Fsyl_all
mod_senesduree_Fsyl_multivar = glmer(senesduree ~ P_summer + P_autumn + P_juin + P_juillet + Tmoy_GS + Tmoy_AS + 
                                       Tnight30j_jsenes10 +
                                       GDDinv25_jsenes10 + 
                                       Tnight30j_jsenes50 +
                                       GDDinv25_jsenes50 + 
                                       (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson)
dredge(mod_senesduree_Fsyl_multivar)
mod_senesduree_Fsyl_multivar = glmer(senesduree ~ GDDinv25_jsenes50 + Tmoy_AS + Tnight30j_jsenes50 +
                                       (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson) # AICc = 426.8


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
AICc(glmer(senesduree ~ GDDinv25_jsenes50 * Tmoy_AS + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson))
AICc(glmer(senesduree ~ GDDinv25_jsenes50 + Tmoy_AS * Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson))
AICc(glmer(senesduree ~ GDDinv25_jsenes50 * Tnight30j_jsenes50 + Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson))
# les interactions n'améliorent pas l'AICc

bestmod_senesduree_Fsyl = glmer(senesduree ~ GDDinv25_jsenes50 + Tmoy_AS + Tnight30j_jsenes50 +
                                  (1|ID_zone) + (1|yearQ), senesduree_Fsyl, family=poisson) # AICc = 426.8
bestmods = c(list("Hetre"=bestmod_senesduree_Fsyl), bestmods)


summary(bestmod_senesduree_Fsyl)
# la sénescence est accélérée lorque le débourrement a été plus précoce
# les précipitations de début de saison plus importantes ralentissent la sénescence, mais celles d'automne l'accélèrent
# les températures plus chaudes pendant la saison de croissance accélèrent la sénescence
# des températures nocturnes juste avant et pendant la sénescence plus froides accélèrent la sénescence

# qqnorm(resid(bestmod_senesduree_Fsyl))
# qqline(resid(bestmod_senesduree_Fsyl))

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
mod_train <- glmer(formula(bestmod_senesduree_Fsyl), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Hetre", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Hetre",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senesduree_Fsyl))),
                                        coef = coef(summary(bestmod_senesduree_Fsyl))[,1],
                                        std = coef(summary(bestmod_senesduree_Fsyl))[,2],
                                        pval = coef(summary(bestmod_senesduree_Fsyl))[,4],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Fsyl)[rownames(coef(summary(bestmod_senesduree_Fsyl)))], error = function(e) return(NA))))




##############################################-
#*---- Sorbier ----

senesduree_Sacu = senesduree_Alps[senesduree_Alps$species == "Sorbier",]
# # on retire les lignes où on n'a pas de données de température (2005)
senesduree_Sacu = senesduree_Sacu[!is.na(senesduree_Sacu$Tmoy_GS),]
ggplot(senesduree_Sacu, aes(x=senesduree)) + geom_histogram(binwidth=5)



#*-------- 1) On regarde un modèle simple comme dans Bison et al. 2019, avec seulement altitude et année ---- 
mod_senesduree_Sacu_Altyear <- glmer(senesduree ~ altitude + year + (year|ID_zone), senesduree_Sacu, family=poisson)
summary(mod_senesduree_Sacu_Altyear)
altyearmods = c(list("Sorbier"=mod_senesduree_Sacu_Altyear), altyearmods)

# On note aussi le R2 de ce modèle, pour comparer avec les résultats des modèles avec température
R2_models[R2_models$species == "Sorbier", c("R2_altyear_fixef","R2_altyear_allef")] = r.squaredGLMM(mod_senesduree_Sacu_Altyear)
# + Validation croisée
n = nrow(senesduree_Sacu)
trainIndex <- sample(1:n, size = round(0.8*n), replace=FALSE)
train <- senesduree_Sacu[trainIndex ,]
test <- senesduree_Sacu[-trainIndex ,]
test1 <- test[!test$ID_zone%in%(unique(senesduree_Sacu$ID_zone[drop=T])[!unique(senesduree_Sacu$ID_zone[drop=T])%in%unique(train$ID_zone[drop=T])]),]
mod_train <- glmer(formula(mod_senesduree_Sacu_Altyear), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
R2_models[R2_models$species == "Sorbier", "R2_altyear_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]


#*-------- 2) On teste l'effet de différentes variables climatiques à la place des variables année et altitude ---- 
# On garde tout de même l'effet année en effet aléatoire (cf autres param climatiques que la température)

# # /!\ on ne peut pas comparer le modèle avec débourrement avec les autres puisqu'il n'y a pas le même nombre de données en entrée (cf NA dans debou10)
# #     => on compare tous les modèles avec la DB sans NA dans debou10, si c'est moins bon on prend la DB globale sans considérer la variable debou10
senesduree_Sacu_all = senesduree_Sacu
senesduree_Sacu = senesduree_Sacu[!is.na(senesduree_Sacu$debou10),] #senesduree_Sacu_all

mod_senesduree_Sacu_debou10 <- glmer(senesduree ~ debou10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson) 
mod_senesduree_Sacu_P_summer <- glmer(senesduree ~ P_summer + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_juin <- glmer(senesduree ~ P_juin + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_juillet <- glmer(senesduree ~ P_juillet + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_aout <- glmer(senesduree ~ P_aout + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_septembre <- glmer(senesduree ~ P_septembre+ (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_octobre <- glmer(senesduree ~ P_octobre + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_P_autumn <- glmer(senesduree ~ P_autumn + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tmoy_GS <- glmer(senesduree ~ Tmoy_GS + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tmoy_MJ <- glmer(senesduree ~ Tmoy_MJ + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tmoy_AS <- glmer(senesduree ~ Tmoy_AS + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tmoy_ASO <- glmer(senesduree ~ Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight21j_jsenes10 <- glmer(senesduree ~ Tnight21j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight30j_jsenes10 <- glmer(senesduree ~ Tnight30j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight40j_jsenes10 <- glmer(senesduree ~ Tnight40j_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv25_jsenes10 <- glmer(senesduree ~ GDDinv25_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv20_jsenes10 <- glmer(senesduree ~ GDDinv20_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv15_jsenes10 <- glmer(senesduree ~ GDDinv15_jsenes10 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight21j_jsenes50 <- glmer(senesduree ~ Tnight21j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight30j_jsenes50 <- glmer(senesduree ~ Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_Tnight40j_jsenes50 <- glmer(senesduree ~ Tnight40j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv25_jsenes50 <- glmer(senesduree ~ GDDinv25_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv20_jsenes50 <- glmer(senesduree ~ GDDinv20_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
mod_senesduree_Sacu_GDDinv15_jsenes50 <- glmer(senesduree ~ GDDinv15_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)

AIC(mod_senesduree_Sacu_debou10, 
    mod_senesduree_Sacu_P_summer,mod_senesduree_Sacu_P_juin,mod_senesduree_Sacu_P_juillet,mod_senesduree_Sacu_P_aout,mod_senesduree_Sacu_P_septembre,mod_senesduree_Sacu_P_octobre,mod_senesduree_Sacu_P_autumn,
    mod_senesduree_Sacu_Tmoy_GS,mod_senesduree_Sacu_Tmoy_MJ,mod_senesduree_Sacu_Tmoy_AS,mod_senesduree_Sacu_Tmoy_ASO,
    mod_senesduree_Sacu_Tnight21j_jsenes10,mod_senesduree_Sacu_Tnight30j_jsenes10,mod_senesduree_Sacu_Tnight40j_jsenes10,
    mod_senesduree_Sacu_GDDinv25_jsenes10,mod_senesduree_Sacu_GDDinv20_jsenes10,mod_senesduree_Sacu_GDDinv15_jsenes10,
    mod_senesduree_Sacu_Tnight21j_jsenes50,mod_senesduree_Sacu_Tnight30j_jsenes50,mod_senesduree_Sacu_Tnight40j_jsenes50,
    mod_senesduree_Sacu_GDDinv25_jsenes50,mod_senesduree_Sacu_GDDinv20_jsenes50,mod_senesduree_Sacu_GDDinv15_jsenes50)
# la date de débourrement donne qqch d'intéressant à cette étape !
# => fonctionner avec le jeu de données tronqué



# On peut aussi tester d'intégrer plusieurs variables différentes, et faire un choix backward/forward
summary(glmer(senesduree ~ debou10 + P_summer + P_autumn + P_juin + P_juillet + P_aout + P_septembre + P_octobre + Tmoy_GS + Tmoy_MJ + Tmoy_AS + Tmoy_ASO + 
                Tnight21j_jsenes10 + Tnight30j_jsenes10 + Tnight40j_jsenes10 +
                GDDinv25_jsenes10 + GDDinv20_jsenes10 + GDDinv15_jsenes10 + 
                Tnight21j_jsenes50 + Tnight30j_jsenes50 + Tnight40j_jsenes50 +
                GDDinv25_jsenes50 + GDDinv20_jsenes50 + GDDinv15_jsenes50 +
                (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson), corr=F)
# # On enlève les variables une par une (en fonction des pvalues), et on regarde également l'AIC et la corrélation entre variables
# # En combinant ces 3 critères (AIC, corrélations, significativité des effets), on retient donc :
# mod_senesduree_Sacu_multivar = glmer(senesduree ~ P_summer + GDDinv25_jsenes50 + 
#                                (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson) # AIC = 13422.9
# # visreg(mod_senesduree_Sacu_multivar, "GDDinv25_jsenes50")

# Plutôt que de faire "a la mano", on peut utiliser la fonction "dredge" de MuMin pour faire la sélection de variable en testant toutes les combinaisons possibles
# Une première sélection basée sur le résumé du modèle global (summary) permet de limiter la lourdeur/durée de la sélection dredge
# /!\ trop de variables incluses et trop de corrélations dans leurs effets
mod_senesduree_Sacu_multivar = glmer(senesduree ~ debou10 + P_autumn + P_septembre + Tmoy_GS + Tmoy_MJ + Tmoy_ASO +
                                       GDDinv25_jsenes10 + 
                                       Tnight30j_jsenes50 +
                                       GDDinv25_jsenes50 + 
                                       (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson)
dredge(mod_senesduree_Sacu_multivar)
mod_senesduree_Sacu_multivar = glmer(senesduree ~ debou10 + GDDinv25_jsenes10 + Tmoy_ASO + Tnight30j_jsenes50 +
                                       (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson) # AIC = 4302.7


#*-------- 3) On complexifie le modèle ---- 
# On peut aussi complexifier en ajoutant des interactions :
AICc(glmer(senesduree ~ debou10 * GDDinv25_jsenes10 + Tmoy_ASO + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
AICc(glmer(senesduree ~ debou10 * Tmoy_ASO + GDDinv25_jsenes10 + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
AICc(glmer(senesduree ~ debou10 * Tnight30j_jsenes50 + GDDinv25_jsenes10 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
AICc(glmer(senesduree ~ debou10 + GDDinv25_jsenes10 * Tmoy_ASO + Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
AICc(glmer(senesduree ~ debou10 + GDDinv25_jsenes10 * Tnight30j_jsenes50 + Tmoy_ASO + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
AICc(glmer(senesduree ~ debou10 + GDDinv25_jsenes10 + Tmoy_ASO * Tnight30j_jsenes50 + (1|ID_zone) + (1|yearQ), senesduree_Sacu, family=poisson))
# les interactions qui améliorent l'AIC font que le modèle ne converge pas

bestmod_senesduree_Sacu = glmer(senesduree ~ debou10 + GDDinv25_jsenes10 + Tmoy_ASO + Tnight30j_jsenes50 +
                                  (1|ID_zone) + (1|yearQ), senesduree_Sacu %>% mutate(debou10=debou10/100), family=poisson) # AICc = 13008.8 #on rescale la variable de débourrement pour pouvoir calculer les indices
bestmods = c(list("Sorbier"=bestmod_senesduree_Sacu), bestmods)


summary(bestmod_senesduree_Sacu)
# la sénescence est accélérée lorque le débourrement a été plus précoce
# les précipitations de début de saison plus importantes ralentissent la sénescence, mais celles d'automne l'accélèrent
# les températures plus chaudes pendant la saison de croissance accélèrent la sénescence
# des températures nocturnes juste avant et pendant la sénescence plus froides accélèrent la sénescence

# qqnorm(resid(bestmod_senesduree_Sacu))
# qqline(resid(bestmod_senesduree_Sacu))

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
mod_train <- glmer(formula(bestmod_senesduree_Sacu), train, family=poisson)
predictions <- mod_train %>% predict(test1)
test2 <- cbind(test1, predictions)
# summary(lm(predictions~julian_day, test2))
# # R2 ajusté = 0.7489
R2_models[R2_models$species == "Sorbier", "R2_bestmod_calibval"] = summary(lm(predictions~senesduree, test2))[["adj.r.squared"]]



resultats = rbind(resultats, data.frame(species = "Sorbier",
                                        periode = "2006-2025",
                                        variable = rownames(coef(summary(bestmod_senesduree_Sacu))),
                                        coef = coef(summary(bestmod_senesduree_Sacu))[,1],
                                        std = coef(summary(bestmod_senesduree_Sacu))[,2],
                                        pval = coef(summary(bestmod_senesduree_Sacu))[,4],
                                        varexpli = tryCatch(calcVarPart(bestmod_senesduree_Sacu)[rownames(coef(summary(bestmod_senesduree_Sacu)))], error = function(e) return(NA))))






##############################################-
#*---- Enregistrement des résultats ----


write.csv(resultats[-1,], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree_20062025__coeff.csv", row.names = F)
write.csv(R2_models[!is.na(R2_models$R2_bestmod_fixef),], "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree_20062025__R2.csv", row.names = F)
save(bestmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senesduree_20062025.Rdata")
save(altyearmods, file="/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/mods_senesduree_20062025_v1altyear.Rdata")









############################################################################################-
# INDICE PHÉNOCLIM BASÉ SUR CES MODÈLES                                                  ----
############################################################################################-

##############################################-
#*----- Par espace protégé ----

# L'idée serait de calculer pour chaque année x espèce x espprot x altitude (ou classe d'altitude) la date de début de sénescence et sa durée, en fonction des 
# conditions climatiques moyennes mesurées une année donnée dans la classe d'altitude x espprot.


phenoclim = read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/Phenoclim/data/_CLEANED_data_pheno.csv")

# Espaces protégés considérés : PN Vanoise, PN Écrins, RNs + Parcs italiens du Gran Paradisio (ID = 555580376) et Mont Avic (ID = 555528148) + Espace Mt Blanc
# => data France : https://data.naturefrance.fr/geonetwork/srv/fre/catalog.search#/metadata/cc1bfe07-625e-4df6-86a2-49fd07d193d7
# => data Italie : UNEP-WCMC and IUCN (2025), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], May 2025, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net. 
# => Espace Mont Blanc : https://www.espace-mont-blanc.com/asset/img/carte_emb_dl.jpg

# *---- 1) On récupère le tableau d'espaces protégés déjà créé pour l'indice de débourrement
tab_recap_Phenoclim_x_espprot5km = read.csv("data/_lim_admin/espaces_proteges/tab_recap_Phenoclim_x_espprot5km.csv", row.names = 1)
colnames(tab_recap_Phenoclim_x_espprot5km) = tab_recap_Phenoclim_x_espprot5km["colnames",]

phenoclim = merge(phenoclim, tab_recap_Phenoclim_x_espprot5km %>% dplyr::select(-c(coord_x_2154,coord_y_2154)), by="id_base_site", all.x=T, all.y=F)


# *---- 2) On crée un tableau avec les informations "espace protégé" et toutes les données météo (T + Précip)
# TEMPERATURES
climperiodes = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/meteo_reconstruc/varTautomne_allphenosites.csv", row.names = 1)
climperiodes = merge(climperiodes, 
                  phenoclim[!duplicated(phenoclim$id_base_site),c(colnames(tab_recap_Phenoclim_x_espprot5km),"altitude","cl_alt","cl_alt2","cl_alt3")],
                  by.x="sitePheno", by.y="id_base_site", all.x=T)

# PRECIPITATIONS
PRECIP<-brick("/Users/ninonfontaine/Desktop/projetsR/TEST/data/_meteo/ERA5land_precip.grib") 
# /!\ téléchargement au 27 novembre --> on n'a pas novembre et décembre 2025 ! Donc on retire ces deux derniers mois
names(PRECIP) = paste("P",rep(2004:2025, each=12), rep(1:12,length(2004:2025)),sep="_")[1:length(names(PRECIP))]

climperiodes[,c("P_juin","P_juillet","P_aout", "P_septembre","P_octobre")] = 
  t(apply(cbind(climperiodes$annee, 
                extract(PRECIP[[grep("_6", names(PRECIP))]], crds(project(vect(climperiodes,geom=c("coord_x_2154", "coord_y_2154"), "epsg:2154"), crs(PRECIP)))),
                extract(PRECIP[[grep("_7", names(PRECIP))]], crds(project(vect(climperiodes,geom=c("coord_x_2154", "coord_y_2154"), "epsg:2154"), crs(PRECIP)))),
                extract(PRECIP[[grep("_8", names(PRECIP))]], crds(project(vect(climperiodes,geom=c("coord_x_2154", "coord_y_2154"), "epsg:2154"), crs(PRECIP)))),
                extract(PRECIP[[grep("_9", names(PRECIP))]], crds(project(vect(climperiodes,geom=c("coord_x_2154", "coord_y_2154"), "epsg:2154"), crs(PRECIP)))),
                extract(PRECIP[[grep("_10", names(PRECIP))]], crds(project(vect(climperiodes,geom=c("coord_x_2154", "coord_y_2154"), "epsg:2154"), crs(PRECIP))))),
          1, function(x){x[x[1]-c(2002,1980,1958,1936,1914)]}))

climperiodes$P_summer = climperiodes$P_juin + climperiodes$P_juillet + climperiodes$P_aout
climperiodes$P_autumn = climperiodes$P_septembre + climperiodes$P_octobre



# # Aperçu des zones (espace protégé x classe d'altitude) où on a des données
# table(climperiodes$espprot, climperiodes$cl_alt)
# table(phenoclim$espprot, phenoclim$year)
# table(climperiodes$espprot, climperiodes$annee)

clim_pourpred = climperiodes %>% filter(!is.na(EMB)) %>% group_by(EMB, annee, esp, cl_alt, cl_alt2, cl_alt3) %>% 
  summarise(across(c(starts_with("P_"),ends_with("ann")), mean)) %>% rename(espprot="EMB")
for(espaceprot in colnames(tab_recap_Phenoclim_x_espprot5km)[-c(1:5)]){
  clim_pourpred = rbind(clim_pourpred,
                     climperiodes[,colnames(climperiodes) != "espprot"] %>% rename(espprot=espaceprot) %>% filter(!is.na(espprot)) %>% group_by(espprot, annee, esp, cl_alt, cl_alt2, cl_alt3) %>% 
                       summarise(across(c(starts_with("P_"), ends_with("ann")), mean)) )
  
}



##########################################################################################################################################-
# /!\ il n'y a pas des données pour tous les espaces et toutes les années --> ajuster les prédictions ?
##########################################################################################################################################-
# pres_obs_espprot = apply(as.matrix(table(phenoclim$espprot, phenoclim$year)),1,function(X){c(2005:2025)[X!=0]}) 
pres_obs_espprot = apply(phenoclim[,colnames(tab_recap_Phenoclim_x_espprot5km)[-c(1:5)]], 2, function(X){unique(phenoclim$year[!is.na(X)])}) 
# (pour l'Espace Mont-Blanc on a bien des obs toutes les années, de 2005 à 2025)
pres_obs_espprot = c(pres_obs_espprot, list("EMB"=2005:2025))

clim_pourpred = clim_pourpred %>% filter(annee %in% pres_obs_espprot[[espprot]])


# On peut donc utiliser les modèles de chaque espèce pour faire les prédictions de début de sénescence et de "vitesse" de sénescence (temps entre le passage de 10 à 50%)

indice_senes = clim_pourpred %>% rename_with(~str_replace(.,"_ann","")) %>% rename(year = annee)
indice_senes$cl_alt = factor(indice_senes$cl_alt, levels=c("150-450","450-750" ,  "750-1050"   ,"1050-1350" , "1350-1650" ,"1650-1950", "1950-2250" ), 
                           ordered = T)
indice_senes$yearQ = factor(indice_senes$year)
# On retire les cas qui n'ont pas été utilisés dans la calibration des modèles parce que pas de donnée climatique 
indice_senes = indice_senes %>% filter(!is.na(Tmoy_GS))


# Pour certains modèles, il faut aussi la date de débourrement --> on la récupère via l'indice de débourrement
indice_deb = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_deb_espprotbuff5km___20062025_varTmedann__allspecies.csv")
indice_senes = merge(indice_senes, indice_deb[,c("espprot","year","esp","cl_alt","pred_deb_fixef", "pred_deb_allef")],
                     by=c("espprot","year","esp","cl_alt"), all.x=T, all.y=F)

# PREDICTIONS
# RQ : on peut choisir de prédire en prenant en compte les effets aléatoires ou sans. Ici les effets aléatoires sont ID_zone (effet local, qu'on
#      cherche à lisser en agrégeant par département) et yearQ (effet annuel lié aux spécificités climatiques non prises en compte dans le modèle,
#      ce qu'on peut vouloir garder). Dans les prédictions V1, l'effet ID_zone n'est pas pris en compte, et l'effet year est pris en compte comme 
#      effet fixe dans le modèle. Dans la V2, on choisit de faire avec et sans l'effet aléatoire yearQ (en factor).

# *---- 1) Début de sénescence
# Chargement des best models
load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes10_20062025.Rdata") 

for (esp in c("Bouleau_verruqueux","Bouleau_pubescent","Hetre","Meleze","Sorbier")){
  print(esp)
  # Prédiction

  indice_senes$pred_senes10_fixef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_fixef) %>% mutate(debou10=debou10/100), 
                                                                     re.form=NA, allow.new.levels=T)
  indice_senes$pred_senes10_allef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_allef) %>% mutate(debou10=debou10/100), 
                                                                     re.form=~(1|yearQ), allow.new.levels=T)

}

# *---- 2) Milieu de sénescence
# Chargement des best models
load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senes50_20062025.Rdata") 

for (esp in c("Bouleau_verruqueux","Bouleau_pubescent","Hetre","Meleze","Sorbier")){
  print(esp)
  # Prédiction
  
  indice_senes$pred_senes50_fixef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_fixef) %>% mutate(debou10=debou10/100), 
                                                                     re.form=NA, allow.new.levels=T)
  indice_senes$pred_senes50_allef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_allef) %>% mutate(debou10=debou10/100), 
                                                                     re.form=~(1|yearQ), allow.new.levels=T)
  
}

# *---- 3) "Vitesse" de sénescence
# Chargement des best models
load("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/bestmods_senesduree_20062025.Rdata") 

for (esp in c("Bouleau_verruqueux","Bouleau_pubescent","Hetre","Meleze","Sorbier")){
  print(esp)
  # Prédiction
  indice_senes$pred_senesduree_fixef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_fixef) %>% mutate(debou10=debou10/100), 
                                                                        re.form=NA, allow.new.levels=T, type="response")
  indice_senes$pred_senesduree_allef[indice_senes$esp == esp] = predict(bestmods[[esp]], indice_senes[indice_senes$esp == esp,] %>% rename(debou10 = pred_deb_allef) %>% mutate(debou10=debou10/100), 
                                                                        re.form=~(1|yearQ), allow.new.levels=T, type="response")
  
}

#==============================================================================================================-
# /!\ il y a des indices de durée de sénescence aberrants pour le hêtre !! (plus de 6 mois de sénescence !?)
#     => on filtre les valeurs qui s'écartent de la moyenne de plus de 2.5 fois la median absolute deviation (cf méthode Vitasse et al 2017 pour filtrer les données phéno)
valaberrantes = indice_senes %>% group_by(esp) %>% 
  summarise(MADsenesduree = mad(pred_senesduree_allef, na.rm=T),
            nb_pred = length(pred_senesduree_allef),
            mean = mean(pred_senesduree_allef, na.rm=T),
            nb_pred_in = length(pred_senesduree_allef[abs(pred_senesduree_allef - mean(pred_senesduree_allef, na.rm=T))<= 2.5*mad(pred_senesduree_allef, na.rm=T)]),
            mean_in = mean(pred_senesduree_allef[abs(pred_senesduree_allef - mean(pred_senesduree_allef, na.rm=T))<= 2.5*mad(pred_senesduree_allef, na.rm=T)], na.rm=T))

indice_senes = indice_senes %>% group_by(esp) %>% mutate(valaberrante_senesduree = ifelse(abs(pred_senesduree_allef - mean(pred_senesduree_allef, na.rm=T))<= 2.5*mad(pred_senesduree_allef, na.rm=T), "non","oui"),
                                                         valaberrante_senes10 = ifelse(abs(pred_senes10_allef - mean(pred_senes10_allef, na.rm=T))<= 2.5*mad(pred_senes10_allef, na.rm=T), "non","oui"),
                                                         valaberrante_senes50 = ifelse(abs(pred_senes50_allef - mean(pred_senes50_allef, na.rm=T))<= 2.5*mad(pred_senes50_allef, na.rm=T), "non","oui")) %>% 
  ungroup()


write.csv(indice_senes, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_espprotbuff5km___20062025_varTmedann__allspecies.csv", row.names = F)





indice_senes_espprot = indice_senes %>% filter(valaberrante_senesduree=="non" & valaberrante_senes10=="non" & valaberrante_senes50=="non") %>% group_by(year, espprot, cl_alt, cl_alt2) %>% 
  summarise(nb_esp = length(unique(esp)),
            pred_senes10_fixef = mean(pred_senes10_fixef, na.rm=T),
            pred_senes10_allef = mean(pred_senes10_allef, na.rm=T),
            pred_senes50_fixef = mean(pred_senes50_fixef, na.rm=T),
            pred_senes50_allef = mean(pred_senes50_allef, na.rm=T),
            pred_senesduree_fixef = mean(pred_senesduree_fixef, na.rm=T),
            pred_senesduree_allef = mean(pred_senesduree_allef, na.rm=T)
            )

write.csv(indice_senes_espprot, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_espprotbuff5km___20062025_varTmedann.csv", row.names = F)


indice_senes_Alps = indice_senes %>% filter(valaberrante_senesduree=="non" & valaberrante_senes10=="non") %>% group_by(year, cl_alt, cl_alt2) %>% 
  summarise(nb_esp = length(unique(esp)),
            pred_senes10_fixef = mean(pred_senes10_fixef, na.rm=T),
            pred_senes10_allef = mean(pred_senes10_allef, na.rm=T),
            pred_senes50_fixef = mean(pred_senes50_fixef, na.rm=T),
            pred_senes50_allef = mean(pred_senes50_allef, na.rm=T),
            pred_senesduree_fixef = mean(pred_senesduree_fixef, na.rm=T),
            pred_senesduree_allef = mean(pred_senesduree_allef, na.rm=T)
  )


#----- Visualisation des indices par espace protégé

# INFO sur le nombre d'observations réelles pour l'espace en question
phenoclim_senes = phenoclim %>% filter(pheno_stade_value=="Changement couleur")
INFO = phenoclim_senes %>% filter(species %in% c("Bouleau_verruqueux","Bouleau_pubescent","Hetre","Meleze","Sorbier")) %>% group_by(espprot) %>% summarise(
  nb_obs = length(unique(id_base_visit)),
  nb_arbres = length(unique(id_base_site)),
  min_annee = min(year),
  max_annee = max(year)
)
# INFO = INFO[INFO$espprot %in% indice_senes_espprot$espprot,]

INFO = rbind(INFO,
             phenoclim %>% filter(species %in% c("Bouleau_verruqueux","Bouleau_pubescent","Hetre","Meleze","Sorbier")) %>% filter(!is.na(EMB)) %>% summarise(
               espprot = "EMB",
               nb_obs = length(unique(id_base_visit)),
               nb_arbres = length(unique(id_base_site)),
               min_annee = min(year),
               max_annee = max(year)
             ))

for (espaceprot in colnames(tab_recap_Phenoclim_x_espprot5km)[-c(1:5)]){
  # for (espaceprot in INFO$espprot){
  print(espaceprot)
  if(espaceprot %in% INFO$espprot){
    INFO$nb_obs_buff5kminclus[INFO$espprot == espaceprot] = length(unique(phenoclim_senes$id_base_visit[!is.na(phenoclim_senes[,espaceprot])]))
    INFO$nb_arbres_buff5kminclus[INFO$espprot == espaceprot] = length(unique(phenoclim_senes$id_base_site[!is.na(phenoclim_senes[,espaceprot])]))
  } else {
    INFO = rbind(INFO, data.frame(espprot = espaceprot, nb_obs = 0, nb_arbres=0, min_annee=NA, max_annee=NA, 
                                  nb_obs_buff5kminclus = length(unique(phenoclim_senes$id_base_visit[!is.na(phenoclim_senes[,espaceprot])])), 
                                  nb_arbres_buff5kminclus = length(unique(phenoclim_senes$id_base_site[!is.na(phenoclim_senes[,espaceprot])]))))
  }
}
INFO = INFO[!is.na(INFO$espprot),]

write.csv(INFO, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_espprotbuff5km___20062025_varTmedann__NBOBS.csv", row.names = F)




#*---------- Visualisation des indices PAR ESPACE PROTEGE ----
indice_senes_espprot$cl_2alt = ifelse(indice_senes_espprot$cl_alt2 == "<1050", "Inf1050", "Sup1050")

indice_senes_espprot_V2 = indice_senes_espprot[!is.na(indice_senes_espprot$espprot),] %>% group_by(espprot,cl_2alt,year) %>% 
  summarise(pred_senes10_fixef = mean(pred_senes10_fixef, na.rm=T),
            pred_senes10_allef = mean(pred_senes10_allef, na.rm=T),
            pred_senes50_fixef = mean(pred_senes50_fixef, na.rm=T),
            pred_senes50_allef = mean(pred_senes50_allef, na.rm=T),
            pred_senesduree_fixef = mean(pred_senesduree_fixef, na.rm=T),
            pred_senesduree_allef = mean(pred_senesduree_allef, na.rm=T)
  )
indice_senes_espprot_V2 = merge(indice_senes_espprot_V2,
                              indice_senes_espprot %>% group_by(espprot,cl_2alt) %>% summarise(pred_senes10_fixef_ref = mean(pred_senes10_fixef),
                                                                                             pred_senes10_allef_ref = mean(pred_senes10_allef),
                                                                                             pred_senes50_fixef_ref = mean(pred_senes50_fixef),
                                                                                             pred_senes50_allef_ref = mean(pred_senes50_allef)),
                              by=c("espprot", "cl_2alt"))
indice_senes_espprot_V2$diffsenes10_allef = indice_senes_espprot_V2$pred_senes10_allef - indice_senes_espprot_V2$pred_senes10_allef_ref
indice_senes_espprot_V2$diffsenes10_fixef = indice_senes_espprot_V2$pred_senes10_fixef - indice_senes_espprot_V2$pred_senes10_fixef_ref
indice_senes_espprot_V2$diffsenes50_allef = indice_senes_espprot_V2$pred_senes50_allef - indice_senes_espprot_V2$pred_senes50_allef_ref
indice_senes_espprot_V2$diffsenes50_fixef = indice_senes_espprot_V2$pred_senes50_fixef - indice_senes_espprot_V2$pred_senes50_fixef_ref


# Différentes représentations sont possibles, pour visualiser à la fois le début, les 50% et/ou une vitesse de sénescence
# - point pour 10% de sénescence + flèche à droite du point pour la vitesse
# - point pour 50% de sénescence + flèche à gauche du point pour la vitesse
# - point pour 10% et point pour 50%
# On choisit une option ou l'autre en fonction de la qualité des modèles pour 10%, 50% et la durée de sénescence.

#-----------------------------------------------------------------------------------------------------------------------*
# Choix des deux meilleurs estimateurs en fonction des R2
quali_senesduree = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senesduree_20062025__R2.csv")
quali_senes10 = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes10_20062025__R2.csv")
quali_senes50 = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/resultats_modeles_senes50_20062025__R2.csv")
# on compare les R2_bestmod_allef
quali = as.data.frame.array(t(data.frame(quali_senes10=quali_senes10$R2_bestmod_allef,
                   quali_senes50=quali_senes50$R2_bestmod_allef,
                   quali_senesduree=quali_senesduree$R2_bestmod_allef)))
colnames(quali)=quali_senes10$species
quali$R2_allef_mean = apply(quali,1,mean)
# => on choisit donc les prédictions de senes10 et senes50, qui ont les meilleurs R2
#-----------------------------------------------------------------------------------------------------------------------*

# plot_indice_espprot_fixef = ggplot(indice_senes_espprot_V2, 
#                                    aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
#                                        x=diffsenes10_fixef, col=cl_2alt)) + 
#   geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
#                 xmin = -Inf, xmax = Inf, color=NULL, 
#                 fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
#   geom_hline(yintercept = c(2006:2026)-0.5, col="white", lty=2, lwd=0) +
#   geom_vline(xintercept = 0, col="black", lty=2) +
#   geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
#   geom_point(shape=15, size=2) + 
#   geom_segment(aes(x=diffsenes10_fixef, xend=-Inf, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=3, lwd=0.7)+
#   scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
#   # scale_color_discrete(type=terrain.colors(7))+
#   scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
#   scale_x_continuous(breaks=seq(-10,10, by=5), labels=seq(-10,10, by=5), limits=c(-12,12)) +
#   labs(x="", y="", title = "Indice d'automne - V2 fixed effects")   + facet_wrap(~espprot, scales="free") + #xlim(50,170)+
#   theme(legend.position = "none", 
#         panel.border = element_blank(),  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank())
# plot_indice_espprot_allef = ggplot(indice_senes_espprot_V2, 
#                                    aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
#                                        x=diffsenes10_allef, col=cl_2alt)) + 
#   geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
#                 xmin = -Inf, xmax = Inf, color=NULL, 
#                 fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
#   geom_hline(yintercept = c(2006:2026)-0.5, col="white", lty=2, lwd=0) +
#   geom_vline(xintercept = 0, col="black", lty=2) +
#   geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
#   geom_point(shape=15, size=2) + 
#   geom_segment(aes(x=diffsenes10_allef, xend=-Inf, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=3, lwd=0.7)+
#   scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
#   # scale_color_discrete(type=terrain.colors(7))+
#   scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
#   scale_x_continuous(breaks=seq(-10,10, by=5), labels=seq(-10,10, by=5), limits=c(-12,12)) +
#   labs(x="", y="", title = "Indice d'automne - V2 all effects")   + facet_wrap(~espprot, scales="free") + #xlim(50,170)+
#   theme(legend.position = "none", 
#         panel.border = element_blank(),  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank())

plot_indice_espprot_fixef = ggplot(indice_senes_espprot_V2, 
                                aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                    x=diffsenes10_fixef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diffsenes10_fixef, xend=-Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  geom_segment(aes(x=diffsenes10_fixef, xend=diffsenes10_fixef + pred_senes50_fixef - pred_senes10_fixef, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
  scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
  scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
  labs(x="", y="", title = "Indice d'automne - V2 fixed effects")   + facet_wrap(~espprot) +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
plot_indice_espprot_allef = ggplot(indice_senes_espprot_V2, 
                                aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                    x=diffsenes10_allef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diffsenes10_allef, xend=-Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  geom_segment(aes(x=diffsenes10_allef, xend=diffsenes10_allef + pred_senes50_allef - pred_senes10_allef, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
  scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
  scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
  labs(x="", y="", title = "Indice d'automne - V2 all effects")   + facet_wrap(~espprot) +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())




pdf("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_espprotbuff5km___20062025_varTmedann.pdf", width=25, height = 35)
# pdf("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_espprot___20062025_varTmedann.pdf", width=15, height = 17)
gridExtra::grid.table(INFO, rows=NULL)
plot_indice_espprot_fixef ; plot_indice_espprot_allef
dev.off()


#*---------- Visualisation des indices POUR TOUTES LES ALPES ----

indice_senes_Alps$cl_2alt = ifelse(indice_senes_Alps$cl_alt2 == "<1050", "Inf1050", "Sup1050")

# indice_senes_Alps_V2 = indice_senes_Alps %>% group_by(cl_2alt,year) %>% summarise(
#             pred_senes10_fixef = mean(pred_senes10_fixef, na.rm=T),
#             pred_senes10_allef = mean(pred_senes10_allef, na.rm=T),
#             pred_senesduree_fixef = mean(pred_senesduree_fixef, na.rm=T),
#             pred_senesduree_allef = mean(pred_senesduree_allef, na.rm=T)
#             
# )
# Pour le hêtre, les durées de sénescence sont particulièrement étalées !! MODELE TOUT NUL ??
indice_senes_Alps_V2 = indice_senes_Alps %>% filter() %>% group_by(cl_2alt,year) %>% summarise(
  pred_senes10_fixef = mean(pred_senes10_fixef, na.rm=T),
  pred_senes10_allef = mean(pred_senes10_allef, na.rm=T),
  pred_senes50_fixef = mean(pred_senes50_fixef, na.rm=T),
  pred_senes50_allef = mean(pred_senes50_allef, na.rm=T),
  pred_senesduree_fixef = mean(pred_senesduree_fixef, na.rm=T),
  pred_senesduree_allef = mean(pred_senesduree_allef, na.rm=T)
  
)
indice_senes_Alps_V2 = merge(indice_senes_Alps_V2,
                           indice_senes_Alps %>% group_by(cl_2alt) %>% summarise(pred_senes10_fixef_ref = mean(pred_senes10_fixef),
                                                                               pred_senes10_allef_ref = mean(pred_senes10_allef),
                                                                               pred_senes50_fixef_ref = mean(pred_senes50_fixef),
                                                                               pred_senes50_allef_ref = mean(pred_senes50_allef)),
                           by=c("cl_2alt"))
indice_senes_Alps_V2$diffsenes10_allef = indice_senes_Alps_V2$pred_senes10_allef - indice_senes_Alps_V2$pred_senes10_allef_ref
indice_senes_Alps_V2$diffsenes10_fixef = indice_senes_Alps_V2$pred_senes10_fixef - indice_senes_Alps_V2$pred_senes10_fixef_ref
indice_senes_Alps_V2$diffsenes50_allef = indice_senes_Alps_V2$pred_senes50_allef - indice_senes_Alps_V2$pred_senes50_allef_ref
indice_senes_Alps_V2$diffsenes50_fixef = indice_senes_Alps_V2$pred_senes50_fixef - indice_senes_Alps_V2$pred_senes50_fixef_ref



# plot_indice_Alps_fixef = ggplot(indice_senes_Alps_V2, 
#                                 aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
#                                     x=diffsenes10_fixef, col=cl_2alt)) + 
#   geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
#                 xmin = -Inf, xmax = Inf, color=NULL, 
#                 fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
#   geom_vline(xintercept = 0, col="black", lty=2) +
#   geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
#   geom_point(shape=15, size=2) + 
#   geom_segment(aes(x=diffsenes10_fixef, xend=-Inf, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=3, lwd=0.7)+
#   geom_segment(aes(x=diffsenes10_fixef, xend=diffsenes10_fixef + pred_senesduree_fixef, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
#   scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
#   # scale_color_discrete(type=terrain.colors(7))+
#   scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
#   scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
#   labs(x="", y="", title = "Indice d'automne - V2 fixed effects")  +
#   theme(legend.position = "none", 
#         panel.border = element_blank(),  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank())
# plot_indice_Alps_allef = ggplot(indice_senes_Alps_V2, 
#                                 aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
#                                     x=diffsenes10_allef, col=cl_2alt)) + 
#   geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
#                 xmin = -Inf, xmax = Inf, color=NULL, 
#                 fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
#   geom_vline(xintercept = 0, col="black", lty=2) +
#   geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
#   geom_point(shape=15, size=2) + 
#   geom_segment(aes(x=diffsenes10_allef, xend=-Inf, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=3, lwd=0.7)+
#   geom_segment(aes(x=diffsenes10_allef, xend=diffsenes10_allef + pred_senesduree_allef, 
#                    y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
#                    yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
#                    col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
#   scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
#   # scale_color_discrete(type=terrain.colors(7))+
#   scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
#   scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
#   labs(x="", y="", title = "Indice d'automne - V2 all effects")  +
#   theme(legend.position = "none", 
#         panel.border = element_blank(),  
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.ticks = element_blank())

plot_indice_Alps_fixef = ggplot(indice_senes_Alps_V2, 
                                aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                    x=diffsenes10_fixef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diffsenes10_fixef, xend=-Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  geom_segment(aes(x=diffsenes10_fixef, xend=diffsenes10_fixef + pred_senes50_fixef - pred_senes10_fixef, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
  scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
  scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
  labs(x="", y="", title = "Indice d'automne - V2 fixed effects")  +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())
plot_indice_Alps_allef = ggplot(indice_senes_Alps_V2, 
                                aes(y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))), 
                                    x=diffsenes10_allef, col=cl_2alt)) + 
  geom_rect(aes(ymax = year + 0.5, ymin = year - 0.5, 
                xmin = -Inf, xmax = Inf, color=NULL, 
                fill = factor(ifelse(year %% 2 == 0, 1,0))), alpha = 0.8) + scale_fill_manual(values=c("gray90","white"))+
  geom_vline(xintercept = 0, col="black", lty=2) +
  geom_vline(xintercept = c(-10,-5,5,10), col="gray40", lty=1, lwd=0.05) +
  geom_point(shape=15, size=2) + 
  geom_segment(aes(x=diffsenes10_allef, xend=-Inf, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=3, lwd=0.7)+
  geom_segment(aes(x=diffsenes10_allef, xend=diffsenes10_allef + pred_senes50_allef - pred_senes10_allef, 
                   y=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))) , 
                   yend=year + 0.1*as.numeric(as.character(factor(cl_2alt, levels=unique(cl_2alt), labels=c(-1,1)))),
                   col=cl_2alt ), lty=1, lwd=0.5, arrow = arrow(length=unit(0.1, "cm"))) +
  scale_color_manual(values=c("darkorange2","gold2"), breaks=c("Inf1050","Sup1050")) +
  # scale_color_discrete(type=terrain.colors(7))+
  scale_y_continuous(breaks=seq(2006,2025, by=1), labels=seq(2006,2025, by=1))+
  scale_x_continuous(breaks=seq(-10,30, by=5), labels=seq(-10,30, by=5), limits=c(-15,35)) +
  labs(x="", y="", title = "Indice d'automne - V2 all effects")  +
  theme(legend.position = "none", 
        panel.border = element_blank(),  
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank())



pdf("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/indice_senes_Alps___20062025_varTmedann.pdf", width=4, height = 5)
plot_indice_Alps_fixef ; plot_indice_Alps_allef
dev.off()







