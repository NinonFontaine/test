#########################################################################################*
# Projet Floraison d'Altitude : premières explorations ----
#########################################################################################*


library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)

library(brms)
# library(tidybayes)
# library(ciTools)

library(rinat)


#########################################################################################*
# MISE EN FORME DES DONNEES ----
#########################################################################################*

# (rebasculé dans le script 1_mise_en_forme_BDD.R)


data_long <- read.csv("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/06. ANALYSES/FloraisonAltitude/data_long.csv")

# Aperçu de la localisation des points
cle_ggmap ="" # à récupérer
ggmap::register_google(key=cle_ggmap)
map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 8,
                     maptype = "terrain", scale = 2)
ggmap(map_base) +
  geom_point(data=data_long, aes(x=coord_x_4326, y=coord_y_4326))



# Transformation du tableau avec une ligne = une date un quadrat et un stade, 
#             en un tableau avec une ligne = une date un quadrat tous les stades
data_parQ <- data_long %>%
  pivot_wider(names_from = etape, values_from = value, names_prefix = "stage")

data_parQ$counts = rowSums(data_parQ[grep("stage", colnames(data_parQ))], na.rm = T)
data_parQ$visit_date = as.Date(data_parQ$visit_date, format="%d/%m/%Y")


# Calcul des proportions de chaque stade pour chaque quadrat-date
data_parQ <- data_parQ %>% mutate(across(starts_with("stage"), ~ .x/counts, .names = "prop_{.col}"))
data_parQ <- data_parQ %>% mutate(across(starts_with("prop_stage"), ~ replace_na(.x, 0)))


# Simplification du tableau en gardant les informations principales
resume_parQ <-  data_parQ[,c("id_site","quadrat","visit_date","counts", paste0("prop_stage",1:4))]

ggplot(resume_parQ, aes(x=visit_date)) + facet_wrap(. ~ id_site) + 
  geom_point(aes(y=prop_stage1), col="green")+ 
  geom_point(aes(y=prop_stage2), col="pink")+ 
  geom_point(aes(y=prop_stage3), col="lightblue")+ 
  geom_point(aes(y=prop_stage4), col="blue")


# Résumé par site (= moyennes entre les quadrats)
# /!\ QUESTION : moyenne des proportions ou moyenne totale ? 
#     Sommer les comptages de tous les quadrats me semble plus juste, puisqu'on réfléchit à la phéno du site --> moyenne totale
#     (si on voulait la moyenne des proportions : across(starts_with("prop_stage"), ~ mean(.x), .names = "{.col}") )
resume <- data_parQ[data_parQ$counts != 0,] %>% group_by(id_site, visit_date) %>%
  summarise(tot_counts = sum(counts),
            sd_counts = sd(counts, na.rm=T),
            across(starts_with("prop_stage"), ~ sd(.x, na.rm=T), .names = "sd_{.col}"),
            across(starts_with("stage"), ~ sum(.x, na.rm=T), .names = "{.col}"))
resume <- resume %>% mutate(across(starts_with("stage"), ~ .x/tot_counts, .names = "prop_{.col}"))

resume$date_no <- yday(resume$visit_date)

#########################################################################################*
# VISUALISATION DES DONNEES BRUTES ----
#########################################################################################*

# Aperçu des données phéno sur les différents sites

colors= c("stage1"="green", "stage2"="pink", "stage3"="lightblue", "stage4"="blue")
ggplot(na.omit(resume), aes(x=date_no)) + facet_wrap(. ~ id_site) + 
  geom_point(aes(y=prop_stage1, col="stage1"))+ 
  geom_errorbar(aes(ymin=prop_stage1 - sd_prop_stage1, ymax=prop_stage1 + sd_prop_stage1, col="stage1"))+
  geom_line(aes(y=prop_stage1, col="stage1"))+ 
  geom_point(aes(y=prop_stage2, col="stage2"))+ 
  geom_errorbar(aes(ymin=prop_stage2 - sd_prop_stage2, ymax=prop_stage2 + sd_prop_stage2, col="stage2"))+
  geom_line(aes(y=prop_stage2, col="stage2"))+ 
  geom_point(aes(y=prop_stage3, col="stage3"))+ 
  geom_errorbar(aes(ymin=prop_stage3 - sd_prop_stage3, ymax=prop_stage3 + sd_prop_stage3, col="stage3"))+
  geom_line(aes(y=prop_stage3, col="stage3"))+ 
  geom_point(aes(y=prop_stage4, col="stage4"))+ 
  geom_errorbar(aes(ymin=prop_stage4 - sd_prop_stage4, ymax=prop_stage4 + sd_prop_stage4, col="stage4"))+
  geom_line(aes(y=prop_stage4, col="stage4"))+
  labs(y="Proportion of counts per stage", x="Date", col="")+
  scale_color_manual(values=colors)

#########################################################################################*
# EXPLORATIONS ----
#########################################################################################*

# Aperçu de la variabilité de production de fruits intersites et intrasites 
# (pour réfléchir au nombre de Q et nombre de sites pour la caractérisation
# de la fructification en fonction des conditions environnementales)
#
# --> on se concentre sur la fin de la fructification, quand il n'y a plus que des fruits mûrs (stage 4)
data_estimvarfruit = data_parQ[data_parQ$counts != 0 & data_parQ$stage4 == data_parQ$counts & !is.na(data_parQ$stage4),]

var_inter = var(data_estimvarfruit$stage4) 
testvar = data_estimvarfruit %>% group_by(id_site, visit_date) %>% summarise(
  var_intra = var(stage4),
  moy = mean(stage4)
)
var_inter2 = var(testvar$moy)

# => var_inter >> var_intra : mieux vaut se concentrer sur un plus grand nombre de sites avec un plus petit nombre de quadrats par site

# --------------------------------------------------------------------------*

#########################################################################################*
# TESTS DIVERS POUR FITTER DES COURBES PHENO AVEC LES DONNEES FLO D'ALTI ----
#########################################################################################*


# # --------------------------------------------------------------------------------------*
# # Tests pour avoir les phenophases, en utilisant les fonctions de base du package **phenor**
# # if(!require(remotes)){install.packages("remotes")}
# # remotes::install_github("bluegreen-labs/phenor")
# # install.packages("phenocamr")
# library(phenocamr)
# library(phenor)
# 
# phenocamr::download_phenocam(
#   frequency = 3,
#   veg_type = "DB",
#   roi_id = 1000,
#   site = "bartlettir",
#   phenophase = TRUE,
#   out_dir = "."
# )
# df <- read.table("bartlettir_DB_1000_3day.csv", header = TRUE, sep = ",")
# 
# # --------------------------------------------------------------------------------------*
# # Tests pour avoir les phenophases, en utilisant les fonctions de base du package **phenology**
# library(phenology)
# 
# resume$visit_date = as.Date(resume$visit_date, format="%Y/%m/%d")
# datatest_fruc24010 = add_phenology(resume[resume$id_site == "24010",c("visit_date", "prop_stage4")],
#                                    reference = as.Date("2024/05/01"), format="%Y/%m/%d")
# 
# # Generate initial points for the optimisation
# parg <- par_init(datatest_fruc24010, fixed.parameters=NULL)
# # Run the optimisation
# result_fruc24010 <- fit_phenology(data=datatest_fruc24010, fitted.parameters=parg,
#                                 fixed.parameters=NULL)
# # data(result_Gratiot)
# # Plot the phenology and get some stats
# output <- plot(result_fruc24010)


# --------------------------------------------------------------------------------------*
# Tests pour avoir les phenophases, en utilisant des fonctions de base de smoothing

# Focalisation sur un site (24010 = myrtille Mont Lachat) et une transition (stade 3 à stade 4)
sitetest = "24010"
stagetest = "stage4"
datatest = resume[resume$id_site == sitetest,c("date_no", stagetest, paste0("prop_",stagetest), "counts")]
colnames(datatest) = c("date_no","stage","prop_stage","counts")


# Pour des proportions comme on a ici (proportion de trucs comptés qui sont au stade fleur ou fruit pas mûr ou fruit),
# il y a plusieurs options / types de réponse possibles : logit, betarégression, zero-inflated beta reg, ...
# https://rcompanion.org/handbook/J_02.html
# - une régression logistique serait valable si on avait toujours le même nombre total de trucs comptés, ce qui n'est pas le cas
# - une régression de poisson fonctionnerait puisqu'on part de comptes, mais il faudrait mettre un offset correspondant au nombre
#   total de trucs comptés
# ...
# Au final, avec le type de données qu'on a sur Flo d'Alti, i.e. "Proportion / Ratio (including zero and one)", on peut faire
# du zero-one-inflated beta (soit en bayésien : brm(family = zero_one_inflated_beta()))
# https://strengejacke.github.io/regressionmodels/

#########################################################################################*
#*---- Quelques tests avec différentes fonctions de smoothing ----
# --------------------------------------------------------------------------------------*
# - avec glm (fonction logit)
sigmfit_logit <- glm(data = datatest,
                         prop_stage~date_no, family="binomial")
# # - avec loess
# sigmfit_loess <- loess(data = datatest,
#                          prop_stage~date_no) # à voir pour ajuster le span
# - en bayésien (cf tests Marjo : le bayésien évite les problèmes de convergence) - package brms
#           pour la version en bayésien, Marjo avait mis des poids aux obs, pour que chaque jour ait un poids total de 1
#               (i.e. s'il y a x obs différentes pour un jour, chacune a un poids de 1/x)
# #       -> logit (proportion de succès/échec).  BIF BOF
# sigmfit_baylogit <- brm(data = datatest, 
#                       family = binomial(),
#                       bf(stage | trials(counts) ~ date_no))
#       -> betareg    /!\ il ne faut pas de 0 ou de 1 !
datatest_sans01 = datatest
datatest_sans01$prop_stage[datatest_sans01$prop_stage == 0] = 0.0001
datatest_sans01$prop_stage[datatest_sans01$prop_stage == 1] = 0.9999
sigmfit_baybeta <- brm(data = datatest_sans01, 
                         family = Beta(),
                         prop_stage ~ date_no) 
#       -> zero inflated (parce qu'on a beaucoup de 0, notamment si on s'intéresse au début de la courbe de transition)
datatest_sans1 = datatest
datatest_sans1$prop_stage[datatest_sans1$prop_stage == 1] = 0.9999 # on en peut pas avoir de 1
sigmfit_bay0infl <- brm(data = datatest_sans1, 
                         family = zero_inflated_beta(),
                         # prop_stage~date_no)
                         bf(prop_stage ~ date_no, # la proportion de chaque stade dépend de la date
                            phi ~ 1, # la précision ne dépend pas de la date
                            zi ~ date_no)) # la proba d'avoir 0 dépend de la date
#       -> zero-one inflated (parce qu'à part à la transition, on a beaucoup de 0 et 1)
sigmfit_bay01infl <- brm(data = datatest, 
                         family = zero_one_inflated_beta(),
                         # prop_stage~date_no)
                         bf(prop_stage ~ date_no, # la proportion de chaque stade dépend de la date
                            phi ~ 1, # la précision ne dépend pas de la date
                            zoi ~ date_no, # la proba d'avoir des 0 ou des 1 dépend de la date
                            coi ~ date_no)) # la proba d'avoir 1 sachant qu'on avait 0 ou 1 dépend de la date



# smoothing
timeserie <- seq(min(datatest$date_no), max(datatest$date_no), length.out=100)
smooth_logit <- predict(sigmfit_logit, data.frame(date_no=timeserie), type="response")
# smooth_loess <- predict(sigmfit_loess, data.frame(date_no=timeserie), type="response")
smooth_baylogit <- predict(sigmfit_baylogit, data.frame(date_no=timeserie, counts=100), type="response")
smooth_baybeta <- predict(sigmfit_baybeta, data.frame(date_no=timeserie), type="response")
smooth_bay0infl <- predict(sigmfit_bay0infl, data.frame(date_no=timeserie), type="response")
smooth_bay01infl <- predict(sigmfit_bay01infl, data.frame(date_no=timeserie), type="response")

# # # Pour avoir le smoothing prédit avec intervalle de confiance 
# pred_logit = predict(sigmfit_logit, data.frame(date_no=timeserie), type="link", se.fit=T)
# pred_logit = data.frame(date_no=timeserie, 
#                         pred=sigmfit_logit$family$linkinv(pred_logit$fit), 
#                         uprCI=sigmfit_logit$family$linkinv(pred_logit$fit + (qnorm(0.975) * pred_logit$se.fit)), 
#                         lwrCI=sigmfit_logit$family$linkinv(pred_logit$fit - (qnorm(0.975) * pred_logit$se.fit)))
# # /!\ bizarre : très grands intervalles de confiance quand on est à 0% et 100% de stade d'intérêt... 
# #               -> du coup je ne peux pas extraire un intervalle de confiance autour des dates de transition...!


# récupération des dates seuil
fun__dates_seuils <-  function(data_smooth, seuils){
  dates = vector(,length=length(seuils))
  for (i in 1:length(seuils)){
    dates[i] = min(data_smooth$timeserie[data_smooth$smooth > seuils[i]])}
  return(dates)
}

# Choix de seuils
seuils = c(0.1, 0.5, 0.8)
dates_logit = fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth_logit),
                     seuils = seuils)
dates_baylogit = fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth_baylogit[,1]),
                                seuils = seuils)
dates_baybeta = fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth_baybeta[,1]),
                                  seuils = seuils)
dates_bay0infl = fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth_bay0infl[,1]),
                                    seuils = seuils)
dates_bay01infl = fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth_bay01infl[,1]),
                                    seuils = seuils)

# Visualisation du fitting pheno et des dates
plot(x=datatest$date_no, y=datatest$prop_stage, bg="black", pch=21, ylim = c(0,1),
     xlab = "Day of the year", ylab = "Transition to stage 4")
lines(x=timeserie, y=smooth_logit, col="red")
# lines(x=timeserie, y=smooth_baylogit[,1], col="pink")
lines(x=timeserie, y=smooth_baybeta[,1], col="purple")
lines(x=timeserie, y=smooth_bay0infl[,1], col="hotpink")
lines(x=timeserie, y=smooth_bay01infl[,1], col="blue")
# lines(x=timeserie, y=smooth_loess, col="blue")
for (i in 1:length(seuils)){rect(0, -1, dates_logit[i], seuils[i], lty = 2, border ="red")}
# for (i in 1:length(seuils)){rect(0, -1, dates_baylogit[i], seuils[i], lty = 2, border ="pink")}
for (i in 1:length(seuils)){rect(0, -1, dates_baybeta[i], seuils[i], lty = 2, border ="purple")}
for (i in 1:length(seuils)){rect(0, -1, dates_bay0inf[i], seuils[i], lty = 2, border ="hotpink")}
for (i in 1:length(seuils)){rect(0, -1, dates_bay01infl[i], seuils[i], lty = 2, border="blue")}

# Finalement pour ce cas-là (maturité des myrtilles sur le site 24010) :
# - le fitting bayésien logit (rose) ne colle pas du tout aux données
# - le fitting beta (violet) ne colle pas trop aux données (démarre trop tôt et pas une sigmoïde)
# - les fittings logit (rouge) ou betareg01infl (bleu) sont équivalents pour la fin de la transition... 
# - ...MAIS le logit est mieux pour le début de la transition = l'onset (la betareg01infl inflated démarre trop tôt)
#
# => cela suggère de garder le smoothing en logit


#########################################################################################*
#*---- Généralisation du fitting des courbes phéno pour tous les sites et toutes les transitions ----
# (rebasculé dans le script 2_identification_dates_transition.R)
# ----------------------------------------------------------------------------------------*

# Fonctions pour le smoothing et la récupération des dates
fun__dates_seuils <-  function(data_smooth, seuils){
  dates = vector(,length=length(seuils))
  for (i in 1:length(seuils)){
    dates[i] = min(data_smooth$timeserie[data_smooth$smooth > seuils[i]])}
  return(dates)
}
smoothingfun_list = list("glm - logit link" = "glmlogit", 
                         "bayesian - zero-inflated betaregression" = "bay0infl", 
                         "bayesian - zero-one-inflated betaregression" = "bay01infl")
fun__id_dates <- function(data, thresholds, smoothfun = c("glmlogit","bay0infl","bay01infl")){ # , direction="rising" # à voir si je rajoute ce paramètre
  
  timeserie <- seq(min(data$date_no), max(data$date_no), length.out=100)
  
  if(smoothfun=="glmlogit"){
    # Fitting
    mod_fit <- glm(data = data,
                   prop_stage~date_no, family="binomial")
    pred <- predict(mod_fit, data.frame(date_no=data$date_no), type="response")
    smooth <- as.data.frame(predict(mod_fit, data.frame(date_no=timeserie), type="response"))
    # Identification of transition dates 
    dates = data.frame("Estimate"=fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth = smooth[,1]),
                                                    seuils = thresholds))
    row.names(dates) = thresholds
  }
  else if(smoothfun=="bay0infl"){
    data$prop_stage[data$prop_stage == 1] = 0.9999 # on ne peut pas avoir de 1 dans une betaregression zero inflated
    # Fitting
    mod_fit <- brm(data = data, 
                   family = zero_inflated_beta(),
                   bf(prop_stage ~ date_no, # la proportion de chaque stade dépend de la date
                      phi ~ 1, # la précision ne dépend pas de la date
                      zi ~ date_no)) # la proba d'avoir des 0 dépend de la date
    pred <- predict(mod_fit, data.frame(date_no=data$date_no), type="response")[,1]
    smooth <- as.data.frame(predict(mod_fit, data.frame(date_no=timeserie), type="response"))
    # Identification of transition dates 
    dates = as.data.frame(apply(smooth[,c("Estimate","Q2.5","Q97.5")], 2, 
                                FUN = function(X){fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth=X), 
                                                                    seuils=thresholds)}),
                          row.names = thresholds)
  }
  else if(smoothfun=="bay01infl"){
    # Fitting
    mod_fit <- brm(data = data, 
                   family = zero_one_inflated_beta(),
                   bf(prop_stage ~ date_no, # la proportion de chaque stade dépend de la date
                      phi ~ 1, # la précision ne dépend pas de la date
                      zoi ~ date_no, # la proba d'avoir des 0 ou des 1 dépend de la date
                      coi ~ date_no)) # la proba d'avoir 1 sachant qu'on avait 0 ou 1 dépend de la date
    pred <- predict(mod_fit, data.frame(date_no=data$date_no), type="response")[,1]
    smooth <- as.data.frame(predict(mod_fit, data.frame(date_no=timeserie), type="response"))
    # Identification of transition dates 
    dates = as.data.frame(apply(smooth[,c("Estimate","Q2.5","Q97.5")], 2, 
                                FUN = function(X){fun__dates_seuils(data_smooth = data.frame(timeserie = timeserie, smooth=X), 
                                                                    seuils=thresholds)}),
                          row.names = thresholds)
  }
  pred_output=cbind(data.frame("date_no"=timeserie),smooth)
  colnames(pred_output)[1:2] = c("date_no","Estimate")
  # Diagnostic of the smoothing
  RMSE = sqrt(mean((data$prop_stage - pred)^2, na.rm=T))
  
  return(list(dates=dates, 
              pred=pred_output, 
              diag=RMSE))
  
}

# INITIALISATION
# Thresholds of interest (to be adjusted)
seuils = c(0.1, 0.5, 0.8)
# Smoothing function
smoothingfun = smoothingfun_list[["bayesian - zero-inflated betaregression"]]
# smoothingfun = smoothingfun_list[["glm - logit link"]]
# Calculations are made using day number instead of date
resume$date_no <- yday(resume$visit_date)
# Initialisation of the result tab
tab_dates_transition = data.frame(id_site = unique(resume$id_site),
                                  year = 2024 # adapter pour quand il y aura plusieurs années
)

# COMPUTATION for each site and transition
for (stage in 2:4){
  # tab_dates_transition[,paste0("tr_date_stage",stage-1, "to",stage,"_",seuils*100,"pourc")] = NA
  # print(stage)
  for (site in unique(tab_dates_transition$id_site)){
    # print(site)
    data = resume[resume$id_site == site, c("date_no", paste0("prop_stage",stage))]
    colnames(data) = c("date_no","prop_stage")
    # On ne peut pas fitter une courbe de transition dans tous les cas : on considère que le critère qui doit être vérifié c'est
    #   un écart d'au moins 50% entre les proportions minimales et maximales évaluées au cours du temps
    # On se concentre aussi sur les transitions 'rising' (dans un premier temps ?)
    min_prop = min(data$prop_stage, na.rm = T) ; max_prop = max(data$prop_stage, na.rm = T)
    if(max_prop - min_prop >= 0.5 & 
       mean(data$date_no[data$prop_stage == min_prop], na.rm=T) < mean(data$date_no[data$prop_stage == max_prop], na.rm=T)){
      data = data[data$date_no <= max(data$date_no[data$prop_stage > 0.9 * max_prop], na.rm=T),] # On sélectionne le timeframe pour lequel les proportions comptées sont en rising (en se basant sur la proportion maximale atteinte)
      # print(data)
      
      # Smoothing
      smooth_res = fun__id_dates(data = data, thresholds = seuils, smoothfun = smoothingfun)
      # print(dates)
      
      # Diagnostic of the smoothing
      tab_dates_transition[tab_dates_transition$id_site == site,
                           paste0("trD_st",stage-1, "to",stage,"_RMSE")] = smooth_res[["diag"]]
      plot(x=data$date_no, y=data$prop_stage, bg="black", pch=21, 
           ylim = c(0,1.1),#xlim=c(min(resume$date_no),max(resume$date_no)),
           xlab = "Day of the year", ylab="Transition")
      lines(x=smooth_res[["pred"]][,"date_no"], y=smooth_res[["pred"]][,"Estimate"], col="red")
      text(x=min(data$date_no), y=1, pos=4,
           labels = paste0("Site ",site," - Transition from stage ",stage-1," to ", stage, 
                           "\n(",names(smoothingfun_list)[unlist(smoothingfun_list)==smoothingfun],")"))
      
      # Dates from the smoothing
      tab_dates_transition[tab_dates_transition$id_site == site,
                           paste0("trD_st",stage-1, "to",stage,"_",seuils*100,"pourc_Est")] = smooth_res[["dates"]][,"Estimate"]
      if(smoothingfun %in% c("bay01infl","bay0infl")){ # Si on est en bayésien, on a les intervalles de confiance autour de la date estimée
        tab_dates_transition[tab_dates_transition$id_site == site,
                             paste0("trD_st",stage-1, "to",stage,"_",seuils*100,"pourc_Q2.5")] = smooth_res[["dates"]][,"Q97.5"]
        tab_dates_transition[tab_dates_transition$id_site == site,
                             paste0("trD_st",stage-1, "to",stage,"_",seuils*100,"pourc_Q97.5")] = smooth_res[["dates"]][,"Q2.5"]
        
      }
    }
  }
}


info_sites <- data_parQ[,c("id_site", "base_site_name", "altitude", "coord_x_4326", "coord_y_4326", "nom_cite")] %>% 
  distinct()

tab_dates_transition = merge(info_sites, tab_dates_transition, by="id_site")



# write.csv(tab_dates_transition, "TransitionDates_per_site__logit.csv", row.names = FALSE)
write.csv(tab_dates_transition, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/TransitionDates_per_site__bay0infl.csv", row.names=F) #../../../../../../Desktop/projetsR/FloraisonAltitude/TransitionDates_per_site__bay0infl.csv", row.names = FALSE) 
# write.csv(tab_dates_transition, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/TransitionDates_per_site__glmlogit.csv", row.names=F) #../../../../../../Desktop/projetsR/FloraisonAltitude/TransitionDates_per_site__glmlogit.csv", row.names = FALSE)



# # si on veut comparer deux fonctions de smoothing en se basant sur les RMSE :
# # on crée 2 tableaux de résultats en changeant la fonction de smoothing, et on compare les RMSE par un t.test apparié
# t.test(na.omit(as.vector(unlist(tab_dates_transition1[,grep("RMSE", colnames(tab_dates_transition1))]))),
#        na.omit(as.vector(unlist(tab_dates_transition2[,grep("RMSE", colnames(tab_dates_transition2))]))),
#        paired = T)#, alternative = "less")


#########################################################################################*
#*----   Notes - to improve ----
# --------------------------------------------------------------------------------------*

# - Je n'ai fait que du smoothing basé sur les transitions "rising". En théorie quand on change de stade on peut avoir 
#   la courbe descendante pour la fin du stade i et la courbe ascendante pour le début du stade i+1 : on pourrait faire 
#   les calculs en croisant ces deux informations
#
# - Je n'ai utilisé que du smoothing avec la fonction logit, alors que d'autres fonctions de smoothing pourraient  
#   être plus adaptées : loess ? Weibull (Rzanny et al 2024) ?
#         => le test en bayésien avec des fonctions de type betareg (ou betareg en zero-one-inflated au vu des données qu'on a)
#            n'est pas concluant : la transition démarre trop tôt / n'est pas suffisamment abrupte
#         => en bayésien ça semble mieux avec du zero-inflated, même si ce n'est pas optimal quand il y a peu de points
#         => si on compare les approches en calculant la RMSE, celle-ci apparaît significativement plus grande pour la 
#            zero-one-inflated-betaregression en bayésien (en plus le temps de calcul est beaucoup plus long !) par rapport
#            à la fonction logit en glm (différence moyenne de 0.069 entre les RMSE, p-value < 0.001)
#
# - J'ai pris les seuils de transition à 10%, 50% et 80%, mais d'autres métriques phénologiques sont possibles, par 
#   exemple le moment où le maximum est atteint
#         => pour rester cohérent avec ce qui a été fait précédemment pour les arbres, ce sont les seuils à 10% qu'on calcule
#            (la date où le maximum est atteint peut être calculée pour la floraison, mais pas pour les autres stades. Ça peut
#            être intéressant par rapport à la pollinisation ?)
#
# - Je n'ai pas intégré la variabilité intrasite (entre quadrats)
#
# - Je n'ai pas calculé les intervalles de confiance autour des dates de transition estimées
#         => l'erreur standard avec le modèle logit est énorme autour des valeurs à 0 et 100% du stade d'intérêt, du coup 
#             je ne peux pas extraire un intervalle de confiance autour des dates de transition...!
#
# (- Il faudrait adapter le script pour quand il y aura plusieurs années de données)

# --------------------------------------------------------------------------------------*



#########################################################################################*
# VISUALISATION DES TENDANCES PHENO ----
#########################################################################################*


tab_dates_transition = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/TransitionDates_per_site__bay0infl.csv")


colors= c("stage 1 to 2"="pink", "stage 2 to 3"="lightblue", "stage 3 to 4"="blue")

ggplot(tab_dates_transition, aes(x=altitude, label=id_site)) +
  geom_point(aes(y=trD_st1to2_10pourc_Est, colour = "stage 1 to 2"))+ # geom_text(aes(y=tr_date_stage1to2_10pourc))+
  geom_point(aes(y=trD_st2to3_10pourc_Est, colour = "stage 2 to 3")) +
  geom_point(aes(y=trD_st3to4_10pourc_Est, colour = "stage 3 to 4")) +
  facet_wrap(~nom_cite) +
  labs(y="Day of the year for transition between stages", colour="Transition")+
  theme(legend.position="bottom", legend.direction = "vertical")+
  scale_color_manual(values = colors) + 
  geom_smooth(method="lm", aes(y=trD_st3to4_10pourc_Est, x=altitude, colour = "stage 3 to 4"), se=F)+ 
  geom_smooth(method="lm", aes(y=trD_st2to3_10pourc_Est, x=altitude, colour = "stage 2 to 3"), se=F)


#########################################################################################*
#*---- Lien entre altitude et phénologie ----
# --------------------------------------------------------------------------------------*

# (pour pouvoir caler les dates d'observation sur SPOT)
modFruc50 = lm(data=tab_dates_transition[tab_dates_transition$nom_cite == "Vaccinium myrtillus L., 1753",],
               formula= trD_st3to4_50pourc_Est ~ altitude)
# => 6 jours de décalage par 100m d'altitude pour 50% de fruits mûrs (R2 adjusted = 0.94, pval=0.020)
# => 5 jours de décalage par 100m d'altitude pour 10% de fruits mûrs (mais R2 beaucoup moins bon : R2 adjusted = 0.48 et pval>0.05)

# (pas assez de données pour la floraison !! (= passage du stade 1 au stade 2))
# modFlow50 = lm(data=tab_dates_transition[tab_dates_transition$nom_cite == "Vaccinium myrtillus L., 1753",],
#                formula= trD_st1to2_50pourc_Est ~ altitude)




#########################################################################################*
# EXPLORATIONS BASÉES SUR LES DONNÉES PHÉNO EXISTANTES DANS iNATURALIST ----
#########################################################################################*

iNat_myrt = rbind(cbind(get_inat_obs(taxon_name = "Vaccinium myrtillus", annotation = c(12, 13), bounds = c(43.6, 4.8, 46.6, 8.2)), data.frame("PHENO" = "Flowering")),
                  cbind(get_inat_obs(taxon_name = "Vaccinium myrtillus", annotation = c(12, 14), bounds = c(43.6, 4.8, 46.6, 8.2)), data.frame("PHENO" = "Fruiting")),
                  cbind(get_inat_obs(taxon_name = "Vaccinium myrtillus", annotation = c(12, 15), bounds = c(43.6, 4.8, 46.6, 8.2)), data.frame("PHENO" = "Flower budding")),
                  cbind(get_inat_obs(taxon_name = "Vaccinium myrtillus", annotation = c(12, 21), bounds = c(43.6, 4.8, 46.6, 8.2)), data.frame("PHENO" = "No Evidence of Flowering")))
iNat_myrt$observed_on = as.Date(iNat_myrt$observed_on, format="%Y-%m-%d")
iNat_myrt$JULDAY = yday(iNat_myrt$observed_on)
  
# Localisation des observations (dans l'extraction, on s'est restreint aux Alpes occidentales)
map_base <- get_map(location = c(lon = 6.5, lat = 45.55), zoom = 7,
                    maptype = "terrain", scale = 2)
ggmap(map_base) +
  geom_point(data=iNat_myrt, aes(x=longitude, y=latitude))

# Aperçu des observations, par stade phénologique
ggplot(iNat_myrt, aes(x=JULDAY)) + geom_histogram() 
ggplot(iNat_myrt, aes(x=JULDAY, y=PHENO)) + geom_violin() + geom_point()

# Lien avec l'altitude au niveau des points observés
DEM_Europe = rast("/Users/ninonfontaine/Google Drive/Drive partagés/05. RECHERCHE/05. DONNEES/Topographie/MNT/DEM_COP89_EU_W_100m.tif")
iNat_myrt$ALTITUDE = extract(DEM_Europe, project(vect(iNat_myrt, geom=c("longitude", "latitude"), "epsg:4326"), crs(DEM_Europe)))[,2]
write.csv(iNat_myrt, "/Users/ninonfontaine/Desktop/projetsR/TEST/output/FlodAlti/apercu_dataiNat_pheno.csv", row.names = F)

ggplot(iNat_myrt[iNat_myrt$PHENO %in% c("Flowering","Fruiting"), ], 
       aes(x=JULDAY, y=factor(ifelse(ALTITUDE < 1500, "< 1500 m", ifelse(ALTITUDE < 2000, "1500 - 2000 m", "> 2000 m")),
                              levels = c("> 2000 m","1500 - 2000 m", "< 1500 m")))) + 
  geom_violin() + geom_point() + facet_wrap(PHENO~., nrow=2) + labs(x="Day of year",y="")
# ggplot(iNat_myrt, aes(x=ALTITUDE, y=JULDAY, colour = PHENO)) + geom_point() + geom_smooth(method="lm")


modFruc_iNat = lm(data=iNat_myrt[iNat_myrt$PHENO == "Fruiting",],
               formula= JULDAY ~ ALTITUDE)
modFlow_iNat = lm(data=iNat_myrt[iNat_myrt$PHENO == "Flowering",],
                  formula= JULDAY ~ ALTITUDE)
# le modèle de flowering en fonction de l'altitude marche plutôt bien (R2 ajusté = 0.646 - 3 jours de retard par 100m d'altitude), mais celui de 
# fructification moins (R2 ajusté = 0.050 - 2 jours de retard par 100m d'altitude)
# HYPOTHÈSES : 
# - c'est un stade qui dure plus longtemps, donc les observations sont plus disparates en terme de dates (on ne distingue pas le 10-20-50 % de 
#   fruits mûrs, et donc un certain moment de cette phénophase)
# - d'autres paramètres que l'altitude entrent en compte, par exemple l'effet sécheresse ?

