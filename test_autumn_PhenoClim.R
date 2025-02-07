library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)
library(terra)
library(brms)
library(sjPlot)


phenoclim = read.csv("/Users/ninonfontaine/Desktop/projetsR/TEST/output/PhenoClim/Phenoclim_data_cleaned.csv") 

pheno_aut10 = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 10%",]
ggplot(pheno_aut10, aes(x=yearQ, y=julian_day, col=altitude)) +
  geom_point(size=0.5)+
  geom_smooth(method="lm", formula=y~x) +
  facet_wrap(.~ nom_cite)


pheno_aut50 = phenoclim[phenoclim$pheno_etape_value == "Changement de couleur - Ok 50%",]
ggplot(pheno_aut50, aes(x=yearQ, y=julian_day)) +
  geom_point(size=0.5)+
  geom_smooth(method="lm", formula=y~x) +
  facet_wrap(.~ nom_cite)


# TESTS SUR LE MÉLÈZE
pheno_aut10_meleze = pheno_aut10[pheno_aut10$nom_cite == "Larix decidua Mill., 1768",]
model_meleze <- brm(
  bf(julian_day ~ altitude+yearQ+(yearQ|nom_zone)),
  data = pheno_aut10_meleze,
  init = "0",
  chains = 4, iter = 5000, warmup = 1000,
  cores = 4
) # /!\ PB convergence !!!

# plot(model_meleze) # to look at chain mixing
# summary(model_meleze)
summary(model_meleze)$fixed
## + 2.8 j / décennie

# Visualisation des effets aléatoires
# plot_model(model_meleze, type="re")
# plot_model(model_meleze, type="re", axis.lim=c(-2,2))
plot_model(model_meleze, type="emm", terms="yearQ")



#____________________________________
#
# L'idée serait d'intégrer d'autres variables que l'année, des variables climatiques notamment, pour mieux comprendre ce qui affecte
# la phénologie automnale.
# Parmi les variables potentielles (biblio, notamment Zani et al 2020) :
# - date de débourrement OU date de feuillaison
# - température estivale (en lien avec la productivité estivale et la photosynthèse)
# - précipitations estivales (en lien avec la sécheresse, qui limite les transferts de sucres entre organes)
# (- température automnale (ce n'est pas un facteur explicatif d'après l'étude de Zani et al 2020, mais c'est couramment utilisé dans 
#    d'autres études))
# (- photopériode (ce n'est pas un facteur explicatif d'après l'étude de Zani et al 2020, mais c'est couramment utilisé dans d'autres 
#    études))
#
#____________________________________


