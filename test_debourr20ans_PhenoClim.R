library(ggmap)
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)


# Concernant les données PhenoClim et leur analyse, certains éléments ont déjà été explorés (notamment par Marjorie Bison).
# Les analyses se sont principalement concentrées sur les Alpes, la phase de débourrement, et sur une première série de données, de 
# 2004 à 2016. Les phases de floraison et feuillaison restent à explorer, tout comme le cas des autres massifs (Pyrénées, Massif Central,
# Vosges...).
# 
# Par rapport aux analyses de débourrement dans les Alpes, il y a désormais 8 années supplémentaires de données, mais les premières 
# explorations de ces données supplémentaires tendent à montrer que les tendances d'avancement du débourrement ne se retrouvent plus 
# avec les 8 ans supplémentaires, peut-être à cause d'une grande variabilité interzone avec des outliers.
# Cela reste à explorer, avec plusieurs pistes : 
# - Piste 1 : élimiter les outliers (Vitasse et al 2017)
# - Piste 2 : contraindre la distribution en jouant sur les priors dans les méthodes bayésiennes
# - Piste 3 : intégrer des données de température comme variable explicative, plus directement que le paramètre 'altitude'


# Pour la piste 3... ----
#*---- Données climatiques réanalysées ----
# - Copernicus CERRA : données météo Europe, jusqu'à 2022, à une résoluton de 5.5 km 
#                     (https://climate.copernicus.eu/copernicus-regional-reanalysis-europe-cerra)
# - ERA 5-Land : données météo Monde, en temps réel j-5, à une résolution de 9 km
#                     (https://climate.copernicus.eu/climate-reanalysis)
# - CHELSA : données météo Europe, 1979-2013, à une résolution de 1 km
#                     (https://chelsa-climate.org/timeseries/)
# - SAFRAN - MeteoFrance :        ???
#
# => PB : il n'y a pas de données récentes avec une bonne résolution


