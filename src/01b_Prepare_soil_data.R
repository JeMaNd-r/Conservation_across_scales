#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Prepare soil data                #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v1.csv"))
data_glob

#data_cont <- read_csv()

#data_regi <- read_csv()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add missing environmental covariates ####

# Elevation
elev <- terra::rast(paste0(here::here(), "/data_raw/Elevation_Worldclim_Global.tif"))
data_glob$Elevation <- terra::extract(elev, data_glob[,c("Longitude_c", "Latitude_c")])[,2]

# Annual temperature
annual_temp <- terra::rast(paste0(here::here(), "/data_raw/CHELSA_bio1_1981-2010_V.2.1.tif"))
data_glob$AnnualTemp <- terra::extract(annual_temp, data_glob[,c("Longitude_c", "Latitude_c")])[,2]
data_glob$AnnualTemp <- data_glob$AnnualTemp / 10

# Annual precipitation
annual_prec <- terra::rast(paste0(here::here(), "/data_raw/CHELSA_bio12_1981-2010_V.2.1.tif"))
data_glob$AnnualPrec <- terra::extract(annual_prec, data_glob[,c("Longitude_c", "Latitude_c")])[,2]

head(data_glob)




#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_glob <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

# add column with information about protected and non-protected sites
data_glob <- data_glob %>% 
  mutate("PA"=protect_glob %>% dplyr::select(id.y, PA) %>% unique() %>% dplyr::select(PA) %>% unlist())
rm(protect_glob)

write_csv(data_glob, file=paste0(here::here(), "/intermediates/Data_global.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Rename land cover types in data ####
data_glob$LC <- data_glob$Eco_c
unique(data_glob$LC)
data_glob[data_glob$LC=="Forest", "LC"] <- "Woodland"
data_glob[data_glob$LC=="Moss_heath", "LC"] <- "Other"
unique(data_glob$LC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Check for colinearity between environmental variables ####





