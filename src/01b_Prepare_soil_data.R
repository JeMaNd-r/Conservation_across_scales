#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Prepare soil data                #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(usdm) # for colinearity check (VIF)

source("src/00_Parameters_functions.R")

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

# remove double entries (i.e., 2+ protected area types for 1 PA)
protect_type <- data.frame("PA_type" = c("Ia", "Ib", "II", "III", "IV", "V", "VI",
                                         "Not Reported", "Not Assigned", "Not Applicable"),
                           "PA_rank" = 1:10,
                           "PA_protected" = c(rep(1,7), rep(0,3)))
protect_type

protect_glob <- protect_glob %>% full_join(protect_type, by=c("IUCN_CAT" = "PA_type"))

# take the highest = min (or lowest = max) protection level
protect_glob <- protect_glob %>% 
  group_by(id.y) %>% 
  dplyr::select(id.y, PA_rank) %>%
  summarize(across(everything(), list("min"=min, "max"=max), na.rm=T)) %>%
  mutate("PA_rank_min" = ifelse(PA_rank_min == Inf, NA, PA_rank_min),
         "PA_rank_max" = ifelse(PA_rank_max == -Inf, NA, PA_rank_max)) %>%
  full_join(protect_type %>% rename("PA_type_min" = PA_type,
                                    "PA_rank_min" = PA_rank,
                                    "PA_protected_min" = PA_protected), 
            by="PA_rank_min") %>%
  full_join(protect_type %>% rename("PA_type_max" = PA_type,
                                    "PA_rank_max" = PA_rank,
                                    "PA_protected_max" = PA_protected), 
            by="PA_rank_max")

# We decided to take the minimum (i.e., highest) level of protection because
# of sites such as ID 49 (Order_ID 302) that has min=1 (Ia) and max=8 (Not Reported)

# add column with information about protected and non-protected sites
data_glob <- data_glob %>% 
  mutate("id.x" = 1:nrow(data_glob)) %>%
  full_join(protect_glob %>% dplyr::select(id.y, PA_type_min, PA_protected_min, PA_rank_min),
            by=c("id.x"="id.y")) %>%
  rename("PA" = PA_protected_min,
         "PA_type" = PA_type_min,
         "PA_rank" = PA_rank_min) %>%
  mutate("PA" = ifelse(is.na(PA), 0, PA))

nrow(data_glob %>% filter(PA==1)) #81
nrow(data_glob %>% filter(!is.na(PA_type))) #136
nrow(data_glob %>% filter(PA==0)) #305

rm(protect_glob)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Keep complete cases only ####
data_glob <- data_glob[complete.cases(data_glob[,c(mahal_vars, fns)]),]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Rename land cover types in data ####
data_glob$LC <- data_glob$Eco_c
unique(data_glob$LC)
data_glob[data_glob$LC=="Forest", "LC"] <- "Woodland"
data_glob[data_glob$LC=="Moss_heath", "LC"] <- "Other"
unique(data_glob$LC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save ####
write_csv(data_glob, file=paste0(here::here(), "/intermediates/Data_global.csv"))


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Check for colinearity between environmental variables ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Calculate variable inflation factor (VIF)
# VIF is the extent of correlation between one predictor and all others.
# The lower VIF, the better we can tell what predictor contributed (most) to the model
env_glob <- as.data.frame(data_glob %>% 
                dplyr::select(Latitude_c, Longitude_c, Soil_pH, 
                              Soil_salinity, Clay_silt_c, Elevation,
                              AnnualTemp, AnnualPrec))

## VIF basen on raw data (explanatory raster stack)
env_vif <- usdm::vif(env_glob)

# which predictors should be excluded?
vif_cor <- usdm::vifcor(env_glob, th=0.8)  #th = threshold correlation for exclusion
# how: first find a pair of variables which has the maximum linear correlation 
# (greater than th), and exclude one of them which has greater VIF. The 
# procedure is repeated until no variable with a high correlation coefficient 
# (grater than threshold) with other variables remains.

vif_step <- usdm::vifstep(env_glob, th=10) #VIF >10

# merge both data.frames
env_vif <- env_vif %>% rename("VIF_raw" = VIF) %>% full_join(vif_step@results) %>%
  full_join(as.data.frame(vif_cor@corMatrix) %>% mutate("Variables"=rownames(vif_cor@corMatrix)))

env_vif

write.csv(env_vif, file=paste0(here::here(), "/results/VIF_envVars_global.csv"), row.names = F)

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate correlations between predictors
# https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses

corMatSpearman <- cor(env_glob, use="complete.obs", method="spearman") %>% round(2)
corMatPearson <- cor(env_glob, use="complete.obs", method="pearson") %>% round(2)

write.csv(corMatSpearman, paste0("./results/corMatSpearman_envVars_global.csv"), row.names = F)
write.csv(corMatPearson, paste0("./results/corMatPearson_envVars_global.csv"), row.names = F)

rm(env_vif, vif_cor, corMatPearson, corMatSpearman, env_glob, vif_step)




