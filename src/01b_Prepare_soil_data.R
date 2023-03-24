#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Prepare soil data                #
#          author: Romy Zeiss               #
#            date: 2023-03-23               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(usdm) # for colinearity check (VIF)

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v1.csv"))
data_glob

#data_cont <- read_csv()

#data_regi <- read_csv()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add missing environmental covariates ####

data_glob <- f_extract_env(data = data_glob,
                           col_lon = "Longitude_c", col_lat = "Latitude_c")
head(data_glob)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_glob <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

data_glob <- f_add_protect(data = data_glob, data_pa = protect_glob, col_id="Order_ID")
rm(protect_glob)

nrow(data_glob %>% filter(PA==1)) #81
nrow(data_glob %>% filter(!is.na(PA_type))) #136
nrow(data_glob %>% filter(PA==0)) #305

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

list_colinear <- f_colinearity(data = data_glob, 
                              col_lon = "Longitude_c", 
                              col_lat = "Latitude_c", 
                              vars_env = c("Soil_pH", "Soil_salinity", 
                                          "Clay_silt_c", "Elevation", 
                                          "AnnualTemp", "AnnualPrec"))

# Save
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/VIF_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatSpearman, paste0("./results/corMatSpearman_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatPearson, paste0("./results/corMatPearson_envVars_global.csv"), row.names = F)






