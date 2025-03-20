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
## 1. GLOBAL ----------------------------------------- ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####
coord_glob <- readxl::read_xlsx(paste0(here::here(), "/data_raw/Global_Atlas_drylands_V2.xlsx"),
                               sheet = "Database")
coord_glob #551 
coord_glob %>% filter(!is.na(ID_sequencing_16S)) #338

data_glob <- read_csv(here::here("data_raw/Diversity_global.csv")) #342 

sort(unique(data_glob$SampleID)) #342
sort(unique(coord_glob$ID_sequencing_16S)) #338

data_glob <- coord_glob %>%
  dplyr::select("Plot_ID", "ID_sequencing_16S", "Latitude_c", "Longitude_c", 
                "Vegetation", "Soil_pH", "Soil_salinity", "Soil_fine_texture",
                "ORC", "BGL veg", "FOS veg", "MIN veg",
                #"WHC veg", 
                "Plant_cover",             
                "TON", "Soil_total_P", "AVP", "DIN", "DON") %>%
  mutate(ID_numb = str_replace(ID_sequencing_16S,"FM", "1")) %>% #ID_sequencing_16s goes to 249, FMxxx
  mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
  mutate(ID_numb = as.numeric(ID_numb)) %>%
  filter(!is.na(ID_numb)) %>%
  full_join(data_glob %>%
              mutate(ID_numb = str_replace(SampleID, "FM", "1")) %>%
              mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
              mutate(ID_numb = as.numeric(ID_numb)), 
            by = "ID_numb") %>%
  filter(!is.na(Latitude_c)) #missing coordinates are samples for other projects
data_glob #338
rm(coord_glob)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename & calculate variables (if needed) ####
data_glob <- data_glob %>%
  mutate("Pathogen_control" = 1 - (Plant_pathogen_richness / Fungi_richness))

data_glob <- data_glob %>%
  rename(Latitude = Latitude_c,
         Longitude = Longitude_c,
         Soil_texture = Soil_fine_texture)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add missing environmental covariates (elevation & annual climate) ####

data_glob <- f_extract_env(data = data_glob,
                           col_lon = "Longitude", col_lat = "Latitude")
head(data_glob)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_glob <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

data_glob <- f_add_protect(data = data_glob, 
                           data_pa = protect_glob %>%
                             mutate(ID_numb = str_replace(ID_sequencing_16S, "FM", "1")) %>%
                             mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
                             mutate(ID_numb = as.numeric(ID_numb)), 
                           col_id="ID_numb")
rm(protect_glob)

nrow(data_glob %>% filter(PA==1)) #54
nrow(data_glob %>% filter(!is.na(PA_type))) #106
nrow(data_glob %>% filter(PA==0)) #285

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename land cover types in data ####
# treat all global data as Drylands 
data_glob$LC <- "Dryland"
# data_glob$LC <- data_glob$Vegetation
# unique(data_glob$LC)
# data_glob[data_glob$LC=="Forest" & !is.na(data_glob$LC), "LC"] <- "Woodland"
unique(data_glob$LC)

table(data_glob$LC)

# rename ID
data_glob <- data_glob %>% mutate(SampleID = ID_numb)

# functions
data_glob <- data_glob %>%
  filter(!is.na(SampleID)) %>%
  rename(Soil_carbon_service = ORC,
         Soil_stability_service = Plant_cover) %>% #,
         #Water_regulation_service = `WHC veg`) %>%
  mutate(OM_decomposition_service = rowSums(cbind(scale(`BGL veg`)[,1], 
                                          scale(`FOS veg`)[,1],
                                          scale(`MIN veg`)[,1]),
                                        na.rm=TRUE)/3,
         Nutrient_service = (scale(TON)[,1]+
                               scale(Soil_total_P)[,1]+ 
                               scale(AVP)[,1]+
                               scale(DIN)[,1]+
                               scale(DON)[,1])/6)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Keep complete cases only ####
data_glob <- data_glob[complete.cases(data_glob[,c(mahal_vars, fns, "LC")]),] #248 (73%)   
data_glob

# subset relevant columns
data_glob <- data_glob %>% 
  dplyr::select(all_of(c("SampleID", "LC", mahal_vars, fns, "PA", "PA_type", "PA_rank")))
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save ####
write_csv(data_glob, file=paste0(here::here(), "/intermediates/Data_clean_global.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for colinearity between environmental variables ####

list_colinear <- f_colinearity(data = data_glob, 
                              col_lon = "Longitude", 
                              col_lat = "Latitude", 
                              vars_env = c("Soil_pH", "Soil_salinity", 
                                          "Soil_texture", "Elevation", 
                                          "AnnualTemp", "AnnualPrec"))
list_colinear
list_colinear$env_vif %>% filter(is.na(VIF))
#-> none excluded based on VIF

# Save
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/VIF_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatSpearman, paste0(here::here(), "/results/corMatSpearman_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatPearson, paste0(here::here(), "/results/corMatPearson_envVars_global.csv"), row.names = F)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 2. CONTINENTAL ------------------------------------ ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####
coord_cont <- read_csv(paste0(here::here(), "/data_raw/LUCAS_2018_iDiv_20221018.csv"))
coord_cont <- coord_cont %>% dplyr::filter(Soil_data_survey == "2018")
coord_cont #885 

data_cont <- read_csv(here::here("data_raw/Diversity_continental.csv")) #885 
data_cont

data_cont <- coord_cont %>% 
  dplyr::select("BARCODE_ID", "Longitude", "Latitude",  
                "LC_3",  "pH_H2O", "Electrical_conductivity", 
                "Clay_content", "Silt_content", # need to be summed up
                "Organic_carbon", 
                "N_actylglucosaminidase", "Basal_respiration", # "OM_decomposition_service",
                #"Water_regulation_service", #not available 
                "Mean_width_diameter", "Water_stable_aggregates",             
                "Phosphorus_content", "Total_nitrogen_content", "Extractable_potassium_content" #"Nutrient_service"
  ) %>%
  full_join(data_cont, by=c("BARCODE_ID" = "SampleID")) %>%
  arrange(BARCODE_ID) %>%
  rename("SampleID" = BARCODE_ID)
data_cont
rm(coord_cont)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename & calculate variables (if needed) ####
data_cont <- data_cont %>%
  mutate("Pathogen_control" = 1 - (Plant_pathogen_richness / Fungi_richness))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add missing environmental covariates (elevation & annual climate) ####

data_cont <- f_extract_env(data = data_cont,
                           col_lon = "Longitude", col_lat = "Latitude")
head(data_cont)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_cont <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_continental.csv"))

data_cont <- f_add_protect(data = data_cont, data_pa = protect_cont, col_id="SampleID")
rm(protect_cont)

nrow(data_cont %>% filter(PA==1)) #68
nrow(data_cont %>% filter(!is.na(PA_type))) #149
nrow(data_cont %>% filter(PA==0)) #818

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename columns in data ####
# land cover types 
data_cont$LC <- data_cont$LC_3
unique(data_cont$LC)
data_cont[data_cont$LC=="Bare land and lichens/moss" & !is.na(data_cont$LC), "LC"] <- "Other"
data_cont[data_cont$LC=="Wetlands" & !is.na(data_cont$LC), "LC"] <- "Other"
data_cont[data_cont$LC=="Artificial Land" & !is.na(data_cont$LC), "LC"] <- "Other"
data_cont[data_cont$LC=="Moss_heath" & !is.na(data_cont$LC), "LC"] <- "Other"
data_cont <- data_cont %>% filter(!is.na(LC))

table(data_cont$LC) # Shrubland and other only 30 and 34... probably to be removed

# functions
data_cont <- data_cont %>%
  rename(Soil_pH = pH_H2O,
         Soil_salinity = Electrical_conductivity,
         Soil_carbon_service = Organic_carbon) %>%
  mutate(Soil_texture = Clay_content + Silt_content,
         OM_decomposition_service = (scale(N_actylglucosaminidase)[,1]+ 
                                       scale(Basal_respiration)[,1])/2,
         Water_regulation_service = NA,
         Nutrient_service = (scale(Phosphorus_content)[,1]+
                               scale(Total_nitrogen_content)[,1]+ 
                               scale(Extractable_potassium_content)[,1])/3,
         Soil_stability_service = (scale(Mean_width_diameter)[,1]+
                                    scale(Water_stable_aggregates)[,1])/2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Keep complete cases only ####
data_cont <- data_cont[complete.cases(data_cont[,c(mahal_vars, fns[fns != "Water_regulation_service"], "LC")]),] #818 (all)
data_cont

# subset relevant columns
data_cont <- data_cont %>% 
  dplyr::select(all_of(c("SampleID", "LC", mahal_vars, fns, "PA", "PA_type", "PA_rank")))
data_cont

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save ####
write_csv(data_cont, file=paste0(here::here(), "/intermediates/Data_clean_continental.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for colinearity between environmental variables ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

list_colinear <- f_colinearity(data = data_cont,
                               col_lon = "Longitude",
                               col_lat = "Latitude",
                               vars_env = c("Soil_pH", "Soil_salinity",
                                            "Soil_texture", "Elevation",
                                            "AnnualTemp", "AnnualPrec"))
list_colinear
list_colinear$env_vif %>% filter(is.na(VIF))
#-> Latitude VIF_raw 36.6

# Save
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/VIF_envVars_continental.csv"), row.names=F)
write.csv(list_colinear$corMatSpearman, paste0("./results/corMatSpearman_envVars_continental.csv"), row.names=T)
write.csv(list_colinear$corMatPearson, paste0("./results/corMatPearson_envVars_continental.csv"), row.names=T)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 3. REGIONAL --------------------------------------- ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####
coord_regi <- read_csv(paste0(here::here(), "/data_raw/SoilReCon_Data_211123.csv"))
coord_regi #407 

data_regi <- read_csv(here::here("data_raw/Diversity_regional.csv")) #407 
data_regi

data_regi <- coord_regi %>% 
  dplyr::select(SampleID, Longitude, Latitude, Landuse, 
                pH_H2O,
                `Electr_conductivity uS/cm`,
                ClaySilt,
                Org_M, # no carbon available (Organic carbon)
                Water_stable_aggregates_mean, #soil aggregate stablity: no mean weight diameter? 
                Phosphorous, Potassium, Ca, Mg, # no N measured
                Microbial_Respiration, # respiration & N-Acetylglucosaminidase(?) = OM decomposition
                SWR, # Soil water repelency... binned (water regulation service)  
                GD_Elevation, GD_MAP, GD_MAT
                ) %>%
  full_join(data_regi,
            by = "SampleID") %>%
  filter(!is.na(SampleID))
data_regi #407
rm(coord_regi)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename & calculate variables (if needed) ####
data_regi <- data_regi %>%
  mutate("Pathogen_control" = 1 - (Plant_pathogen_richness / Fungi_richness))

data_regi <- data_regi %>%
  mutate("P_z" = scale(Phosphorous)[,1], "K_z" = scale(Potassium)[,1], 
         "Ca_z" = scale(Ca)[,1], "Mg_z" = scale(Mg)[,1],
         "Soil_stability_service" = scale(Water_stable_aggregates_mean)[,1], 
         "OM_decomposition_service" = scale(Microbial_Respiration)[,1]) %>%
  rename("Soil_texture" = ClaySilt,
         "Soil_pH" = pH_H2O,
         "Soil_carbon_service" = Org_M, 
         "Water_regulation_service" = SWR,
         "Soil_salinity" = 'Electr_conductivity uS/cm')
# calculate Nutrient_service separately (mutate and mean did not work)
data_regi$Nutrient_service <- rowMeans(data_regi[,c("P_z", "K_z", "Ca_z", "Mg_z")], na.rm = T)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add missing environmental covariates (elevation & annual climate) ####

data_regi <- f_extract_env(data = data_regi,
                           col_lon = "Longitude", col_lat = "Latitude")
data_regi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_regi <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_regional.csv"))

data_regi <- f_add_protect(data = data_regi, data_pa = protect_regi, col_id="SampleID")
rm(protect_regi)

nrow(data_regi %>% filter(PA==1)) #60
nrow(data_regi %>% filter(!is.na(PA_type))) #107
nrow(data_regi %>% filter(PA==0)) #348

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename land cover types in data ####

data_regi <- data_regi %>%
  mutate("LC" = ifelse(Landuse=="CROP", "Cropland",
                       ifelse(Landuse=="PAST", "Grassland",
                              ifelse(Landuse=="FOR", "Woodland",
                                     ifelse(Landuse=="EXO", "Woodland",
                                            ifelse(Landuse=="URB", "Other",
                                                   ifelse(Landuse=="PERM", "Cropland",
                                                          NA)))))))
unique(data_regi$LC)
table(data_regi$LC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Keep complete cases only ####
data_regi <- data_regi[complete.cases(data_regi[,c(mahal_vars, fns[fns != "Water_regulation_service"], "LC")]),] #354 (not all)
data_regi

# subset relevant columns
data_regi <- data_regi %>% 
  dplyr::select(all_of(c("SampleID", "LC", mahal_vars, fns, "PA", "PA_type", "PA_rank")))
data_regi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save ####
write_csv(data_regi, file=paste0(here::here(), "/intermediates/Data_clean_regional.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for colinearity between environmental variables ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

list_colinear <- f_colinearity(data = data_regi,
                               col_lon = "Longitude",
                               col_lat = "Latitude",
                               vars_env = c("Soil_pH", "Soil_salinity",
                                            "Soil_texture", "Elevation",
                                            "AnnualTemp", "AnnualPrec"))
list_colinear
list_colinear$env_vif %>% filter(is.na(VIF))
#-> Elevation VIF_raw 253.9

# Save
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/VIF_envVars_regional.csv"), row.names=F)
write.csv(list_colinear$corMatSpearman, paste0("./results/corMatSpearman_envVars_regional.csv"), row.names=T)
write.csv(list_colinear$corMatPearson, paste0("./results/corMatPearson_envVars_regional.csv"), row.names=T)

