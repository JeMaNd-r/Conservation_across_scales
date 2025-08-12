#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Sensitivity analysis             #
#   more similar indicator attributes       #
#          author: Romy Zeiss               #
#            date: 2025-08-12               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# remove in cont.: "N_actylglucosaminidase" for OM decomposition, "Total_nitrogen_content" for Nutrient cycling, "Mean_width_diameter" for Soil stability
# remove in regi:  Ca and Mg for Nutrient cycling

library(here)
library(tidyverse)
library(terra)
library(psych)

# create directory for intermediate results/sensitivity_similarAttributes
if(!dir.exists(paste0(here::here(), "/results/sensitivity_similarAttributes/"))){
  dir.create(paste0(here::here(), "/results/sensitivity_similarAttributes/"))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare data ####

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
                "Basal_respiration", # "OM_decomposition_service",
                #"Water_regulation_service", #not available 
                "Water_stable_aggregates",             
                "Phosphorus_content", "Extractable_potassium_content" #"Nutrient_service"
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
         OM_decomposition_service = scale(Basal_respiration)[,1],
         Water_regulation_service = NA,
         Nutrient_service = (scale(Phosphorus_content)[,1]+
                               scale(Extractable_potassium_content)[,1])/2,
         Soil_stability_service = scale(Water_stable_aggregates)[,1])

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
write_csv(data_cont, file=paste0(here::here(), "/results/sensitivity_similarAttributes/Data_clean_continental.csv"))

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
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/sensitivity_similarAttributes/VIF_envVars_continental.csv"), row.names=F)
write.csv(list_colinear$corMatSpearman, paste0("./results/sensitivity_similarAttributes/corMatSpearman_envVars_continental.csv"), row.names=T)
write.csv(list_colinear$corMatPearson, paste0("./results/sensitivity_similarAttributes/corMatPearson_envVars_continental.csv"), row.names=T)

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
                Phosphorous, Potassium, # no N measured
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
         "Soil_stability_service" = scale(Water_stable_aggregates_mean)[,1], 
         "OM_decomposition_service" = scale(Microbial_Respiration)[,1]) %>%
  rename("Soil_texture" = ClaySilt,
         "Soil_pH" = pH_H2O,
         "Soil_carbon_service" = Org_M, 
         "Water_regulation_service" = SWR,
         "Soil_salinity" = 'Electr_conductivity uS/cm')
# calculate Nutrient_service separately (mutate and mean did not work)
data_regi$Nutrient_service <- rowMeans(data_regi[,c("P_z", "K_z")], na.rm = T)

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
write_csv(data_regi, file=paste0(here::here(), "/results/sensitivity_similarAttributes/Data_clean_regional.csv"))

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
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/sensitivity_similarAttributes/VIF_envVars_regional.csv"), row.names=F)
write.csv(list_colinear$corMatSpearman, paste0("./results/sensitivity_similarAttributes/corMatSpearman_envVars_regional.csv"), row.names=T)
write.csv(list_colinear$corMatPearson, paste0("./results/sensitivity_similarAttributes/corMatPearson_envVars_regional.csv"), row.names=T)



#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Pairing & effect sizes ####

#### Function to check for pairing ####
# same as function in 00_Functions.R, except for 2 changes:
# 1. changed directory to sensitivity subfolder in results/sensitivity_similarAttributes

f_check_pairs_dist <- function(data, col_id, col_lc, vars_z){
  
  # data:   dataframe containing following columns
  # col_id: name of column with site IDs
  # col_lc: column name with land-cover types
  # vars_z: vector with column names for mahalanobis variables
  
  nonpa <- data[data$PA==0,c(col_id, col_lc, "PA", vars_z)]
  pa <- data[data$PA==1,c(col_id, col_lc, "PA", vars_z)]
  pa_noPair <- c()
  all_nonPA <- nonpa[0,]
  
  for(i in 1:length(lc_names)){
    temp_nonPA  <- nonpa[nonpa[,col_lc]==lc_names[i],]
    temp_PA <- pa[pa[,col_lc]==lc_names[i],]
    
    sigma <- cov(temp_nonPA[,vars_z]) 
    for(j in 1:nrow(temp_PA)){
      mu <- as.numeric(temp_PA[j,vars_z])
      temp_nonPA[,as.character(temp_PA[j,col_id])] <- 
        mahalanobis(temp_nonPA[,vars_z], mu, sigma, tol=1e-30)
      #print(j)
    }
    # if continental data, remove nonPA with geographical distance >500km
    if(temp_scale == "continental"){
      for(k in as.character(pull(temp_PA[,col_id]))){
        temp_nonPA[!(as.character(pull(temp_nonPA[,col_id])) %in% list_dist[[k]]), k] <- NA
      }
    }
    
    min_mahal <- apply(X = temp_nonPA[, as.character(temp_PA %>% pull(col_id))], MARGIN = 2, FUN = min, na.rm = TRUE)
    pa_noPair <- rbind(pa_noPair, cbind(names(min_mahal[min_mahal>mahal_thres]), min_mahal[min_mahal>mahal_thres]))
    
    all_nonPA <- full_join(all_nonPA, temp_nonPA, by=c(col_id, col_lc, "PA", vars_z))
  }
  pa_noPair
  nrow(pa_noPair) #nrow=12
  
  unpaired_pa <- data[data[,col_id] %in% pa_noPair,] 
  write.csv(unpaired_pa, file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale,  "/Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)
  
  cat("#---------------------------------------------------", sep="\n")
  cat(paste0("Saved: csv file with unpaired sites is saved under: "), 
      paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Unpaired_protected_sites_", Sys.Date(), ".csv"),
      sep="\n")
  cat("#---------------------------------------------------", sep="\n")
  
  # look at Mahalanoubis distance values for each nonPA (Order_ID) and PA (columns)
  #all_nonPA
  
  # count how many "options" exist for one PA
  count_nonPA <- data.frame("SampleID"=NA, "No_nonPA"=NA)
  for(i in colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))]){
    temp_column <- all_nonPA[,i]
    
    count_nonPA <- rbind(count_nonPA,
                         c(i, nrow(all_nonPA[which(temp_column<=mahal_thres),col_id]) ))
    
  }
  
  count_nonPA <- count_nonPA %>% 
    mutate("No_nonPA"=as.numeric(No_nonPA)) %>% 
    arrange(No_nonPA)
  
  print("Number of nonPA sites per PA (col_ID)")
  print(head(count_nonPA))
  
  rm(temp_PA, temp_nonPA, unpaired_pa, pa_noPair, min_mahal, nonpa, pa)
  
  return(list(count_nonPA, all_nonPA)) #number of nonPA sites per PA (Order_ID)
}


df_summary <- data.frame()

for(temp_scale in c("continental", "regional")){
  
  source(paste0(here::here(), "/src/00_Parameters.R"))
  source(paste0(here::here(), "/src/00_Functions.R"))

    if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  }
  # define each land cover type
  lc_names_all <- c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "Other")

  if(!dir.exists(paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale))){
    dir.create(paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale))
  }
  
  
 #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Load data ####
  
  data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/Data_clean_", temp_scale, ".csv"))
  data_clean
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Exclude LC types if needed ####
  
  data_clean <- data_clean %>% filter(LC %in% lc_names)
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Scale variables for mahalanobis distance ####
  
  data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Calculate distance between all samples ####
  
  nonpa <- data_clean[data_clean$PA==0,c("SampleID", "LC", "PA", mahal_vars)]
  pa <- data_clean[data_clean$PA==1,c("SampleID", "LC", "PA", mahal_vars)]
  
  temp_dist <- terra::distance(terra::vect(pa %>% dplyr::select(Longitude, Latitude), 
                                           geom=c("Longitude", "Latitude"), 
                                           crs = "+proj=longlat +datum=WGS84"),
                               terra::vect(nonpa %>% dplyr::select(Longitude, Latitude), 
                                           geom=c("Longitude", "Latitude"), 
                                           crs = "+proj=longlat +datum=WGS84"),
                               unit = "m") %>%
    as_tibble()
  
  colnames(temp_dist) <- nonpa$SampleID
  temp_dist <- temp_dist %>% mutate(ID = pa$SampleID) %>% dplyr::select(ID, everything())
  temp_dist
  
  # extract nonPA with distance <= radius_thres
  list_dist <- lapply(1:nrow(temp_dist), function(i){
    colnames(temp_dist[,-1])[temp_dist[i, -1] <= radius_thres]
  })
  names(list_dist) <- temp_dist$ID
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Build pairs ####
  
  count_nonPA <- f_check_pairs_dist(data = data_clean, 
                                    col_id = "SampleID", 
                                    col_lc = "LC", 
                                    vars_z = mahal_vars_z)
  all_nonPA <- count_nonPA[[2]]
  count_nonPA <- count_nonPA[[1]] # >= min_nonPA (e.g. 10) to build pairs and keep PA in dataset
  
  # sites that have lower number of possible pairing options than min_nonPA
  count_nonPA %>% filter(No_nonPA < min_nonPA)
  count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Remove sites ####
  unpaired_pa <- read.csv(here::here(paste0("results/sensitivity_similarAttributes/", temp_scale, "/Unpaired_protected_sites_", Sys.Date(), ".csv")))
  head(unpaired_pa) 
  
  # Note: If there were nonPA sites with mahalanobis distance below threshold
  data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
  nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
  
  # Remove sites that can only be paired less than min_nonPA times
  data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA < min_nonPA, "SampleID"]),]
  nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
  data_clean %>% group_by(LC, PA) %>% count()
  
  # save data
  write_csv(data_clean, paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Data_paired_", temp_scale, ".csv"))
  data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Data_paired_", temp_scale, ".csv"))
  
  # check LC types and number of (protected) sites
  table(data_clean$LC, data_clean$PA)
  
  # based on number of sites per LC, exclude LC
  # e.g. global: only 3 unprotected on 0 protected for Other & onyl 1 PA for Woodland & 0 for Cropland
  # continental: no Shrubland protected, 28 unprotected; and 0 PA in Other
  # regional: PA only min. 7 -> decrease minimum size number to 7, 
  #           exclude Shrublands & Others to get it running (otherwise no complete pairing achieved) 
  min_size <- min(table(data_clean$LC, 
                        data_clean$PA)[table(data_clean$LC, 
                                             data_clean$PA)
                                       >0])
  min_size
  
  # The following function will print the number of times it successfully 
  # paired sites. It will show the same number multiple times if it didn't 
  # reach a successful pairing (i.e., min_size pairs per lc_names).
  list_pairs <- f_pairing(data = data_clean, 
                          col_id = "SampleID", col_lc = "LC",
                          vars_z = mahal_vars_z)
  
  
  # show total count of unpaired (and removed) PAs and compare with number of paired sites
  table(list_pairs$missing_pa[,2])  # can be larger than 0, 0 is perfect
  table(list_pairs$pa_pairs$nonPA) # counts should be lower or equal to number of runs (i.e. times)
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #### Check & Save pairing ####
  
  pa_pairs <- list_pairs$pa_pairs
  
  # check for runs that failed (i.e. count < number of PA sites), and remove the respective pairs
  # note: other result objects are not effected as they were overwritten
  nrow(pa_pairs)  # should be length(lc_names) * min_size * 1000
  pa_pairs <- pa_pairs %>% add_count(times.with.error) %>%
    filter(n==length(lc_names)*min_size)
  nrow(pa_pairs)  # exactly length(lc_names)* min_size * 1000
  
  write_csv(pa_pairs, file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Pairs_paNonpa_1000trails_", Sys.Date(), ".csv"))


  
  # load data as above
  data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/Data_clean_", temp_scale, ".csv"))
  data_clean
  
  data_clean <- data_clean %>% filter(LC %in% lc_names)
  
  data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
  
  # load pairs
  pa_pairs <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Pairs_paNonpa_1000trails_", Sys.Date(), ".csv"))
  
  # subset sites that were actually paired
  data_locations <- data_clean %>% 
    filter(SampleID %in% unique(pa_pairs$ID) | 
             SampleID %in% unique(pa_pairs$nonPA)) %>%
    dplyr::select(Longitude,Latitude,SampleID, PA, LC)
  data_locations #G: nrow=131, C: 316, R: 161
  write_csv(data_locations, file = paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale ,"/Locations_", temp_scale, ".csv"))
  
  # add information to output list
  df_summary <- rbind(df_summary, data.frame("scale" = temp_scale, 
                                             "min_size" = min_size,
                                             "sites_lower_min_nonPA" = pull(count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()),
                                             "n_pa" = nrow(data_locations %>% filter(PA==1)),
                                             "n_nonpa" = nrow(data_locations %>% filter(PA==0))))
  
  df_summary
}  
  
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Calculate d values & count number of significances ####
for(temp_scale in c("continental", "regional")){
    
  df_summary$n_d_sign <- NA
  df_summary$n_d_sign_pos <- NA
  df_summary$n_d_sign_neg <- NA
  
        
  # load data as above
  data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/Data_clean_", temp_scale, ".csv"))
  data_clean
  
  data_clean <- data_clean %>% filter(LC %in% lc_names)
  
  data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
  
  # load PA-nonPA pairs
  pa_pairs <- read_csv(file=here::here("results", "sensitivity_similarAttributes", temp_scale, paste0("Pairs_paNonpa_1000trails_", Sys.Date(), ".csv")))
  head(pa_pairs)


  # using Cohen's D
  try({
    d_list <- f_compare_pa_nonpa(data = data_clean,
                                 data_pairs = pa_pairs,
                                 col_id = "SampleID",
                                 col_fns = fns)
    head(d_list)
    
    # save total df with effect sizes 
    save(d_list,  file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/d_1000_trails.RData"))
    
    # summarize
    d_df <- do.call(rbind, d_list)
    str(d_df)
    
    d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
      mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
             "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
    
    # save data for plot
    write.csv(d_df, file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Data_d-value.csv"), row.names = FALSE)
    
    d_df <- read_csv(file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Data_d-value.csv"))
    
    d_summary <- d_df %>% 
      dplyr::select(-run) %>%
      #pivot_longer(cols = c(p_value, ci_lower, ci_upper, t_stats),
      #             names_to = "metric") %>%
      group_by(lc, fns, Label, Group_function, Organism) %>%
      summarize(across(effect, .fns = list("mean"=mean, "SD"=sd, "median"=median,
                                           "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
                                           "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
                                           "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
                                           "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE))))
    d_summary
    write.csv(d_summary, file=paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Results_d-value_summary.csv"), row.names = FALSE)
    
    # significance
    d_summary <- d_summary %>%
      mutate(effect_direction = as.factor(sign(effect_mean)))%>%
      mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                         ifelse(effect_direction=="1", "positive", "0"))) %>%
      mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
             effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
      mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
    
    # add to summary table
    df_summary[df_summary$scale==temp_scale,]$n_d_sign <- d_summary %>% 
      filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% nrow()
    df_summary[df_summary$scale==temp_scale,]$n_d_sign_pos <- d_summary %>% 
      filter(effect_significance=="positive" & !is.na(effect_significance)) %>% nrow()
    df_summary[df_summary$scale==temp_scale,]$n_d_sign_neg <- d_summary %>% 
      filter(effect_significance=="negative" & !is.na(effect_significance)) %>% nrow()
  })
  
  df_summary
  
  write_csv(df_summary, file = paste0(here::here(), "/results/sensitivity_similarAttributes/Summary_allRadius_", temp_scale, ".csv"))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## check differences in significances #### 
d_all <- list()
for(temp_scale in c("continental", "regional")){
  d_summary <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/", temp_scale, "/Results_d-value_summary.csv"))
  
  # significance
  d_summary <- d_summary %>%
    mutate(effect_direction = as.factor(sign(effect_mean)))%>%
    mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                       ifelse(effect_direction=="1", "positive", "0"))) %>%
    mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
           effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
    mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
  
  # add radius_thres
  d_summary$scale <- temp_scale
  
  # add to overall table
  d_all <- c(d_all, list(d_summary))
}  

d_all <- do.call(rbind, d_all)

d_significant <- d_all %>% 
  unique() %>%
  filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% 
  count(scale, fns, effect_significance, lc) %>%
  pivot_wider(names_from = scale, values_from = n) 
d_significant

write_csv(d_significant, paste0(here::here(), "/results/sensitivity_similarAttributes/Results_d-value_summary_allScales.csv"))



# effect sizes
d_significant <- read_csv(paste0(here::here(), "/results/sensitivity_similarAttributes/Results_d-value_summary_allScales.csv"))
d_significant

# compare to main analysis





