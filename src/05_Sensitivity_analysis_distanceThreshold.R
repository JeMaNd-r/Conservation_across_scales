#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Sensitivity analysis             #
#  distance threshold for continental pairs #
#          author: Romy Zeiss               #
#            date: 2025-02-14               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(terra)

temp_scale <- "continental"

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]

# create directory for intermediate results
if(!dir.exists(paste0(here::here(), "/results/sensitivity_distanceThreshold"))){
  dir.create(paste0(here::here(), "/results/sensitivity_distanceThreshold"))
}

if(!dir.exists(paste0(here::here(), "/results/sensitivity_distanceThreshold/", temp_scale))){
  dir.create(paste0(here::here(), "/results/sensitivity_distanceThreshold/", temp_scale))
}

#### Function to check for pairing ####
# same as function in 00_Functions.R, except for 2 changes:
# 1. changed directory to sensitivity subfolder in results
# 2. file naming changed from Sys.Date() to radius_thres

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
  write.csv(unpaired_pa, file=paste0(here::here(), "/results/sensitivity_distanceThreshold/", temp_scale,  "/Unpaired_protected_sites_", radius_thres, ".csv"), row.names = F)
  
  cat("#---------------------------------------------------", sep="\n")
  cat(paste0("Saved: csv file with unpaired sites is saved under: "), 
      paste0(here::here(), "/intermediates/", temp_scale, "/Unpaired_protected_sites_", Sys.Date(), ".csv"),
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


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Loop over radius thresholds ####

for(radius_thres in c(500000, 200000, 300000, 400000, 600000, 700000, 1000000)){

  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Load data ####
  
  data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  data_clean
  
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
                               col_id = "SampleID", col_lc = "LC", 
                               vars_z = mahal_vars_z)
  all_nonPA <- count_nonPA[[2]]
  count_nonPA <- count_nonPA[[1]] # >= min_nonPA (e.g. 10) to build pairs and keep PA in dataset
  
  # sites that have lower number of possible pairing options than min_nonPA
  count_nonPA %>% filter(No_nonPA < min_nonPA)
  count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ### Remove sites ####
  unpaired_pa <- read.csv(here::here(paste0("results/sensitivity_distanceThreshold/", temp_scale, "/Unpaired_protected_sites_", radius_thres, ".csv")))
  head(unpaired_pa) 
  
  # Note: If there were nonPA sites with mahalanobis distance below threshold
  data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
  nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
  
  # Remove sites that can only be paired less than min_nonPA times
  data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA < min_nonPA, "SampleID"]),]
  nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
  data_clean %>% group_by(LC, PA) %>% count()
  
  # save data
  write_csv(data_clean, paste0(here::here(), "/results/sensitivity_distanceThreshold", temp_scale, "Data_paired_", radius_thres, ".csv"))
  data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_distanceThreshold", temp_scale, "Data_paired_", radius_thres, ".csv"))
  
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
  
  write_csv(pa_pairs, file=paste0(here::here(), "/results/sensitivity_distanceThreshold/", temp_scale, "/Pairs_paNonpa_1000trails_", radius_thres, ".csv"))
}


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Numbers ####

#...


