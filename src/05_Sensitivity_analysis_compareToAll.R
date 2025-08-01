#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Sensitivity analysis             #
#        comparison of all sites            #
#          author: Romy Zeiss               #
#            date: 2025-07-31               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

options(scipen = 999) # to keep numbers as they are

# Script to compare all protected areas to all unprotected areas based on 
# no (very large) and picked threshold (500km), as well as with and without
# selecting sites with minimum Mahalanobis distance

library(here)
library(tidyverse)
library(terra)
library(psych)

#### Function to check for pairing ####
# same as function in 00_Functions.R, except for 2 changes:
# 1. changed directory to sensitivity subfolder in results
# 2. file naming changed from Sys.Date() to radius_thres

#temp_scale <- "continental"
for(temp_scale in c("continental", "regional", "global")){  

  source(paste0(here::here(), "/src/00_Parameters.R"))
  source(paste0(here::here(), "/src/00_Functions.R"))
  rm(mahal_thres, radius_thres)
  
  
  f_check_pairs_dist <- function(data, col_id, col_lc, vars_z, mahal_thres){
    
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
    write.csv(unpaired_pa, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale,  "/Unpaired_protected_sites_r", radius_thres, "-m", mahal_thres, ".csv"), row.names = F)
    
    cat("#---------------------------------------------------", sep="\n")
    cat(paste0("Saved: csv file with unpaired sites is saved under: "), 
        paste0(here::here(), "/intermediates/", temp_scale, "/Unpaired_protected_sites_", Sys.Date(), "_r", radius_thres, "-m", mahal_thres, ".csv"),
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
  
  
  f_pairing <- function(data, col_id, col_lc, vars_z, mahal_thres){
    
    # data:         dataframe containing following columns
    # col_id:       name of column with site IDs
    # col_lc:       column name with land-cover types
    # vars_z:       vector with column names for mahalanobis variables
    
    # split data
    nonpa <- data[data$PA==0,c(col_id, col_lc, "PA", vars_z)]
    pa <- data[data$PA==1,c(col_id, col_lc, "PA", vars_z)]
    
    pa_pairs <- data.frame(ID=NULL, nonPA=NULL, mahal.min=NULL, LC=NULL, run=NULL)
    missing_pa <- data.frame(run=NA, pa.site=NA)[0,]
    times <- 0; times.with.error <- 0; set.seed(1) 
    
    repeat {
      times <- times+1; if(times > number_times) {break} #stop loop if reached 1000 trails
      times.with.error <- times.with.error + 1
      
      stop_loop <- FALSE #to stop loop if needed (see further below)
      
      # add columns to store pairing temporary
      pa[,c("nonPA", "mahal.min")] <- NA  
      
      for(i in 1:length(lc_names)){
        
        temp_PA <- pa[pa[,col_lc]==lc_names[i],] 
        temp_nonPA  <- nonpa[nonpa[,col_lc]==lc_names[i],]
        
        # select environmental data only as matrix, remove rows with NAs
        temp_PA <- temp_PA[complete.cases(temp_PA[,vars_z]),c(col_id, vars_z)]
        temp_nonPA <- temp_nonPA[complete.cases(temp_nonPA[,vars_z]),c(col_id, vars_z)]
        
        sigma <- cov(temp_nonPA[,vars_z]) # covariance/correlation between variables
        
        # random order of PA sites
        temp_PA <- temp_PA[order(f_resample(temp_PA %>% pull(all_of(col_id)), size = nrow(temp_PA))),]
        
        # add empty columns
        temp_nonPA[,as.character(temp_PA %>% pull(all_of(col_id)))] <- NA
        
        # calculate Mahalanobis distance
        for(j in 1:nrow(temp_PA)){
          mu = as.numeric(temp_PA[j,vars_z])
          temp_nonPA[,as.character(temp_PA[j,col_id])] <- 
            mahalanobis(temp_nonPA[,vars_z], mu, sigma, tol=1e-30)
          #print(j)
        }
        
        # if continental data, remove nonPA with geographical distance >radius_thres
        if(temp_scale == "continental"){
          for(k in as.character(pull(temp_PA[,col_id]))){
            temp_nonPA[!(as.character(pull(temp_nonPA[,col_id])) %in% list_dist[[k]]), k] <- NA
          }
        }
        
        # add column to PA data with respective nonPA sites 
        # based on minimal mahalanobis distance
        temp_PA[, c("nonPA", "mahal.min")]  <- NA
        #temp.col <- 0
        
        for(k in f_resample(temp_PA %>% pull(all_of(col_id)), size = nrow(temp_PA))){
          
          # if Mahal. distances compared to all relevant nonPA are above threshold...
          if(min(temp_nonPA[,as.character(k)], na.rm = TRUE)>=mahal_thres){
            #... stop and re-do run
            missing_pa <- rbind(missing_pa,c(times, k))  # to know why we've stopped
            print(paste0("Not all PA sites paired. Check PA ", k, "."))
            
            pa[pa %>% pull(all_of(col_id))==k, c("nonPA", "mahal.min")] <-
              cbind(NA, min(temp_nonPA[,as.character(k)]))
            temp_PA <- temp_PA[temp_PA %>% pull(all_of(col_id))!=k,] # remove respective PA site
            
          }else{
            # select (max. 10) nonPA sites with mahalanobis distance below threshold
            nonPA.pair <- temp_nonPA[temp_nonPA[,as.character(k)]<=mahal_thres & !is.na(temp_nonPA[,as.character(k)]),c(col_id, as.character(k))]
            nonPA.pair[nrow(nonPA.pair)+1:10,] <- NA #add empty rows for when there are less than 10 sites
            nonPA.pair <- nonPA.pair %>% arrange(2)
            nonPA.pair <- nonPA.pair[c(1:10),1]
            
            # sample one of the top 10 nonPA site that isn't NA 
            temp_PA[temp_PA %>% pull(all_of(col_id))==k,"nonPA"] <- f_resample(as.numeric(nonPA.pair[!is.na(nonPA.pair)]),1)
            
            # add value of distance
            temp_PA[temp_PA %>% pull(all_of(col_id))==k,"mahal.min"] <- unique(temp_nonPA[temp_nonPA %>% pull(all_of(col_id))==as.numeric(temp_PA[temp_PA %>% pull(all_of(col_id))==k,"nonPA"]),as.character(k)])
            temp_nonPA <- temp_nonPA[temp_nonPA %>% pull(all_of(col_id))!=as.character(temp_PA[temp_PA %>% pull(all_of(col_id))==k,"nonPA"]),]
          }
        }
        # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
        if(nrow(temp_PA) < min_size) {
          times <- times-1
          print("Not enough PA sites paired.")
          
          if(times < -1){ # if running into negative numbers with times
            print("No pairing possible. Try decreasing min_size, radius_thresh or col_lc (data grouping).")
            if(exists("pairing_possible")){ # relevant for sensitivity analysis of distance threshold
              # add information to output list
              pairing_possible <- c(pairing_possible, list(cbind("radius_thres" = radius_thres, 
                                                                 "min_size" = min_size,
                                                                 "sites_lower_min_nonPA" = pull(count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()),
                                                                 "pairing_possible" = 0)))
            }
            stop_loop <- TRUE
            break() #stop looping into negative numbers
          }
          
        } else {
          temp_PA <- temp_PA %>% sample_n(size = min_size, replace = FALSE)
          
          # add to result table
          pa[(pa %>% pull(all_of(col_id))) %in% (temp_PA %>% pull(all_of(col_id))), c("nonPA", "mahal.min")] <-
            cbind(temp_PA[order(as.numeric(rownames(temp_PA))),c("nonPA", "mahal.min")])
          
          # add to result table to analyse all at once below
          pa_pairs <- rbind(pa_pairs, cbind(temp_PA %>% dplyr::select(all_of(c(col_id, "nonPA", "mahal.min"))) %>% rename(ID=col_id),
                                            "LC"=lc_names[i], times.with.error, times))
          
        }
        print(times)
      }
      
      if(stop_loop) break()
    }
    
    rm(temp_PA, temp_nonPA, pa, nonpa, nonPA.pair)
    return(list(missing_pa = missing_pa, pa_pairs = pa_pairs))
    
  }
  
  # define combinations to be tested
  mahal_thres_vect <- c(qchisq(.975, df=length(mahal_vars)), qchisq(1, df=length(mahal_vars))) #17.53455, Inf
  radius_thres_vect <- c(500000, Inf)
  
  if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
    #min_size <- 39 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
    #min_size <- 14 # number of samples/ sites that should be paired per LC type
    fns <- fns[fns != "Water_regulation_service"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
    #min_size <- 7 # number of samples/ sites that should be paired per LC type
    fns <- fns[fns != "Water_regulation_service"]
  }
  
  # create directory for intermediate results
  if(!dir.exists(paste0(here::here(), "/results/sensitivity_compareToAll"))){
    dir.create(paste0(here::here(), "/results/sensitivity_compareToAll"))
  }
  
  if(!dir.exists(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale))){
    dir.create(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale))
  }
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Loop over radius thresholds ####
  
  pairing_possible <- list()
  
  for(radius_thres in radius_thres_vect){
    for(mahal_thres in mahal_thres_vect){
      print("#####################################################")
      print(radius_thres)
      print(mahal_thres)
      
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ### Load data ####
      
      data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
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
                                        col_id = "SampleID", col_lc = "LC", 
                                        vars_z = mahal_vars_z,
                                        mahal_thres = mahal_thres)
      all_nonPA <- count_nonPA[[2]]
      count_nonPA <- count_nonPA[[1]] # >= min_nonPA (e.g. 10) to build pairs and keep PA in dataset
      
      # sites that have lower number of possible pairing options than min_nonPA
      count_nonPA %>% filter(No_nonPA < min_nonPA)
      count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()
      
      #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ### Remove sites ####
      unpaired_pa <- read.csv(here::here(paste0("results/sensitivity_compareToAll/", temp_scale, "/Unpaired_protected_sites_r", radius_thres, "-m", mahal_thres,  ".csv")))
      head(unpaired_pa) 
      
      # Note: If there were nonPA sites with mahalanobis distance below threshold
      data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
      nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
      
      # Remove sites that can only be paired less than min_nonPA times
      data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA < min_nonPA, "SampleID"]),]
      nrow(data_clean); nrow(data_clean[data_clean$PA,]) 
      data_clean %>% group_by(LC, PA) %>% count()
      
      # save data
      write_csv(data_clean, paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_paired_r", radius_thres, "-m", mahal_thres,  ".csv"))
      data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_paired_r", radius_thres, "-m", mahal_thres,  ".csv"))
      
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
                              vars_z = mahal_vars_z,
                              mahal_thres = mahal_thres)
      
      
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
      
      write_csv(pa_pairs, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Pairs_paNonpa_1000trails_r", radius_thres, "-m", mahal_thres,  ".csv"))
      
      # add information to output list
      if(!(radius_thres %in% pairing_possible$radius_thres)){
        pairing_possible <- c(pairing_possible, list(cbind("radius_thres" = radius_thres, 
                                                           "mahal_thres" = mahal_thres,
                                                           "min_size" = min_size,
                                                           "sites_lower_min_nonPA" = pull(count_nonPA %>% filter(No_nonPA < min_nonPA) %>% count()),
                                                           "pairing_possible" = 1)))
    }}
  }
  
  pairing_possible <- do.call(rbind, pairing_possible) %>% 
    as_tibble() %>% 
    arrange(radius_thres, mahal_thres)
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Numbers ####
  
  pairing_possible
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Selected samples ####
  
  pairing_possible$n_pa <- NA
  pairing_possible$n_nonpa <- NA
  
  for(radius_thres in radius_thres_vect){
    for(mahal_thres in mahal_thres_vect){
      
      # load data as above
      data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
      data_clean
      
      data_clean <- data_clean %>% filter(LC %in% lc_names)
      
      data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
      
      # load pairs
      pa_pairs <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Pairs_paNonpa_1000trails_r", radius_thres, "-m", mahal_thres,  ".csv"))
      
      # subset sites that were actually paired
      data_locations <- data_clean %>% 
        filter(SampleID %in% unique(pa_pairs$ID) | 
                 SampleID %in% unique(pa_pairs$nonPA)) %>%
        dplyr::select(Longitude,Latitude,SampleID, PA, LC)
      data_locations #G: nrow=131, C: 316, R: 161
      write_csv(data_locations, file = paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale ,"/Locations_r", radius_thres, "-m", mahal_thres,  ".csv"))
      pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_pa <- nrow(data_locations %>% filter(PA==1)) #G: 28 PAs, G-together: 39; C: 48, R: 36
      pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_nonpa <- nrow(data_locations %>% filter(PA==0)) #G: 93 PAs, G-together: 92; C: 269, R: 125
      
  }}
  
  pairing_possible
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Calculate d values & count number of significances ####
  
  pairing_possible$n_d_sign <- NA
  pairing_possible$n_d_sign_pos <- NA
  pairing_possible$n_d_sign_neg <- NA
  
  for(radius_thres in radius_thres_vect){
    for(mahal_thres in mahal_thres_vect){
      
      if(pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,"n_pa"] == 0) { 
        print("No pairing, no effect sizes.")
      } else {
        
        # load data as above
        data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
        data_clean
        
        data_clean <- data_clean %>% filter(LC %in% lc_names)
        
        data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
        
        # load PA-nonPA pairs
        pa_pairs <- read_csv(file=here::here("results", "sensitivity_compareToAll", temp_scale, paste0("Pairs_paNonpa_1000trails_r", radius_thres, "-m", mahal_thres,  ".csv")))
        head(pa_pairs)
        
        # using Cohen's D
        try({
          d_list <- f_compare_pa_nonpa(data = data_clean,
                                       data_pairs = pa_pairs,
                                       col_id = "SampleID",
                                       col_fns = fns)
          head(d_list)
          
          # save total df with effect sizes 
          save(d_list,  file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/d_1000_trails_r", radius_thres, "-m", mahal_thres,  ".RData"))
          
          # summarize
          d_df <- do.call(rbind, d_list)
          str(d_df)
          
          d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
            mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
                   "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
          
          # save data for plot
          write.csv(d_df, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_d-value_r", radius_thres, "-m", mahal_thres,  ".csv"), row.names = FALSE)
          
          d_df <- read_csv(file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_d-value_r", radius_thres, "-m", mahal_thres,  ".csv"))
          
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
          write.csv(d_summary, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Results_d-value_summary_r", radius_thres, "-m", mahal_thres,  ".csv"), row.names = FALSE)
          
          # significance
          d_summary <- d_summary %>%
            mutate(effect_direction = as.factor(sign(effect_mean)))%>%
            mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                               ifelse(effect_direction=="1", "positive", "0"))) %>%
            mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
                   effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
            mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
          
          # add to summary table
          pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign <- d_summary %>% 
            filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% nrow()
          pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign_pos <- d_summary %>% 
            filter(effect_significance=="positive" & !is.na(effect_significance)) %>% nrow()
          pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign_neg <- d_summary %>% 
            filter(effect_significance=="negative" & !is.na(effect_significance)) %>% nrow()
        })
    }
  }}
  
  pairing_possible
  
  write_csv(pairing_possible, file = paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_", temp_scale, ".csv"))
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## mean distance of pairs & all samples #### 
  
  # nonpa <- data_clean[data_clean$PA==0,c("SampleID", "LC", "PA", mahal_vars)]
  # pa <- data_clean[data_clean$PA==1,c("SampleID", "LC", "PA", mahal_vars)]
  # 
  # temp_dist <- terra::distance(terra::vect(pa %>% dplyr::select(Longitude, Latitude), 
  #                                          geom=c("Longitude", "Latitude"), 
  #                                          crs = "+proj=longlat +datum=WGS84"),
  #                              terra::vect(nonpa %>% dplyr::select(Longitude, Latitude), 
  #                                          geom=c("Longitude", "Latitude"), 
  #                                          crs = "+proj=longlat +datum=WGS84"),
  #                              unit = "m") %>%
  #   as_tibble()
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## check differences in significances #### 
  
  d_all <- list()
  
  for(radius_thres in radius_thres_vect){
    for(mahal_thres in mahal_thres_vect){ try({
      d_summary <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Results_d-value_summary_r", radius_thres, "-m", mahal_thres,  ".csv"))
      
      # significance
      d_summary <- d_summary %>%
        mutate(effect_direction = as.factor(sign(effect_mean)))%>%
        mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                           ifelse(effect_direction=="1", "positive", "0"))) %>%
        mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
               effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
        mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
      
      # add radius_thres
      d_summary$radius_thres <- radius_thres
      d_summary$mahal_thres <- mahal_thres
      
      # add to overall table
      d_all <- c(d_all, list(d_summary))
  })}}
  
  d_all <- do.call(rbind, d_all)
  
  d_significant <- d_all %>% 
    unique() %>%
    filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% 
    count(radius_thres, mahal_thres, fns, effect_significance, lc) %>%
    pivot_wider(names_from = c(radius_thres, mahal_thres), values_from = n) 
  d_significant
  
  write_csv(d_significant, paste0(here::here(), "/results/sensitivity_compareToAll/Results_d-value_summary_allRadius_", temp_scale,".csv"))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pairing without any Mahalanobis selection #### 
# Note that with Mahalnobis threshold == Inf, only few unprotected sites are "used" for pairing

#temp_scale <- "continental"
for(temp_scale in c("continental", "regional", "global")){  
  
  source(paste0(here::here(), "/src/00_Parameters.R"))
  source(paste0(here::here(), "/src/00_Functions.R"))
  rm(mahal_thres, radius_thres)
  mahal_thres <- "NA"
  
  # define combinations to be tested
  radius_thres_vect <- c(500000, Inf)
  
  if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
    #min_size <- 39 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
    #min_size <- 14 # number of samples/ sites that should be paired per LC type
    fns <- fns[fns != "Water_regulation_service"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
    #min_size <- 7 # number of samples/ sites that should be paired per LC type
    fns <- fns[fns != "Water_regulation_service"]
  }
  
  # create directory for intermediate results
  if(!dir.exists(paste0(here::here(), "/results/sensitivity_compareToAll"))){
    dir.create(paste0(here::here(), "/results/sensitivity_compareToAll"))
  }
  
  if(!dir.exists(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale))){
    dir.create(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale))
  }
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Loop over radius thresholds ####
  
  pairing_possible <- list()
  
  for(radius_thres in radius_thres_vect){
    pairing_possible <- c(pairing_possible, list(cbind("radius_thres" = radius_thres, 
                                                       "mahal_thres" = mahal_thres,
                                                       "min_size" = NA,
                                                       "sites_lower_min_nonPA" = NA,
                                                       "pairing_possible" = 1)))
  }
  
  pairing_possible <- do.call(rbind, pairing_possible) %>% 
    as_tibble() %>% 
    arrange(radius_thres, mahal_thres)
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Numbers ####
  
  pairing_possible
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Calculate d values & count number of significances ####
  
  pairing_possible$n_pa <- NA
  pairing_possible$n_nonpa <- NA
  
  pairing_possible$n_d_sign <- NA
  pairing_possible$n_d_sign_pos <- NA
  pairing_possible$n_d_sign_neg <- NA
  
  for(radius_thres in radius_thres_vect){

      # load data as above
      data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
      data_clean
      
      data_clean <- data_clean %>% filter(LC %in% lc_names)
      
      data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)
    
      pa_pairs <- data.frame()
      temp_pairs$nonPA <- NA
      
      for(times in 1:number_times){
        print(times)
        for(lc in lc_names){
          temp_pairs <- data_clean %>%
            filter(PA == 1) %>%
            dplyr::select(SampleID, LC) %>%
            filter(LC == lc) %>%
            rename("ID" = SampleID)
  
          # if continental data, remove nonPA with geographical distance >radius_thres
          temp_nonPA <- data_clean %>% filter(PA == 0) %>% 
            filter(LC == lc) %>%
            pull(SampleID)
          
          if(temp_scale == "continental"){
            temp_nonPA <- temp_nonPA[!(temp_nonPA %in% unlist(list_dist[temp_PA$SampleID]))]
          }
          
          temp_pairs$nonPA <- temp_nonPA %>% f_resample(nrow(temp_pairs), replace = FALSE)
          
          temp_pairs <- temp_pairs %>% 
            mutate(mahal.min = 1,
                   times.with.error = 0,
                   times = times,
                   n = 0) %>%
            dplyr::select(ID, nonPA, mahal.min, LC, everything())
          pa_pairs <- rbind(pa_pairs, temp_pairs)
        }
      }
      
      # load PA-nonPA pairs
      write_csv(pa_pairs, file=here::here("results", "sensitivity_compareToAll", temp_scale, paste0("Pairs_paNonpa_1000trails_r", radius_thres, "-m", mahal_thres,  ".csv")))
      head(pa_pairs)
      
      # subset sites that were actually paired
      data_locations <- data_clean %>% 
        filter(SampleID %in% unique(pa_pairs$ID) | 
                 SampleID %in% unique(pa_pairs$nonPA)) %>%
        dplyr::select(Longitude,Latitude,SampleID, PA, LC)
      data_locations 
      
      write_csv(data_locations, file = paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale ,"/Locations_r", radius_thres, "-m", mahal_thres,  ".csv"))
      pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_pa <- nrow(data_locations %>% filter(PA==1)) #G: 28 PAs, G-together: 39; C: 48, R: 36
      pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_nonpa <- nrow(data_locations %>% filter(PA==0)) #G: 93 PAs, G-together: 92; C: 269, R: 125
      
      # using Cohen's D
      try({
        d_list <- f_compare_pa_nonpa(data = data_clean,
                                     data_pairs = pa_pairs,
                                     col_id = "SampleID",
                                     col_fns = fns)
        head(d_list)
        
        # save total df with effect sizes 
        save(d_list,  file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/d_1000_trails_r", radius_thres, "-m", mahal_thres,  ".RData"))
        
        # summarize
        d_df <- do.call(rbind, d_list)
        str(d_df)
        
        d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
          mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
                 "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
        
        # save data for plot
        write.csv(d_df, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_d-value_r", radius_thres, "-m", mahal_thres,  ".csv"), row.names = FALSE)
        
        d_df <- read_csv(file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Data_d-value_r", radius_thres, "-m", mahal_thres,  ".csv"))
        
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
        write.csv(d_summary, file=paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Results_d-value_summary_r", radius_thres, "-m", mahal_thres,  ".csv"), row.names = FALSE)
        
        # significance
        d_summary <- d_summary %>%
          mutate(effect_direction = as.factor(sign(effect_mean)))%>%
          mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                             ifelse(effect_direction=="1", "positive", "0"))) %>%
          mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
                 effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
          mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
        
        # add to summary table
        pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign <- d_summary %>% 
          filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% nrow()
        pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign_pos <- d_summary %>% 
          filter(effect_significance=="positive" & !is.na(effect_significance)) %>% nrow()
        pairing_possible[pairing_possible$radius_thres==radius_thres & pairing_possible$mahal_thres==mahal_thres,]$n_d_sign_neg <- d_summary %>% 
          filter(effect_significance=="negative" & !is.na(effect_significance)) %>% nrow()
      })
    }
  
  pairing_possible
  
  write_csv(pairing_possible, file = paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_", temp_scale, "_mNA.csv"))
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## mean distance of pairs & all samples #### 
  
  # nonpa <- data_clean[data_clean$PA==0,c("SampleID", "LC", "PA", mahal_vars)]
  # pa <- data_clean[data_clean$PA==1,c("SampleID", "LC", "PA", mahal_vars)]
  # 
  # temp_dist <- terra::distance(terra::vect(pa %>% dplyr::select(Longitude, Latitude), 
  #                                          geom=c("Longitude", "Latitude"), 
  #                                          crs = "+proj=longlat +datum=WGS84"),
  #                              terra::vect(nonpa %>% dplyr::select(Longitude, Latitude), 
  #                                          geom=c("Longitude", "Latitude"), 
  #                                          crs = "+proj=longlat +datum=WGS84"),
  #                              unit = "m") %>%
  #   as_tibble()
  
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## check differences in significances #### 
  
  d_all <- list()
  
  for(radius_thres in radius_thres_vect){try({
      d_summary <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Results_d-value_summary_r", radius_thres, "-m", mahal_thres,  ".csv"))
      
      # significance
      d_summary <- d_summary %>%
        mutate(effect_direction = as.factor(sign(effect_mean)))%>%
        mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                           ifelse(effect_direction=="1", "positive", "0"))) %>%
        mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
               effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
        mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))
      
      # add radius_thres
      d_summary$radius_thres <- radius_thres
      d_summary$mahal_thres <- mahal_thres
      
      # add to overall table
      d_all <- c(d_all, list(d_summary))
    })}
  
  d_all <- do.call(rbind, d_all)
  
  d_significant <- d_all %>% 
    unique() %>%
    filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% 
    count(radius_thres, mahal_thres, fns, effect_significance, lc) %>%
    pivot_wider(names_from = c(radius_thres, mahal_thres), values_from = n) 
  d_significant
  
  write_csv(d_significant, paste0(here::here(), "/results/sensitivity_compareToAll/Results_d-value_summary_allRadius_", temp_scale,"_mNA.csv"))
}


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## some numbers #### 

radius_thres <- 500000
mahal_thres <- Inf

# number of PA sites per habitat 
temp_data <- read_csv(file = paste0(here::here(), "/results/sensitivity_compareToAll/", temp_scale, "/Locations_r", radius_thres, "-m", mahal_thres, ".csv"))
temp_data %>% count(LC, PA)

# effect sizes
d_significant <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/Results_d-value_summary_allRadius_continental.csv"))
d_significant

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Summary across scales #### 

d_cont <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_continental.csv"))
d_regi <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_regional.csv"))
d_glob <- read_csv(paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_global.csv"))

d_all <- d_cont %>% mutate(scale = "continental") %>%
  full_join(d_regi %>% mutate(scale = "regional") %>%
              filter(radius_thres != 500000)) %>%
  full_join(d_glob %>% mutate(scale = "global") %>%
              filter(radius_thres != 500000)) %>%
  dplyr::select(scale, radius_thres, mahal_thres, everything())
d_all

write_csv(d_all %>%
            dplyr::select(-pairing_possible) %>%
            mutate(mahal_thres=round(mahal_thres,3)) %>%
            rename("Spatial scale"= scale, "Spatial distance threshold [km]"=radius_thres,
                   "Mahalanobis distance threshold" = mahal_thres,
                   "Number of paired protected sites per habitat type" = min_size,
                   "Excluded protected sites"=sites_lower_min_nonPA,
                   "Number of protected sites" = n_pa,
                   "Number of unprotected sites" = n_nonpa,
                   "Number of significant effect sizes"=n_d_sign,
                   "Number of significant positive effect sizes"=n_d_sign_pos,
                   "Number of significant negative effect sizes"=n_d_sign_neg,
            ),
          paste0(here::here(), "/results/sensitivity_compareToAll/Summary_allRadius_all.csv"))
