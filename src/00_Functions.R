#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#               Functions                   #
#          author: Romy Zeiss               #
#            date: 2023-03-22               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#- - - - - - - - - - - - - - - - - - - - - -
## PREPARE DATA ####
#- - - - - - - - - - - - - - - - - - - - - -

### Extract PA status ####
f_extract_pa <- function(data, col_lon, col_lat, col_id, shp, col_pa){
  
  # data:     data frame with one column for longitude and latitude each
  # col_lon:  column name of longitude
  # col_lat:  column name of latitude
  # col_id:   column name with Sample ID
  # shp:      shapefile containing polygons of protected areas
  # col_pa:   name of layer in shp that contains information of PA 
  #           (if !is.na(col_pa), PA=1)
  
  temp_vect <- terra::vect(data[,c(col_lon, col_lat)], geom=c(col_lon, col_lat))
  #plot(temp_vect)
  
  # extract protection status
  temp_pa <- terra::extract(shp, temp_vect)
  temp_pa
  
  temp_pa$PA <- 0
  temp_pa[!is.na(temp_pa[,col_pa]), "PA"] <- 1
  
  temp_pa <- temp_pa %>% 
    full_join(data %>% 
                dplyr::select(all_of(col_id, col_lon, col_lat)) %>% 
                mutate(id.y=1:nrow(data))) %>%
    dplyr::select(-id.y) %>%
    dplyr::select(all_of(col_id, col_lon, col_lat, PA), everything())
  
  return(temp_pa)
}

#- - - - - - - - - - - - - - - - - - - - - -
### Extract elevation & climate data ####
f_extract_env <- function(data, col_lon, col_lat){
  
  # data:     data frame with one column for longitude and latitude each
  # col_lon:  column name of longitude
  # col_lat:  column name of latitude
  
  # Elevation
  elev <- terra::rast(paste0(here::here(), "/data_raw/Elevation_Worldclim_Global.tif"))
  data$Elevation <- terra::extract(elev, data[,c(col_lon, col_lat)])[,2]
  
  # Annual temperature
  annual_temp <- terra::rast(paste0(here::here(), "/data_raw/CHELSA_bio1_1981-2010_V.2.1.tif"))
  data$AnnualTemp <- terra::extract(annual_temp, data[,c(col_lon, col_lat)])[,2]
  data$AnnualTemp <- data$AnnualTemp / 10
  
  # Annual precipitation
  annual_prec <- terra::rast(paste0(here::here(), "/data_raw/CHELSA_bio12_1981-2010_V.2.1.tif"))
  data$AnnualPrec <- terra::extract(annual_prec, data[,c(col_lon, col_lat)])[,2]
  
  return(data)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add location in or outside of protected area ####
f_add_protect <- function(data, data_pa, col_id){
  
  # data:    focal dataframe 
  # data_pa: dataframe containing information about protection status
  # col_id:  name of column with site IDs (in both dataframes)
  
  ## Remove double entries (i.e., 2+ protected area types for 1 PA)
  data_pa <- data_pa %>% full_join(protect_type, by=c("IUCN_CAT" = "PA_type"))
  
  # take the highest = min (or lowest = max) protection level
  data_pa <- data_pa[, c(col_id, "PA_rank")] %>%
    group_by(data_pa[,col_id]) %>%
    summarize(across(PA_rank, list("min"=~ min(.x, na.rm = TRUE), 
                                   "max"=~ max(.x, na.rm = TRUE)))) %>%
    mutate("PA_rank_min" = ifelse(PA_rank_min == Inf, NA, PA_rank_min),
           "PA_rank_max" = ifelse(PA_rank_max == -Inf, NA, PA_rank_max)) %>%
    left_join(protect_type %>% rename("PA_type_min" = PA_type,
                                      "PA_rank_min" = PA_rank,
                                      "PA_protected_min" = PA_protected), 
              by="PA_rank_min") %>%
    left_join(protect_type %>% rename("PA_type_max" = PA_type,
                                      "PA_rank_max" = PA_rank,
                                      "PA_protected_max" = PA_protected), 
              by="PA_rank_max")
  
  # We decided to take the minimum (i.e., highest) level of protection because
  # of sites in the global dataset such as ID 49 (Order_ID 302) that has 
  # min=1 (Ia) and max=8 (Not Reported)
  
  # add column with information about protected and non-protected sites
  data <- data %>% 
    full_join(data_pa %>% dplyr::select(all_of(col_id), PA_type_min, PA_protected_min, PA_rank_min),
              by=as.character(col_id)) %>%
    rename("PA" = PA_protected_min,
           "PA_type" = PA_type_min,
           "PA_rank" = PA_rank_min) %>%
    mutate("PA" = ifelse(is.na(PA), 0, PA))
  
  rm(data_pa)
  return(data)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for colinearity between environmental variables ####

f_colinearity <- function(data, col_lon, col_lat, vars_env){
  
  # data:  focal dataframe containing the following columns
  # col_lon:  column name of longitude
  # col_lat:  column name of latitude
  # vars_env: column names of environmental variables that are used for pairing 
  #           based on environmental similarity (mahalanobis distance) 
  
  ## Calculate variable inflation factor (VIF)
  # VIF is the extent of correlation between one predictor and all others.
  # The lower VIF, the better we can tell what predictor contributed (most) to the model
  temp_env <- as.data.frame( data[,c(col_lon, col_lat, vars_env)] )
  
  ## VIF basen on raw data (explanatory raster stack)
  env_vif <- usdm::vif(temp_env)
  
  # which predictors should be excluded?
  vif_cor <- usdm::vifcor(temp_env, th=0.8)  #th = threshold correlation for exclusion
  # how: first find a pair of variables which has the maximum linear correlation 
  # (greater than th), and exclude one of them which has greater VIF. The 
  # procedure is repeated until no variable with a high correlation coefficient 
  # (grater than threshold) with other variables remains.
  
  vif_step <- usdm::vifstep(temp_env, th=10) #VIF >10
  
  # merge both data.frames
  env_vif <- env_vif %>% rename("VIF_raw" = VIF) %>% full_join(vif_step@results) %>%
    full_join(as.data.frame(vif_cor@corMatrix) %>% mutate("Variables"=rownames(vif_cor@corMatrix)))
  
  #env_vif
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Calculate correlations between predictors
  # https://github.com/joaofgoncalves/GoncalvesAna_et_al_2021/tree/master/RCODE/PostModelAnalyses
  
  corMatSpearman <- cor(temp_env, use="complete.obs", method="spearman") %>% round(2)
  corMatPearson <- cor(temp_env, use="complete.obs", method="pearson") %>% round(2)
  
  rm(vif_cor, temp_env, vif_step)
  
  return(list("env_vif" = env_vif,
              "corMatSpearman" = corMatSpearman,
              "corMatPearson" = corMatPearson))
}


#- - - - - - - - - - - - - - - - - - - - - -
## PAIRING ####
#- - - - - - - - - - - - - - - - - - - - - -

### Scale variables for mahalanobis distance ####
f_scale_vars <- function(data, vars){
  
  # data: data frame containing vars
  # vars: vector of column names that are going to be scaled
  
  if(sum(vars %in% colnames(data)) != length(vars)){ 
    cat("Error: Please check your data. Not all 'vars' are present in your dataset.")
  } else {
    # scale each variable individually
    for(i in vars){
      data[,paste0(i,".z")] <- as.numeric(scale(data[,i]))
    }
    
    cat("Done: Scaling of all variables completed.")
    return(data)
  }
}

#- - - - - - - - - - - - - - - - - - - - - -
### Check for pairing ####
f_check_pairs <- function(data, col_id, col_lc, vars_z){
  
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
    min_mahal <- apply(X = temp_nonPA[,as.character(temp_PA %>% pull(col_id))], MARGIN = 2, FUN = min, na.rm = TRUE)
    pa_noPair <- rbind(pa_noPair, cbind(names(min_mahal[min_mahal>mahal_thres]), min_mahal[min_mahal>mahal_thres]))
    
    all_nonPA <- full_join(all_nonPA, temp_nonPA, by=c(col_id, col_lc, "PA", vars_z))
  }
  pa_noPair
  nrow(pa_noPair) #nrow=12
  
  unpaired_pa <- data[data[,col_id] %in% pa_noPair,] #12 PA sites can't be paired
  write.csv(unpaired_pa, file=paste0(here::here(), "/intermediates/", temp_scale,  "/Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)
  
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
  
  rm(temp_PA, temp_nonPA, temp_column, unpaired_pa, pa_noPair, min_mahal, nonpa, pa)
  
  return(list(count_nonPA, all_nonPA)) #number of nonPA sites per PA (Order_ID)
}

#- - - - - - - - - - - - - - - - - - - - - -
### Re-define sample function ####
# if only 1 value in sample(), it considers it as vector
# used in f_pairing
f_resample <- function(x, ...) x[sample.int(length(x), ...)]

#- - - - - - - - - - - - - - - - - - - - - -
### Build pairs ####
f_pairing <- function(data, col_id, col_lc, vars_z){
  
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
      temp_PA <- temp_PA[order(sample(temp_PA %>% pull(col_id))),]
      
      # add empty columns
      temp_nonPA[,as.character(temp_PA %>% pull(col_id))] <- NA
      
      # calculate Mahalanobis distance
      for(j in 1:nrow(temp_PA)){
        mu = as.numeric(temp_PA[j,vars_z])
        temp_nonPA[,as.character(temp_PA[j,col_id])] <- 
          mahalanobis(temp_nonPA[,vars_z], mu, sigma, tol=1e-30)
        #print(j)
      }
      
      # add column to PA data with respective nonPA sites 
      # based on minimal mahalanobis distance
      temp_PA[, c("nonPA", "mahal.min")]  <- NA
      #temp.col <- 0
      
      for(k in temp_PA %>% pull(col_id) %>% sample(size = nrow(temp_PA))){
        
        # if Mahal. distances compared to all relevant nonPA are above threshold...
        if(min(temp_nonPA[,as.character(k)])>=mahal_thres){
          #... stop and re-do run
          missing_pa <- rbind(missing_pa,c(times, k))  # to know why we've stopped
          print(paste0("Not all PA sites paired. Check PA ", k, "."))
          
          pa[pa %>% pull(col_id)==k, c("nonPA", "mahal.min")] <-
            cbind(NA, min(temp_nonPA[,as.character(k)]))
          temp_PA <- temp_PA[temp_PA %>% pull(col_id)!=k,] # remove respective PA site
          
        }else{
          # select (max. 10) nonPA sites with mahalanobis distance below threshold
          nonPA.pair <- temp_nonPA[temp_nonPA[,as.character(k)]<=mahal_thres,c(col_id, as.character(k))]
          nonPA.pair[nrow(nonPA.pair)+1:10,] <- NA #add empty rows for when there are less than 10 sites
          nonPA.pair <- nonPA.pair %>% arrange(2)
          nonPA.pair <- nonPA.pair[c(1:10),1]
          
          # sample one of the top 10 nonPA site that isn't NA 
          temp_PA[temp_PA %>% pull(col_id)==k,"nonPA"] <- f_resample(as.numeric(nonPA.pair[!is.na(nonPA.pair)]),1)
          
          # add value of distance
          temp_PA[temp_PA %>% pull(col_id)==k,"mahal.min"] <- unique(temp_nonPA[temp_nonPA %>% pull(col_id)==as.numeric(temp_PA[temp_PA %>% pull(col_id)==k,"nonPA"]),as.character(k)])
          temp_nonPA <- temp_nonPA[temp_nonPA %>% pull(col_id)!=as.character(temp_PA[temp_PA %>% pull(col_id)==k,"nonPA"]),]
        }
      }
      # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
      if(nrow(temp_PA) < min_size) {
        times <- times-1
        print("Not enough PA sites paired.")
        
      } else {
        temp_PA <- temp_PA %>% sample_n(size = min_size, replace = FALSE)
        
        # add to result table
        pa[(pa %>% pull(col_id)) %in% (temp_PA %>% pull(col_id)), c("nonPA", "mahal.min")] <-
          cbind(temp_PA[order(as.numeric(rownames(temp_PA))),c("nonPA", "mahal.min")])
        
        # add to result table to analyse all at once below
        pa_pairs <- rbind(pa_pairs, cbind(temp_PA %>% dplyr::select(all_of(col_id), "nonPA", "mahal.min") %>% rename(ID=col_id),
                                          "LC"=lc_names[i], times.with.error, times))
        
      }
      print(times)
    }
  }
  
  rm(temp_PA, temp_nonPA, pa, nonpa, nonPA.pair)
  return(list(missing_pa = missing_pa, pa_pairs = pa_pairs))
  
}

#- - - - - - - - - - - - - - - - - - - - - -
## COMPARISON ####
#- - - - - - - - - - - - - - - - - - - - - -
# Compare difference between PA and nonPA sites

#- - - - - - - - - - - - - - - - - - - - - -
### Calculate effect size (Cohen's d) ####
f_compare_pa_nonpa <- function(data, data_pairs, col_id, col_fns){
  
  # data:       focal dataframe
  # data_pairs: list of length = number_times with pairings    
  # col_id:     name of column with site IDs
  # col_fns:    vector with column names for response variables
  
  # take one run (1:number_times), estimate p and effect size (and Bayesian), take next
  d_list <- vector("list", length = number_times)
  
  # progress bar
  progress_bar <- txtProgressBar(min = 0, max = length(d_list), initial = 0,
                                 style=3, title = "Processing") 
  
  for(x in unique(data_pairs$times)){
    
    setTxtProgressBar(progress_bar,x)
    #print("#########################")
    #print(paste0("Run number ", x))
    
    # subset data
    temp_pairs <- data_pairs[data_pairs$times==x,]
    
    # Perform tests
    d_list[[x]] <- list()
    
    # individual tests per LC type
    for(lc in lc_names){
      
      # subset data based on land cover type
      temp_PA <- temp_pairs %>% 
        left_join(data[,c(col_id, "PA", col_fns)], by=c("ID"=as.character(col_id))) %>% 
        filter(LC==lc) %>% 
        dplyr::select(-nonPA, -mahal.min)
      temp_nonPA <- temp_pairs %>% 
        left_join(data[,c(col_id, "PA", col_fns)], by=c("nonPA"=as.character(col_id))) %>% 
        filter(LC==lc) %>% 
        mutate("ID"=nonPA) %>% 
        dplyr::select(-nonPA, -mahal.min)
      
      temp_d <- psych::cohen.d(rbind(temp_PA, temp_nonPA)[,c("PA",col_fns)], group="PA")
      d_list[[x]][[lc]] <- data.frame(lc=lc, 
                                      data.frame(temp_d$cohen.d),
                                      fns = rownames(data.frame(temp_d$cohen.d)),
                                      run=x,
                                      row.names = NULL)
      
      #print(sum(temp_pairs[temp_pairs$LC==lc,"nonPA"] == temp_nonPA[,"ID"])==min_size) # make sure that the sites are pairing properly
      # print should give min_size (TRUE = fitting pairs) for each LC types * 1000 runs
      
    }
    
    # combine individual list elements (df) into one df
    d_list[[x]] <- do.call(rbind, d_list[[x]])
    
    # remove rownames
    rownames(d_list[[x]]) <- NULL
    
    close(progress_bar)
  }
  
  return(d_list)
}

#- - - - - - - - - - - - - - - - - - - - - -
### Compare PA types (incl. nonPA) ####

f_compare_pa_types <- function(data, protect_levels, col_fns) {
  
  # data:           focal dataframe
  # protect_levels: dataframe containing translation of PA type levels 
  # col_fns:        names of columns of response variables
  
  # create empty list
  pars_list <- vector("list", length = length(lc_names))
  names(pars_list) <- lc_names
  pars_list <- lapply(pars_list, 
                      function(x) {
                        x <- vector("list", length = length(col_fns))
                        names(x) <- col_fns
                        
                        lapply(x, 
                               function(y) {
                                 y <- vector("list", length = 1)
                               }
                        )
                      }
  )
  
  stan_data <- pars_list
  stan_fit <- pars_list
  
  # individual tests per LC type
  for(lc in lc_names){
    
    print("#####################################")
    print(paste0("Run land cover ", lc))
    
    # subset data based on land cover type
    temp_data <- data %>%
      full_join(protect_levels, by=c("LC", "PA_rank", "PA_type")) %>%
      filter(LC==lc) %>%
      arrange(PA_rank)
    
    for(no_fns in 1:(length(col_fns))){
      
      temp_n_groups <- pull(protect_levels %>% filter(LC==lc) %>% count())
      
      stan_data[[lc]][[no_fns]] <- list(n = nrow(temp_data),
                                        y = pull(temp_data[,col_fns[no_fns]]),
                                        group = temp_data$PA_rank_new,
                                        n_group = temp_n_groups)
      
      stan_fit[[lc]][[no_fns]] <- sampling(stan_model, data = stan_data[[lc]][[no_fns]],
                                           chains = 4, iter = 10000, warmup = 2000,
                                           show_messages = TRUE)
      #print(stan_fit, digits = 3, probs = c(0.025, 0.975))
      
      temp_pars <- rstan::extract(stan_fit[[lc]][[no_fns]], pars=c(paste0("a[", 1:temp_n_groups, "]"), "sigma", "mu_a", "sigma_a"))
      names(temp_pars)[1:temp_n_groups] <- unique(temp_data$PA_rank)
      
      pars_list[[lc]][[no_fns]] <- temp_pars
      
      stan_fit[[lc]][[no_fns]] <- NULL
      stan_data[[lc]][[no_fns]] <- NULL
      
    }
    
    gc()
    #unlink(file.path("tmp", "Rtmp*"), recursive = T)
    
    # remove temporary files
    dso_filenames <- dir(tempdir(), pattern=.Platform$dynlib.ext)
    filenames  <- dir(tempdir())
    for (i in seq(dso_filenames))
      try(dyn.unload(file.path(tempdir(), dso_filenames[i])))
    for (i in seq(filenames))
      if (file.exists(file.path(tempdir(), filenames[i])) & nchar(filenames[i]) < 42) # some files w/ long filenames that didn't like to be removeed
        file.remove(file.path(tempdir(), filenames[i]))
    
  }
  return(pars_list)
  
}  

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine Bayesian list into df ####
# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values

f_combine_pars_list <- function(pars_list){
  
  # pars_list: output of f_compare_pa_types
  
  pars_sample <- lapply(pars_list, function (x){ #across lc types
    lapply(x, function(y){ #across functions
      dplyr::bind_rows(y)
    })
  })
  
  rm(pars_list)
  gc()
  
  pars_sample <- lapply(pars_sample, function (x){ #across lc types
    dplyr::bind_rows(x, .id = "fns")
  })
  
  pars_sample <- dplyr::bind_rows(pars_sample, .id="lc")
  
  return(pars_sample)
}


#- - - - - - - - - - - - - - - - - - - - - -
## ... ####
#- - - - - - - - - - - - - - - - - - - - - -