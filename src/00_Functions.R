#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#               Functions                   #
#          author: Romy Zeiss               #
#            date: 2023-03-22               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#- - - - - - - - - - - - - - - - - - - - - -
## PAIRING ####
#- - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - -
### Re-define sample function ####
# if only 1 value in sample(), it considers it as vector
#resample <- function(x, ...) x[sample.int(length(x), ...)]

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
    return(data)
    cat("Done: Scaling of all variables completed.")
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
  
  unpaired_pa <- data_glob[data_glob$Order_ID %in% pa_noPair,] #12 PA sites can't be paired
  write.csv(unpaired_pa, file=paste0(here::here(), "/intermediates/", "Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)
  
  cat("#---------------------------------------------------", sep="\n")
  cat(paste0("Saved: csv file with unpaired sites is saved under: "), 
      paste0(here::here(), "/intermediates/", "Unpaired_protected_sites_", Sys.Date(), ".csv"),
      sep="\n")
  cat("#---------------------------------------------------", sep="\n")
  
  # look at Mahalanoubis distance values for each nonPA (Order_ID) and PA (columns)
  #all_nonPA
  
  # count how many "options" exist for one PA
  count_nonPA <- data.frame("Order_ID"=NA, "No_nonPA"=NA)
  for(i in colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))]){
    temp_column <- all_nonPA[,i]
    
    count_nonPA <- rbind(count_nonPA,
                         c(i, nrow(all_nonPA[which(temp_column<=mahal_thres),"Order_ID"]) ))
    
  }
  
  count_nonPA <- count_nonPA %>% 
    mutate("No_nonPA"=as.numeric(No_nonPA)) %>% 
    arrange(No_nonPA)
  
  print("Number of nonPA sites per PA (col_ID)")
  print(head(count_nonPA))
  return(count_nonPA) #number of nonPA sites per PA (Order_ID)
  
  rm(temp_PA, temp_nonPA, temp_column, unpaired_pa, pa_noPair, min_mahal, nonpa, pa)
}

