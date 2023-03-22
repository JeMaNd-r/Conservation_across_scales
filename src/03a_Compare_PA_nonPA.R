#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#        (Frequentist approach)             #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(psych)   # to calculate CI of Cohens d effect size

source("src/00_Parameters_functions.R")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load PA-nonPA pairs ####
pa_pairs <- read.csv(file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_2023-02-22.csv"))

# # add columns of functions for each PA and nonPA (wide format)
# pa_pairs <- merge(pa_pairs, data_glob[,c("Order_ID", fns)], by="Order_ID")
# pa_pairs <- merge(pa_pairs, data_glob[,c("Order_ID", fns)], by.x="nonPA", by.y="Order_ID",
#                   suffixes = c(".pa", ".nonpa"))

head(pa_pairs)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA sites ####
# take one run (1:1000), estimate p and effect size (and Bayesian), take next
d_list <- vector("list", length = 1000)

# progress bar
progress_bar <- txtProgressBar(min = 0, max = length(p_list), initial = 0,
                               style=3, title = "Processing") 

for(x in unique(pa_pairs$times)){
  
  setTxtProgressBar(progress_bar,x)
  #print("#########################")
  #print(paste0("Run number ", x))
  
  # subset data
  temp_pairs <- pa_pairs[pa_pairs$times==x,]
  
  # Perform tests
  d_list[[x]] <- list()
  
  # individual tests per LC type
  for(lc in lc_names){
    
    # subset data based on land cover type
    temp_PA <- temp_pairs %>% 
      left_join(data_glob[,c("Order_ID", "PA", fns)], by=c("Order_ID"="Order_ID")) %>% 
      filter(LC==lc) %>% 
      dplyr::select(-nonPA, -mahal.min)
    temp_nonPA <- temp_pairs %>% 
      left_join(data_glob[,c("Order_ID", "PA", fns)], by=c("nonPA"="Order_ID")) %>% 
      filter(LC==lc) %>% 
      mutate("Order_ID"=nonPA) %>% 
      dplyr::select(-nonPA, -mahal.min)
    
    temp_d <- psych::cohen.d(rbind(temp_PA, temp_nonPA)[,c("PA",fns)], "PA")
    d_list[[x]][[lc]] <- data.frame(lc=lc, 
                                    data.frame(temp_d$cohen.d),
                                    fns = rownames(data.frame(temp_d$cohen.d)),
                                    run=x,
                                    row.names = NULL)
    
    #print(sum(temp_pairs[temp_pairs$LC==lc,"nonPA"] == temp_nonPA[,"Order_ID"])==min_size) # make sure that the sites are pairing properly
    # print should give min_size (TRUE = fitting pairs) for each LC types * 1000 runs
    
  }
  
  # combine individual list elements (df) into one df
  d_list[[x]] <- do.call(rbind, d_list[[x]])
  
  # remove rownames
  rownames(d_list[[x]]) <- NULL
  
  close(progress_bar)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables & effect sizes ####
save(d_list,  file=paste0(here::here(), "/results/d_1000_trails.RData"))




