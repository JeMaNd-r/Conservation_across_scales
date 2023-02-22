#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#     Pair protected with unprotected       #
#          author: Romy Zeiss               #
#            date: 2023-02-21               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# pair one PA with most similar nonPA sample point
# based on physical properties and distance

library(here)
library(tidyverse)

library(igraph)  # to map network (pairing)

source("src/00_Parameters_functions.R")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

## Explore data
summary(as.factor(data_glob$PA)) #248 nonPA and 135 PAs
# number of observations (raw)
nrow(data_glob); nrow(data_glob[data_glob$PA,])  #383 with 135 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Scale variables for mahalanobis distance ####
for(i in mahal_vars){
  data_glob[,paste0(i,".z")] <- as.numeric(scale(data_glob[,i]))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check for pairing
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars_z)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars_z)]
pa_noPair <- c()
all_nonPA <- nonpa[0,]

for(i in 1:length(lc_names)){
    data_nonPA  <- nonpa[nonpa$LC==lc_names[i],]
    data_PA <- pa[pa$LC==lc_names[i],]
    sigma <- cov(data_nonPA[,mahal_vars_z]) 
    for(j in 1:nrow(data_PA)){
      mu <- as.numeric(data_PA[j,mahal_vars_z])
      data_nonPA[,as.character(data_PA[j,"Order_ID"])] <- 
        mahalanobis(data_nonPA[,mahal_vars_z], mu, sigma, tol=1e-30)
      #print(j)
    }
    min_mahal <- apply(X = data_nonPA[,as.character(data_PA$Order_ID)], MARGIN = 2, FUN = min, na.rm = TRUE)
    pa_noPair <- rbind(pa_noPair, cbind(names(min_mahal[min_mahal>mahal_thres]), min_mahal[min_mahal>mahal_thres]))
    
    all_nonPA <- full_join(all_nonPA, data_nonPA)
}
pa_noPair
nrow(pa_noPair) #nrow=12

unpaired_pa <- data_glob[data_glob$Order_ID %in% pa_noPair,] #12 PA sites can't be paired
write.csv(unpaired_pa, file=paste0(here::here(), "/intermediates/", "Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)

# look at Mahalanoubis distance values for each nonPA (Order_ID) and PA (columns)
#all_nonPA

# count how many "options" exist for one PA
count_nonPA <- data.frame("Order_ID"=NA, "n"=NA)
for(i in colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))]){
  temp_column <- all_nonPA[,i]
  
  count_nonPA <- rbind(count_nonPA,
                       c(i, nrow(all_nonPA[which(temp_column<=mahal_thres),"Order_ID"]) ))
  
}

count_nonPA <- count_nonPA %>% mutate("n"=as.numeric(n)) %>% arrange(n)
count_nonPA #number of nonPA sites per PA (Order_ID)

#View(all_nonPA %>% dplyr::select(Order_ID, count_nonPA[count_nonPA$n<10 & !is.na(count_nonPA$Order_ID), "Order_ID"]))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Network of possible combinations ####
# l <- lapply(colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))],
#        function(x) tibble("PA_ID"=rep(x, nrow(all_nonPA[which(all_nonPA[,x]<=mahal_thres),])),
#                          all_nonPA[which(all_nonPA[,x]<=mahal_thres),"Order_ID"]))
# l <- do.call(rbind, l)
# l
# network <- igraph::graph_from_data_frame(d=l, directed=T) 
# plot(network,
#      vertex.size=0.5, vertex.label.cex=0.2,
#      edge.arrow.size = 0.2, edge.arrow.width = 0.2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Remove sites ####

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with Order_ID. They have to been removed.
data_glob <- data_glob[!(data_glob$Order_ID %in% unpaired_pa$Order_ID),] 
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=371 with 123 PAs

# Remove sites that can only be paired once
data_glob <- data_glob[!(data_glob$Order_ID %in% count_nonPA[count_nonPA$n<3, "Order_ID"]),]
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=364 with 116 PAs
data_glob %>% group_by(LC, PA) %>% count()

#data_glob <- data_glob[data_glob$Order_ID!=28302282,]
#data_glob <- data_glob[data_glob$Order_ID!=29542352,]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Split data ####
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars_z)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars_z)]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Randomization ####
pa_pairs <- data.frame(Order_ID=NULL, nonPA=NULL, mahal.min=NULL, LC=NULL, run=NULL)
missing_pa <- data.frame(run=NA, pa.site=NA)[0,]
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
  times <- times+1; if(times > 1000) {break} #stop loop if reached 1000 trails
  times.with.error <- times.with.error + 1
  
  # add columns to store pairing temporary
  pa[,c("nonPA", "mahal.min")] <- NA  
  
  for(i in 1:length(lc_names)){
    data_PA <- pa[pa$LC==lc_names[i],] 
    data_nonPA  <- nonpa[nonpa$LC==lc_names[i],]
    
    # select environmental data only as matrix, remove rows with NAs
    data_PA <- data_PA[complete.cases(data_PA[,mahal_vars_z]),c("Order_ID", mahal_vars_z)]
    data_nonPA <- data_nonPA[complete.cases(data_nonPA[,mahal_vars_z]),c("Order_ID", mahal_vars_z)]
    
    sigma <- cov(data_nonPA[,mahal_vars_z]) # covariance/correlation between variables
    
    # random order of PA sites
    data_PA <- data_PA[order(sample(data_PA$Order_ID)),]
    
    # add empty columns
    data_nonPA[,as.character(data_PA$Order_ID)] <- NA
    
    # calculate Mahalanobis distance
    for(j in 1:nrow(data_PA)){
      mu = as.numeric(data_PA[j,mahal_vars_z])
      data_nonPA[,as.character(data_PA[j,"Order_ID"])] <- 
        mahalanobis(data_nonPA[,mahal_vars_z], mu, sigma, tol=1e-30)
      #print(j)
    }
    
    # add column to PA data with respective nonPA sites 
    # based on minimal mahalanobis distance
    data_PA[, c("nonPA", "mahal.min")]  <- NA
    #temp.col <- 0
    
    for(k in data_PA$Order_ID){
      
      # if Mahal. distances compared to all relevant nonPA are above threshold...
      if(min(data_nonPA[,as.character(k)])>=mahal_thres){
        #... stop and re-do run
        missing_pa <- rbind(missing_pa,c(times, k))  # to know why we've stopped
        print(paste0("Not all PA sites paired. Check PA ", k, "."))
        
        pa[pa$Order_ID==k, c("nonPA", "mahal.min")] <-
          cbind(NA, min(data_nonPA[,as.character(k)]))
        data_PA <- data_PA[data_PA$Order_ID!=k,] # remove respective PA site
      
      }else{
        # select (max. 10) nonPA sites with mahalanobis distance below threshold
        nonPA.pair <- data_nonPA[data_nonPA[,as.character(k)]<=mahal_thres,c("Order_ID", as.character(k))]
        nonPA.pair[nrow(nonPA.pair)+1:10,] <- NA #add empty rows for when there are less than 10 sites
        nonPA.pair <- nonPA.pair %>% arrange(2)
        nonPA.pair <- nonPA.pair[c(1:10),1]
        
        # sample one of the top 10 nonPA site that isn't NA 
        data_PA[data_PA$Order_ID==k,"nonPA"] <- resample(as.numeric(nonPA.pair[!is.na(nonPA.pair)]),1)
        
        # add value of distance
        data_PA[data_PA$Order_ID==k,"mahal.min"] <- unique(data_nonPA[data_nonPA$Order_ID==as.numeric(data_PA[data_PA$Order_ID==k,"nonPA"]),as.character(k)])
        data_nonPA <- data_nonPA[data_nonPA$Order_ID!=as.character(data_PA[data_PA$Order_ID==k,"nonPA"]),]
      }
    }
    # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
    if(nrow(data_PA) < min_size) {
      times <- times-1
      stop("Not enough PA sites paired.")
      
    } else {
      data_PA <- data_PA %>% sample_n(size = min_size, replace = FALSE)
      
      # add to result table
      pa[pa$Order_ID %in% data_PA$Order_ID, c("nonPA", "mahal.min")] <-
        cbind(data_PA[order(as.numeric(rownames(data_PA))),c("nonPA", "mahal.min")])
      
      # add to result table to analyse all at once below
      pa_pairs <- rbind(pa_pairs, cbind(data_PA[,c("Order_ID", "nonPA", "mahal.min")],
                                        "LC"=lc_names[i], times.with.error, times))
      
    }
    print(times)
  }
}
# show total count of unpaired (and removed) PAs and compare with number of paired sites
table(missing_pa[,2])  # can be larger than 0
#table(pa_pairs$nonPA) # counts should be lower or equal to number of runs (i.e. times)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Check & Save ####

# check for runs that failed (i.e. count < number of PA sites), and remove the respective pairs
# note: other result objects are not effected as they were overwritten
nrow(pa_pairs)  # should be 3 * min_size * 1000
pa_pairs <- pa_pairs %>% add_count(times.with.error) %>%
  filter(n==3*min_size)
nrow(pa_pairs)  # exactly 3 * min_size * 1000

# look what non-protected sites have (not) been paired to any PA
hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(data_glob[data_glob$PA==0,]$Order_ID, pa_pairs$nonPA))  # nonPA sites never used

#save(pa_pairs, file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_", Sys.Date(),".RData"))
write.csv(pa_pairs, file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_", Sys.Date(),".csv"), row.names=F)

