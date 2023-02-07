#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(psych)   # to calculate CI of Cohens d effect size

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

## Explore data
summary(as.factor(data_glob$PA))

# number of observations (raw)
nrow(data_glob); nrow(data_glob[data_glob$PA,])  #383 with 135 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## pair one PA with most similar nonPA sample point ####
# based on physical properties and distance

# define functions to be compared
fns = c("Soil_carbon_service", "OM_decomposition_service", "Water_regulation_service", 
        "Soil_stability_service", "Nutrient_service", "Pathogen_control", 
        "Richness_bacteria", "Fungi_18s_richness", "Invertebrate_18s_richness", 
        "Protist_18s_richness", "Ectomycorrhizal_18s_richness", 
        "Arbuscular_mycorrhizal_18s_richness", "Decomposers_18s_richness", 
        "Diss_Bacteria_std", "Diss_Fungi_std", "Diss_Protists_std", "Diss_invert_std")

## define variables to be compared between PA and nonPA, & threshold
mahal_vars <- c("Latitude_c", "Longitude_c", "Elevation", "AnnualPrec", "AnnualTemp", 
                #"MonthlyPrecipSum","MonthlyMeanTemp", 
                "Soil_pH", "Soil_salinity", "Clay_silt_c")

# keep complete cases only
data_glob <- data_glob[complete.cases(data_glob[,c(mahal_vars, fns)]),]
nrow(data_glob); nrow(data_glob[data_glob$PA==1,])  # n = 383 including 135 PAs

mahal_thres <- qchisq(.975, df=length(mahal_vars)) #21.92005

# scale variables for mahalanobis distance
for(i in mahal_vars){
  data_glob[,paste0(i,".z")] <- as.numeric(scale(data_glob[,i]))
}

mahal_vars <- paste0(mahal_vars, ".z")

# define each land cover type
lc.names <- c("Grassland", "Shrubland", "Woodland", "Other")

# rename land cover types in data
data_glob$LC <- data_glob$Eco_c
unique(data_glob$LC)
data_glob[data_glob$LC=="Forest", "LC"] <- "Woodland"
data_glob[data_glob$LC=="Moss_heath", "LC"] <- "Other"
unique(data_glob$LC)

# check for pairing
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars, fns)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars, fns)]
pa_noPair <- c()
all_nonPA <- nonpa[0,]

for(i in 1:length(lc.names)){
    data_nonPA  <- nonpa[nonpa$LC==lc.names[i],]
    data_PA <- pa[pa$LC==lc.names[i],]
    sigma <- cov(data_nonPA[,mahal_vars]) 
    for(j in 1:nrow(data_PA)){
      mu <- as.numeric(data_PA[j,mahal_vars])
      data_nonPA[,as.character(data_PA[j,"Order_ID"])] <- 
        mahalanobis(data_nonPA[,mahal_vars], mu, sigma, tol=1e-30)
      #print(j)
    }
    min_mahal <- apply(X = data_nonPA[,as.character(data_PA$Order_ID)], MARGIN = 2, FUN = min, na.rm = TRUE)
    pa_noPair <- rbind(pa_noPair, cbind(names(min_mahal[min_mahal>mahal_thres]), min_mahal[min_mahal>mahal_thres]))
    
    all_nonPA <- full_join(all_nonPA, data_nonPA)
}
pa_noPair
nrow(pa_noPair) #nrow=12

unpaired_pa <- data_glob[data_glob$Order_ID %in% pa_noPair,] #12 sites can't be paired
write.csv(unpaired_pa, file=paste0(here::here(), "/intermediates/", "Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)


all_nonPA

# count how many "options" exist for one PA
count_nonPA <- data.frame("Order_ID"=NA, "n"=NA)
for(i in colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))]){
  temp_column <- all_nonPA[,i]
  
  count_nonPA <- rbind(count_nonPA,
                       c(i, nrow(all_nonPA[which(temp_column<=mahal_thres),"Order_ID"]) ))
  
}

count_nonPA <- count_nonPA %>% mutate("n"=as.numeric(n)) %>% arrange(n)
count_nonPA

View(all_nonPA %>% dplyr::select(Order_ID, count_nonPA[count_nonPA$n<10 & !is.na(count_nonPA$Order_ID), "Order_ID"]))


# Note: There were no nonPA sites at all with mahalanobis distance below threshold
# for three protected sites with Order_ID. They have to been removed.
data_glob <- data_glob[!(data_glob$Order_ID %in% unpaired_pa$Order_ID),] 
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=371 with 123 PAs

# Remove sites that can only be paired once
data_glob <- data_glob[!(data_glob$Order_ID %in% count_nonPA[count_nonPA$n<5, "Order_ID"]),]

#data_glob <- data_glob[data_glob$Order_ID!=28302282,]
#data_glob <- data_glob[data_glob$Order_ID!=29542352,]

# split data
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars, fns)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars, fns)]

# randomization
p_list_total <- vector("list", length = 1000)
effect_size_d <- vector("list", length=1000)
pa_pairs <- data.frame(Order_ID=NULL, nonPA=NULL, mahal.min=NULL, LC=NULL, run=NULL)
missing_pa <- data.frame(run=NA, pa.site=NA)[0,]
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
  times <- times+1; if(times > 100) {break} #stop loop if reached 1000 trails
  times.with.error <- times.with.error + 1
  
  check.error <- try({  # for doing run again if error occurs (e.g. in pairing)
    
    # add columns to store pairing temporary
    pa[,c("nonPA", "mahal.min")] <- NA  
    
    for(i in 1:length(lc.names)){
      data_PA <- pa[pa$LC==lc.names[i],] 
      data_nonPA  <- nonpa[nonpa$LC==lc.names[i],]
      
      # select environmental data only as matrix, remove rows with NAs
      data_PA <- data_PA[complete.cases(data_PA[,mahal_vars]),c("Order_ID", mahal_vars)]
      data_nonPA <- data_nonPA[complete.cases(data_nonPA[,mahal_vars]),c("Order_ID", mahal_vars)]
      
      sigma <- cov(data_nonPA[,mahal_vars]) # covariance/correlation between variables
      
      # random order of PA sites
      data_PA <- data_PA[order(sample(data_PA$Order_ID)),]
      
      # add empty columns
      data_nonPA[,as.character(data_PA$Order_ID)] <- NA
      
      # calculate Mahalanobis distance
      for(j in 1:nrow(data_PA)){
        mu = as.numeric(data_PA[j,mahal_vars])
        data_nonPA[,as.character(data_PA[j,"Order_ID"])] <- 
          mahalanobis(data_nonPA[,mahal_vars], mu, sigma, tol=1e-30)
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
            c(NA, min(data_nonPA[,as.character(k)]))
          data_PA <- data_PA[data_PA$Order_ID!=k,] # remove respective PA site
        }
        
        # select (max. 10) nonPA sites with mahalanobis distance below threshold
        nonPA.pair <- data_nonPA[data_nonPA[,as.character(k)]<=mahal_thres,c("Order_ID", as.character(k))]
        nonPA.pair[nrow(nonPA.pair)+1:10,] <- NA #add empty rows for when there are less than 10 sites
        nonPA.pair <- nonPA.pair %>% arrange(2)
        nonPA.pair <- nonPA.pair[c(1:10),1]
        
        # sample one of the top 10 nonPA site that isn't NA 
        data_PA[data_PA$Order_ID==k,"nonPA"] <- sample(as.character(nonPA.pair[!is.na(nonPA.pair)]),1)
        
        # add value of distance
        data_PA[data_PA$Order_ID==k,"mahal.min"] <- unique(data_nonPA[data_nonPA$Order_ID==as.character(data_PA[data_PA$Order_ID==k,"nonPA"]),as.character(k)])
        data_nonPA <- data_nonPA[data_nonPA$Order_ID!=data_PA[data_PA$Order_ID==k,"nonPA"],]
      }
      
      
      # add to result table
      pa[pa$Order_ID %in% data_PA$Order_ID, c("nonPA", "mahal.min")] <-
        data_PA[order(as.numeric(rownames(data_PA))),c("nonPA", "mahal.min")]
      
      # add to result table to analyse all at once below
      pa_pairs <- rbind(pa_pairs, cbind(data_PA[,c("Order_ID", "nonPA", "mahal.min")],
                                        lc.names[i], times.with.error, times))
    }
    
#})}  # stop here if only calculating p-value for all runs together...
#table(missing_pa[,2])    


    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Compare difference between PA and nonPA sites ####
    
    # merge pairs of PA and nonPA in one table ####
    nonpa.subset <- nonpa[unique(nonpa$Order_ID) %in% pa$nonPA,]
    data_glob.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$Order_ID),])
    #head(data_glob.paired)
    
    # Perform t tests
    p.list <- list()
    effect_size_d[[times]] <- list()
    for(l in lc.names){
      p.list[[l]] <- vector("list", length(fns))
      names(p.list[[l]]) <- fns
      
      # subset data based on land cover type
      temp.PA <- data_glob.paired[data_glob.paired$LC==l & data_glob.paired$PA==1,]
      temp.nonPA <- data_glob.paired[data_glob.paired$LC==l & data_glob.paired$PA==0,]
      temp.nonPA <- temp.nonPA[order(match(temp.nonPA$Order_ID, temp.PA$nonPA)), ]
      
      temp.cohens <- psych::cohen.d(rbind(temp.PA, temp.nonPA)[,c("PA",fns)], "PA")
      effect_size_d[[times]][[l]] <- cbind(lc=l, data.frame(temp.cohens$cohen.d), run=times)
      
      for(no.fns in 1:(length(fns))){
        # unpaired t test
        p.list[[l]][[no.fns]] <- t.test(temp.PA[,fns[no.fns]],temp.nonPA[,fns[no.fns]])
      }
      print(summary(temp.PA[,"nonPA"] == temp.nonPA[temp.nonPA$Order_ID %in% temp.PA[,"nonPA"],"Order_ID"])) # make sure that the sites are pairing properly
      
    }
    effect_size_d[[times]] <- do.call(rbind, effect_size_d[[times]])
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Make table of all p values for all functions and LC types
    p_table = as.data.frame(matrix(ncol=length(fns),nrow=5))
    colnames(p_table) = fns
    rownames(p_table) = c("Other","Shrubland","Grassland", "Woodland","total")
    
    for(n in 1:ncol(p_table)){
      p_table[1,n] = p.list[["Shrubland"]][[n]]["p.value"]
      p_table[2,n] = p.list[["Grassland"]][[n]]["p.value"]
      p_table[3,n] = p.list[["Woodland"]][[n]]["p.value"]
      p_table[4,n] = p.list[["Other"]][[n]]["p.value"]
      #... total missing
    }
    
    # add to overall result list
    p_list_total[[times]] <- p_table  #p-values from Chi-squared test
  })
  
  # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
  if(is(check.error,"try-error")) {
    times <- times-1; print(times)
  } else {
      print(check.error); print(times)}
  
} #end of randomization

table(missing_pa[,2])


save(effect_size_d,  file="d_1000_trails.RData")
