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

# save sample function (if only 1 value in sample(), it considers it as vector)
resample <- function(x, ...) x[sample.int(length(x), ...)]

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
lc.names <- c("Grassland", "Shrubland", "Woodland")

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

# split data
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars, fns)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars, fns)]

# randomization
p_list_total <- vector("list", length = 1000)
effect_size_d <- vector("list", length=1000)
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
 
}
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Compare difference between PA and nonPA sites ####
    
    # merge pairs of PA and nonPA in one table ####
    nonpa.subset <- nonpa[unique(nonpa$Order_ID) %in% pa$nonPA,]
    #data_glob.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$Order_ID),])
    #head(data_glob.paired)
    
    # Perform t tests
    p.list <- list()
    effect_size_d[[times]] <- list()
    for(l in lc.names){
      p.list[[l]] <- vector("list", length(fns))
      names(p.list[[l]]) <- fns
      
      # subset data based on land cover type
      temp.PA <- pa[pa$LC==l & pa$PA==1,]
      temp.nonPA <- nonpa.subset[nonpa.subset$LC==l & nonpa.subset$PA==0,]
      #NEEDS TO BE FIXED: temp.nonPA <- full_join(temp.PA, temp.nonPA, by=c("nonPA"="Order_ID"), suffix)
      
      temp.cohens <- psych::cohen.d(rbind(temp.PA, temp.nonPA)[,c("PA",fns)], "PA")
      effect_size_d[[times]][[l]] <- cbind(lc=l, data.frame(temp.cohens$cohen.d), run=times)
      
      for(no.fns in 1:(length(fns))){
        # unpaired t test
        p.list[[l]][[no.fns]] <- t.test(temp.PA[,fns[no.fns]],temp.nonPA[,fns[no.fns]])
      }
      print(summary(temp.PA[,"nonPA"] == temp.nonPA[temp.nonPA$Order_ID %in% as.numeric(unlist(temp.PA[,"nonPA"])),"Order_ID"])) # make sure that the sites are pairing properly
      
    }
    effect_size_d[[times]] <- do.call(rbind, effect_size_d[[times]])
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Make table of all p values for all functions and LC types
    p_table <- as.data.frame(matrix(ncol=length(fns),nrow=5))
    colnames(p_table) <- fns
    rownames(p_table) <- c("Other","Shrubland","Grassland", "Woodland","total")
    
    for(n in 1:ncol(p_table)){
      for(lc in 1:length(lc.names)){
        p_table[which(rownames(p_table)==lc.names[lc]),n] <- p.list[[lc.names[lc]]][[n]]["p.value"]
      }
    }
    
    # add to overall result list
    p_list_total[[times]] <- p_table  #p-values from Chi-squared test
  

save(effect_size_d,  file="d_1000_trails.RData")
