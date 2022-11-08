#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Prepare data for compare function       #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(terra)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v1.csv"))
data_glob

library(psych)   # to calculate CI of Cohens d effect size

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## load information about protected areas ####
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)

protect_df <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Add location in or outside of protected area ####
# after checking for location with ArcMap
lucas

# load table created in ArcMap with points in PA only
lucas.pa <- read.table("LUCAS_points_withPA.txt", sep=",", header=T)
colnames(lucas.pa) <- c("FID", "LUCAS_ID", "PA")

# add column with information about protected and non-protected sites
lucas <- full_join(lucas, lucas.pa[,2:3], by="LUCAS_ID") #....
#lucas[is.na(lucas$PA), "PA"] <- 0  #replace NA by 0
summary(as.factor(lucas$PA))

# number of observations (raw)
nrow(lucas); nrow(lucas[lucas$PA,])  #893 with 102 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## pair one PA with most similar nonPA sample point ####
# based on physical properties and distance

# define functions to be compared
fns = c("Basal_respiration","Cmic","Xylosidase","Cellulase","beta.glucosidase",
        "N.actylglucosaminidase","Acid.phosphatase","MWD","WSA")

## define variables to be compared between PA and nonPA, & threshold
mahal.vars <- c("TH_LAT", "TH_LONG", "Elevation", "AnnualPrecip", "AnnualTemp", 
                "MonthlyPrecipSum","MonthlyMeanTemp", "Organic_carbon",
                "Clay_content", "Silt_content",#"sand_yearCombine", #remove one because colinear
                "pH_H2O")

# keep complete cases only
lucas <- lucas[complete.cases(lucas[,c(mahal.vars, fns)]),]
nrow(lucas); nrow(lucas[lucas$PA==1,])  # n = 876 including 101 PAs

mahal.thres <- qchisq(.975, df=length(mahal.vars)) #21.92005

# scale variables for mahalanobis distance
for(i in mahal.vars){
  lucas[,paste0(i,".z")] <- as.numeric(scale(lucas[,i]))
}

mahal.vars <- paste0(mahal.vars, ".z")

# define each land cover type
lc.names <- c("Cropland", "Woodland", "Grassland", "Other")

# Note: There were no nonPA sites at all with mahalanobis distance below threshold
# for three protected sites with LUCAS_ID. They have to been removed.
lucas <- lucas[lucas$LUCAS_ID!=31783686,]
lucas <- lucas[lucas$LUCAS_ID!=28302282,]
lucas <- lucas[lucas$LUCAS_ID!=29542352,]

# split data
nonpa <- lucas[lucas$PA==0,c("LUCAS_ID", "LC_5", "PA", mahal.vars, fns)]
pa <- lucas[lucas$PA==1,c("LUCAS_ID", "LC_5", "PA", mahal.vars, fns)]

# randomization
p.list.total <- vector("list", length = 1000)
effect.size.d <- vector("list", length=1000)
pa.pairs <- data.frame(LUCAS_ID=NULL, nonPA=NULL, mahal.min=NULL, LC=NULL, run=NULL)
missing.pa <- data.frame(run=NA, pa.site=NA)[0,]
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
  times <- times+1; if(times > 1000) {break} #stop loop if reached 1000 trails
  times.with.error <- times.with.error + 1
  
  check.error <- try({  # for doing run again if error occurs (e.g. in pairing)
    
    # add columns to store pairing temporary
    pa[,c("nonPA", "mahal.min")] <- NA  
    
    for(i in 1:length(lc.names)){
      data_PA <- pa[pa$LC_5==lc.names[i],] 
      data_nonPA  <- nonpa[nonpa$LC_5==lc.names[i],]
      
      # select environmental data only as matrix, remove rows with NAs
      data_PA <- data_PA[complete.cases(data_PA[,mahal.vars]),c("LUCAS_ID", mahal.vars)]
      data_nonPA <- data_nonPA[complete.cases(data_nonPA[,mahal.vars]),c("LUCAS_ID", mahal.vars)]
      
      sigma <- cov(data_nonPA[,mahal.vars]) # covariance/correlation between variables
      
      # random order of PA sites
      data_PA <- data_PA[order(sample(data_PA$LUCAS_ID)),]
      
      # add empty columns
      data_nonPA[,as.character(data_PA$LUCAS_ID)] <- NA
      
      # calculate Mahalanobis distance
      for(j in 1:nrow(data_PA)){
        mu = as.numeric(data_PA[j,mahal.vars])
        data_nonPA[,as.character(data_PA[j,"LUCAS_ID"])] <- 
          mahalanobis(data_nonPA[,mahal.vars], mu, sigma, tol=1e-30)
        #print(j)
      }
      
      # add column to PA data with respective nonPA sites 
      # based on minimal mahalanobis distance
      data_PA[, c("nonPA", "mahal.min")]  <- NA
      #temp.col <- 0
      
      for(k in data_PA$LUCAS_ID){
        
        # if Mahal. distances compared to all relevant nonPA are above threshold...
        if(min(data_nonPA[,as.character(k)])>=mahal.thres){
          #... stop and re-do run
          missing.pa <- rbind(missing.pa,c(times, k))  # to know why we've stopped
          print(paste0("Not all PA sites paired. Check PA ", k, "."))
          
          pa[pa$LUCAS_ID==k, c("nonPA", "mahal.min")] <-
            c(NA, min(data_nonPA[,as.character(k)]))
          data_PA <- data_PA[data_PA$LUCAS_ID!=k,] # remove respective PA site
        }
        
        # select (max. 10) nonPA sites with mahalanobis distance below threshold
        nonPA.pair <- data_nonPA[data_nonPA[,as.character(k)]<=mahal.thres,c("LUCAS_ID", as.character(k))]
        nonPA.pair[nrow(nonPA.pair)+1:nrow(nonPA.pair)+10,] <- NA #add empty rows for when there are less than 10 sites
        nonPA.pair <- nonPA.pair[order(nonPA.pair[,2]),]
        nonPA.pair <- nonPA.pair[c(1:10),1]
        
        # sample one of the top 10 nonPA site that isn't NA 
        data_PA[data_PA$LUCAS_ID==k,"nonPA"] <- sample(as.character(nonPA.pair[!is.na(nonPA.pair)]),1)
        
        # add value of distance
        data_PA[data_PA$LUCAS_ID==k,"mahal.min"] <- unique(data_nonPA[data_nonPA$LUCAS_ID==data_PA[data_PA$LUCAS_ID==k,"nonPA"],as.character(k)])
        data_nonPA <- data_nonPA[data_nonPA$LUCAS_ID!=data_PA[data_PA$LUCAS_ID==k,"nonPA"],]
      }
      
      
      # add to result table
      pa[pa$LUCAS_ID %in% data_PA$LUCAS_ID, c("nonPA", "mahal.min")] <-
        data_PA[order(as.numeric(rownames(data_PA))),c("nonPA", "mahal.min")]
      
      # add to result table to analyse all at once below
      pa.pairs <- rbind(pa.pairs, cbind(data_PA[,c("LUCAS_ID", "nonPA", "mahal.min")],
                                        lc.names[i], times.with.error, times))
    }
    
    #})}  # stop here if only calculating p-value for all runs together...
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Compare difference between PA and nonPA sites ####
    
    # merge pairs of PA and nonPA in one table ####
    nonpa.subset <- nonpa[unique(nonpa$LUCAS_ID) %in% pa$nonPA,]
    lucas.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$LUCAS_ID),])
    #head(lucas.paired)
    
    # Perform t tests
    p.list <- list()
    effect.size.d[[times]] <- list()
    for(l in lc.names){
      p.list[[l]] <- vector("list", length(fns))
      names(p.list[[l]]) <- fns
      
      # subset data based on land cover type
      temp.PA <- lucas.paired[lucas.paired$LC_5==l & lucas.paired$PA==1,]
      temp.nonPA <- lucas.paired[lucas.paired$LC_5==l & lucas.paired$PA==0,]
      temp.nonPA <- temp.nonPA[order(match(temp.nonPA$LUCAS_ID, temp.PA$nonPA)), ]
      
      temp.cohens <- psych::cohen.d(rbind(temp.PA, temp.nonPA)[,c("PA",fns)], "PA")
      effect.size.d[[times]][[l]] <- cbind(lc=l, data.frame(temp.cohens$cohen.d), run=times)
      
      for(no.fns in 1:(length(fns))){
        # unpaired t test
        p.list[[l]][[no.fns]] <- t.test(temp.PA[,fns[no.fns]],temp.nonPA[,fns[no.fns]])
      }
      print(summary(temp.PA[,"nonPA"] == temp.nonPA[temp.nonPA$LUCAS_ID %in% temp.PA[,"nonPA"],"LUCAS_ID"])) # make sure that the sites are pairing properly
      
    }
    effect.size.d[[times]] <- do.call(rbind, effect.size.d[[times]])
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Make table of all p values for all functions and LC types
    p.table = as.data.frame(matrix(ncol=length(fns),nrow=5))
    colnames(p.table) = fns
    rownames(p.table) = c("Other","Cropland","Grassland", "Woodland","total")
    
    for(n in 1:ncol(p.table)){
      p.table[1,n] = p.list[["Cropland"]][[n]]["p.value"]
      p.table[2,n] = p.list[["Grassland"]][[n]]["p.value"]
      p.table[3,n] = p.list[["Woodland"]][[n]]["p.value"]
      p.table[4,n] = p.list[["Other"]][[n]]["p.value"]
      #... total missing
    }
    
    # add to overall result list
    p.list.total[[times]] <- p.table  #p-values from Chi-squared test
  })
  
  # do run again if there is an error (i.e. if no nonPA site with distance lower than threshold)
  if(is(check.error,"try-error")) {times <- times-1} else {print(check.error)}
  
} #end of randomization

save(effect.size.d,  file="d_1000_trails.RData")
