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
library(terra) # to get distance of points

library(igraph)  # to map network (pairing)

source("src/00_Parameters.R")
source("src/00_Functions.R")


#temp_scale <- "global"
#temp_scale <- "continental"
temp_scale <- "regional"

# create directory for intermediate results
if(!dir.exists(paste0(here::here(), "/intermediates/", temp_scale))){
  dir.create(paste0(here::here(), "/intermediates/", temp_scale))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####

data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

## Explore data
summary(as.factor(data_clean$PA)) #G: 206 nonPA and 42 PAs, C: 743 vs. 61, R: 275 vs. 49
# number of observations (raw)
nrow(data_clean); nrow(data_clean[data_clean$PA,])  #G: 248 with 42 PAs, C: 807 vs. 64, R: 324 vs. 49

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Scale variables for mahalanobis distance ####

data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for pairing ####

## if continental data, remove nonPA with geographical distance >500km
if(temp_scale == "continental"){
  # calculate geographical distance
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
  
  # save for later
  save(list_dist, file = paste0(here::here(), "/intermediates/", temp_scale, "/Geographically_close_sites.RData"))
}

if(temp_scale == "continental"){
  load(paste0(here::here(), "/intermediates/", temp_scale, "/Geographically_close_sites.RData")) #list_dist

  # check number of sites remaining (e.g. min_size)
  do.call(rbind, lapply(list_dist, length)) %>% sort() #2 lower than min_size=10 (840 and 368), may be more when considering LC types
  
}

count_nonPA <- f_check_pairs(data = data_clean, 
                             col_id = "SampleID", col_lc = "LC", 
                             vars_z = mahal_vars_z)
all_nonPA <- count_nonPA[[2]]
count_nonPA <- count_nonPA[[1]] #G: 14 <5, 19 <7, 23 <10; G-together: 3 <22; C: 16 without enough (10) nonPAs for pairing; R: 8 without enough
count_nonPA

#View(all_nonPA %>% dplyr::select(SampleID, count_nonPA[count_nonPA$n<10 & !is.na(count_nonPA$SampleID), "SampleID"]))

# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ### Network of possible combinations ####
# l <- lapply(colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))],
#        function(x) tibble("PA_ID"=rep(x, nrow(all_nonPA[which(all_nonPA[,x]<=mahal_thres),])),
#                          all_nonPA[which(all_nonPA[,x]<=mahal_thres),"SampleID"]))
# l <- do.call(rbind, l)
# l
# network <- igraph::graph_from_data_frame(d=l, directed=T)
# plot(network,
#      vertex.size=0.5, vertex.label.cex=0.2,
#      edge.arrow.size = 0.2, edge.arrow.width = 0.2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Remove sites ####

unpaired_pa <- read.csv(sort(list.files(here::here("intermediates", temp_scale), pattern = "Unpaired", full.names = TRUE), decreasing = TRUE)[1])
head(unpaired_pa) 

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with SampleID. They have to been removed.
data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #G/G-t: nrow=248 with 42 PAs, C: 807 vs. 64; R. 324 with 49

# Remove sites that can only be paired less than min_nonPA times
# start with something small, then check how many possible;
#if(temp_scale == "global") min_nonPA <- 5 #not needed when analyzing all LC types together
data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA < min_nonPA, "SampleID"]),]
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #nrow=234 with 28 PAs; G-together: 245 with 39; C: 791 wit 48; R: 316 with 41
data_clean %>% group_by(LC, PA) %>% count()

# save data
write_csv(data_clean, paste0(here::here(), "/intermediates/Data_paired_", temp_scale, ".csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Randomization ####

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

if(temp_scale == "global"){
  lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
  #min_size <- 39 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
} 
if(temp_scale == "continental"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  #min_size <- 14 # number of samples/ sites that should be paired per LC type
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  #min_size <- 7 # number of samples/ sites that should be paired per LC type
}
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
### Check & Save ####

pa_pairs <- list_pairs$pa_pairs

# check for runs that failed (i.e. count < number of PA sites), and remove the respective pairs
# note: other result objects are not effected as they were overwritten
nrow(pa_pairs)  # should be length(lc_names) * min_size * 1000
pa_pairs <- pa_pairs %>% add_count(times.with.error) %>%
  filter(n==length(lc_names)*min_size)
nrow(pa_pairs)  # exactly length(lc_names)* min_size * 1000

#save(pa_pairs, file=paste0(here::here(), "/intermediates/", temp_scale, "/Pairs_paNonpa_1000trails_", Sys.Date(),".RData"))
write_csv(pa_pairs, file=paste0(here::here(), "/intermediates/", temp_scale, "/Pairs_paNonpa_1000trails_", Sys.Date(),".csv"))

