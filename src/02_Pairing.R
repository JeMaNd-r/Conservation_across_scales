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

source("src/00_Parameters.R")
source("src/00_Functions.R")


# temp_scale <- "global"
temp_scale <- "continental"
#temp_scale <- "regional"

# create directory for intermediate results
if(!dir.exists(paste0(here::here(), "/intermediates/", temp_scale))){
  dir.create(paste0(here::here(), "/intermediates/", temp_scale))
}

# set date of latest analysis
if(temp_scale == "global") temp_date <- "2023-12-01"
if(temp_scale == "continental") temp_date <- "2023-12-14"
if(temp_scale == "regional") temp_date <- "2023-12-14"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####

data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

## Explore data
summary(as.factor(data_clean$PA)) #G: 309 nonPA and 74 PAs, C: 753 vs. 65, R: 302 vs. 52
# number of observations (raw)
nrow(data_clean); nrow(data_clean[data_clean$PA,])  #G: 383 with 74 PAs, C: 818 vs. 65, R: 354 vs. 52

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Scale variables for mahalanobis distance ####

data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for pairing ####

count_nonPA <- f_check_pairs(data = data_clean, 
                             col_id = "SampleID", col_lc = "LC", 
                             vars_z = mahal_vars_z)
all_nonPA <- count_nonPA[[2]]
count_nonPA <- count_nonPA[[1]]
#head(count_nonPA)

#View(all_nonPA %>% dplyr::select(SampleID, count_nonPA[count_nonPA$n<10 & !is.na(count_nonPA$SampleID), "SampleID"]))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Network of possible combinations ####
l <- lapply(colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))],
       function(x) tibble("PA_ID"=rep(x, nrow(all_nonPA[which(all_nonPA[,x]<=mahal_thres),])),
                         all_nonPA[which(all_nonPA[,x]<=mahal_thres),"SampleID"]))
l <- do.call(rbind, l)
l
network <- igraph::graph_from_data_frame(d=l, directed=T)
plot(network,
     vertex.size=0.5, vertex.label.cex=0.2,
     edge.arrow.size = 0.2, edge.arrow.width = 0.2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Remove sites ####

unpaired_pa <- read.csv(paste0(here::here(), "/intermediates/", temp_scale, "/Unpaired_protected_sites_", temp_date, ".csv"))
head(unpaired_pa) 

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with SampleID. They have to been removed.
data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #G: nrow=383 with 74 PAs, C: 818 vs. 65

# Remove sites that can only be paired less than 10 times
data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA<10, "SampleID"]),]
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #nrow=365 with 56 PAs; C: 814 vs. 61
data_clean %>% group_by(LC, PA) %>% count()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Randomization ####

# check LC types and number of (protected) sites
table(data_clean$LC, data_clean$PA)

# based on number of sites per LC, exclude LC
# e.g. global: only 3 unprotected on 0 protected for Other
# continental: only 1 Shrubland protected, 28 unprotected; and 0 PA in Other
# regional: PA only min. 7 -> decrease minimum size number to 7, 
#           exclude Shrublands & Others to get it running (otherwise no complete pairing achieved) 
if(temp_scale == "global") lc_names <- lc_names[lc_names != "Other"]
if(temp_scale == "continental") lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
  min_size <- 7
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

