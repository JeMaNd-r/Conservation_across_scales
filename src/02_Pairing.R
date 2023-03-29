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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## GLOBAL ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

## Explore data
summary(as.factor(data_glob$PA)) #248 nonPA and 135 PAs
# number of observations (raw)
nrow(data_glob); nrow(data_glob[data_glob$PA,])  #383 with 135 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Scale variables for mahalanobis distance ####

data_glob <- f_scale_vars(data = data_glob, vars = mahal_vars)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for pairing ####

count_nonPA <- f_check_pairs(data = data_glob, 
                             col_id = "Order_ID", col_lc = "LC", 
                             vars_z = mahal_vars_z)
head(count_nonPA)

#View(all_nonPA %>% dplyr::select(Order_ID, count_nonPA[count_nonPA$n<10 & !is.na(count_nonPA$Order_ID), "Order_ID"]))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Network of possible combinations ####
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
### Remove sites ####

unpaired_pa <- read.csv(paste0(here::here(), "/intermediates/", "Unpaired_protected_sites_global.csv"))
head(unpaired_pa)

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with Order_ID. They have to been removed.
data_glob <- data_glob[!(data_glob$Order_ID %in% unpaired_pa$Order_ID),] 
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=371 with 123 PAs

# # Remove sites that can only be paired once
# data_glob <- data_glob[!(data_glob$Order_ID %in% count_nonPA[count_nonPA$n<3, "Order_ID"]),]
# nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=364 with 116 PAs
data_glob %>% group_by(LC, PA) %>% count()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Randomization ####

list_pairs <- f_pairing(data = data_glob, 
                        col_id = "Order_ID", col_lc = "LC",
                        vars_z = mahal_vars_z)


# show total count of unpaired (and removed) PAs and compare with number of paired sites
table(list_pairs$missing_pa[,2])  # can be larger than 0, 0 is perfect
table(list_pairs$pa_pairs$nonPA) # counts should be lower or equal to number of runs (i.e. times)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check & Save ####

pa_pairs <- list_pairs$pa_pairs

# check for runs that failed (i.e. count < number of PA sites), and remove the respective pairs
# note: other result objects are not effected as they were overwritten
nrow(pa_pairs)  # should be 3 * min_size * 1000
pa_pairs <- pa_pairs %>% add_count(times.with.error) %>%
  filter(n==3*min_size)
nrow(pa_pairs)  # exactly 3 * min_size * 1000

# look what non-protected sites have (not) been paired to any PA
hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(data_glob[data_glob$PA==0,]$Order_ID, pa_pairs$nonPA))  # nonPA sites never used

mean(table(pa_pairs$nonPA))
mean(table(pa_pairs$ID))

sd(table(pa_pairs$nonPA))
sd(table(pa_pairs$ID))

min(table(pa_pairs$nonPA))
min(table(pa_pairs$ID))
max(table(pa_pairs$nonPA))
max(table(pa_pairs$ID))

pa_pairs %>% group_by(LC) %>% dplyr::select(ID) %>% unique() %>% count()
pa_pairs %>% group_by(LC) %>% dplyr::select(nonPA) %>% unique() %>% count()

#save(pa_pairs, file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_", Sys.Date(),".RData"))
write.csv(pa_pairs, file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_", Sys.Date(),".csv"), row.names=F)

