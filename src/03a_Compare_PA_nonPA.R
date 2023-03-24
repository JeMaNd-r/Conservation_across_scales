#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#          author: Romy Zeiss               #
#            date: 2023-03-23               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(psych) # to calculate CI of Cohens d effect size

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load PA-nonPA pairs ####
pa_pairs <- read.csv(file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_global.csv"))
head(pa_pairs)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA sites ####
# using Cohen's D
d_list <- f_compare_pa_nonpa(data = data_glob,
                             data_pairs = pa_pairs,
                             col_id = "Order_ID",
                             col_fns = fns)
head(d_list)

## Save total df with effect sizes 
save(d_list,  file=paste0(here::here(), "/results/d_1000_trails_global.RData"))


