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
 temp_scale <- "global"
# temp_scale <- "continental"
# temp_scale <- "regional"

# set date of latest analysis
if(temp_scale == "global") temp_date <- "2024-09-12"
if(temp_scale == "continental") temp_date <- "2024-08-01"
if(temp_scale == "regional") temp_date <- "2024-08-01"

if(temp_scale == "global"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
  min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
} 
if(temp_scale == "continental"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
  min_size <- 10 # number of samples/ sites that should be paired per LC type
  fns <- fns[fns != "Water_regulation_service"]
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
  min_size <- 7 # number of samples/ sites that should be paired per LC type
  fns <- fns[fns != "Water_regulation_service"]
}
data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load PA-nonPA pairs ####
pa_pairs <- read_csv(file=paste0(here::here(), "/intermediates/", temp_scale, "/Pairs_paNonpa_1000trails_", temp_date,".csv"))
head(pa_pairs)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA sites ####
# using Cohen's D
d_list <- f_compare_pa_nonpa(data = data_clean,
                             data_pairs = pa_pairs,
                             col_id = "SampleID",
                             col_fns = fns)
head(d_list)

## Save total df with effect sizes 
save(d_list,  file=paste0(here::here(), "/results/d_1000_trails_", temp_scale, ".RData"))


