#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#           Sensitivity analysis            #
#    extended global dataset (Bacteria)     #
#          author: Romy Zeiss               #
#            date: 2024-08-05               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# This script is mostly similar to the scripts named in the headers.
# If relevant, only parts of scripts on temp_scale == "global" was copied.
# Changes to these scripts are the following:
# - subset fns vector to only Bacteria fns
# - change directories to save output in 1 new folder
# - change one function where directory is written inside

library(here)
library(tidyverse)

library(rstan)
options(mc.cores = 4) # number of CPU cores

library(emmeans) # to estimate contrast i.e. EM means

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))


# subset functions
fns <- c("Bac_richness", "Bac_shannonDiv", "Bac_JaccDist_av")

# create directory for intermediate results
if(!dir.exists(paste0(here::here(), "/results/sensitivity_globalBacteria"))){
  dir.create(paste0(here::here(), "/results/sensitivity_globalBacteria"))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 01b_Prepare_soil_data.R ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

### Load soil biodiversity data ####
coord_glob <- readxl::read_xlsx(paste0(here::here(), "/data_raw/Global_Atlas_drylands_V2.xlsx"),
                                sheet = "Database")
coord_glob #551 
coord_glob %>% filter(!is.na(ID_sequencing_16S)) #338

data_glob <- read_csv(here::here("data_raw/Diversity_global.csv")) #342 

sort(unique(data_glob$SampleID)) #342
sort(unique(coord_glob$ID_sequencing_16S)) #338

data_glob <- coord_glob %>%
  dplyr::select("Plot_ID", "ID_sequencing_16S", "Latitude_c", "Longitude_c", 
                "Vegetation", "Soil_pH", "Soil_salinity", "Soil_fine_texture",
                "ORC", "BGL veg", "FOS veg", "MIN veg",
                #"WHC veg", 
                "Plant_cover",             
                "TON", "Soil_total_P", "AVP", "DIN", "DON") %>%
  mutate(ID_numb = str_replace(ID_sequencing_16S,"FM", "1")) %>% #ID_sequencing_16s goes to 249, FMxxx
  mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
  mutate(ID_numb = as.numeric(ID_numb)) %>%
  filter(!is.na(ID_numb)) %>%
  full_join(data_glob %>%
              mutate(ID_numb = str_replace(SampleID, "FM", "1")) %>%
              mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
              mutate(ID_numb = as.numeric(ID_numb)), 
            by = "ID_numb") %>%
  filter(!is.na(Latitude_c)) #missing coordinates are samples for other projects
data_glob #338
rm(coord_glob)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename & calculate variables (if needed) ####
data_glob <- data_glob %>%
  mutate("Pathogen_control" = 1 - (Plant_pathogen_richness / Fungi_richness))

data_glob <- data_glob %>%
  rename(Latitude = Latitude_c,
         Longitude = Longitude_c,
         Soil_texture = Soil_fine_texture)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add missing environmental covariates (elevation & annual climate) ####

data_glob <- f_extract_env(data = data_glob,
                           col_lon = "Longitude", col_lat = "Latitude")
head(data_glob)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Add location in or outside of protected area ####

## load information about protected areas
# based on World Database on Protected Areas (WDPA) (protectedplanet.net)
protect_glob <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

data_glob <- f_add_protect(data = data_glob, 
                           data_pa = protect_glob %>%
                             mutate(ID_numb = str_replace(ID_sequencing_16S, "FM", "1")) %>%
                             mutate(ID_numb = str_replace(ID_numb, "X", "")) %>%
                             mutate(ID_numb = as.numeric(ID_numb)), 
                           col_id="ID_numb")
rm(protect_glob)

nrow(data_glob %>% filter(PA==1)) #54
nrow(data_glob %>% filter(!is.na(PA_type))) #106
nrow(data_glob %>% filter(PA==0)) #285

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Rename land cover types in data ####
data_glob$LC <- data_glob$Vegetation
unique(data_glob$LC)
data_glob[data_glob$LC=="Forest" & !is.na(data_glob$LC), "LC"] <- "Woodland"
unique(data_glob$LC)

table(data_glob$LC)

# rename ID
data_glob <- data_glob %>% mutate(SampleID = ID_numb)

# functions
data_glob <- data_glob %>%
  filter(!is.na(SampleID)) %>%
  rename(Soil_carbon_service = ORC,
         Soil_stability_service = Plant_cover) %>% #,
  #Water_regulation_service = `WHC veg`) %>%
  mutate(OM_decomposition_service = rowSums(cbind(scale(`BGL veg`)[,1], 
                                                  scale(`FOS veg`)[,1],
                                                  scale(`MIN veg`)[,1]),
                                            na.rm=TRUE)/3,
         Nutrient_service = (scale(TON)[,1]+
                               scale(Soil_total_P)[,1]+ 
                               scale(AVP)[,1]+
                               scale(DIN)[,1]+
                               scale(DON)[,1])/6)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Keep complete cases only ####
data_glob <- data_glob[complete.cases(data_glob[,c(mahal_vars, fns, "LC")]),] #330 (97.6%)   
data_glob

# subset relevant columns
data_glob <- data_glob %>% 
  dplyr::select(all_of(c("SampleID", "LC", mahal_vars, fns, "PA", "PA_type", "PA_rank")))
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Save ####
write_csv(data_glob, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_global.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for colinearity between environmental variables ####

list_colinear <- f_colinearity(data = data_glob, 
                               col_lon = "Longitude", 
                               col_lat = "Latitude", 
                               vars_env = c("Soil_pH", "Soil_salinity", 
                                            "Soil_texture", "Elevation", 
                                            "AnnualTemp", "AnnualPrec"))
list_colinear
list_colinear$env_vif %>% filter(is.na(VIF))
#-> none excluded based on VIF

# Save
write.csv(list_colinear$env_vif, file=paste0(here::here(), "/results/sensitivity_globalBacteria/VIF_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatSpearman, paste0(here::here(), "/results/sensitivity_globalBacteria/corMatSpearman_envVars_global.csv"), row.names = F)
write.csv(list_colinear$corMatPearson, paste0(here::here(), "/results/sensitivity_globalBacteria/corMatPearson_envVars_global.csv"), row.names = F)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 02_Pairing.R ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####

data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_", temp_scale, ".csv"))
data_clean

## Explore data
summary(as.factor(data_clean$PA)) #G: 277 nonPA and 53 PAs, C: 743 vs. 61, R: 275 vs. 49
# number of observations (raw)
nrow(data_clean); nrow(data_clean[data_clean$PA,])  #G: 330 with 53 PAs, C: 807 vs. 64, R: 324 vs. 49

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Scale variables for mahalanobis distance ####

data_clean <- f_scale_vars(data = data_clean, vars = mahal_vars)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check for pairing ####
### CHANGED function (directory) ####
f_check_pairs <- function(data, col_id, col_lc, vars_z){
  
  # data:   dataframe containing following columns
  # col_id: name of column with site IDs
  # col_lc: column name with land-cover types
  # vars_z: vector with column names for mahalanobis variables
  
  nonpa <- data[data$PA==0,c(col_id, col_lc, "PA", vars_z)]
  pa <- data[data$PA==1,c(col_id, col_lc, "PA", vars_z)]
  pa_noPair <- c()
  all_nonPA <- nonpa[0,]
  
  for(i in 1:length(lc_names)){
    temp_nonPA  <- nonpa[nonpa[,col_lc]==lc_names[i],]
    temp_PA <- pa[pa[,col_lc]==lc_names[i],]
    
    sigma <- cov(temp_nonPA[,vars_z]) 
    for(j in 1:nrow(temp_PA)){
      mu <- as.numeric(temp_PA[j,vars_z])
      temp_nonPA[,as.character(temp_PA[j,col_id])] <- 
        mahalanobis(temp_nonPA[,vars_z], mu, sigma, tol=1e-30)
      #print(j)
    }
    # if continental data, remove nonPA with geographical distance >500km
    if(temp_scale == "continental"){
      for(k in as.character(pull(temp_PA[,col_id]))){
        temp_nonPA[!(as.character(pull(temp_nonPA[,col_id])) %in% list_dist[[k]]), k] <- NA
      }
    }
    
    min_mahal <- apply(X = temp_nonPA[, as.character(temp_PA %>% pull(col_id))], MARGIN = 2, FUN = min, na.rm = TRUE)
    pa_noPair <- rbind(pa_noPair, cbind(names(min_mahal[min_mahal>mahal_thres]), min_mahal[min_mahal>mahal_thres]))
    
    all_nonPA <- full_join(all_nonPA, temp_nonPA, by=c(col_id, col_lc, "PA", vars_z))
  }
  pa_noPair
  nrow(pa_noPair) #nrow=12
  
  unpaired_pa <- data[data[,col_id] %in% pa_noPair,] 
  write.csv(unpaired_pa, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Unpaired_protected_sites_", Sys.Date(), ".csv"), row.names = F)
  
  cat("#---------------------------------------------------", sep="\n")
  cat(paste0("Saved: csv file with unpaired sites is saved under: "), 
      paste0(here::here(), "/intermediates/", temp_scale, "/Unpaired_protected_sites_", Sys.Date(), ".csv"),
      sep="\n")
  cat("#---------------------------------------------------", sep="\n")
  
  # look at Mahalanoubis distance values for each nonPA (Order_ID) and PA (columns)
  #all_nonPA
  
  # count how many "options" exist for one PA
  count_nonPA <- data.frame("SampleID"=NA, "No_nonPA"=NA)
  for(i in colnames(all_nonPA)[!is.na(as.numeric(colnames(all_nonPA)))]){
    temp_column <- all_nonPA[,i]
    
    count_nonPA <- rbind(count_nonPA,
                         c(i, nrow(all_nonPA[which(temp_column<=mahal_thres),col_id]) ))
    
  }
  
  count_nonPA <- count_nonPA %>% 
    mutate("No_nonPA"=as.numeric(No_nonPA)) %>% 
    arrange(No_nonPA)
  
  print("Number of nonPA sites per PA (col_ID)")
  print(head(count_nonPA))
  
  rm(temp_PA, temp_nonPA, temp_column, unpaired_pa, pa_noPair, min_mahal, nonpa, pa)
  
  return(list(count_nonPA, all_nonPA)) #number of nonPA sites per PA (Order_ID)
}
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

count_nonPA <- f_check_pairs(data = data_clean, 
                             col_id = "SampleID", col_lc = "LC", 
                             vars_z = mahal_vars_z)
all_nonPA <- count_nonPA[[2]]
count_nonPA <- count_nonPA[[1]] #G: 7 <5, 7 <7, 14 <10; C: 16 without enough (10) nonPAs for pairing
#head(count_nonPA)

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

unpaired_pa <- read.csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Unpaired_protected_sites_", temp_date, ".csv"))
head(unpaired_pa) 

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with SampleID. They have to been removed.
data_clean <- data_clean[!(data_clean$SampleID %in% unpaired_pa$SampleID),] 
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #G: nrow=330 with 53 PAs, C: 807 vs. 64; R. 324 with 49

# Remove sites that can only be paired less than min_nonPA times
# start with something small, then check how many possible;
if(temp_scale == "global") min_nonPA <- 5
data_clean <- data_clean[!(data_clean$SampleID %in% count_nonPA[count_nonPA$No_nonPA < min_nonPA, "SampleID"]),]
nrow(data_clean); nrow(data_clean[data_clean$PA,]) #nrow=323 with 46 PAs; C: 791 vs. 48; R: 318 with 43
data_clean %>% group_by(LC, PA) %>% count()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Randomization ####

# check LC types and number of (protected) sites
table(data_clean$LC, data_clean$PA)

# based on number of sites per LC, exclude LC
# e.g. global: only 3 unprotected on 0 protected for Other & onyl 1 PA for Woodland & 0 for Cropland
# continental: no Shrubland protected, 28 unprotected; and 0 PA in Other
# regional: PA only min. 7 -> decrease minimum size number to 7, 
#           exclude Shrublands & Others to get it running (otherwise no complete pairing achieved) 
lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
# Note: could be 10 as minimum number of samples per LC = 10

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

#save(pa_pairs, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Pairs_paNonpa_1000trails_", Sys.Date(),".RData"))
write_csv(pa_pairs, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Pairs_paNonpa_1000trails_", Sys.Date(),".csv"))


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 03a_Compare_PA_nonPA.R ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

### Load soil biodiversity data ####
temp_scale <- "global"

# set date of latest analysis
temp_date <- "2024-09-12"
lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC

data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_", temp_scale, ".csv"))
data_clean

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load PA-nonPA pairs ####
pa_pairs <- read_csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Pairs_paNonpa_1000trails_", temp_date,".csv"))
head(pa_pairs)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare difference between PA and nonPA sites ####
# using Cohen's D
d_list <- f_compare_pa_nonpa(data = data_clean,
                             data_pairs = pa_pairs,
                             col_id = "SampleID",
                             col_fns = fns)
head(d_list)

### Save total df with effect sizes 
save(d_list,  file=paste0(here::here(), "/results/sensitivity_globalBacteria/d_1000_trails_", temp_scale, ".RData"))


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 03b_Compare_PA_types_Bayesian.R ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load soil biodiversity data ####
temp_scale <- "global"

# set date of latest analysis
temp_date <- "2024-09-12"
lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC

data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_", temp_scale, ".csv"))
data_clean

data_clean <- data_clean %>% 
  mutate("PA_rank" = ifelse(is.na(PA_rank), 11, PA_rank)) %>%
  filter(LC %in% lc_names) %>%
  arrange(LC, PA_rank)
data_clean

# extract sample size and number and type of PA_ranks for each LC
protect_legend <- data_clean  %>%
  dplyr::select(LC, PA_rank, PA_type) %>% unique() %>% 
  group_by(LC) %>% arrange(LC, PA_rank) %>% #count() %>% 
  ungroup() 
protect_legend <- split(protect_legend, f = protect_legend$LC)
protect_legend <- lapply(protect_legend, function(x) mutate(x, "PA_rank_new" = 1:nrow(x)))
protect_legend <- do.call(rbind, protect_legend)
protect_legend
#unique(data_clean$PA_rank)

# Reverse the PA_rank column
protect_legend$PA_rank_rev <- max(protect_legend$PA_rank) - protect_legend$PA_rank + 1

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare difference between PA types (incl. nonPA) ####
# using Bayesian model

# define stan code for linear model
stan_code = '
  data {
    int n;
    int n_group;
    real y[n];
    int group[n];
  }
  parameters {
    real a[n_group];
    real<lower=0> sigma;
    real mu_a;
    real<lower=0> sigma_a;
  }
  model {
    // priors
    mu_a ~ normal(0,10);
    sigma_a ~ normal(0,10);
    
    for (j in 1:n_group){
     a[j] ~ normal(mu_a,sigma_a);
    }
    
    sigma ~ normal(0,10);
    
    // likelihood
    for(i in 1:n){
      y[i] ~ normal(a[ group[i] ], sigma);
    }
  }
'

stan_model <- stan_model(model_code = stan_code)

pars_list <- f_compare_pa_types(data = data_clean, 
                                protect_levels = protect_legend,
                                col_fns = fns)

save(pars_list, file=paste0(here::here(), "/results/sensitivity_globalBacteria/pars_PAtypes_Bayesian_", temp_scale, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values

load(file=paste0(here::here(), "/results/sensitivity_globalBacteria/pars_PAtypes_Bayesian_", temp_scale, ".RData")) #pars_list

pars_sample <- f_combine_pars_list(pars_list = pars_list)
str(pars_sample)

rm(pars_list)
gc()

#### Save total list with p tables & effect sizes ####
save(pars_sample, file=paste0(here::here(), "/results/sensitivity_globalBacteria/pars_PAtypes_Bayesian_df_", temp_scale, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Calculate contrast between protection types ####
# 
# load(file=paste0(here::here(), "/results/sensitivity_globalBacteria/pars_PAtypes_Bayesian_df_", temp_scale, ".RData")) #pars_sample
# 
# # group protection type Ia-VI (1-7) together, and group Not... and Unprotected (8-11) together
# pars_sample %>% 
#   pivot_longer(cols="1":"2") %>%
#   mutate("PA" = ifelse(name %in% as.character(1:7), 1, 
#                        ifelse(name %in% as.character(8:11), 0, NA))) %>%
#   group_by(lc, fns, PA) %>%
#   summarize(across(value, list("mean"  = function(x) mean(x, na.rm=TRUE), 
#                                 "median" = function(x) median(x, na.rm=TRUE), 
#                                 "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
#                                 "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
#                                 "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
#                                 "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE))))
# 
# # Estimated marginal means
# emmeans::emmeans()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare difference using brms & linear regression model ####

library(brms)
library(modelr)
library(tidybayes)

source(paste0(here::here(), "/src/00_Parameters.R")) 

# subset functions
fns <- c("Bac_richness", "Bac_shannonDiv", "Bac_JaccDist_av")

# set date of latest analysis
temp_scale <- "global"
temp_date <- "2024-09-12"
lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC

data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_", temp_scale, ".csv"))
data_clean

data_clean <- data_clean %>% 
  mutate("PA_rank_rev" = ifelse(is.na(PA_rank), 1, 11-PA_rank+1)) %>%
  filter(LC %in% lc_names) %>%
  arrange(LC, PA_rank_rev)

# store results
fixed_effects <- vector("list")
pred_list <- vector("list")

for(temp_fns in fns){
  
  data_temp <- data_clean %>% dplyr::select(all_of(c("LC", temp_fns, "PA_rank_rev")))
  model_output <- brm(brmsformula(paste(temp_fns, "~ LC * PA_rank_rev")), data = data_temp,
                      chains = 4, iter = 10000, warmup = 2000)
  
  # sink(paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_", temp_fns, ".txt"))
  # model_output
  # model_output$fit
  # hypothesis(model_output, "PA_rank_rev<0")
  # hypothesis(model_output, "PA_rank_rev>0")
  # sink()
  
  fixed_effects[[temp_fns]] <- vector("list")
  fixed_effects[[temp_fns]][["fixef"]] <- brms::fixef(model_output)
  fixed_effects[[temp_fns]][["emtrends"]] <- as_tibble(emmeans::emtrends(model_output, specs = "LC", var = "PA_rank_rev"))
  fixed_effects[[temp_fns]][["emmeans"]] <- as_tibble(emmeans::emmeans(model_output, specs = c("PA_rank_rev", "LC")))
  
  # Extract estimates and credible intervals for PA_rank_rev
  pred_list[[temp_fns]] <- data_temp %>%
    group_by(LC) %>%
    modelr::data_grid(PA_rank_rev = modelr::seq_range(PA_rank_rev, n = 51)) %>%
    tidybayes::add_epred_draws(model_output) %>%
    mutate(scale = temp_scale,
           fns = temp_fns)
  
}

pred_list <- do.call(rbind, pred_list)

save(pred_list, file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, ".RData"))
pred_sample <- pred_list %>% group_by(LC) %>% slice_sample(n = 10000)

save(pred_sample, file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, "_sample10k.RData"))

save(fixed_effects, file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, "_summary.RData"))
sink(paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, ".txt"))
fixed_effects
sink()

# extract emtrends
load(file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, "_summary.RData")) #fixed_effects
emtrends <- sapply(fixed_effects,function(x) x[2])
for(i in 1:length(emtrends)){
  emtrends[[i]] <- emtrends[[i]] %>% 
    mutate("fns" = gsub(".emtrends", "", names(emtrends)[i]),
           "scale" = temp_scale)
}
emtrends <- do.call(rbind, emtrends)
write_csv(emtrends, file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_", temp_scale, "_emtrends.csv"))


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 04_Plotting.R ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(here)
library(tidyverse)
library(terra)
library(ggdist) #to plot distributions (Bayesian)
library(gridExtra) #to add table in plot
library(ggrepel) #to add text not overlapping (geom_text_repel)
library(ggtext) #to add icons as axis labels
library(ggh4x) # for free axes in facet_grid
library(corrplot)

temp_scale <- "global"

# load background map
world.inp <- map_data("world")
if(temp_scale != "global"){
  world.inp <- subset(world.inp, region %in% c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan",
                                               "Belarus", "Belgium", "Bosnia and Herzegovina", "Bulgaria",
                                               "Croatia", "Cyprus", "Czechia","Denmark","Estonia","Finland", 
                                               "France","Georgia", "Germany", "Greece","Hungary","Iceland", 
                                               "Ireland", "Italy","Kazakhstan", "Kosovo", "Latvia","Liechtenstein", 
                                               "Lithuania", "Luxembourg","Malta","Moldova","Monaco","Montenegro",
                                               "Macedonia", "Netherlands","Norway","Poland","Portugal","Romania",
                                               "Russia","San Marino","Serbia","Slovakia","Slovenia","Spain",
                                               "Sweden","Switzerland","Turkey","Ukraine","UK","Vatican"))
}


source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

# subset functions
fns <- c("Bac_richness", "Bac_shannonDiv", "Bac_JaccDist_av")

# set date of latest analysis
temp_date <- "2024-09-12"

lc_names <- lc_names[lc_names != "Other" & lc_names != "Cropland"]
min_size <- 5 # number of samples/ sites that should be paired per LC type = min. number of PA per LC
 
# define each land cover type
lc_names_all <- c( "Cropland", "Grassland", "Shrubland", "Woodland", "Other")

# define order of functions
labels_order <- c(
  "Decomposition (OM)", "Nutrient cycling", "Pathogen control", "Soil carbon", "Soil stability", "Water regulation",
  "Bacterial Richness", "Fungal Richness", "Invertebrate Richness", "Protist Richness",
  "AM fungi Richness", "EM fungi Richness",
  "Bacterial Shannon", "Fungal Shannon", "Invertebrate Shannon", "Protist Shannon",
  "Nematode Richness", "Decomposer Richness",
  "Bacterial Dissimilarity", "Fungal Dissimilarity", "Invertebrate Dissimilarity", "Protist Dissimilarity"
)

data_clean <- read_csv(paste0(here::here(), "/results/sensitivity_globalBacteria/Data_clean_", temp_scale, ".csv"))
data_clean

# load pairs of PA and nonPA
pa_pairs <- read_csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Pairs_paNonpa_1000trails_", temp_date, ".csv"))
head(pa_pairs)

# load effect sizes
load(file=paste0(here::here(), "/results/sensitivity_globalBacteria/d_1000_trails_", temp_scale, ".RData")) #d_list

# load Bayesian results from PA_type comparison
load(file=paste0(here::here(), "/results/sensitivity_globalBacteria/pars_PAtypes_Bayesian_df_", temp_scale, ".RData")) #pars_sample

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sampling locations ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FIGURE 1 - Maps ####
# extract list of sampling locations actually used in comparison
data_locations <- data_clean %>% 
  filter(SampleID %in% unique(pa_pairs$ID) | 
           SampleID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Longitude,Latitude,SampleID, PA, LC)
data_locations #G: nrow=190, C: 316, R: 161
write_csv(data_locations, file = paste0(here::here(), "/results/sensitivity_globalBacteria/Locations_", temp_scale, ".csv"))
nrow(data_locations %>% filter(PA==1)) #G: 46 PAs, C: 48, R: 36
nrow(data_locations %>% filter(PA==0)) #G: 144 PAs, C: 268, R: 125

# set limits for point maps
if(temp_scale == "global") temp_limits <- c(-180, 180, -180, 180)

ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey90", linewidth = 0.15) + #fill = "grey80", color="grey75"
  #xlim(-180, 180) +  ylim(-180, 180) + #global
  #xlim(-10, 35) +  ylim(35, 70) + #continental
  #xlim(-9.5, -6) +  ylim(40.5, 42.5) + #regional
  xlim(temp_limits[1], temp_limits[2])+
  ylim(temp_limits[3], temp_limits[4])+
  
  geom_point(data=data_locations, aes(x=Longitude, y=Latitude, 
                                      shape = as.character(PA), color=LC, 
                                      size = as.character(PA)),
             stroke = 1)+ #increase circle line width; G: 0.9, C+R:3
  scale_shape_manual(values = c("0" = 19, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" =1.5, "1" = 2.5))+ #G: 0.4,1.2, C+R:3,8
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000"))+ 
  coord_map()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="right",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill= "grey80"))
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_locations_", temp_scale,".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pairing ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Save some numbers
sink(paste0(here::here(), "/results/sensitivity_globalBacteria/Numbers_", temp_scale, ".txt"))

print(temp_scale)
# look what non-protected sites have (not) been paired to any PA
cat("Number of non-protected sites have never been paired to any PA:")
hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(data_clean[data_clean$PA==0,]$SampleID, pa_pairs$nonPA))  
# nonPA sites never used: G ..., C: 475

cat(paste0("On average, each protected and unprotected site was used in ",
           round(mean(table(pa_pairs$ID)),0), 
           " (SD=", round(sd(table(pa_pairs$ID)),0),
           ", min=", round(min(table(pa_pairs$ID)),0),
           ", max=", round(max(table(pa_pairs$ID)),0), ") and ",
           round(mean(table(pa_pairs$nonPA)),0),
           " (SD=", round(sd(table(pa_pairs$nonPA)),0),
           ", min=", round(min(table(pa_pairs$nonPA)),0),
           ", max=", round(max(table(pa_pairs$nonPA)),0), 
           ") of the 1000 randomizations, respectively."))

# cat("Mean number of nonPA sites per PA.")
# mean(table(pa_pairs$nonPA))
# 
# cat("SD number of nonPA sites per PA.")
# sd(table(pa_pairs$nonPA))
# 
# cat("Min/max number of nonPA sites per PA.")
# min(table(pa_pairs$nonPA))
# max(table(pa_pairs$nonPA))
# 
# cat("Mean/SD/min/max number of times individual PA were used.")
# mean(table(pa_pairs$ID))
# sd(table(pa_pairs$ID))
# min(table(pa_pairs$ID))
# max(table(pa_pairs$ID))

table(pa_pairs$ID)

cat("Number of sites (PA/ nonPA) per LC.")
pa_pairs %>% group_by(LC) %>% dplyr::select(ID) %>% unique() %>% count()
pa_pairs %>% group_by(LC) %>% dplyr::select(nonPA) %>% unique() %>% count()

sink()

## Plot pairs on map
data_pairs <- pa_pairs %>%
  full_join(data_clean %>%
              filter(SampleID %in% unique(pa_pairs$ID) | 
                       SampleID %in% unique(pa_pairs$nonPA)) %>%
              dplyr::select(Longitude,Latitude,SampleID, PA, LC),
            by = c("ID" = "SampleID", "LC")) 

# Merge data to get coordinates for each ID
data_merged <- pa_pairs %>%
  inner_join(data_clean, by = c("ID" = "SampleID")) %>%
  inner_join(data_clean, by = c("nonPA" = "SampleID"), suffix = c("_1", "_2"))

# Calculate frequency of appearance of pairs
pair_freq <- data_merged %>%
  group_by(ID, nonPA) %>%
  summarise(freq = n())

# Create plot with ggplot2
ggplot() +
  geom_segment(data = data_merged %>% unique(), aes(x = Longitude_1, y = Latitude_1, xend = Longitude_2, yend = Latitude_2), alpha = 0.5, color = "blue") +
  geom_point(data = data_clean, aes(x = Longitude, y = Latitude)) +
  #scale_size_continuous(name = "Frequency", guide = "legend") +
  labs(x = "Longitude", y = "Latitude", title = "Paths Between Points with Width Based on Frequency") +
  theme_minimal()
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Locations_connection_", temp_scale, ".pdf"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boxplot mahalanobis distance ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

ggplot(pa_pairs, aes(x = LC, y = mahal.min, fill = LC))+
  geom_boxplot()+
  geom_violin(alpha = 0.3, adjust = 0.3)+
  theme_bw() +
  labs(x="Land-cover type",y="Mahalanobis distance") +
  theme(axis.text.x=element_text(size=15),text = element_text(size=20),  
        legend.position = "none", axis.text.y = element_text(size=15), legend.title = element_blank())+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000")) #"gold3", "limegreen", "forestgreen"
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_boxplot_mahal.distance_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Effect size per LC type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
d_df <- do.call(rbind, d_list)
str(d_df)

d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Violin plot effect size per LC type ####
# 
# ggplot(data=d_df, aes(x=round(effect,2), y=Label))+
#   geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
#   geom_vline(aes(xintercept=-0.8, linetype = "-0.8"), color="grey60")+ #large effect
#   geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
#   geom_vline(aes(xintercept=-0.5, linetype = "-0.5"), color="grey60")+ #medium effect
#   geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+ #small effect
#   geom_vline(aes(xintercept=-0.2, linetype = "-0.2"), color="grey60")+
#   #geom_vline(aes(xintercept=0, linetype = "0"))+
#   scale_linetype_manual(values = c("-0.8" = "dotted",
#                                    "-0.5" = "dashed", "-0.2" = "solid",
#                                    "0.2" = "solid","0.5" = "dashed",
#                                    "0.8" = "dotted"))+
#   
#   geom_violin(aes(fill = Organism), scale="width", width=0.5)+ #scale by width of violin bounding box
#   stat_summary(fun = "mean", geom = "pointrange", color = "black")+
#   facet_wrap(vars(lc))+
#   scale_x_continuous(labels = function(x) format(x, nsmall = 0))+#round c-value axis label
#   scale_fill_viridis_d()+
#   xlab("Effect size")+
#   theme_bw() + # use a white background
#   theme(legend.position = "right", axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),  
#         axis.text.x = element_text(size=7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.spacing = unit(1, "lines"))
# ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_violin_d-value_global.png"),
#        plot = last_plot())

# save data for plot
write.csv(d_df, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_pointrange_d-value_", temp_scale, ".csv"), row.names = FALSE)

d_df <- read_csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_pointrange_d-value_", temp_scale, ".csv"))

d_summary <- d_df %>% 
  dplyr::select(-run) %>%
  #pivot_longer(cols = c(p_value, ci_lower, ci_upper, t_stats),
  #             names_to = "metric") %>%
  group_by(lc, fns, Label, Group_function, Organism) %>%
  summarize(across(effect, .fns = list("mean"=mean, "SD"=sd, "median"=median,
                                       "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
                                       "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
                                       "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
                                       "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE))))
d_summary
write.csv(d_summary, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_summary_", temp_scale, ".csv"))

# mean per lc type
d_summary <- read.csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_summary_", temp_scale, ".csv"))

d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median, effect_ci_2.5:effect_ci_97.5), 
                   function(x) mean(x)))
d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# mean per Group_function
d_summary %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x), na.rm=TRUE),
                        sd=function(x) sd(abs(x), na.rm=TRUE))))

# check for significant mean p-values
#d_summary %>% arrange(p_value_mean)

sink(file=paste0(here::here(), "/results/sensitivity_globalBacteria/D-values_", temp_scale, ".txt"))
cat("#################################################", sep="\n")
cat("##  Significant d-values from global analysis  ##", sep="\n")
cat(paste0("##  Sys.Date() ", Sys.Date(), "  ##"), sep="\n")
cat("#################################################", sep="\n")
print(d_summary %>% ungroup() %>%
        filter(abs(effect_mean) >= 0.2) %>%
        dplyr::select(lc, fns, starts_with("effect")) %>%
        arrange(desc(abs(effect_median))),
      #n=nrow(d_summary %>% filter(abs(effect_mean) >= 0.2))
) #for print() command
sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pointrange effect size per LC type ####
# with confidence intervals 
ggplot(data = d_df %>% 
         filter(lc %in% lc_names) %>%
         mutate(Label=factor(Label, levels = rev(fns_labels$Label))), 
       aes(y = effect, x = Label))+
  
  geom_hline(aes(yintercept=0.8, linetype = "-0.8 / 0.8"), color="grey60")+
  geom_hline(aes(yintercept=-0.8, linetype = "-0.8 / 0.8"), color="grey60")+ #large effect
  geom_hline(aes(yintercept=0.5, linetype = "-0.5 / 0.5"), color="grey60")+
  geom_hline(aes(yintercept=-0.5, linetype = "-0.5 / 0.5"), color="grey60")+ #medium effect
  geom_hline(aes(yintercept=0.2, linetype = "-0.2 / 0.2"), color="grey90")+ #small effect
  geom_hline(aes(yintercept=-0.2, linetype = "-0.2 / 0.2"), color="grey90")+
  #geom_hline(aes(yintercept=0, linetype = "0"))+
  scale_linetype_manual(values = c("-0.8 / 0.8" = "dotted",
                                   "-0.5 / 0.5" = "dashed", 
                                   "-0.2 / 0.2" = "solid"), name="Strength of effect")+
  annotate("rect", ymin = -0.2, ymax = 0.2, xmin=-Inf, xmax=Inf, fill = "grey95")+
  
  # fill background of Group_function (code order matters)
  #annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=Inf, fill = "grey", alpha=0.1)+ #chocolate4
  #annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=15+0.5, fill = "grey", alpha=0.1)+
  #annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=8+0.5, fill = "grey", alpha=0.1)+
  #annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=4+0.5, fill = "grey", alpha=0.1)+
  
  ggdist::stat_pointinterval(fatten_point=1.2, shape=21) +
  coord_flip()+
  ylab("Effect size")+
  facet_wrap(vars(lc))+
  theme_bw() + # use a white background
  theme(legend.position = "bottom", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_", temp_scale, ".png"),
       plot = last_plot(),
       width=5, height=4)

## one-sided
ggplot(data = d_summary %>%
         filter(lc %in% lc_names) %>%
         mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
                                         abs(effect_ci_17), 
                                         abs(effect_ci_83))),
       aes(y = abs(effect_ci_min66), x = abs(effect_median),
           color=as.factor(sign(effect_median))
       ))+
  
  # geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
  # geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
  # geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+
  # scale_linetype_manual(values = c("0.8" = "dotted",
  #                                  "0.5" = "dashed",
  #                                  "0.2" = "solid"),
  #                       name="Strength of effect")+
  #annotate("rect", ymin = 0, ymax = 0.2, xmin=0, xmax=0.2, fill = "grey", alpha=0.1)+ #chocolate4
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.2, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.5, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.8, fill = "grey", alpha=0.3)+
  
  # geom_pointrange(aes(xmin = abs(effect_ci_2.5), xmax = abs(effect_ci_97.5),
  #                     shape = Group_function))+
  # geom_linerange(aes(xmin = ifelse(sign(effect_ci_17) != sign(effect_ci_83),
  #                                  abs(effect_ci_17)), xmax = abs(effect_ci_83),
  #                     #shape = Group_function
  #                    ), linewidth=1, alpha = 0.2)+
  geom_point(aes(shape = Group_function), size = 5)+
  ggrepel::geom_text_repel(aes(label = Label, shape = Group_function), size = 3)+
  #geom_text(aes(label = Label), size = 1.5, nudge_y = 0.005)+
  
  xlab("Median effect")+ ylab("Lower CI (17%)")+
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8))+
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))+
  # scale_y_continuous(limits = c(0, 0.7),
  #                    expand = c(0,0), breaks = c(0, 0.25, 0.5))+
  #scale_color_manual(c(""))+
  
  facet_wrap(vars(lc), scales = "free")+
  theme_bw() + # use a white background
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        #axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        #panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_medianSD_", temp_scale, ".png"),
       plot = last_plot())

## heatmap
ggplot(data = d_summary %>%
         filter(lc %in% lc_names) %>%
         mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
                                         abs(effect_ci_17), 
                                         abs(effect_ci_83))) %>%
         mutate(effect_ci_min66f = cut(effect_ci_min66,
                                       breaks=c(0, 0.2, 0.5, 0.8, Inf),
                                       labels=c("ns", "small", "medium", "large"))) %>%
         filter(!is.na(effect_ci_min66)),
       aes(x = lc, y = Label, alpha=effect_ci_min66f, 
           fill=as.factor(sign(effect_median))))+
  
  geom_tile()+
  scale_fill_manual(values = c("-1" = "#fc8d59", "0" = "#ffffbf", "1" = "#91bfdb"),
                    name = "Direction of effect")+
  scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
                     name = "Minimum effect size (66% CI)")+
  theme_bw() + # use a white background
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        #axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_medianSD_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Heatmap all 3 scales ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### FIGURE 2 - Heatmap ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mean per lc type and all 3 scales
d_sum_glob <- read.csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_summary_global.csv"))
d_sum_cont <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_continental.csv"))
d_sum_regi <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_regional.csv"))

d_sum_all <- rbind(d_sum_glob %>% mutate("scale" = "global"), 
                   d_sum_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_sum_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_sum_glob, d_sum_cont, d_sum_regi)

d_sum_all %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median, effect_ci_2.5:effect_ci_97.5), 
                   function(x) mean(x)))
d_sum_all %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# mean per Group_function
d_sum_all %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x), na.rm=TRUE),
                        sd=function(x) sd(abs(x), na.rm=TRUE))))

## switch lc and scales
d_plot_all <- d_sum_all %>%
  filter(lc %in% lc_names_all) %>%
  dplyr::select(-Label) %>%
  right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
             by = c("fns" = "Function")) %>%
  # mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
  #                                 abs(effect_ci_17), 
  #                                 abs(effect_ci_83))) %>%
  # mutate(effect_ci_min66 = ifelse(sign(effect_ci_17)!= sign(effect_ci_83), 0.01, effect_ci_min66)) %>%
  # mutate(effect_ci_min66f = cut(effect_ci_min66,
  #                               breaks=c(0, 0.2, 0.5, 0.8, Inf),
  #                               labels=c("marginal", "small", "medium", "large"))) %>%
  mutate(effect_mean_f = cut(abs(effect_mean),
                             breaks=c(0, 0.2, 0.5, 0.8, Inf),
                             labels=c("marginal", "small", "medium", "large"))) %>%
  filter(!is.na(effect_mean)) %>% #!is.na(effect_ci_min66) & 
  full_join(expand.grid(scale = c("global", "continental", "regional"), 
                        lc = lc_names_all, 
                        Label = fns_labels$Label)) %>%
  mutate(scale = factor(scale, levels = rev(c("global", "continental", "regional"))),
         Label = factor(Label, levels = rev(fns_labels %>% arrange(Group_function, Label) %>% pull(Label))),
         lc = factor(lc, levels = c("Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
  filter(lc != "Other" & !is.na(Label)) %>% 
  mutate(effect_direction = as.factor(sign(effect_mean)))%>%
  mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                     ifelse(effect_direction=="1", "positive", "0"))) %>%
  #mutate(effect_direction_c = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c)) %>%
  mutate(Label = factor(Label, levels = labels_order)) %>%
  mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c),
         effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
  mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "ns")))

ggplot(data = d_plot_all,
       aes(x = lc, y = scale))+
  
  # geom_tile(aes(alpha=effect_ci_min66f, 
  #                fill=as.factor(sign(effect_median))))+ 
  geom_point(aes(size = effect_mean_f,
                 color = effect_direction_c, 
                 fill= effect_significance,
                 shape = effect_na))+
  facet_wrap(vars(Label), ncol=6, drop=FALSE)+
  scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "ns" = "white"),
                    name = "Direction of effect",
                    na.value = "black", drop = FALSE)+
  scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
                     name = "Direction of effect",
                     na.value = "black")+
  # scale_color_manual(values = c("negative" = "#fc8d59", "ns" = "black", "positive" = "#91bfdb"),
  #                   name = "Direction of effect",
  #                   na.value = "grey60")+
  scale_size_manual(values = c("marginal" = 2, "ns" = 5, "small" = 5, "medium" = 10, "large" = 15),
                    name = "Effect size",
                    na.value = 5)+
  # scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
  #                    name = "Effect size")+
  scale_shape_manual(values = c("not available" = 4),
                     name = "Missing data",
                     na.value = 21)+
  scale_x_discrete(labels = c(
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
  ))+
  
  scale_y_discrete(labels = c(
    "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
    "continental" = "<img src='figures/icon_location-black.png' width='30'>",
    "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
  ))+
  xlab("")+ylab("")+
  theme_bw() + # use a white background
  
  guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black"), 
                                                 shape = 21, size = 5)), #tell legend to use different point shape
         color = "none", #don't show legend
         shape = guide_legend(override.aes = list(size = 5)))+
  theme(#legend.position = c(0.96, -0.01),
    legend.position = "bottom", #c(0.96, -0.05),
    legend.justification = c(1, 0),
    legend.box = "horizontal",
    legend.direction = "vertical",
    #legend.key.size = unit(20, "pt"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    #axis.title.y =element_blank(),
    axis.text.y = ggtext::element_markdown(hjust = 0),
    axis.ticks = element_blank(),
    axis.text.x = ggtext::element_markdown(vjust = 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(fill="white", color = "white"), #chocolate4
    strip.text = element_text(color="black", size = 15, hjust = 0)) #white
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_meanCI_allScales_fns.png"),
       plot = last_plot(), 
       width = 4400, height = 3800, units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### Summarizing stats ####
# number significant effects
d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% nrow() #23
# number ns
d_plot_all %>% filter(effect_significance=="ns" & !is.na(effect_significance)) %>% nrow() #112

# number significant per lc
table(d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% dplyr::select(scale, lc))
table(d_plot_all %>% filter(!is.na(effect_significance)) %>% dplyr::select(scale, lc))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX S2 Pointrange plot grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

d_df_glob <- read_csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/Data_pointrange_d-value_global.csv"))
d_df_cont <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_continental.csv"))
d_df_regi <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_regional.csv"))

d_df_all <- rbind(d_df_glob %>% mutate("scale" = "global"), 
                  d_df_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_df_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_df_glob, d_df_cont, d_df_regi)

d_df_grouped <- d_df_all %>% filter(!is.na(lc)) %>%
  mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='30'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='30'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>", NA)))) %>%
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>", 
                                                    "<img src='figures/icon_location-black.png' width='30'>", 
                                                    "<img src='figures/icon_flag-Portugal.png' width='30'>"))) %>%
  mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                 ifelse(Group_function=="Diversity", "Shannon", Group_function))) %>%
  mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity")))
write_csv(d_df_grouped %>%
            group_by(Group_function, scale, lc) %>%
            summarize(across(effect, list(median = median, 
                                          ci2.5 = function(x) quantile(x, 0.025, na.rm = TRUE), 
                                          ci92.5 = function(x) quantile(x, 0.925, na.rm = TRUE)))) %>%
            mutate(across(c(effect_median, effect_ci2.5, effect_ci92.5), function(x) round(x, 3))) %>%
            mutate(ci_95 = paste0("[", effect_ci2.5, "; ", effect_ci92.5, "]")) %>%
            dplyr::select(-effect_ci2.5, -effect_ci92.5), 
          paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_medianCI_allScales_grouped.csv"))

ggplot(data = d_df_grouped,
       aes(fill = lc, color = lc, 
           y = scale_icon, x = effect))+
  
  geom_vline(aes(xintercept=0), color="black")+
  ggdist::stat_pointinterval(fatten_point=1.2, shape=21,
                             position=position_dodgejust(width=0.5)) +
  coord_flip()+
  facet_wrap(vars(Group_function), drop=FALSE)+
  #ylab("Effect size")+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000"))+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000"))+
  # scale_x_discrete(labels = c(
  #   "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
  #   "continental" = "<img src='figures/icon_location-black.png' width='30'>",
  #   "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
  # ))+
  scale_x_continuous(breaks = c(-2, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 2))+
  theme_void()+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.title.y =element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.text.y = element_text(size = 13, hjust = 1),  
        axis.text.x = ggtext::element_markdown(vjust = 1, margin = unit(c(1, 1, 20, 1), "pt")),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        strip.text = element_text(size = 20, hjust = 0),
        plot.background = element_rect(fill = "white", color = "white"))

ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_d-value_medianCI_allScales_grouped.png"),
       plot = last_plot(), 
       width = 2700, height = 2200,
       units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boxplots values ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# extract list of sampling locations actually used in comparison
data_values <- data_clean %>% 
  filter(SampleID %in% unique(pa_pairs$ID) | 
           SampleID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(SampleID, PA, LC, all_of(fns)) %>%
  mutate(across(all_of(fns), .fns = function(x) { (x - mean(x)) / sd(x)})) %>%
  
  pivot_longer(cols = c(fns), names_to = "fns") %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
data_values

data_values <- data_values %>% filter(LC %in% lc_names)

ggplot(data_values, aes(x = Label, y = value)) +
  geom_hline(yintercept = 0, lty="dashed", color="grey60") +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  
  xlab("")+ ylab("")+
  facet_wrap(vars(LC), ncol=1, nrow=4)+
  theme_bw() +
  theme(axis.text.x=element_text(size=5, angle=45, hjust=1),
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", legend.text = element_text(size=15), axis.text.y = element_text(size=15))
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_boxplot_estimates_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bayesian results (PA types) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### FIGURE (OLD) - Bayesian PA types ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# transform parameters to long format and assign labels
pars_long <- pars_sample %>% 
  dplyr::select(-sigma, -sigma_a, -mu_a) %>% 
  pivot_longer(cols = colnames(.)[!(colnames(.) %in% c("lc", "fns"))])%>%
  full_join(protect_type %>% 
              mutate(PA_rank = as.character(PA_rank)), 
            by=c("name" = "PA_rank"))  %>%
  mutate("Label_pa" = ifelse(name=="11", "Unprotected", PA_type)) %>%
  mutate("Label_pa" = factor(Label_pa, levels = rev(c(protect_type$PA_type, "Unprotected")))) %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label_fns" = factor(Label, fns_labels$Label))
head(pars_long)

# save summary (i.e., data from figure)
pars_summary <- pars_long %>% group_by(lc, fns) %>% 
  summarize(across(value, list("mean"  = function(x) mean(x, na.rm=TRUE), 
                               "median" = function(x) median(x, na.rm=TRUE), 
                               "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
                               "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
                               "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
                               "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE)))) %>%
  arrange(lc, fns)
pars_summary
write_csv(pars_summary, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_parsBayesian_summary_", temp_scale, ".csv"))

# extract sample size
n_table <- data_clean %>% filter(LC!="Other") %>%
  group_by(LC, PA, PA_type, PA_rank) %>% count() %>%
  pivot_wider(id_cols = c("PA", "PA_type", "PA_rank"), 
              names_from = LC, values_from = n) %>% 
  arrange(PA_rank) %>% ungroup() %>%
  dplyr::select(-PA_rank)
n_table
write_csv(n_table, file=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_parsBayesian_nTable_", temp_scale, ".csv"))

ggplot(pars_long %>% filter(!is.na(Label)) %>% #filter(!is.na(PA_type)) %>%
         # add number of sizes to plot
         rbind(data_clean %>% filter(LC!="Other") %>%
                 group_by(LC, PA_type) %>% count() %>%
                 rename(lc=LC, value=n) %>%
                 mutate(fns="Number of sites",
                        name=NA,
                        PA_protected=NA,
                        Label_pa=ifelse(is.na(PA_type), "Unprotected", PA_type),
                        Label="Number of sites",
                        Label_fns = Label,
                        Label_short = Label,
                        Group_function=NA,
                        Organism=NA) %>%
                 dplyr::select(colnames(pars_long))) %>%
         #filter(Label_fns != "Water regulation") %>%%>%
         filter(!is.na(Label_pa)) %>%
         mutate(Label_fns = factor(Label_fns, levels = c(labels_order, "Number of sites"))),
       aes(y=Label_pa, x=value, color=lc))+
  
  ## adapt for scale
  annotate("rect", ymin = -Inf, ymax = 4+0.5, xmin=-Inf, xmax=Inf, fill = "grey90", alpha=0.5)+ #global: 4, C: 3, R: 2
  annotate("rect", ymin = -Inf, ymax = 1+0.5, xmin=-Inf, xmax=Inf, fill = "grey85", alpha=0.5)+
  
  ggdist::stat_pointinterval(fatten_point=1, shape=3, 
                             position=position_dodgejust(width=0.5))+ 
  facet_wrap(vars(Label_fns), scales = "free_x", ncol=6, drop=FALSE)+
  
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000"), name="Habitat type")+
  ylab("")+ xlab("")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 13))
ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_parsBayesian_", temp_scale, ".png"),
       plot = last_plot(),
       width=12, height=10)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bayesian results (PA ranks/ levels) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

pred_list <- rbind(get(load(paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_global_sample10k.RData"))), 
                   get(load(paste0(here::here(), "/results/PAranks_Bayesian_continental_sample10k.RData")))) %>% 
  rbind(get(load(paste0(here::here(), "/results/PAranks_Bayesian_regional_sample10k.RData"))))

for(temp_scale in c("global", "continental", "regional")){
  ggplot(data = pred_list %>% filter(scale == temp_scale) %>%
           right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
                      by = c("fns" = "Function")) %>%
           mutate(Label = factor(Label, levels = labels_order)), 
         
         aes(x = PA_rank_rev, y = .epred, color = ordered(LC))) +
    stat_lineribbon() +
    #geom_line(aes(y = .epred, group = paste(LC, .draw)), alpha = 0.2)+
    #geom_point(data = pred_list) +
    facet_wrap(vars(Label), scales = "free_y", ncol=6)+
    scale_fill_brewer(palette = "Greys") +
    scale_color_manual(values=c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000"), name="Habitat type")+
    scale_x_continuous(limits = c(1, 10), breaks = c(2, 10), minor_breaks = c(2,4,6,8, 10))+
    theme_void()+
    theme(axis.text = element_text(),
          panel.grid.major.y = element_line(color = "grey"),
          panel.grid.minor.x =  element_line(color = "grey"),
          strip.text = element_text(size = 15, hjust=0),
          legend.position = c(0.8, 0.1),
          legend.box = "horizontal")
  ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_regressions_parsBayesian_", temp_scale,".png"),
         plot = last_plot(),
         width=15, height=10)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bayesian pointrange grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
pars_glob <- read_csv(file=paste0(here::here(), "/results/sensitivity_globalBacteria/PAranks_Bayesian_global_emtrends.csv"))
pars_cont <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_continental_emtrends.csv"))
pars_regi <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_regional_emtrends.csv"))

pars_all <- rbind(pars_glob, pars_cont) %>%
  rbind(pars_regi)  %>%
  as_tibble()
rm(pars_glob, pars_cont, pars_regi)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### Bayesian summary ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combined diversity and function estimates (x axis)
# scale on y axis
# habitat types as colors
# pointrange plot

# summarize into groups
pars_sum <- pars_all %>%
  full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
  dplyr::select(-fns, -lower.HPD, -upper.HPD) %>%
  group_by(Group_function, LC, scale) %>%
  summarize(across(everything(), list("mean" = function(x) mean(x, na.rm=T),
                                      "SE" = function(x) sd(x, na.rm=T) / sqrt(length(x))))) %>%
  mutate(# Calculate confidence intervals
    "PA_rank_rev.trend_CI_lower" = PA_rank_rev.trend_mean - (1.96 * PA_rank_rev.trend_SE),  # 1.96 is the Z-value for a 95% confidence interval
    "PA_rank_rev.trend_CI_upper" = PA_rank_rev.trend_mean + (1.96 * PA_rank_rev.trend_SE)
  )

# plot
ggplot(pars_all %>%
         full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
         mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                    ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                           ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
         mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                        ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                           "<img src='figures/icon_location-black.png' width='70'>",
                                                           "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
         mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
         mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
         mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
         filter(!is.na(LC) & !is.na(scale)),
       aes(x = PA_rank_rev.trend, y = LC,
           xmin = lower.HPD, xmax = upper.HPD,
           fill = significance, color = LC,
           linetype = significance))+
  
  geom_vline(aes(xintercept=0), color="black")+
  geom_pointrange(position = position_jitter(height = 0.3),
                  size = 0.5, shape = 21) +
  
  geom_pointrange(data = pars_sum %>% filter(!is.na(LC) & !is.na(scale)) %>%
                    mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
                    mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                               ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                                      ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
                    mutate(Group_function = ifelse(Group_function=="Service", "Function",
                                                   ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                           scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                                      "<img src='figures/icon_location-black.png' width='70'>",
                                                                      "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
                    mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
                    mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                           significance =  ifelse(sign(PA_rank_rev.trend_CI_lower)!= sign(PA_rank_rev.trend_CI_upper), "ns", as.character(as.factor(LC)))) %>%
                    mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns"))),
                  
                  aes(fill = significance, color = LC,
                      y = LC, x = PA_rank_rev.trend_mean, 
                      xmin = PA_rank_rev.trend_CI_lower, xmax = PA_rank_rev.trend_CI_upper),
                  
                  position = position_dodge(width = 0.1),
                  size = 1.5, shape = 21, linewidth = 1.5) +
  
  #coord_flip()+
  coord_cartesian(clip = "off")+
  ggh4x::facet_grid2(scale_icon ~ Group_function, drop=FALSE, 
                     scales = "free", independent = "all", switch = "y",
                     shrink = FALSE)+
  
  #ylab("Effect size")+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000",
                             "ns" = "white"),
                    drop = FALSE)+
  scale_linetype_manual(values=c("Cropland" = "solid",
                                 "Grassland" = "solid",
                                 "Shrubland" = "solid", 
                                 "Woodland" = "solid", 
                                 "Other" = "solid",
                                 "ns" = "longdash"),
                        drop = FALSE)+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000",
                              "ns" = "#000000"),
                     drop = FALSE)+
  theme_void()+
  
  guides(fill = guide_legend(reverse = F, override.aes = list(shape = 21, 
                                                              color = c("Cropland" = "#4A2040",
                                                                        "Grassland" = "#E69F00",
                                                                        "Shrubland" = "#0072B2", 
                                                                        "Woodland" = "#009E73", 
                                                                        "ns" = "#000000"))), 
         color = "none")+
  theme(legend.position = "bottom", 
        axis.title.y =element_blank(),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(0.5, "cm"),
        legend.box.spacing = unit(c(1,0,0,0), "cm"),
        legend.title = element_blank(),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        axis.text.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.spacing = unit(0.8, "cm"),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        #strip.text.x.bottom = ggtext::element_markdown(vjust = 1),
        strip.text.x = element_text(size = 30, hjust = 0, vjust = 1),
        strip.text.y = ggtext::element_markdown(vjust = 0.5))

ggsave(filename=paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_BayesianTrends_allScales_grouped.png"),
       plot = last_plot(), 
       width = 5000, height = 4000,
       units = "px")

# table all
pars_all %>%
  full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
  mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                 ifelse(Group_function=="Diversity", "Shannon", Group_function)),
         scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                    "<img src='figures/icon_location-black.png' width='70'>",
                                                    "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
  mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
  mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
         significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
  mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
  filter(!is.na(LC) & !is.na(scale)) %>%
  filter(significance != "ns") %>%
  arrange(abs(PA_rank_rev.trend)) %>%
  print(n=100)

### APPENDIX S3 - Table with stats ####
options("scipen"=100, "digits"=4)
write_csv(pars_sum %>% 
            mutate("trend" = round(PA_rank_rev.trend_mean, 4),
                   "SE" = round(PA_rank_rev.trend_SE, 4),
                   "95% CI" = paste0("[", round(PA_rank_rev.trend_CI_lower, 4),"; ", 
                                     round(PA_rank_rev.trend_CI_upper, 4), "]"),
                   "Scale" = factor(scale, levels = c("global", "continental", "regional"))) %>%
            dplyr::select(Group_function, Scale, LC, trend, SE, "95% CI") %>%
            arrange(Group_function, Scale, LC),
          paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_BayesianTrends_allScales_grouped.csv"))

write_csv(pars_all %>%
            full_join(fns_labels %>% dplyr::select(Function, Group_function, Label_short), by = c("fns" = "Function")) %>%
            mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                       ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                              ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
            mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                           ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                   scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                              "<img src='figures/icon_location-black.png' width='70'>",
                                                              "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
            mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
            mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                   significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
            mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
            filter(!is.na(LC) & !is.na(scale)) %>%
            arrange(abs(PA_rank_rev.trend)) %>%
            mutate("Habitat" = LC,
                   "Group" = Group_function,
                   "Slope [HPD]" = paste0(round(PA_rank_rev.trend, 3), 
                                          " [", round(lower.HPD, 3), "; ", round(upper.HPD, 3), "] ", 
                                          ifelse(significance == "ns", "ns", "*")),
                   "Variable" = Label_short,
                   "Scale_long" = factor(scale, levels = c("global", "continental", "regional")))  %>%
            mutate("Scale" = factor(substr(scale, 1, 4), levels = c("glob", "cont", "regi"))) %>%
            dplyr::select(Group, Variable, Scale, Habitat, "Slope [HPD]") %>%
            pivot_wider(names_from = "Habitat", values_from = "Slope [HPD]") %>%
            arrange(Group, Variable, Scale),
          paste0(here::here(), "/results/sensitivity_globalBacteria/Results_pointrange_BayesianTrends_allScales.csv"))

