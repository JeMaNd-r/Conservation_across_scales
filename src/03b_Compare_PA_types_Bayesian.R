#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#    incl. PA type (Bayesian approach)      #
#          author: Romy Zeiss               #
#            date: 2023-03-24               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(rstan)
options(mc.cores = 4) # number of CPU cores

library(emmeans) # to estimate contrast i.e. EM means

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
# temp_scale <- "global"
# temp_scale <- "continental"
temp_scale <- "regional"

if(temp_scale == "global") lc_names <- lc_names[lc_names != "Other"]

# set date of latest analysis
if(temp_scale == "global") temp_date <- "2023-12-01"
if(temp_scale == "continental") temp_date <- "2023-12-04"
if(temp_scale == "regional") temp_date <- "2023-12-05"

data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA types (incl. nonPA) ####
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

save(pars_list, file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian_", temp_scale, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values

load(file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian_", temp_scale, ".RData")) #pars_list

pars_sample <- f_combine_pars_list(pars_list = pars_list)
str(pars_sample)

rm(pars_list)
gc()

### Save total list with p tables & effect sizes ####
save(pars_sample, file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_", temp_scale, ".RData"))

# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ### Calculate contrast between protection types ####
# 
# load(file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_", temp_scale, ".RData")) #pars_sample
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
## Compare difference using brms ####

# library(brms)

