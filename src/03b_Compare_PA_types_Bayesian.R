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

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob <- data_glob %>% 
  mutate("PA_rank" = ifelse(is.na(PA_rank), 11, PA_rank)) %>%
  filter(LC %in% lc_names) %>%
  arrange(LC, PA_rank)
data_glob

# extract sample size and number and type of PA_ranks for each LC
protect_legend <- data_glob %>% 
  dplyr::select(LC, PA_rank, PA_type) %>% unique() %>% 
  group_by(LC) %>% arrange(LC, PA_rank) %>% #count() %>% #3 times per LC
  ungroup() %>%
  mutate("PA_rank_new" = rep(1:8, 3))
protect_legend
#unique(data_glob$PA_rank)

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

pars_list <- f_compare_pa_types(data = data_glob, 
                                protect_levels = protect_legend,
                                col_fns = fns)

save(pars_list, file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian_global.RData"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values

pars_sample <- f_combine_pars_list(pars_list = pars_list)
head(pars_sample)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables & effect sizes 
save(pars_sample, file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_global.RData"))

