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
## Compare difference between PA types (incl. nonPA) ####
# using Bayesian model - random intercept model

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
## Compare difference using brms & linear regression model ####
# random-slope model

library(brms)
library(modelr)
library(tidybayes)

for(temp_scale in c("global")){ #, "continental", "regional"
  source(paste0(here::here(), "/src/00_Parameters.R")) 

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
  
  save(pred_list, file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, ".RData"))
  pred_sample <- pred_list %>% group_by(LC) %>% slice_sample(n = 10000)
  
  save(pred_sample, file=paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_sample10k.RData"))
  
  save(fixed_effects, file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, "_summary.RData"))
  sink(paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, ".txt"))
  fixed_effects
  sink()
}

# extract emtrends
for(temp_scale in c("global")){ #, "continental", "regional"
  load(file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, "_summary.RData")) #fixed_effects
  emtrends <- sapply(fixed_effects,function(x) x[2])
  for(i in 1:length(emtrends)){
    emtrends[[i]] <- emtrends[[i]] %>% 
      mutate("fns" = gsub(".emtrends", "", names(emtrends)[i]),
             "scale" = temp_scale)
  }
  emtrends <- do.call(rbind, emtrends)
  write_csv(emtrends, file=paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_emtrends.csv"))
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA (Bayesian) ####
# random intercept model but with 2 groups only

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

pars_list <- f_compare_pa_dummy(data = data_clean, 
                                col_fns = fns)
# 2 levels (0 = nonPA = coded as 1 and PA = 1 = coded as 2)

save(pars_list, file=paste0(here::here(), "/intermediates/pars_PAdummy_Bayesian_", temp_scale, ".RData"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values

load(file=paste0(here::here(), "/intermediates/pars_PAdummy_Bayesian_", temp_scale, ".RData")) #pars_list

pars_sample <- f_combine_pars_list(pars_list = pars_list)
str(pars_sample)

rm(pars_list)
gc()

### Save total list with p tables & effect sizes ####
save(pars_sample, file=paste0(here::here(), "/results/pars_PAdummy_Bayesian_df_", temp_scale, ".RData"))
