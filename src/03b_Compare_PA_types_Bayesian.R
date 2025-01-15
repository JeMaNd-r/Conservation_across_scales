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

#library(rstan)
#options(mc.cores = 4) # number of CPU cores

library(emmeans) # to estimate contrast i.e. EM means
library(brms)
library(modelr)
library(tidybayes)

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
#temp_scale <- "global"
# temp_scale <- "continental"
 temp_scale <- "regional"

data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

data_clean <- data_clean %>% 
  mutate("PA_rank" = ifelse(is.na(PA_rank), 11, PA_rank)) %>%
  filter(LC %in% lc_names) %>%
  arrange(LC, PA_rank)
data_clean

min_size <- min(table(data_clean$LC, 
                      data_clean$PA)[table(data_clean$LC, 
                                           data_clean$PA)
                                     >0])

if(temp_scale == "global"){
  lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
  } 
if(temp_scale == "continental"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  fns <- fns[fns != "Water_regulation_service"]
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  fns <- fns[fns != "Water_regulation_service"]
}

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

# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## Compare difference between PA types (incl. nonPA) ####
# # using Bayesian model - random intercept model
# 
# # define stan code for linear model
# stan_code = '
#   data {
#     int n;
#     int n_group;
#     real y[n];
#     int group[n];
#   }
#   parameters {
#     real a[n_group];
#     real<lower=0> sigma;
#     real mu_a;
#     real<lower=0> sigma_a;
#   }
#   model {
#     // priors
#     mu_a ~ normal(0,10);
#     sigma_a ~ normal(0,10);
#     
#     for (j in 1:n_group){
#      a[j] ~ normal(mu_a,sigma_a);
#     }
#     
#     sigma ~ normal(0,10);
#     
#     // likelihood
#     for(i in 1:n){
#       y[i] ~ normal(a[ group[i] ], sigma);
#     }
#   }
# '
# 
# stan_model <- stan_model(model_code = stan_code)
# 
# pars_list <- f_compare_pa_types(data = data_clean, 
#                                 protect_levels = protect_legend,
#                                 col_fns = fns)
# 
# save(pars_list, file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian_", temp_scale, ".RData"))
# 
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # combine individual list elements (c) per fns & lc into one vector
# # that is, we have one list element per fns and lc containing the values
# 
# load(file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian_", temp_scale, ".RData")) #pars_list
# 
# pars_sample <- f_combine_pars_list(pars_list = pars_list)
# str(pars_sample)
# 
# rm(pars_list)
# gc()
# 
# ### Save total list with p tables & effect sizes ####
# save(pars_sample, file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_", temp_scale, ".RData"))
# 
# Calculating emmeans and contrasts more complicated then in brms, therefore
# switch to brms (2025-01-13)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference using brms & linear regression model ####
# random-slope model & random-intercept, both using brms

source(paste0(here::here(), "/src/00_Parameters.R")) 

min_size <- min(table(data_clean$LC, 
                      data_clean$PA)[table(data_clean$LC, 
                                           data_clean$PA)
                                     >0])

if(temp_scale == "global"){
  lc_names <- "Dryland"
} 
if(temp_scale == "continental"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  fns <- fns[fns != "Water_regulation_service"]
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
  fns <- fns[fns != "Water_regulation_service"]
}

data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

data_clean <- data_clean %>% 
  mutate("PA_rank_rev" = ifelse(is.na(PA_rank), 1, 11-PA_rank+1)) %>%
  filter(LC %in% lc_names) %>%
  arrange(LC, PA_rank_rev)

# store results
fixed_effects_slope <- vector("list")
fixed_effects_intercept <- vector("list")
pred_list_slope <- vector("list")
pred_list_intercept <- vector("list")

for(temp_fns in fns){
  
  data_temp_slope <- data_clean %>% dplyr::select(all_of(c("LC", temp_fns, "PA_rank_rev")))
  data_temp_intercept <- data_clean %>% dplyr::select(all_of(c("LC", temp_fns, "PA_rank"))) %>%
    mutate(PA_rank = ifelse(is.na(PA_rank), "Unprotected", PA_rank))
  
  if(temp_scale == "global"){
    output_slope <- brm(brmsformula(paste(temp_fns, "~ PA_rank_rev")), data = data_temp_slope,
                        chains = 4, iter = 10000, warmup = 2000)
    
    output_intercept <- brm(brmsformula(paste(temp_fns, "~ PA_rank")), data = data_temp_intercept,
                        chains = 4, iter = 10000, warmup = 2000)
    
  }else{
    output_slope <- brm(brmsformula(paste(temp_fns, "~ LC * PA_rank_rev")), data = data_temp_slope,
                        chains = 4, iter = 10000, warmup = 2000)
    
    output_intercept <- brm(brmsformula(paste(temp_fns, "~ LC + PA_rank")), data = data_temp_intercept,
                        chains = 4, iter = 10000, warmup = 2000)
  }
  
  # sink(paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_", temp_fns, ".txt"))
  # output_slope
  # output_slope$fit
  # hypothesis(output_slope, "PA_rank_rev<0")
  # hypothesis(output_slope, "PA_rank_rev>0")
  # sink()
  
  fixed_effects_slope[[temp_fns]] <- vector("list")
  fixed_effects_intercept[[temp_fns]] <- vector("list")
  fixed_effects_slope[[temp_fns]][["fixef"]] <- brms::fixef(output_slope)
  fixed_effects_intercept[[temp_fns]][["fixef"]] <- brms::fixef(output_intercept)
  
  if(temp_scale == "global"){
    fixed_effects_slope[[temp_fns]][["emtrends"]] <- as_tibble(emmeans::emtrends(output_slope, var = "PA_rank_rev"))
    fixed_effects_slope[[temp_fns]][["emmeans"]] <- as_tibble(emmeans::emmeans(output_slope, specs = c("PA_rank_rev")))

    fixed_effects_intercept[[temp_fns]][["emmeans"]] <- as_tibble(emmeans::emmeans(output_intercept, specs = "PA_rank"))
  }else{
    fixed_effects_slope[[temp_fns]][["emtrends"]] <- as_tibble(emmeans::emtrends(output_slope, specs = "LC", var = "PA_rank_rev"))
    fixed_effects_slope[[temp_fns]][["emmeans"]] <- as_tibble(emmeans::emmeans(output_slope, specs = c("PA_rank_rev", "LC")))
    
    fixed_effects_intercept[[temp_fns]][["emmeans"]] <- as_tibble(emmeans::emmeans(output_intercept, specs = c("PA_rank", "LC")))
  }
  
  # Extract estimates and credible intervals for PA_rank_rev
  pred_list_slope[[temp_fns]] <- data_temp_slope %>%
    group_by(LC) %>%
    modelr::data_grid(PA_rank_rev = modelr::seq_range(PA_rank_rev, n = 51)) %>%
    tidybayes::add_epred_draws(output_slope) %>%
    mutate(scale = temp_scale,
           fns = temp_fns)

  pred_list_intercept[[temp_fns]] <- data_temp_intercept %>%
    group_by(LC) %>%
    tidybayes::add_epred_draws(output_intercept) %>%
    mutate(scale = temp_scale,
           fns = temp_fns)

}

pred_list_slope <- do.call(rbind, pred_list_slope)
pred_list_intercept <- do.call(rbind, pred_list_intercept)

save(pred_list_slope, file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, ".RData"))
save(pred_list_intercept, file=paste0(here::here(), "/intermediates/PAtypes_Bayesian_", temp_scale, ".RData"))

# subset
pred_sample_slope <- pred_list_slope %>% group_by(LC) %>% slice_sample(n = 10000)
pred_sample_intercept <- pred_list_intercept %>% group_by(LC) %>% slice_sample(n = 10000)

save(pred_sample_slope, file=paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_sample10k.RData"))
save(pred_sample_intercept, file=paste0(here::here(), "/results/PAtypes_Bayesian_", temp_scale, "_sample10k.RData"))

# save fixed effects & emmeans
save(fixed_effects_slope, file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, "_summary.RData"))
save(fixed_effects_intercept, file=paste0(here::here(), "/intermediates/PAtypes_Bayesian_", temp_scale, "_summary.RData"))

sink(paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, ".txt"))
fixed_effects_slope
sink()

sink(paste0(here::here(), "/intermediates/PAtypes_Bayesian_", temp_scale, ".txt"))
fixed_effects_intercept
sink()


# extract emtrends (slope model)
load(file=paste0(here::here(), "/intermediates/PAranks_Bayesian_", temp_scale, "_summary.RData")) #fixed_effects_slope
emtrends <- sapply(fixed_effects_slope,function(x) x[2])
for(i in 1:length(emtrends)){
  emtrends[[i]] <- emtrends[[i]] %>% 
    mutate("fns" = gsub(".emtrends", "", names(emtrends)[i]),
           "scale" = temp_scale)
}
emtrends <- do.call(rbind, emtrends)

if(temp_scale == "global") emtrends <- emtrends %>% mutate(LC = "Dryland", .before = 1) %>% dplyr::select(-PA_rank_rev)

write_csv(emtrends, file=paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_emtrends.csv"))


# extract emmeans (intercept model)
load(file=paste0(here::here(), "/intermediates/PAtypes_Bayesian_", temp_scale, "_summary.RData")) #fixed_effects_intercept
emtrends_mean <- sapply(fixed_effects_intercept,function(x) x[2])
for(i in 1:length(emtrends_mean)){
  emtrends_mean[[i]] <- emtrends_mean[[i]] %>% 
    mutate("fns" = gsub(".emmeans", "", names(emtrends_mean)[i]),
           "scale" = temp_scale)
}
emtrends_mean <- do.call(rbind, emtrends_mean)
emtrends_mean

if(temp_scale == "global") emtrends_mean <- emtrends_mean %>% mutate(LC = "Dryland", .before = 1) 

write_csv(emtrends_mean, file=paste0(here::here(), "/results/PAtypes_Bayesian_", temp_scale, "_emmeans.csv"))

