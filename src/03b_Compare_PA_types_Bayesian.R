#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#    incl. PA type (Bayesian approach)      #
#          author: Romy Zeiss               #
#            date: 2023-03-17               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(rstan)
options(mc.cores = 4) # number of CPU cores

source(paste0(here::here(), "/src/00_Parameters_functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob <- data_glob %>% 
  mutate("PA_rank" = ifelse(is.na(PA_rank), 11, PA_rank))
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA sites ####
# but include effect of PA_type

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

pars_list <- vector("list", length = length(lc_names))
names(pars_list) <- lc_names
pars_list <- lapply(pars_list, 
                     function(x) {
                       x <- vector("list", length = length(fns))
                       names(x) <- fns
                       
                       lapply(x, 
                              function(y) {
                                y <- vector("list", length = 1)
                              }
                       )
                     }
)

stan_data <- pars_list
stan_fit <- pars_list

# individual tests per LC type
for(lc in lc_names){
  
  print("#####################################")
  print(paste0("Run land cover ", lc))
  
  # subset data based on land cover type
  temp_data <- data_glob %>% filter(LC==lc)
  
  for(no_fns in 1:(length(fns))){
    
    stan_data[[lc]][[no_fns]] <- list(n = nrow(temp_data),
                                      y = c(unlist(temp_data[,fns[no_fns]])),
                                      group = temp_data$PA_rank,
                                      n_group = 11)
    
    stan_fit[[lc]][[no_fns]] <- sampling(stan_model, data = stan_data[[lc]][[no_fns]],
                              chains = 4, iter = 10000, warmup = 2000,
                              show_messages = TRUE)
    #print(stan_fit, digits = 3, probs = c(0.025, 0.975))
    
    pars_list[[lc]][[no_fns]]<- rstan::extract(stan_fit[[lc]][[no_fns]], pars=c(paste0("a[", 1:11, "]"), "sigma")) 
    
    stan_fit[[lc]][[no_fns]] <- NULL
    stan_data[[lc]][[no_fns]] <- NULL
    
  }
  
  gc()
  #unlink(file.path("tmp", "Rtmp*"), recursive = T)
  
  # remove temporary files
  dso_filenames <- dir(tempdir(), pattern=.Platform$dynlib.ext)
  filenames  <- dir(tempdir())
  for (i in seq(dso_filenames))
    try(dyn.unload(file.path(tempdir(), dso_filenames[i])))
  for (i in seq(filenames))
    if (file.exists(file.path(tempdir(), filenames[i])) & nchar(filenames[i]) < 42) # some files w/ long filenames that didn't like to be removeed
      file.remove(file.path(tempdir(), filenames[i]))
  
}

save(pars_list, file=paste0(here::here(), "/intermediates/pars_PAtypes_Bayesian.RData"))

# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values 

pars_sample <- lapply(pars_list, function (x){ #across lc types
  lapply(x, function(y){ #across functions
      bind_rows(y)
  })
})

pars_sample <- lapply(pars_sample, function (x){ #across lc types
    bind_rows(x, .id = "fns")
})

pars_sample <- bind_rows(pars_sample, .id="lc")
head(pars_sample)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables & effect sizes ####
save(pars_sample, file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df.RData"))

