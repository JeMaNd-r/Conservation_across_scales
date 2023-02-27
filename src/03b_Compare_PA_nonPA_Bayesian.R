#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#         (Bayesian approach)               #
#          author: Romy Zeiss               #
#            date: 2023-02-27               #
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
data_glob

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load PA-nonPA pairs ####
pa_pairs <- read.csv(file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_2023-02-22.csv"))

# # add columns of functions for each PA and nonPA (wide format)
# pa_pairs <- merge(pa_pairs, data_glob[,c("Order_ID", fns)], by="Order_ID")
# pa_pairs <- merge(pa_pairs, data_glob[,c("Order_ID", fns)], by.x="nonPA", by.y="Order_ID",
#                   suffixes = c(".pa", ".nonpa"))

head(pa_pairs)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Compare difference between PA and nonPA sites ####
unique(pa_pairs$times)

# take one run (1:unique(pa_pairs$times)), estimate p and effect size (and Bayesian), take next
delta_list <- vector("list", length = length(lc_names))
names(delta_list) <- lc_names
delta_list <- lapply(delta_list, 
                     function(x) {
                       x <- vector("list", length = length(fns))
                       names(x) <- fns
                       
                       lapply(x, 
                              function(y) {
                                y <- vector("list", length = length(unique(pa_pairs$times)))
                                names(y) <- unique(pa_pairs$times)
                              }
                       )
                       }
                     )

# progress bar
progress_bar <- txtProgressBar(min = 0, max = length(unique(pa_pairs$times)), initial = 0,
                               style=3, title = "Processing") 

for(x in unique(pa_pairs$times)){
  
  setTxtProgressBar(progress_bar,x)
  #print("#########################")
  #print(paste0("Run number ", x))
  
  # subset data
  temp_pairs <- pa_pairs[pa_pairs$times==x,]
  
  # individual tests per LC type
  for(lc in lc_names){
    
    # subset data based on land cover type
    temp_PA <- temp_pairs %>% 
      left_join(data_glob[,c("Order_ID", "PA", fns)], by=c("Order_ID"="Order_ID")) %>% 
      filter(LC==lc) %>% 
      dplyr::select(-nonPA, -mahal.min)
    temp_nonPA <- temp_pairs %>% 
      left_join(data_glob[,c("Order_ID", "PA", fns)], by=c("nonPA"="Order_ID")) %>% 
      filter(LC==lc) %>% 
      mutate("Order_ID"=nonPA) %>% 
      dplyr::select(-nonPA, -mahal.min)
    
    for(no_fns in 1:(length(fns))){
      
      data <- list(n1 = nrow(temp_PA),
                   n2 = nrow(temp_nonPA),
                   y1 = temp_PA[,fns[no_fns]],
                   y2 = temp_nonPA[,fns[no_fns]])
      
      # unpaired t test
      stan_code = '
        data {
          int n1;
          vector[n1] y1;
          int n2;
          vector[n2] y2;
        }
        parameters {
          real mu1;
          real mu2;
          real<lower=0> sigma;
        }
        model {
          // priors
          mu1 ~ normal(0, 10);
          mu2 ~ normal(0, 10);
          sigma ~ normal(0, 10);
          // likelihood
          y1 ~ normal(mu1, sigma);
          y2 ~ normal(mu2, sigma);
        }
        generated quantities{
          real delta;
          delta = mu2-mu1;
        }
      '
      
      stan_model <- stan_model(model_code = stan_code)
      stan_fit <- sampling(stan_model, data = data,
                           chains = 4, iter = 4000, warmup = 1000,
                           show_messages = FALSE)
      #print(stan_fit, digits = 3, probs = c(0.025, 0.975))
      
      delta_list[[lc]][[no_fns]][[x]] <- extract(stan_fit, pars="delta") 
    }
    
  }
  
  close(progress_bar)
}

save(delta_list, file=paste0(here::here(), "/intermediates/delta_1000_trails_Bayesian.RData"))

# combine individual list elements (c) per fns & lc into one vector
# that is, we have one list element per fns and lc containing the values 
# from all 1000 runs = 4000*1000 estimates
delta_list[[lc]][[fns]] <- do.call(c, delta_list[[lc]][[fns]])

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables & effect sizes ####
save(delta_list, file=paste0(here::here(), "/results/delta_1000_trails_Bayesian.RData"))



