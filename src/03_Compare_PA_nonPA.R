#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#   Compare protected with unprotected      #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)

library(psych)   # to calculate CI of Cohens d effect size

source("src/00_Parameters_functions.R")

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
# take one run (1:1000), estimate p and effect size (and Bayesian), take next
p_list <- vector("list", length = 1000)
d_list <- vector("list", length = 1000)
lapply(unique(pa_pairs$times), function(x){
  
  # subset data
  temp_pairs <- pa_pairs[pa_pairs$times==x,]
  
  # Perform tests
  p_list[[x]] <- list()
  d_list[[x]] <- list()
  
  # individual tests per LC type
  for(lc in lc_names){
    temp_p[[lc]] <- vector("list", length(fns))
    names(temp_p[[lc]]) <- fns
    
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
    
    temp_d <- psych::cohen.d(rbind(temp_PA, temp_nonPA)[,c("PA",fns)], "PA")
    d_list[[x]][[lc]] <- data.frame(lc=lc, 
                                    data.frame(temp_d$cohen.d),
                                    fns = rownames(data.frame(temp_d$cohen.d)),
                                    run=x,
                                    row.names = NULL)
    
    for(no_fns in 1:(length(fns))){
      # unpaired t test
      temp_p <- t.test(temp_PA[,fns[no_fns]], 
                       temp_nonPA[,fns[no_fns]])[c("p.value", 
                                                   "conf.int", # 95% confidence intervals
                                                   #"estimate", # mean of each group (1. PA, 2. nonPA)
                                                   "statistic")] # t statistic
      
      p_list[[x]][[lc]][[no_fns]] <- data.frame(lc=lc, 
                                           fns = fns[no_fns],
                                           p_value = temp_p$p.value,
                                           ci_lower = temp_p$conf.int[1],
                                           ci_upper = temp_p$conf.int[2],
                                           t_stats = temp_p$statistic,
                                           run=x,
                                           row.names = NULL) 
    }
    #print(sum(temp_pairs[temp_pairs$LC==lc,"nonPA"] == temp_nonPA[,"Order_ID"])==min_size) # make sure that the sites are pairing properly
    # print should give min_size (TRUE = fitting pairs) for each LC types * 1000 runs
    # You may wanna run this is a for loop to see the output of each print statement.
    
    p_list[[x]][[lc]] <- do.call(rbind, p_list[[x]][[lc]])
    
  }
  
  # combine individual list elements (df) into one df
  d_list[[x]] <- do.call(rbind, d_list[[x]])
  p_list[[x]] <- do.call(rbind, p_list[[x]])
  
  # remove rownames
  rownames(d_list[[x]]) <- NULL
  rownames(p_list[[x]]) <- NULL
})



  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Compare difference between PA and nonPA sites ####
  
  # merge pairs of PA and nonPA in one table ####
  nonpa.subset <- nonpa[unique(nonpa$Order_ID) %in% pa$nonPA,]
  #data_glob.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$Order_ID),])
  #head(data_glob.paired)
  
  # Perform t tests
  p.list <- list()
  effect_size_d[[times]] <- list()
  for(l in lc_names){
    p.list[[l]] <- vector("list", length(fns))
    names(p.list[[l]]) <- fns
    
    # subset data based on land cover type
    temp.PA <- pa[pa$LC==l & pa$PA==1,]
    temp.nonPA <- nonpa.subset[nonpa.subset$LC==l & nonpa.subset$PA==0,]
    #NEEDS TO BE FIXED: temp.nonPA <- full_join(temp.PA, temp.nonPA, by=c("nonPA"="Order_ID"), suffix)
    
    temp.cohens <- psych::cohen.d(rbind(temp.PA, temp.nonPA)[,c("PA",fns)], "PA")
    effect_size_d[[times]][[l]] <- cbind(lc=l, data.frame(temp.cohens$cohen.d), run=times)
    
    for(no.fns in 1:(length(fns))){
      # unpaired t test
      p.list[[l]][[no.fns]] <- t.test(temp.PA[,fns[no.fns]],temp.nonPA[,fns[no.fns]])
    }
    print(summary(temp.PA[,"nonPA"] == temp.nonPA[temp.nonPA$Order_ID %in% as.numeric(unlist(temp.PA[,"nonPA"])),"Order_ID"])) # make sure that the sites are pairing properly
    
  }
  effect_size_d[[times]] <- do.call(rbind, effect_size_d[[times]])
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Make table of all p values for all functions and LC types
  p_table <- as.data.frame(matrix(ncol=length(fns),nrow=5))
  colnames(p_table) <- fns
  rownames(p_table) <- c("Other","Shrubland","Grassland", "Woodland","total")
  
  for(n in 1:ncol(p_table)){
    for(lc in 1:length(lc_names)){
      p_table[which(rownames(p_table)==lc_names[lc]),n] <- p.list[[lc_names[lc]]][[n]]["p.value"]
    }
  }
  
  # add to overall result list
  p_list_total[[times]] <- p_table  #p-values from Chi-squared test
})
  


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables & effect sizes ####
#save(p.list.total, file="p_1000_trails.RData")
#write.csv(p.list.total, file="p_1000_trails.csv")
save(effect_size_d,  file="d_1000_trails.RData")

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Extract Cohen's d effect size for each run ####
d.table = as.data.frame(matrix(ncol=length(fns)+1,nrow=4000))
colnames(d.table) = c("lc", fns)
d.table$lc <- rep(lc_names, each=1000)
d.table$run <- rep(1:1000, 4)
#rownames(d.table) = c("Cropland","Grassland", "Woodland","Other","total")

dat <- effect.size
#dat <- effect.size.t   #for d based on t-value
#dat <- effect.size.g   #for Hedge's g (less biased d)
dat <- effect.size.upp  #lower confidence interval
dat <- effect.size.low  #upper confidence interval

for(n in fns){
  for(i in 1:1000){try({
    d.table[d.table$lc=="Cropland" & d.table$run==i,n] = dat[[i]][1, n]
    d.table[d.table$lc=="Grassland" & d.table$run==i,n] = dat[[i]][3, n]
    d.table[d.table$lc=="Woodland" & d.table$run==i,n] = dat[[i]][2, n]
    d.table[d.table$lc=="Other" & d.table$run==i,n] = dat[[i]][4, n]
    #... total missing
  })}
}
head(d.table)
save(d.table, file="Effect_size_d_table.RData")
#save(d.table, file="Effect_size_dt_table.RData") #for dat = effect.size.t
#save(d.table, file="Effect_size_d_CI_upper_table.RData")
#save(d.table, file="Effect_size_d_CI_lower_table.RData")
#save(d.table, file="Effect_size_gs_table.RData") #if dat = effect.size.g
#write.csv(d.table, file="Effect_size_d_table.csv", row.names = F)

#load(file="Effect_size_d_table.RData") #d.table

## Summarize per land use and function (i.e. get mean and confidence intervals)
d.table.mean <- aggregate(.~lc, d.table, mean)
d.table.min <- aggregate(.~lc, d.table, mean)
d.table.max <- aggregate(.~lc, d.table, mean)
#d.table.min <- aggregate(.~lc, d.table, function(x){quantile(x,0.25)})
#d.table.max <- aggregate(.~lc, d.table, function(x){quantile(x,0.75)})
d.table.mean
#write.csv(d.table.mean, file="Effect_size_d_mean.csv", row.names = F)
#write.csv(d.table.min, file="Effect_size_d_CI_low_mean.csv", row.names = F)
d.table.mean <- read.csv("Effect_size_d_mean.csv")

d.table.mean.long <- pivot_longer(d.table.mean, cols=2:11)
d.table.max.long <- pivot_longer(d.table.max, cols=2:11)
d.table.min.long <- pivot_longer(d.table.min, cols=2:11)

# Merge mean, max and min into one data frame
colnames(d.table.mean.long)[3] <- "mean.D"
colnames(d.table.max.long)[3] <- "CI.upper"
colnames(d.table.min.long)[3] <- "CI.lower"

d.table.summary <- merge(d.table.mean.long, d.table.max.long, 
                         by=c("name", "lc"))
d.table.summary <- merge(d.table.summary, d.table.min.long, by=c("name", "lc"))
d.table.summary <- d.table.summary[d.table.summary$name!="run",]
head(d.table.summary)

#write.csv(d.table.summary, file="Effect_size_d_summary.csv", row.names = F)

# # 3. plot the effect size between PA and NonPA
# #setwd(figu.wd); pdf("LUCAS_Effect_size_CohensD.pdf")
# ggplot() +
#   geom_pointrange(data=effect.size[effect.size$lc!="Other",], aes(x=funct, y=eff.mean, ymin=eff.low, ymax=eff.upp, 
#                    color=lc, shape=lc), position=position_dodge(width=0.7)) + 
#   geom_hline(yintercept=0, lty=1) +  # add a dotted line at x=0 after flip
#   geom_hline(yintercept=c(-0.2,0.2), lty=2) + # add line for significant difference
#   geom_hline(yintercept=c(-0.5,0.5), lty=3) +
#   #geom_text(aes(y=rep(c(0.7,1.1),48), label=c(rbind(d2$df, d3$df))), cex=2) +
#   #geom_text(aes(y=-2.3, label=rep(d1$df, each=2)), cex=2) +
#   #facet_grid(cols=vars(deca)) + 
#   coord_flip() +  # flip coordinates (puts labels on y axis)
#   xlab("Label") + ylab("Cohens d") +
#   #scale_color_manual(values=c()) +
#   #scale_shape_manual(values=c(1:3, 5)) +
#   scale_linetype_manual(values=c("solid", "longdash"))+
#   theme_bw() # use a white background
# dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combine different results in one table
setwd(data.wd); load("p_1000_trails.RData") #p.list.total
head(p.list.total)

#pvals <- vector("list", length=4)
pvals.agri <- data.frame(run=1:1000); pvals.agri[,fns] <- NA
pvals.gras <- data.frame(run=1:1000); pvals.gras[,fns] <- NA
pvals.fore <- data.frame(run=1:1000); pvals.fore[,fns] <- NA
pvals.othe <- data.frame(run=1:1000); pvals.othe[,fns] <- NA

for(i in 1:length(p.list.total)){
  try({
    pvals.agri[i,-1] <- p.list.total[[i]]["Cropland",]
    pvals.gras[i,-1] <- p.list.total[[i]]["Grassland",]
    pvals.fore[i,-1] <- p.list.total[[i]]["Woodland",]
    pvals.othe[i,-1] <- p.list.total[[i]]["Other",]
  })
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Plot p-values of all runs ####
pvals.agri <- melt(pvals.agri, id.vars=1, na.rm=T)
pvals.gras <- melt(pvals.gras, id.vars=1, na.rm=T)
pvals.fore <- melt(pvals.fore, id.vars=1, na.rm=T)
pvals.othe <- melt(pvals.othe, id.vars=1, na.rm=T)

#p.table$lc <- rownames(p.table) 
#pvals.all2 <- melt(p.table, id.vars="lc", na.rm=T)

## Extract mean p-values
p.vals.mean <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.agri, mean), lc="Cropland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.gras, mean), lc="Grassland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.fore, mean), lc="Woodland"))
p.vals.mean <- rbind(p.vals.mean, cbind(aggregate(.~variable, pvals.othe, mean), lc="Other"))
p.vals.mean

#write.csv(p.vals.mean, file="p_1000_trails_mean.csv", row.names = F)

## Extract maximum/upper quartile of p-values
p.vals.max <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.agri, 
                                                function(x) {quantile(x,0.75)}), lc="Cropland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.gras, 
                                                function(x) {quantile(x,0.75)}), lc="Grassland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.fore, 
                                                function(x) {quantile(x,0.75)}), lc="Woodland"))
p.vals.max <- rbind(p.vals.max, cbind(aggregate(.~variable, pvals.othe, 
                                                function(x) {quantile(x,0.75)}), lc="Other"))
p.vals.max

# minimum/ lower quartile values
p.vals.min <- data.frame(variable=NA, run=NA, value=NA, lc=NA)[0,] 
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.agri, 
                                                function(x) {quantile(x,0.25)}), lc="Cropland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.gras, 
                                                function(x) {quantile(x,0.25)}), lc="Grassland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.fore, 
                                                function(x) {quantile(x,0.25)}), lc="Woodland"))
p.vals.min <- rbind(p.vals.min, cbind(aggregate(.~variable, pvals.othe, 
                                                function(x) {quantile(x,0.25)}), lc="Other"))
p.vals.min

# Merge mean, max and min into one data frame
colnames(p.vals.mean)[3] <- "mean.P"
colnames(p.vals.max)[3] <- "quartile75"
colnames(p.vals.min)[3] <- "quartile25"

p.vals.summary <- merge(p.vals.mean[,-2], p.vals.max[,-2], 
                        by=c("variable", "lc"))
p.vals.summary <- merge(p.vals.summary, p.vals.min[,-2], by=c("variable", "lc"))
head(p.vals.summary)

#write.csv(p.vals.summary, file="p_1000_trails_summary.csv", row.names = F)

# merge with d values
statistics <- merge(p.vals.summary, d.table.summary, 
                    by.x=c("variable", "lc"), by.y=c("name", "lc"))
head(statistics)
#write.csv(statistics, file="1000_trails_statistics_p+d.csv", row.names = F)







