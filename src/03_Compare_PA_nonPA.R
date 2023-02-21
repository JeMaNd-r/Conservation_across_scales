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

# save sample function (if only 1 value in sample(), it considers it as vector)
resample <- function(x, ...) x[sample.int(length(x), ...)]

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

## Explore data
summary(as.factor(data_glob$PA))

# number of observations (raw)
nrow(data_glob); nrow(data_glob[data_glob$PA,])  #383 with 135 PAs

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## pair one PA with most similar nonPA sample point ####
# based on physical properties and distance

# define functions to be compared
fns = c("Soil_carbon_service", "OM_decomposition_service", "Water_regulation_service", 
        "Soil_stability_service", "Nutrient_service", "Pathogen_control", 
        "Richness_bacteria", "Fungi_18s_richness", "Invertebrate_18s_richness", 
        "Protist_18s_richness", "Ectomycorrhizal_18s_richness", 
        "Arbuscular_mycorrhizal_18s_richness", "Decomposers_18s_richness", 
        "Diss_Bacteria_std", "Diss_Fungi_std", "Diss_Protists_std", "Diss_invert_std")

## define variables to be compared between PA and nonPA, & threshold
mahal_vars <- c("Latitude_c", "Longitude_c", "Elevation", "AnnualPrec", "AnnualTemp", 
                #"MonthlyPrecipSum","MonthlyMeanTemp", 
                "Soil_pH", "Soil_salinity", "Clay_silt_c")

# keep complete cases only
data_glob <- data_glob[complete.cases(data_glob[,c(mahal_vars, fns)]),]
nrow(data_glob); nrow(data_glob[data_glob$PA==1,])  # n = 383 including 135 PAs

mahal_thres <- qchisq(.975, df=length(mahal_vars)) #21.92005

# scale variables for mahalanobis distance
for(i in mahal_vars){
  data_glob[,paste0(i,".z")] <- as.numeric(scale(data_glob[,i]))
}

mahal_vars <- paste0(mahal_vars, ".z")

# define each land cover type
lc.names <- c("Grassland", "Shrubland", "Woodland")

# Note: There were nonPA sites with mahalanobis distance below threshold
# for three protected sites with Order_ID. They have to been removed.
data_glob <- data_glob[!(data_glob$Order_ID %in% unpaired_pa$Order_ID),] 
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=371 with 123 PAs

# Remove sites that can only be paired once
data_glob <- data_glob[!(data_glob$Order_ID %in% count_nonPA[count_nonPA$n<3, "Order_ID"]),]
nrow(data_glob); nrow(data_glob[data_glob$PA,]) #nrow=364 with 116 PAs
data_glob %>% group_by(LC, PA) %>% count()

#data_glob <- data_glob[data_glob$Order_ID!=28302282,]
#data_glob <- data_glob[data_glob$Order_ID!=29542352,]

# split data
nonpa <- data_glob[data_glob$PA==0,c("Order_ID", "LC", "PA", mahal_vars, fns)]
pa <- data_glob[data_glob$PA==1,c("Order_ID", "LC", "PA", mahal_vars, fns)]

# randomization
p_list_total <- vector("list", length = 1000)
effect_size_d <- vector("list", length=1000)
times <- 0; times.with.error <- 0; set.seed(1) 
repeat {
 
}
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Compare difference between PA and nonPA sites ####
    
    # merge pairs of PA and nonPA in one table ####
    nonpa.subset <- nonpa[unique(nonpa$Order_ID) %in% pa$nonPA,]
    #data_glob.paired <- full_join(pa[order(pa$nonPA),], nonpa.subset[order(nonpa.subset$Order_ID),])
    #head(data_glob.paired)
    
    # Perform t tests
    p.list <- list()
    effect_size_d[[times]] <- list()
    for(l in lc.names){
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
      for(lc in 1:length(lc.names)){
        p_table[which(rownames(p_table)==lc.names[lc]),n] <- p.list[[lc.names[lc]]][[n]]["p.value"]
      }
    }
    
    # add to overall result list
    p_list_total[[times]] <- p_table  #p-values from Chi-squared test
  

save(effect_size_d,  file="d_1000_trails.RData")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save total list with p tables
#save(p.list.total, file="p_1000_trails.RData")
#write.csv(p.list.total, file="p_1000_trails.csv")


# - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Extract Cohen's d effect size for each run ####
d.table = as.data.frame(matrix(ncol=length(fns)+1,nrow=4000))
colnames(d.table) = c("lc", fns)
d.table$lc <- rep(lc.names, each=1000)
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

# ggplot(data=pvals.all2, aes(x=value, fill=variable))+
#   geom_histogram()+
#   geom_vline(xintercept=0.05, linetype="dashed")

## violin plot all together, p values
# ggplot(data=pvals.all2, aes(x=value, y=variable, fill=variable))+
#   geom_violin()+
#   geom_vline(xintercept=0.05, linetype="dashed")+
#   ggtitle("Cropland") +
#   theme(legend.position = "none", axis.title.y =element_blank())

## violin plot per land use type
# agri
aplot <- ggplot(data=pvals.agri, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))

# gras
gplot <- ggplot(data=pvals.gras, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Grassland") +
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

# Woodland
fplot<- ggplot(data=pvals.fore, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Woodland") +
  theme_bw() + # use a white background
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

# others
oplot<- ggplot(data=pvals.othe, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.05, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Others") +
  theme_bw() + # use a white background
  stat_summary(fun = "mean",geom = "point",color = "black")+
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.text.y = element_blank(),  axis.text.x = element_text(size=20))

library(ggpubr)
#setwd(figu.wd); pdf(file="Pvals_violin_1000trails_ttest.pdf", width=15)
ggarrange(aplot, gplot, fplot, oplot, ncol = 4, nrow = 1, align = "none",
          widths=c(2,1,1,1))
dev.off()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Effect size as violin plot ####
setwd(data.wd)
load(file="Effect_size_d_table.RData")
#load(file="Effect_size_dt_table.RData")    #Cohen's d from t value
#load(file="Effect_size_gs_table.RData")    #Hedge's g

dvals.all <- melt(d.table, id.vars=c("lc", "run"))
dvals.agri <- melt(d.table[d.table$lc=="Cropland",], id.vars = c("lc", "run"))
dvals.gras <- melt(d.table[d.table$lc=="Grassland",], id.vars = c("lc", "run"))
dvals.fore <- melt(d.table[d.table$lc=="Woodland",], id.vars = c("lc", "run"))
dvals.othe <- melt(d.table[d.table$lc=="Other",], id.vars = c("lc", "run"))

## violin plot per land use type
#setwd(figu.wd); pdf(file="Gvals_violin_1000trails_in1.pdf", height=15)
ggplot(data=dvals.all[dvals.all$lc!="Other",], aes(x=value, y=variable, fill=lc))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  geom_vline(xintercept=0, linetype="solid")+
  #facet_grid(lc~.) +
  theme(#legend.position = "none", 
    axis.title.y =element_blank(), 
    axis.text = element_text(size=10))
dev.off()

## or as separate plots merged at the end
# agri
aplot <- ggplot(data=dvals.agri, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# gras
gplot <- ggplot(data=dvals.gras, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Grassland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# Woodland
fplot<- ggplot(data=dvals.fore, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Woodland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

# others
oplot<- ggplot(data=dvals.othe, aes(x=value, y=variable))+
  geom_violin()+
  geom_vline(xintercept=0.2, linetype="dashed")+
  geom_vline(xintercept=-0.2, linetype="dashed")+
  stat_summary(fun = "mean",geom = "point",color = "black")+
  geom_vline(xintercept=0, linetype="solid")+
  ggtitle("Others") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20))+
  scale_x_continuous(limits = c(-1.5,0.8))

library(ggpubr)
#setwd(figu.wd); pdf(file="Dvals_violin_1000trails.pdf", height=15)
ggarrange(aplot, gplot, fplot, ncol = 1, nrow = 3, align = "none")
dev.off()

# ## with confidence intervals ####
setwd(data.wd)
d.table.summary <- read.csv("Effect_size_d_summary.csv")
head(d.table.summary)

aplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Cropland",], 
                aes(y=mean.D, x=name, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  geom_hline(yintercept=0, linetype="dashed")+
  ggtitle("Cropland") +
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_blank())+
  scale_y_continuous(limits = c(-1,1))

gplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Grassland",], 
                aes(y=mean.D, x=name, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  ggtitle("Grassland") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_blank())+
  scale_y_continuous(limits = c(-1,1))

fplot <- ggplot(data=d.table.summary[d.table.summary$lc=="Woodland",], 
                aes(y=mean.D, x=name, ymin=CI.lower,ymax=CI.upper))+
  geom_pointrange() +
  ggtitle("Woodland") +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw() + # use a white background
  theme(legend.position = "none", axis.title.y =element_blank(),
        axis.title.x =element_blank(),
        axis.text.y = element_text(size=20),  axis.text.x = element_text(size=20, angle=45, hjust=1))+
  scale_y_continuous(limits = c(-1,1))

library(ggpubr)
#setwd(figu.wd); pdf(file="Dvals_pointrange_1000trails.pdf", height=15)
ggarrange(aplot, gplot, fplot, ncol = 1, nrow = 3, align = "none", heights=c(1,1,1.7))
dev.off()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplots ####

# Standardize
std.data <- lucas[lucas$LUCAS_ID %in% unique(c(pa.pairs$LUCAS_ID, pa.pairs$nonPA)),]
std.data[,fns] <- scale(std.data[,fns])

# Labels for functions
labeling <- c("BAS", "Cmic", "Xylo", "Cellu", "bGluco", "NAG", 
              "Phosp", "MWD", "WSA")

# Split up by land cover type
std.Woodland = std.data[which(std.data$LC_5 == "Woodland"),]
std.grass = std.data[which(std.data$LC_5 == "Grassland"),]
std.agri = std.data[which(std.data$LC_5 == "Cropland"),]
std.other = std.data[which(std.data$LC_5 == "Other"),]

# melt data to get one colume with both PA and nonPA variables
amelt <- melt(std.agri, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)
gmelt <- melt(std.grass, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)
fmelt <- melt(std.Woodland, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns,na.rm=T)
omelt <- melt(std.other, id.vars=c("LUCAS_ID", "PA"), measure.vars = fns, na.rm=T)

## Plot
library(ggplot2)
library(ggpubr)

#Woodland
fplot <- ggplot(fmelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),text = element_text(size=20),
        legend.position = c(0.9,0.9), legend.text = element_text(size=20), axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Woodland") #+
#geom_text(x=2,y=12.5,label="*",size=9,color="red3")

#Cropland
aplot <- ggplot(amelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  scale_x_discrete(labels=labeling) +
  theme(text = element_text(size=20),legend.position = "none", axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Cropland") #+
#geom_text(x=1,y=2.25,label="*",size=9,color="red3")

#grassland
gplot <- ggplot(gmelt,aes(variable,value)) +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  geom_hline(yintercept = 0, lty="dashed") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(),text = element_text(size=20),
        legend.position = "none", axis.text.y = element_text(size=15)) +
  labs(x=NULL,y="Grassland")

#setwd(figu.wd); pdf(file="LUCAS_Boxplot_paNonpa.pdf", width=12)
ggarrange(fplot, gplot, aplot, ncol = 1, nrow = 3, align = "v")
dev.off()





