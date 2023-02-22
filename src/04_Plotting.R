#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Plot data and results             #
#          author: Romy Zeiss               #
#            date: 2022-11-08               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(terra)

# load background map
world.inp <- map_data("world")

source("src/00_Parameters_functions.R")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v1.csv"))
data_glob

# load protected data
protect_df <- read_csv(paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

data_glob <- data_glob %>% 
  mutate("PA"=protect_df %>% dplyr::select(id.y, PA) %>% unique() %>% dplyr::select(PA) %>% unlist())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Sampling locations ####

pdf(paste0(here::here(), "/figures/Data_locations_global.pdf"))
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  xlim(-180, 180) +
  ylim(-180, 180) +
  
  geom_point(data=data_glob, aes(x=Longitude_c, y=Latitude_c, fill=as.factor(PA), shape=Eco_c), alpha=0.5)+
  scale_fill_manual(values = c("black", "olivedrab3"))+
  scale_shape_manual(values = c(21, 22, 23, 24))+ #label = c("Protected", "Unprotected")
  
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="right",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()



#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplot mahalanobis distance ####
pa.pairs <- merge(pa.pairs, lucas[,c("LUCAS_ID", "Country")])

# get sample sizes
no.sample <- pa.pairs %>% count(Country, LC)
#write.csv(no.sample, file="LUCAS_Boxplot_mahal.distance_samplesize.csv", row.names=F)
#no.sample <- as.character(no.sample$n)

#setwd(figu.wd); pdf(file="LUCAS_Boxplot_mahal.distance_wide.pdf", width=15) 
ggplot(pa.pairs[pa.pairs$LC!="Other",],aes(Country,mahal.min, fill=LC)) +
  geom_boxplot()+#,fatten = NULL) + # activate to remove the median lines
  theme_classic() +
  labs(x="Country",y="Mahalanobis distance") +
  theme(axis.text.x=element_text(size=15, angle=45, hjust=1),text = element_text(size=20),  
        legend.position = c(0.6,0.8), axis.text.y = element_text(size=15), legend.title = element_blank())+
  scale_fill_manual(values=c("gold3", "forestgreen", "limegreen"))
dev.off()


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

