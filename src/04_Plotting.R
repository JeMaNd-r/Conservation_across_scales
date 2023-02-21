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


