#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Plot data and results             #
#          author: Romy Zeiss               #
#            date: 2023-03-23               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(terra)
library(ggdist) #to plot distributions (Bayesian)
library(gridExtra) #to add table in plot
library(ggrepel) #to add text not overlapping (geom_text_repel)
library(ggtext) #to add icons as axis labels
library(ggh4x) # for free axes in facet_grid
library(corrplot)
library(magick) # to create plot from multiple pngs

temp_scale <- "global"
#temp_scale <- "continental"
#temp_scale <- "regional"

# load background map
world.inp <- map_data("world")
if(temp_scale != "global"){
  world.inp <- subset(world.inp, region %in% c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan",
                                      "Belarus", "Belgium", "Bosnia and Herzegovina", "Bulgaria",
                                      "Croatia", "Cyprus", "Czechia","Denmark","Estonia","Finland", 
                                      "France","Georgia", "Germany", "Greece","Hungary","Iceland", 
                                      "Ireland", "Italy","Kazakhstan", "Kosovo", "Latvia","Liechtenstein", 
                                      "Lithuania", "Luxembourg","Malta","Moldova","Monaco","Montenegro",
                                      "Macedonia", "Netherlands","Norway","Poland","Portugal","Romania",
                                      "Russia","San Marino","Serbia","Slovakia","Slovenia","Spain",
                                      "Sweden","Switzerland","Turkey","Ukraine","UK","Vatican"))
}


source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

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
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
  min_size <- 7 # number of samples/ sites that should be paired per LC type
}
# define each land cover type
lc_names_all <- c( "Cropland", "Grassland", "Shrubland", "Woodland", "Other")

# define order of functions
labels_order <- c(
  "Decomposition (OM)", "Nutrient cycling", "Pathogen control", "Soil carbon", "Soil stability", "Water regulation",
  "Bacterial Richness", "Fungal Richness", "Invertebrate Richness", "Protist Richness",
  "AM fungi Richness", "EM fungi Richness",
  "Bacterial Shannon", "Fungal Shannon", "Invertebrate Shannon", "Protist Shannon",
  "Nematode Richness", "Decomposer Richness",
  "Bacterial Dissimilarity", "Fungal Dissimilarity", "Invertebrate Dissimilarity", "Protist Dissimilarity"
)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

# load pairs of PA and nonPA
pa_pairs <- read_csv(file=paste0(here::here(), "/intermediates/", temp_scale,"/Pairs_paNonpa_1000trails_", temp_date, ".csv"))
head(pa_pairs)

# load effect sizes
load(file=paste0(here::here(), "/results/d_1000_trails_", temp_scale, ".RData")) #d_list

# load Bayesian results from PA_type comparison
load(file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_", temp_scale, ".RData")) #pars_sample

# load Bayesian results from PA-nonPA comparison
load(file=paste0(here::here(), "/results/pars_PAdummy_Bayesian_df_", temp_scale, ".RData")) #pars_dummy_sample

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Sampling locations ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FIGURE 1 - Maps ####
# extract list of sampling locations actually used in comparison
data_locations <- data_clean %>% 
  filter(SampleID %in% unique(pa_pairs$ID) | 
         SampleID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Longitude,Latitude,SampleID, PA, LC)
data_locations #G: nrow=126, C: 316, R: 161
write_csv(data_locations, file = paste0(here::here(), "/results/Locations_", temp_scale, ".csv"))
nrow(data_locations %>% filter(PA==1)) #G: 28 PAs, C: 48, R: 36
nrow(data_locations %>% filter(PA==0)) #G: 93 PAs, C: 268, R: 125

# set limits for point maps
if(temp_scale == "global") temp_limits <- c(-180, 180, -180, 180)
if(temp_scale == "continental") temp_limits <- c(-10, 35, 35, 70)
if(temp_scale == "regional") temp_limits <- c(-9, -6, 40.5, 42.5)

# plot - FIGURE 1
ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey90", linewidth = 0.15) + #fill = "grey80", color="grey75"
  xlim(temp_limits[1], temp_limits[2])+
  ylim(temp_limits[3], temp_limits[4])+
  
  geom_point(data=data_locations, aes(x=Longitude, y=Latitude, 
                                      shape = as.character(PA), color=LC, 
                                      size = as.character(PA)),
             stroke = 2)+ #increase circle line width; G: 2, C+R:3
  scale_shape_manual(values = c("0" = 19, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" =1.5, "1" = 4.5))+ #G: 1.4,4.5, C+R:3,8
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000"))+ 
  coord_map()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="right",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill= "grey80"))
ggsave(filename=paste0(here::here(), "/figures/Data_locations_", temp_scale,".png"),
       plot = last_plot())


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.1 - Maps for random-slope model ####
# extract list of sampling locations actually used in comparison
data_locations <- data_clean %>%
  filter(LC %in% lc_names)
data_locations #G: nrow=248, C: 745, R: 270
write_csv(data_locations, file = paste0(here::here(), "/results/Locations_PAranks_", temp_scale, ".csv"))
nrow(data_locations %>% filter(PA==1)) #G: 42 PAs, C: 61, R: 40
nrow(data_locations %>% filter(PA==0)) #G: 206 PAs, C: 684, R: 230

# set limits for point maps
if(temp_scale == "global") temp_limits <- c(-180, 180, -180, 180)
if(temp_scale == "continental") temp_limits <- c(-10, 35, 35, 70)
if(temp_scale == "regional") temp_limits <- c(-9, -6, 40.5, 42.5)

# plot - APPENDIX FIGURE 3.1
ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey90", linewidth = 0.15) + #fill = "grey80", color="grey75"
  xlim(temp_limits[1], temp_limits[2])+
  ylim(temp_limits[3], temp_limits[4])+
  
  geom_point(data=data_locations, aes(x=Longitude, y=Latitude, 
                                      shape = as.character(PA), color=LC, 
                                      size = as.character(PA)),
             stroke = 2)+ #increase circle line width; G+C: 2; R:3
  scale_shape_manual(values = c("0" = 19, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" =1.5, "1" = 4.5))+ #G:+C: 2,4; ,R:3,8 
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000"))+ 
  coord_map()+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="right",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill= "grey80"))
ggsave(filename=paste0(here::here(), "/figures/Data_locations_PAranks_", temp_scale,".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pairing ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Save some numbers
sink(paste0(here::here(), "/results/Numbers_d-value_", temp_scale, ".txt"))

print(temp_scale)
# look what non-protected sites have (not) been paired to any PA
cat("Number of non-protected sites have never been paired to any PA:")
hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(data_clean[data_clean$PA==0,]$SampleID, pa_pairs$nonPA))  
# nonPA sites never used: G ..., C: 475

cat(paste0("On average, each protected and unprotected site was used in ",
           round(mean(table(pa_pairs$ID)),0), 
           " (SD=", round(sd(table(pa_pairs$ID)),0),
           ", min=", round(min(table(pa_pairs$ID)),0),
           ", max=", round(max(table(pa_pairs$ID)),0), ") and ",
           round(mean(table(pa_pairs$nonPA)),0),
           " (SD=", round(sd(table(pa_pairs$nonPA)),0),
           ", min=", round(min(table(pa_pairs$nonPA)),0),
           ", max=", round(max(table(pa_pairs$nonPA)),0), 
           ") of the 1000 randomizations, respectively."))

table(pa_pairs$ID)

cat("Number of sites (PA/ nonPA) per LC.")
pa_pairs %>% group_by(LC) %>% dplyr::select(ID) %>% unique() %>% count()
pa_pairs %>% group_by(LC) %>% dplyr::select(nonPA) %>% unique() %>% count()

sink()

## Plot pairs on map
data_pairs <- pa_pairs %>%
  full_join(data_clean %>%
    filter(SampleID %in% unique(pa_pairs$ID) | 
             SampleID %in% unique(pa_pairs$nonPA)) %>%
    dplyr::select(Longitude,Latitude,SampleID, PA, LC),
    by = c("ID" = "SampleID", "LC")) 

# Merge data to get coordinates for each ID
data_merged <- pa_pairs %>%
  inner_join(data_clean, by = c("ID" = "SampleID")) %>%
  inner_join(data_clean, by = c("nonPA" = "SampleID"), suffix = c("_1", "_2"))

# Calculate frequency of appearance of pairs
pair_freq <- data_merged %>%
  group_by(ID, nonPA) %>%
  summarise(freq = n())

# Create plot with ggplot2
ggplot() +
  geom_segment(data = data_merged %>% unique(), aes(x = Longitude_1, y = Latitude_1, xend = Longitude_2, yend = Latitude_2), alpha = 0.5, color = "blue") +
  geom_point(data = data_clean, aes(x = Longitude, y = Latitude)) +
  #scale_size_continuous(name = "Frequency", guide = "legend") +
  labs(x = "Longitude", y = "Latitude", title = "Paths Between Points with Width Based on Frequency") +
  theme_minimal()
ggsave(filename=paste0(here::here(), "/figures/Locations_connection_", temp_scale, ".pdf"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Mahalanobis distance ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boxplot Mahalanobis distance ####

for(temp_scale in c("global", "continental", "regional")){
  # set date of latest analysis
  if(temp_scale == "global") temp_date <- "2024-09-12"
  if(temp_scale == "continental") temp_date <- "2024-08-01"
  if(temp_scale == "regional") temp_date <- "2024-08-01"
  
  # load pairs of PA and nonPA
  pa_pairs <- read_csv(file=paste0(here::here(), "/intermediates/", temp_scale,"/Pairs_paNonpa_1000trails_", temp_date, ".csv"))
  
  # plot
  ggplot(pa_pairs, aes(x = LC, y = mahal.min, fill = LC))+
    geom_boxplot()+
    geom_violin(alpha = 0.3, adjust = 0.3)+
    theme_bw() +
    ggtitle(paste0("Paired samples - ", temp_scale))+
    labs(x="Habitat type",y="Mahalanobis distance") +
    theme(axis.text.x=element_text(size=15),text = element_text(size=20),  
          legend.position = "none", axis.text.y = element_text(size=15), legend.title = element_blank())+
    scale_fill_manual(values=c("Cropland" = "#4A2040",
                               "Grassland" = "#E69F00",
                               "Shrubland" = "#0072B2", 
                               "Woodland" = "#009E73", 
                               "Other" = "#000000")) #"gold3", "limegreen", "forestgreen"
  ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_", temp_scale, ".png"),
         plot = last_plot())
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### APPENDIX FIGURE 2.1 - Boxplot Mahalanobis distance - all vs. pairs ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# calculate & plot Mahalanobis distance for all
for(temp_scale in c("global", "continental", "regional")){
  data <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  
  if(temp_scale == "global"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland"]
  }
  
  # scale variables
  data <- f_scale_vars(data = data, vars = mahal_vars)
  
  nonpa <- data[data$PA==0,c("SampleID", "LC", "PA", mahal_vars_z)]
  pa <- data[data$PA==1,c("SampleID", "LC", "PA", mahal_vars_z)]
  all_nonPA <- data.frame("SampleID" = "a", "PA_ID" = 1, "LC" = "lc", "mahal_dist" = 1)[0,]
  
  # calculate mahalanobis distance per lc type
  for(i in 1:length(lc_names)){try({
    temp_nonPA  <- nonpa[nonpa[,"LC"]==lc_names[i],]
    temp_PA <- pa[pa[,"LC"]==lc_names[i],]
    
    sigma <- cov(temp_nonPA[,mahal_vars_z]) 
    for(j in 1:nrow(temp_PA)){
      mu <- as.numeric(temp_PA[j,mahal_vars_z])
      temp_nonPA[,as.character(temp_PA[j,"SampleID"])] <- 
        mahalanobis(temp_nonPA[,mahal_vars_z], mu, sigma, tol=1e-30)
      #print(j)
    }
    
    temp_nonPA <- temp_nonPA %>% pivot_longer(cols = as.character(temp_PA %>% pull("SampleID")),
                                names_to = "PA_ID", values_to = "mahal_dist") %>%
      dplyr::select(SampleID, PA_ID, LC, mahal_dist) %>% unique()
    
    all_nonPA <- rbind(all_nonPA, temp_nonPA)
  })}
  all_nonPA
  
  write_csv(all_nonPA, paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_", temp_scale, ".csv"))
  
  ggplot(all_nonPA, aes(x = LC, y = mahal_dist, fill = LC))+
    geom_boxplot()+
    geom_violin(alpha = 0.3, adjust = 0.3)+
    theme_bw() +
    ggtitle(paste0("All samples - ", temp_scale))+
    labs(x="Habitat type",y="Mahalanobis distance") +
    theme(axis.text.x=element_text(size=15),text = element_text(size=20),  
          legend.position = "none", axis.text.y = element_text(size=15), legend.title = element_blank())+
    scale_fill_manual(values=c("Cropland" = "#4A2040",
                               "Grassland" = "#E69F00",
                               "Shrubland" = "#0072B2", 
                               "Woodland" = "#009E73", 
                               "Other" = "#000000")) #"gold3", "limegreen", "forestgreen"
  ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_", temp_scale, ".png"),
         plot = last_plot())
}

# put in 1 plot - APPENDIX FIGURE 2.1
pairs_glob <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_global.png"))
pairs_cont <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_continental.png"))
pairs_regi <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_regional.png"))
all_glob <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_global.png"))
all_cont <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_continental.png"))
all_regi <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_regional.png"))

magick::image_write(magick::image_append(c(magick::image_append(c(pairs_glob, pairs_cont, pairs_regi)),
                                           magick::image_append(c(all_glob, all_cont, all_regi))),
                     stack = TRUE),
                    path=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_allScales.png"))
rm(pairs_glob, pairs_cont, pairs_regi, all_glob, all_cont, all_regi)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX TABLE 2.1 - Effect size per LC type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
d_df <- do.call(rbind, d_list)
str(d_df)

d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

# save data for plot
write.csv(d_df, file=paste0(here::here(), "/figures/Data_d-value_", temp_scale, ".csv"), row.names = FALSE)

d_df <- read_csv(file=paste0(here::here(), "/figures/Data_d-value_", temp_scale, ".csv"))

d_summary <- d_df %>% 
  dplyr::select(-run) %>%
  #pivot_longer(cols = c(p_value, ci_lower, ci_upper, t_stats),
  #             names_to = "metric") %>%
  group_by(lc, fns, Label, Group_function, Organism) %>%
  summarize(across(effect, .fns = list("mean"=mean, "SD"=sd, "median"=median,
                                       "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
                                       "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
                                       "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
                                       "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE))))
d_summary
write.csv(d_summary, file=paste0(here::here(), "/figures/Results_d-value_summary_", temp_scale, ".csv"))

# mean per lc type
d_summary <- read.csv(file=paste0(here::here(), "/figures/Results_d-value_summary_", temp_scale, ".csv"))

d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median, effect_ci_2.5:effect_ci_97.5), 
                   function(x) mean(x)))
d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# mean per Group_function
d_summary %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# check for significant mean p-values
#d_summary %>% arrange(p_value_mean)

sink(file=paste0(here::here(), "/results/D-values_", temp_scale, ".txt"))
cat("#################################################", sep="\n")
cat("##  Significant d-values from global analysis  ##", sep="\n")
cat(paste0("##  Sys.Date() ", Sys.Date(), "  ##"), sep="\n")
cat("#################################################", sep="\n")
print(d_summary %>% ungroup() %>%
  filter(abs(effect_mean) >= 0.2) %>%
  dplyr::select(lc, fns, starts_with("effect")) %>%
    arrange(desc(abs(effect_median))),
)
sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 2.2 - Heatmap correlation Mahalanobis & difference ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

all_corr <- data.frame("LC" = "lc", "fns" = "fns", "correlation" = 1, "scale" = "scale", "p_value" = 1)[0,]

for(temp_scale in c("global", "continental", "regional")){
  # set date of latest analysis
  if(temp_scale == "global") temp_date <- "2024-09-12"
  if(temp_scale == "continental") temp_date <- "2024-08-01"
  if(temp_scale == "regional") temp_date <- "2024-08-01"
  
  if(temp_scale == "global"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland"]
  }
  
  data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  
  # load pairs of PA and nonPA
  pa_pairs <- read_csv(file=paste0(here::here(), "/intermediates/", temp_scale,"/Pairs_paNonpa_1000trails_", temp_date, ".csv"))
  
  # add data
  pa_env <- pa_pairs %>% full_join(data_clean, by = c("nonPA" = "SampleID", "LC")) %>% mutate(data = "nonPA") %>%
    rbind(pa_pairs %>% full_join(data_clean, by = c("ID" = "SampleID", "LC")) %>% mutate(data = "PA")) %>%
    filter(!is.na(nonPA)) %>%
    arrange(ID, nonPA, times) 
  
  # calculate difference for fns columns and mean for mahal.min column (= mahal.min because same for both) within each group
  pa_env <- pa_env %>% group_by(ID, nonPA, LC, times, n) %>% 
    summarize(across(all_of(fns), diff),
              across(mahal.min, mean)) %>% 
    dplyr::select(ID, nonPA, mahal.min, LC, all_of(fns))
  
  # Calculate correlations between env and fns for each site
  correlation_matrix <- sapply(fns, function(fn) {
    sapply(lc_names, function(lc) {
      pa_temp <- pa_env %>% filter(LC == lc)
      cor.test(pa_temp[["mahal.min"]], pa_temp[[fn]], 
               use = "na.or.complete", method = "spearman")[["estimate"]][["rho"]] #cor for Pearson, rho for Spearman, tau for Kendall
    })
  })
  
  correlation_p_matrix <- sapply(fns, function(fn) {
    sapply(lc_names, function(lc) {
      pa_temp <- pa_env %>% filter(LC == lc)
      round(cor.test(pa_temp[["mahal.min"]], pa_temp[[fn]], 
                     use = "na.or.complete", method = "spearman")[["p.value"]],3)
    })
  })
  
  # Convert the correlation matrix to a data frame for easier interpretation
  correlation_df <- as.data.frame(correlation_matrix)
  rownames(correlation_df) <- lc_names  
  colnames(correlation_df) <- fns 
  
  correlation_p_df <- as.data.frame(correlation_p_matrix)
  rownames(correlation_p_df) <- lc_names  
  colnames(correlation_p_df) <- fns  
  
  correlation_df <- correlation_df %>% 
    rownames_to_column(var = "LC") %>%  # Convert row names to a column
    pivot_longer(cols = -LC,  # Pivot to long format, excluding the env_variable column
                 names_to = "fns",
                 values_to = "correlation") %>%
    mutate(fns = factor(fns, levels = unique(fns[order(correlation)])),
           LC = factor(LC, levels = lc_names),
           scale = temp_scale)
  
  correlation_p_df <- correlation_p_df %>% 
    rownames_to_column(var = "LC") %>%  # Convert row names to a column
    pivot_longer(cols = -LC,  # Pivot to long format, excluding the env_variable column
                 names_to = "fns",
                 values_to = "p_value") %>%
    mutate(fns = factor(fns, levels = unique(fns[order(p_value)])),
           LC = factor(LC, levels = lc_names),
           scale = temp_scale)
  
  correlation_df <- correlation_df %>%
    full_join(correlation_p_df)
  
  all_corr <- rbind(all_corr, correlation_df)
}
all_corr

# long format
all_corr_plot <- all_corr %>%
  mutate(LC = factor(LC, levels = lc_names_all[lc_names_all != "Other"]),
         scale = factor(scale, levels = c("global", "continental", "regional")),
         scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
  full_join(fns_labels %>% 
              filter(Function %in% fns),
            dplyr::select(Label_short, Function), 
            by = c("fns" = "Function")) %>%
  mutate(Label_short = factor(Label_short, levels = rev(fns_labels$Label_short)),
         scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                    "<img src='figures/icon_location-black.png' width='70'>",
                                                    "<img src='figures/icon_flag-Portugal.png' width='70'>" )))

# add d-values (effect size results)
all_corr_plot <- all_corr_plot %>% 
  full_join(d_plot_all %>%
              filter(!is.na(effect_significance)) %>%
              rename(LC = lc) %>%
              dplyr::select(Label_short, LC, fns, Group_function, effect_significance, scale))

# save
write_csv(all_corr_plot %>%
            mutate("Correlation [p value]" = paste0(round(correlation, 3), 
                                                    " [", ifelse(p_value<0.002, "<0.002", round(p_value, 3)), "] ", 
                                                    ifelse(p_value < 0.05, "*", "ns"))) %>%
            dplyr::select(Group_function, Label_short, scale, LC, `Correlation [p value]`) %>%
            rename(Group = Group_function,
                   Variable = Label_short, 
                   Scale = scale) %>%
            pivot_wider(id_cols = c(Group, Variable, Scale), names_from = LC, values_from = `Correlation [p value]`) %>%
            dplyr::select(Group, Variable, Scale, Cropland, Grassland, Shrubland, Woodland) %>%
            arrange(Group, Variable, Scale),
          paste0(here::here(), "/figures/Correlation_diff_mahal_allScales.csv"))

# plotting
ggplot(data = all_corr_plot)+
  geom_tile(aes(x = LC, y = Label_short, fill = correlation))+
  
  geom_point(data = all_corr_plot %>% filter(p_value >= 0.05),
             aes(x = LC, y = Label_short), shape = 4)+
  
  geom_point(data = all_corr_plot %>% filter(effect_significance != "ns"),
             aes(x = LC, y = Label_short), shape = 0, size = 8)+
  
  scale_shape_manual(values = c(NA, 4))+
  scale_fill_distiller(type = "div", na.value = "white", limits = c(-0.7, 0.7),
                       name = "Correlation")+
  scale_x_discrete(labels = c(
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
  ))+
  facet_grid(cols = vars(scale_icon))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 15),
        strip.text.x = ggtext::element_markdown(vjust = 0.5),
        strip.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.x = ggtext::element_markdown(vjust = 0))
ggsave(paste0(here::here(), "/figures/Correlation_diff_mahal_allScales.png"),
       last_plot(),
       height = 10, width = 8)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Effect sizes (Cohen's D) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Pointrange effect size per LC type ####

# with confidence intervals 
ggplot(data = d_df %>% 
         filter(lc %in% lc_names) %>%
         mutate(Label=factor(Label, levels = rev(fns_labels$Label))), 
       aes(y = effect, x = Label))+
  
  geom_hline(aes(yintercept=0.8, linetype = "-0.8 / 0.8"), color="grey60")+
  geom_hline(aes(yintercept=-0.8, linetype = "-0.8 / 0.8"), color="grey60")+ #large effect
  geom_hline(aes(yintercept=0.5, linetype = "-0.5 / 0.5"), color="grey60")+
  geom_hline(aes(yintercept=-0.5, linetype = "-0.5 / 0.5"), color="grey60")+ #medium effect
  geom_hline(aes(yintercept=0.2, linetype = "-0.2 / 0.2"), color="grey90")+ #small effect
  geom_hline(aes(yintercept=-0.2, linetype = "-0.2 / 0.2"), color="grey90")+
  #geom_hline(aes(yintercept=0, linetype = "0"))+
  scale_linetype_manual(values = c("-0.8 / 0.8" = "dotted",
                                   "-0.5 / 0.5" = "dashed", 
                                   "-0.2 / 0.2" = "solid"), name="Strength of effect")+
  annotate("rect", ymin = -0.2, ymax = 0.2, xmin=-Inf, xmax=Inf, fill = "grey95")+

  # fill background of Group_function (code order matters)
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=Inf, fill = "grey", alpha=0.1)+ #chocolate4
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=15+0.5, fill = "grey", alpha=0.1)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=8+0.5, fill = "grey", alpha=0.1)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=4+0.5, fill = "grey", alpha=0.1)+
  
  ggdist::stat_pointinterval(fatten_point=1.2, shape=21) +
  coord_flip()+
  ylab("Effect size")+
  facet_wrap(vars(lc))+
  theme_bw() + # use a white background
  theme(legend.position = "bottom", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianCI_", temp_scale, ".png"),
       plot = last_plot(),
       width=5, height=4)

## one-sided
ggplot(data = d_summary %>%
         filter(lc %in% lc_names) %>%
         mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
                                         abs(effect_ci_17), 
                                         abs(effect_ci_83))),
       aes(y = abs(effect_ci_min66), x = abs(effect_median),
           color=as.factor(sign(effect_median))
           ))+

  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.2, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.5, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.8, fill = "grey", alpha=0.3)+

  geom_point(aes(shape = Group_function), size = 5)+
  ggrepel::geom_text_repel(aes(label = Label, shape = Group_function), size = 3)+
  
  xlab("Median effect")+ ylab("Lower CI (17%)")+
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8))+
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))+
  
  facet_wrap(vars(lc), scales = "free")+
  theme_bw() + # use a white background
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianSD_", temp_scale, ".png"),
       plot = last_plot())

## heatmap
ggplot(data = d_summary %>%
         filter(lc %in% lc_names) %>%
         mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
                                         abs(effect_ci_17), 
                                         abs(effect_ci_83))) %>%
         mutate(effect_ci_min66f = cut(effect_ci_min66,
                                       breaks=c(0, 0.2, 0.5, 0.8, Inf),
                                       labels=c("ns", "small", "medium", "large"))) %>%
         filter(!is.na(effect_ci_min66)),
       aes(x = lc, y = Label, alpha=effect_ci_min66f, 
           fill=as.factor(sign(effect_median))))+

  geom_tile()+
  scale_fill_manual(values = c("-1" = "#fc8d59", "0" = "#ffffbf", "1" = "#91bfdb"),
                    name = "Direction of effect")+
  scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
                     name = "Minimum effect size (66% CI)")+
  theme_bw() + # use a white background
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        #axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianSD_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### FIGURE 2 & APPENDIX TABLE 2.2 - Heatmap of effect sizes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mean per lc type and all 3 scales
d_sum_glob <- read.csv(file=paste0(here::here(), "/figures/Results_d-value_summary_global.csv"))
d_sum_cont <- read.csv(file=paste0(here::here(), "/figures/Results_d-value_summary_continental.csv"))
d_sum_regi <- read.csv(file=paste0(here::here(), "/figures/Results_d-value_summary_regional.csv"))

d_sum_all <- rbind(d_sum_glob %>% mutate("scale" = "global"), 
                   d_sum_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_sum_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_sum_glob, d_sum_cont, d_sum_regi)

d_sum_all %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median, effect_ci_2.5:effect_ci_97.5), 
                   function(x) mean(x)))
d_sum_all %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# mean per Group_function
d_sum_all %>% ungroup() %>% 
  mutate("Supergroup_function" = ifelse(Group_function == "Function", "Functioning", "Diversity")) %>%
  group_by(Supergroup_function) %>%
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x), na.rm=TRUE),
                        sd=function(x) sd(abs(x), na.rm=TRUE))))

## switch lc and scales
d_plot_all <- d_sum_all %>%
  filter(lc %in% lc_names_all) %>%
  dplyr::select(-Label) %>%
  right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
             by = c("fns" = "Function")) %>%
  mutate(effect_mean_f = cut(abs(effect_mean),
                                breaks=c(0, 0.2, 0.5, 0.8, Inf),
                                labels=c("marginal", "small", "medium", "large"))) %>%
  filter(!is.na(effect_mean)) %>% #!is.na(effect_ci_min66) & 
  full_join(expand.grid(scale = c("global", "continental", "regional"), 
                        lc = lc_names_all, 
                        Label = fns_labels$Label)) %>%
  mutate(scale = factor(scale, levels = rev(c("global", "continental", "regional"))),
         Label = factor(Label, levels = rev(fns_labels %>% arrange(Group_function, Label) %>% pull(Label))),
         lc = factor(lc, levels = c("Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
  filter(lc != "Other" & !is.na(Label)) %>% 
  mutate(effect_direction = as.factor(sign(effect_mean)))%>%
  mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                     ifelse(effect_direction=="1", "positive", "0"))) %>%
  mutate(Label = factor(Label, levels = labels_order)) %>%
  mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c),
         effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
  mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "ns")))

# save table - APPENDIX TABLE 2.2
write_csv(d_sum_all %>% 
            dplyr::select(-X) %>%
            dplyr::select(Group_function, Label, scale, lc, effect_mean:effect_ci_97.5) %>%
            arrange(Group_function, Label, scale, lc) %>%
            mutate("Habitat" = lc,
                   "Group" = Group_function,
                   "Slope [HPD]" = paste0(round(effect_mean, 3), 
                                          " [", round(effect_ci_2.5, 3), "; ", round(effect_ci_97.5, 3), "] "),
                   "Variable" = Label,
                   "Scale_long" = factor(scale, levels = c("global", "continental", "regional")))  %>%
            mutate("Scale" = factor(substr(scale, 1, 4), levels = c("glob", "cont", "regi"))) %>%
            dplyr::select(Group, Variable, Scale, Habitat, "Slope [HPD]") %>%
            pivot_wider(names_from = "Habitat", values_from = "Slope [HPD]") %>%
            arrange(Group, Variable, Scale),
          paste0(here::here(), "/figures/Results_d-value_meanCI_allScales_fns.csv"))

# plot - FIGURE 2
ggplot(data = d_plot_all,
       aes(x = lc, y = scale))+
  
  geom_point(aes(size = effect_mean_f,
                 color = effect_direction_c, 
                 fill= effect_significance,
                 shape = effect_na))+
  facet_wrap(vars(Label), ncol=6, drop=FALSE)+
  scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "ns" = "white"),
                    name = "Direction of effect",
                    na.value = "black", drop = FALSE)+
  scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
                    name = "Direction of effect",
                    na.value = "black")+
  scale_size_manual(values = c("marginal" = 2, "ns" = 5, "small" = 5, "medium" = 10, "large" = 15),
                     name = "Effect size",
                    na.value = 5)+
  scale_shape_manual(values = c("not available" = 4),
                     name = "Missing data",
                     na.value = 21)+
  scale_x_discrete(labels = c(
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
  ))+
  
  scale_y_discrete(labels = c(
    "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
    "continental" = "<img src='figures/icon_location-black.png' width='30'>",
    "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
  ))+
  xlab("")+ylab("")+
  theme_bw() + # use a white background
  
  guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black"), 
                                                 shape = 21, size = 5)), #tell legend to use different point shape
    color = "none", #don't show legend
    shape = guide_legend(override.aes = list(size = 5)))+
  theme(legend.position = "bottom", #c(0.96, -0.05),
        legend.justification = c(1, 0),
        legend.box = "horizontal",
        legend.direction = "vertical",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text.y = ggtext::element_markdown(hjust = 0),
        axis.ticks = element_blank(),
        axis.text.x = ggtext::element_markdown(vjust = 0),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        strip.text = element_text(color="black", size = 15, hjust = 0)) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_meanCI_allScales_fns.png"),
       plot = last_plot(), 
       width = 4400, height = 3800, units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### Summarizing stats ####
# number significant effects
d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% nrow() #33
# number ns
d_plot_all %>% filter(effect_significance=="ns" & !is.na(effect_significance)) %>% nrow() #156

# number significant per lc
table(d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% dplyr::select(scale, lc))
table(d_plot_all %>% filter(!is.na(effect_significance)) %>% dplyr::select(scale, lc))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 2.3 & TABLE 2.3 - Pointrange plot grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

d_df_glob <- read_csv(file=paste0(here::here(), "/figures/Data_d-value_global.csv"))
d_df_cont <- read_csv(file=paste0(here::here(), "/figures/Data_d-value_continental.csv"))
d_df_regi <- read_csv(file=paste0(here::here(), "/figures/Data_d-value_regional.csv"))

d_df_all <- rbind(d_df_glob %>% mutate("scale" = "global"), 
                   d_df_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_df_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_df_glob, d_df_cont, d_df_regi)

d_df_grouped <- d_df_all %>% filter(!is.na(lc)) %>%
  mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='30'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='30'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>", NA)))) %>%
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>", 
                                          "<img src='figures/icon_location-black.png' width='30'>", 
                                          "<img src='figures/icon_flag-Portugal.png' width='30'>"))) %>%
  mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                 ifelse(Group_function=="Diversity", "Shannon", Group_function))) %>%
  mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity")))

write_csv(d_df_grouped %>%
            group_by(Group_function, scale, lc) %>%
            summarize(across(effect, list(median = median, 
                                          ci2.5 = function(x) quantile(x, 0.025, na.rm = TRUE), 
                                          ci92.5 = function(x) quantile(x, 0.925, na.rm = TRUE)))) %>%
            mutate(across(c(effect_median, effect_ci2.5, effect_ci92.5), function(x) round(x, 3))) %>%
            mutate(ci_95 = paste0("[", effect_ci2.5, "; ", effect_ci92.5, "]")) %>%
            dplyr::select(-effect_ci2.5, -effect_ci92.5), 
          paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_grouped.csv"))

# plot
ggplot(data = d_df_grouped,
       aes(fill = lc, color = lc, 
           y = scale_icon, x = effect))+
  
  geom_vline(aes(xintercept=0), color="black")+
  ggdist::stat_pointinterval(fatten_point=1.2, shape=21,
                             position=position_dodgejust(width=0.5)) +
  coord_flip()+
  facet_wrap(vars(Group_function), drop=FALSE)+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000"))+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000"))+
  scale_x_continuous(breaks = c(-2, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 2))+
  theme_void()+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal",
        axis.title.y =element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_blank(),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.text.y = element_text(size = 13, hjust = 1),  
        axis.text.x = ggtext::element_markdown(vjust = 1, margin = unit(c(1, 1, 20, 1), "pt")),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        strip.text = element_text(size = 20, hjust = 0),
        plot.background = element_rect(fill = "white", color = "white"))

#ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_groupedOrganisms.png"), #switch facet_wrap to Organism
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_grouped.png"),
       plot = last_plot(), 
       width = 2700, height = 2200,
       units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Random-intercept models (PA types) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.3 & 3.5 & 3.7 - Bayesian PA types ####

# transform parameters to long format and assign labels
pars_long <- pars_sample %>% 
  dplyr::select(-sigma, -sigma_a, -mu_a) %>% 
  pivot_longer(cols = colnames(.)[!(colnames(.) %in% c("lc", "fns"))])%>%
  full_join(protect_type %>% 
              mutate(PA_rank = as.character(PA_rank)), 
            by=c("name" = "PA_rank"))  %>%
  mutate("Label_pa" = ifelse(name=="11", "Unprotected", PA_type)) %>%
  mutate("Label_pa" = factor(Label_pa, levels = rev(c(protect_type$PA_type, "Unprotected")))) %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label_fns" = factor(Label, fns_labels$Label))
head(pars_long)

# save summary (i.e., data from figure)
pars_summary <- pars_long %>% group_by(lc, fns) %>% 
  summarize(across(value, list("mean"  = function(x) mean(x, na.rm=TRUE), 
                               "median" = function(x) median(x, na.rm=TRUE), 
                               "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
                               "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
                               "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
                               "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE)))) %>%
  arrange(lc, fns)
pars_summary
write_csv(pars_summary, file=paste0(here::here(), "/figures/Results_intercept_parsBayesian_summary_", temp_scale, ".csv"))

# extract sample size
n_table <- data_clean %>% filter(LC!="Other") %>%
  group_by(LC, PA, PA_type, PA_rank) %>% count() %>%
  pivot_wider(id_cols = c("PA", "PA_type", "PA_rank"), 
              names_from = LC, values_from = n) %>% 
  arrange(PA_rank) %>% ungroup() %>%
  dplyr::select(-PA_rank)
n_table
write_csv(n_table, file=paste0(here::here(), "/figures/Results_intercept_parsBayesian_nTable_", temp_scale, ".csv"))

ggplot(pars_long %>% filter(!is.na(Label)) %>% #filter(!is.na(PA_type)) %>%
                     # add number of sizes to plot
                     rbind(data_clean %>% filter(LC!="Other") %>%
                             group_by(LC, PA_type) %>% count() %>%
                             rename(lc=LC, value=n) %>%
                             mutate(fns="Number of sites",
                                    name=NA,
                                    PA_protected=NA,
                                    Label_pa=ifelse(is.na(PA_type), "Unprotected", PA_type),
                                    Label="Number of sites",
                                    Label_fns = Label,
                                    Label_short = Label,
                                    Group_function=NA,
                                    Organism=NA) %>%
                             dplyr::select(colnames(pars_long))) %>%
                     filter(!is.na(Label_pa)) %>%
                     mutate(Label_fns = factor(Label_fns, levels = c(labels_order, "Number of sites"))),
       aes(y=Label_pa, x=value, color=lc))+
  
  ## adapt for scale
  annotate("rect", ymin = -Inf, ymax = 4+0.5, xmin=-Inf, xmax=Inf, fill = "grey90", alpha=0.5)+ #global: 4, C: 3, R: 2
  annotate("rect", ymin = -Inf, ymax = 1+0.5, xmin=-Inf, xmax=Inf, fill = "grey85", alpha=0.5)+
  
  ggdist::stat_pointinterval(fatten_point=1, shape=3, 
                             position=position_dodgejust(width=0.5))+ 
  facet_wrap(vars(Label_fns), scales = "free_x", ncol=6, drop=FALSE)+
  
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000"), name="Habitat type")+
  ylab("")+ xlab("")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 13))
ggsave(filename=paste0(here::here(), "/figures/Results_intercept_parsBayesian_", temp_scale, ".png"),
       plot = last_plot(),
       width=12, height=10)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Random-slope model (PA ranks/ levels) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.3 & 3.6 & 3.8 (PA ranks/ levels) ####

pred_list <- rbind(get(load(paste0(here::here(), "/results/PAranks_Bayesian_global_sample10k.RData"))), 
                   get(load(paste0(here::here(), "/results/PAranks_Bayesian_continental_sample10k.RData")))) %>% 
  rbind(get(load(paste0(here::here(), "/results/PAranks_Bayesian_regional_sample10k.RData"))))

for(temp_scale in c("global", "continental", "regional")){
  ggplot(data = pred_list %>% filter(scale == temp_scale) %>%
           right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
                      by = c("fns" = "Function")) %>%
           mutate(Label = factor(Label, levels = labels_order)), 
         
         aes(x = PA_rank_rev, y = .epred, color = ordered(LC))) +
    stat_lineribbon() +
    facet_wrap(vars(Label), scales = "free_y", ncol=6)+
    scale_fill_brewer(palette = "Greys") +
    scale_color_manual(values=c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000"), name="Habitat type")+
    scale_x_continuous(limits = c(1, 10), breaks = c(2, 10), minor_breaks = c(2,4,6,8, 10))+
    theme_void()+
    theme(axis.text = element_text(),
          panel.grid.major.y = element_line(color = "grey"),
          panel.grid.minor.x =  element_line(color = "grey"),
          strip.text = element_text(size = 15, hjust=0),
          legend.position = c(0.8, 0.1),
          legend.box = "horizontal")
  ggsave(filename=paste0(here::here(), "/figures/Results_slope_parsBayesian_", temp_scale,".png"),
         plot = last_plot(),
         width=15, height=10)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bayesian pointrange grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
pars_glob <- read_csv(file=paste0(here::here(), "/results/PAranks_Bayesian_global_emtrends.csv"))
pars_cont <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_continental_emtrends.csv"))
pars_regi <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_regional_emtrends.csv"))

pars_all <- rbind(pars_glob, pars_cont) %>%
  rbind(pars_regi)  %>%
  as_tibble()
rm(pars_glob, pars_cont, pars_regi)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bayesian summary ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combined diversity and function estimates (x axis)
# scale on y axis
# habitat types as colors
# pointrange plot

# summarize into groups
pars_sum <- pars_all %>%
  full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
  dplyr::select(-fns, -lower.HPD, -upper.HPD) %>%
  group_by(Group_function, LC, scale) %>%
  summarize(across(everything(), list("mean" = function(x) mean(x, na.rm=T),
                                      "SE" = function(x) sd(x, na.rm=T) / sqrt(length(x))))) %>%
  mutate(# Calculate confidence intervals
    "PA_rank_rev.trend_CI_lower" = PA_rank_rev.trend_mean - (1.96 * PA_rank_rev.trend_SE),  # 1.96 is the Z-value for a 95% confidence interval
    "PA_rank_rev.trend_CI_upper" = PA_rank_rev.trend_mean + (1.96 * PA_rank_rev.trend_SE)
  )

# plot
ggplot(pars_all %>%
         full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
         mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                    ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                           ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
         mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                        ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                           "<img src='figures/icon_location-black.png' width='70'>",
                                                           "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
         mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
         mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
         mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
         filter(!is.na(LC) & !is.na(scale)),
       aes(x = PA_rank_rev.trend, y = LC,
           xmin = lower.HPD, xmax = upper.HPD,
           fill = significance, color = LC,
           linetype = significance))+
  
  geom_vline(aes(xintercept=0), color="black")+
  geom_pointrange(position = position_jitter(height = 0.3),
                  size = 0.5, shape = 21) +

  geom_pointrange(data = pars_sum %>% filter(!is.na(LC) & !is.na(scale)) %>%
                    mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
                    mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                               ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                                      ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
                    mutate(Group_function = ifelse(Group_function=="Service", "Function",
                                                   ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                           scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                                      "<img src='figures/icon_location-black.png' width='70'>",
                                                                      "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
                    mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
                    mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                           significance =  ifelse(sign(PA_rank_rev.trend_CI_lower)!= sign(PA_rank_rev.trend_CI_upper), "ns", as.character(as.factor(LC)))) %>%
                    mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns"))),

                  aes(fill = significance, color = LC,
                      y = LC, x = PA_rank_rev.trend_mean, 
                      xmin = PA_rank_rev.trend_CI_lower, xmax = PA_rank_rev.trend_CI_upper),

                  position = position_dodge(width = 0.1),
                  size = 1.5, shape = 21, linewidth = 1.5) +

  #coord_flip()+
  coord_cartesian(clip = "off")+
  ggh4x::facet_grid2(scale_icon ~ Group_function, drop=FALSE, 
             scales = "free", independent = "all", switch = "y",
             shrink = FALSE)+

  #ylab("Effect size")+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000",
                             "ns" = "white"),
                    drop = FALSE)+
  scale_linetype_manual(values=c("Cropland" = "solid",
                             "Grassland" = "solid",
                             "Shrubland" = "solid", 
                             "Woodland" = "solid", 
                             "Other" = "solid",
                             "ns" = "longdash"),
                    drop = FALSE)+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000",
                              "ns" = "#000000"),
                     drop = FALSE)+
  theme_void()+
  
  guides(fill = guide_legend(reverse = F, override.aes = list(shape = 21, 
                                                              color = c("Cropland" = "#4A2040",
                                                                        "Grassland" = "#E69F00",
                                                                        "Shrubland" = "#0072B2", 
                                                                        "Woodland" = "#009E73", 
                                                                        "ns" = "#000000"))), 
         color = "none")+
  theme(legend.position = "bottom", 
        axis.title.y =element_blank(),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(0.5, "cm"),
        legend.box.spacing = unit(c(1,0,0,0), "cm"),
        legend.title = element_blank(),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "cm"),
        axis.text.x = element_text(size = 15),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.spacing = unit(0.8, "cm"),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        #strip.text.x.bottom = ggtext::element_markdown(vjust = 1),
        strip.text.x = element_text(size = 30, hjust = 0, vjust = 1),
        strip.text.y = ggtext::element_markdown(vjust = 0.5))

ggsave(filename=paste0(here::here(), "/figures/Results_slope_BayesianTrends_allScales_grouped.png"),
       plot = last_plot(), 
       width = 5000, height = 4000,
       units = "px")

# table all
pars_stats <- pars_all %>%
  full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
  mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                 ifelse(Group_function=="Diversity", "Shannon", Group_function)),
         scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                    "<img src='figures/icon_location-black.png' width='70'>",
                                                    "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
  mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
  mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
         significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
  mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
  filter(!is.na(LC) & !is.na(scale)) %>%
  arrange(abs(PA_rank_rev.trend))
table(pars_stats$significance, pars_stats$scale)
table(pars_stats$significance, pars_stats$scale, pars_stats$Group_function)
table(pars_stats[pars_stats$significance != "ns",]$scale)
table(pars_stats$significance)
table(pars_stats[pars_stats$significance != "ns",]$Group_function)
pars_stats %>% print(n=100)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX TABLE 3.1 & 3.2 - Tables with stats ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
options("scipen"=100, "digits"=4)
write_csv(pars_sum %>% 
    mutate("trend" = round(PA_rank_rev.trend_mean, 4),
           "SE" = round(PA_rank_rev.trend_SE, 4),
           "95% CI" = paste0("[", round(PA_rank_rev.trend_CI_lower, 4),"; ", 
                             round(PA_rank_rev.trend_CI_upper, 4), "]"),
           "Scale" = factor(scale, levels = c("global", "continental", "regional"))) %>%
    dplyr::select(Group_function, Scale, LC, trend, SE, "95% CI") %>%
    arrange(Group_function, Scale, LC),
  paste0(here::here(), "/figures/Results_slope_BayesianTrends_allScales_grouped.csv"))

write_csv(pars_all %>%
      full_join(fns_labels %>% dplyr::select(Function, Group_function, Label_short), by = c("fns" = "Function")) %>%
      mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                 ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                        ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
      mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                     ifelse(Group_function=="Diversity", "Shannon", Group_function)),
             scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                        "<img src='figures/icon_location-black.png' width='70'>",
                                                        "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
      mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
      mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
             significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
      mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
      filter(!is.na(LC) & !is.na(scale)) %>%
      arrange(abs(PA_rank_rev.trend)) %>%
      mutate("Habitat" = LC,
             "Group" = Group_function,
             "Slope [HPD]" = paste0(round(PA_rank_rev.trend, 3), 
                                    " [", round(lower.HPD, 3), "; ", round(upper.HPD, 3), "] ", 
                                    ifelse(significance == "ns", "ns", "*")),
             "Variable" = Label_short,
             "Scale_long" = factor(scale, levels = c("global", "continental", "regional")))  %>%
      mutate("Scale" = factor(substr(scale, 1, 4), levels = c("glob", "cont", "regi"))) %>%
      dplyr::select(Group, Variable, Scale, Habitat, "Slope [HPD]") %>%
      pivot_wider(names_from = "Habitat", values_from = "Slope [HPD]") %>%
      arrange(Group, Variable, Scale),
  paste0(here::here(), "/figures/Results_slope_BayesianTrends_allScales.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.2 - Heatmap Bayesian slopes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## plot like heatmap
ggplot(data = pars_all %>%
         full_join(fns_labels %>% dplyr::select(Function, Group_function, Label_short), by = c("fns" = "Function")) %>%
         mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
                                    ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
                                           ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
         mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                        ifelse(Group_function=="Diversity", "Shannon", Group_function)),
                scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
                                                           "<img src='figures/icon_location-black.png' width='70'>",
                                                           "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
         mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
         mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","ns")),
                significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
         mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
         filter(!is.na(LC) & !is.na(scale)) %>%
         arrange(abs(PA_rank_rev.trend)) %>%

         filter(LC %in% lc_names_all) %>%
         mutate(PA_trend_f = cut(abs(PA_rank_rev.trend),
                                    breaks=c(-Inf, -10, -1, 0, 1, 10, 100, Inf),
                                    labels=c("<-10", "<-1", "<0", "<1", "<10", "<100", ">100"))) %>%
         right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
                    by = c("fns" = "Function")) %>%
         full_join(expand.grid(scale = c("global", "continental", "regional"), 
                               LC = lc_names_all, 
                               Label = fns_labels$Label)) %>%
         mutate(scale = factor(scale, levels = rev(c("global", "continental", "regional"))),
                Label = factor(Label, levels = rev(fns_labels %>% arrange(Group_function, Label) %>% pull(Label))),
                LC = factor(LC, levels = c("Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
         filter(LC != "Other" & !is.na(Label)) %>% 
         mutate(effect_direction = as.factor(sign(PA_rank_rev.trend)))%>%
         mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                            ifelse(effect_direction=="1", "positive", "0"))) %>%
         mutate(Label = factor(Label, levels = labels_order)) %>%
         mutate(effect_significance = ifelse(significance == "ns", "ns", effect_direction_c),
                effect_na = ifelse(is.na(PA_trend_f), "not available", NA)) %>%
         mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "ns"))),
       aes(x = LC, y = scale))+
  
  geom_point(aes(size = PA_trend_f,
                 color = effect_direction_c, 
                 fill= effect_significance,
                 shape = effect_na))+
  facet_wrap(vars(Label), ncol=6, drop=FALSE)+
  scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "ns" = "white"),
                    name = "Direction of effect",
                    na.value = "black", drop = FALSE)+
  scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
                     name = "Direction of effect",
                     na.value = "black")+
  scale_size_manual(values =  c("<-10" = 10, "<-1" = 5, "<0" = 2, "<1" = 5, "<10" = 10, "<100" = 15, ">100" = 20), #c("marginal" = 2, "ns" = 5, "small" = 5, "medium" = 10, "large" = 15), 
                     name = "Slope estimate",
                     na.value = 5)+
  scale_shape_manual(values = c("not available" = 4),
                     name = "Missing data",
                     na.value = 21)+
  scale_x_discrete(labels = c(
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
  ))+
  
  scale_y_discrete(labels = c(
    "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
    "continental" = "<img src='figures/icon_location-black.png' width='30'>",
    "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
  ))+
  xlab("")+ylab("")+
  theme_bw() + # use a white background
  
  guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black"), 
                                                 shape = 21, size = 5)), #tell legend to use different point shape
         color = "none", #don't show legend
         shape = guide_legend(override.aes = list(size = 5)))+
  theme(
    legend.position = "bottom", #c(0.96, -0.05),
    legend.justification = c(1, 0),
    legend.box = "horizontal",
    legend.direction = "vertical",
    #legend.key.size = unit(20, "pt"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    #axis.title.y =element_blank(),
    axis.text.y = ggtext::element_markdown(hjust = 0),
    axis.ticks = element_blank(),
    axis.text.x = ggtext::element_markdown(vjust = 0),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.background = element_rect(fill="white", color = "white"), #chocolate4
    strip.text = element_text(color="black", size = 15, hjust = 0)) #white
ggsave(filename=paste0(here::here(), "/figures/Results_slope_BayesianTrends_allScales_fns.png"),
       plot = last_plot(), 
       width = 4400, height = 3800, units = "px")


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Raw data plot ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boxplots values ####

# extract list of sampling locations actually used in comparison
data_values <- data_clean %>% 
  filter(SampleID %in% unique(pa_pairs$ID) | 
           SampleID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(SampleID, PA, LC, all_of(fns)) %>%
  mutate(across(all_of(fns), .fns = function(x) { (x - mean(x)) / sd(x)})) %>%
  
  pivot_longer(cols = c(fns), names_to = "fns") %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
data_values

data_values <- data_values %>% filter(LC %in% lc_names)

ggplot(data_values, aes(x = Label, y = value)) +
  geom_hline(yintercept = 0, lty="dashed", color="grey60") +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  
  xlab("")+ ylab("")+
  facet_wrap(vars(LC), ncol=1, nrow=4)+
  theme_bw() +
  theme(axis.text.x=element_text(size=5, angle=45, hjust=1),
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", legend.text = element_text(size=15), axis.text.y = element_text(size=15))
ggsave(filename=paste0(here::here(), "/figures/Results_boxplot_estimates_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bi-plots ####
# protected area categories against functions

# prepare data
data_long <- data_clean %>%
  pivot_longer(cols = all_of(fns), names_to = "fns", values_to = "fns_value") %>%
  filter(LC %in% lc_names) %>%
  right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
            by = c("fns" = "Function")) %>%
  mutate(Label = factor(Label, levels = labels_order))

# plotting
ggplot(data = data_long,
       aes(x = as.factor(PA_rank), y = fns_value, group = LC, color = LC))+
  xlab("Category of protected area (ranked)")+
  ylab("")+
  geom_jitter()+
  geom_smooth(method = "lm", stat = "smooth")+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000"), name="Habitat type")+
  facet_wrap(vars(Label), scales = "free_y", ncol = 6)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15))
ggsave(filename=paste0(here::here(), "/figures/Data_clean_PAranks_fns_", temp_scale, ".png"),
       plot = last_plot(), 
       width = 5500, height = 4000,
       units = "px")

