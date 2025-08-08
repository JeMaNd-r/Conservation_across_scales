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
library(geosphere) # calculate spatial distance
library(ggpubr) #plot correlation in ggplot
library(pwr) #power tests

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

# define each land cover type
lc_names_all <- c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "Other")

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
## Sampling locations ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Locations per habitat type ####
data_locations_all <- data.frame()
data_locations_Bayes <- data.frame()

for(temp_scale in c("global", "continental", "regional")){
  # load locations
  temp_data <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  pa_pairs <- read_csv(file=sort(list.files(here::here("intermediates", temp_scale), pattern = "Pairs_paNonpa_1000trails_", full.names = TRUE), decreasing = TRUE)[1])
  
  temp_subset <- temp_data %>% 
    dplyr::filter(SampleID %in% unique(pa_pairs$ID) | 
             SampleID %in% unique(pa_pairs$nonPA)) %>%
    dplyr::select(Longitude,Latitude,SampleID, PA, PA_type, LC)
  temp_data <- temp_data %>%
    dplyr::select(Longitude,Latitude,SampleID, PA, PA_type, LC)
  
  # add scale
  temp_data$scale <- temp_scale
  temp_subset$scale <- temp_scale 
  
  # add to full dataframe
  data_locations_all <- rbind(data_locations_all, temp_subset)  
  data_locations_Bayes <- rbind(data_locations_Bayes, temp_data)  
}

data_locations_all
data_locations_Bayes

data_locations_all <- data_locations_all %>%
  mutate(LC_f = as.factor(LC)) %>%
  mutate(LC_f = factor(LC_f, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland")))

# plot
ggplot()+
  xlab("")+
  ylab("Number of sampling sites")+
   geom_bar(data = data_locations_all %>%
             mutate(scale = factor(scale, levels = c("global", "continental", "regional"))), 
           aes(x = LC_f, fill = LC_f, alpha = "Unprotected"),
           position = "dodge",
           na.rm=TRUE)+
   scale_fill_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000",
                                "Dryland" = "#000000"))+ 
  
   geom_bar(data = data_locations_all %>%
             filter(PA == 1) %>%
             mutate(scale = factor(scale, levels = c("regional", "continental", "global"))), 
           aes(x = LC_f, alpha = "Protected"), fill = "white", 
           position = "dodge",
           na.rm=TRUE)+
  scale_alpha_manual(values = c("Protected" = 0.5, "Unprotected" = 1))+
  guides(alpha = guide_legend(override.aes = list(fill = "black",
                                                  alpha = c(0.1, 1))))+
  scale_x_discrete(expand = c(0,0.7),
                   drop = TRUE)+
  scale_y_continuous(expand = c(0,0), limits = c(0, 160))+
  facet_grid(cols = vars(scale), drop = TRUE, scales = "free_x", space = "free_x")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "grey80"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey60"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1, "cm"),
        axis.text.y = element_text(size = 19.5, vjust = 0.5),
        axis.title.y = element_text(size = 25.5),
        axis.text.x = element_blank())
        #axis.text.x = ggtext::element_markdown(vjust = 0))
ggsave(paste0(here::here(), "/figures/Data_habitats_allScales.png"),
       width = 10, height = 5,
       plot = last_plot())

# table with number of sites per PA type
write_csv(data_locations_all %>% 
    group_by(scale, LC, PA_type, PA) %>%
    count() %>%
    mutate(PA_type = factor(PA_type, levels = c("Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Applicable", "Not Assigned", "Not Reported", NA))) %>%
    pivot_wider(names_from = scale, values_from = n) %>%
    ungroup() %>%
    dplyr::select(LC, PA_type, regional, continental, global) %>%
    arrange(LC, PA_type) %>%
    rename("Habitat type" = LC, "Protection type" = PA_type),
  paste0(here::here(), "/figures/Data_habitats_allScales.csv")
)

write_csv(data_locations_Bayes %>% 
            group_by(scale, LC, PA_type, PA) %>%
            count() %>%
            mutate(PA_type = factor(PA_type, levels = c("Ia", "Ib", "II", "III", "IV", "V", "VI", "Not Applicable", "Not Assigned", "Not Reported", NA))) %>%
            pivot_wider(names_from = scale, values_from = n) %>%
            ungroup() %>%
            dplyr::select(LC, PA_type, regional, continental, global) %>%
            arrange(LC, PA_type) %>%
            rename("Habitat type" = LC, "Protection type" = PA_type),
          paste0(here::here(), "/figures/Data_habitats_allScales_Bayesian.csv")
)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Mahalanobis distance ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Boxplot Mahalanobis distance ####

for(temp_scale in c("global", "continental", "regional")){
  
  # load pairs of PA and nonPA
  pa_pairs <- read_csv(file=sort(list.files(here::here("intermediates", temp_scale), pattern = "Pairs_paNonpa_1000trails_", full.names = TRUE), decreasing = TRUE)[1])
  
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
                               "Other" = "#000000",
                               "Dryland" = "#000000")) #"gold3", "limegreen", "forestgreen"
  ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_", temp_scale, ".png"),
         plot = last_plot(),
         width = 10, height = 8)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### SUPPLEMENTARY FIGURE S2 - Boxplot Mahalanobis distance - all vs. pairs ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# calculate & plot Mahalanobis distance for all
for(temp_scale in c("global", "continental", "regional")){
  data <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  
  if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names_all[lc_names_all != "Other" & lc_names_all != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
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
                               "Other" = "#000000",
                               "Dryland" = "#000000")) #"gold3", "limegreen", "forestgreen"
  ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_", temp_scale, ".png"),
         plot = last_plot(),
         width = 10, height = 8)
}

# put in 1 plot - APPENDIX FIGURE 2.1
pairs_glob <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_global.png"))
pairs_cont <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_continental.png"))
pairs_regi <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_regional.png"))
all_glob <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_global.png"))
all_cont <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_continental.png"))
all_regi <- magick::image_read(paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_regional.png"))

# Create a spacer (e.g., 50px wide, same height as images)
spacer <- image_blank(width = 0.2*(image_info(all_glob)$width), height = image_info(all_glob)$height, color = "white")

magick::image_write(magick::image_append(c(magick::image_append(c(pairs_regi, pairs_cont, spacer, pairs_glob)),
                                           magick::image_append(c(all_regi, all_cont, spacer, all_glob))),
                     stack = TRUE),
                    path=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_all_allScales.png"))
rm(pairs_glob, pairs_cont, pairs_regi, all_glob, all_cont, all_regi)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
## Effect sizes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for(temp_scale in c("global", "continental", "regional")){
  load(file=paste0(here::here(), "/results/d_1000_trails_", temp_scale, ".RData")) #d_list

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
                     list(mean=function(x) mean(abs(x), na.rm=TRUE),
                          sd=function(x) sd(abs(x), na.rm=TRUE))))
  
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
}


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

# post-hoc averaging... be carefull with interpretation
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
         lc = factor(lc, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
  filter(lc != "Other" & !is.na(Label)) %>% 
  mutate(effect_direction = as.factor(sign(effect_mean)))%>%
  mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                     ifelse(effect_direction=="1", "positive", "0"))) %>%
  mutate(Label = factor(Label, levels = labels_order)) %>%
  mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "not significant", effect_direction_c),
         effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
  mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant")))

# save table - APPENDIX TABLE 2.2
write_csv(d_plot_all %>% 
            dplyr::select(-X) %>%
            dplyr::select(Group_function, Label, scale, lc, effect_mean:effect_ci_97.5, effect_significance) %>%
            arrange(Group_function, Label, scale, lc) %>%
            mutate(effect_significance = ifelse(effect_significance=="not significant", "", "*")) %>%
            mutate("Habitat" = lc,
                   "Group" = Group_function,
                   "Slope [HPD]" = paste0(round(effect_mean, 3), 
                                          " [", round(effect_ci_2.5, 3), "; ", round(effect_ci_97.5, 3), "] ", effect_significance),
                   "Variable" = Label,
                   "Scale_long" = factor(scale, levels = c("global", "continental", "regional")))  %>%
            mutate("Scale" = factor(substr(scale, 1, 4), levels = c("glob", "cont", "regi"))) %>%
            dplyr::select(Group, Variable, Scale, Habitat, "Slope [HPD]") %>%
            pivot_wider(names_from = "Habitat", values_from = "Slope [HPD]") %>%
            arrange(Group, Variable, Scale) %>%
            filter(!is.na(Group)),
          paste0(here::here(), "/figures/Results_d-value_meanCI_allScales_fns.csv"))

# # plot - FIGURE 2
# ggplot(data = d_plot_all %>%
#          filter(lc != "Shrubland"),
#        aes(x = lc, y = scale))+
#   
#   geom_point(aes(size = effect_mean_f,
#                  color = effect_direction_c, 
#                  fill= effect_significance,
#                  shape = effect_na))+
#   facet_wrap(vars(Label), ncol=6, drop=FALSE)+
#   scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "not significant" = "white"),
#                     name = "Direction of effect",
#                     na.value = "black", drop = FALSE)+
#   scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
#                     name = "Direction of effect",
#                     na.value = "black")+
#   scale_size_manual(values = c("marginal" = 2, "not significant" = 2, "small" = 5, "medium" = 10, "large" = 15),
#                      name = "Effect size",
#                     na.value =1)+
#   scale_shape_manual(values = c("not available" = 4),
#                      name = "Missing data",
#                      na.value = 21)+
#   scale_x_discrete(labels = c(
#     "Dryland" = "<img src='figures/icon_land.png' width='30'>",
#     "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
#     "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
#     #"Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
#     "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
#   ))+
#   
#   scale_y_discrete(labels = c(
#     "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
#     "continental" = "<img src='figures/icon_location-black.png' width='30'>",
#     "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
#   ))+
#   xlab("")+ylab("")+
#   theme_bw() + # use a white background
#   
#   guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black", "black"),
#                                                  shape = 21, size = 5)), #tell legend to use different point shape
#     color = "none", #don't show legend
#     shape = guide_legend(override.aes = list(size = 5)))+
#   theme(legend.position = "bottom", #c(0.96, -0.05),
#         legend.justification = c(1, 0),
#         legend.box = "horizontal",
#         legend.direction = "vertical",
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         axis.text.y = ggtext::element_markdown(hjust = 0),
#         axis.ticks = element_blank(),
#         axis.text.x = ggtext::element_markdown(vjust = 0),
#         panel.grid = element_blank(),
#         panel.border = element_blank(),
#         strip.background = element_rect(fill="white", color = "white"), #chocolate4
#         strip.text = element_text(color="black", size = 15, hjust = 0)) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_d-value_meanCI_allScales_fns.png"),
#        plot = last_plot(), 
#        width = 4400, height = 3800, units = "px")

## Pointrange
ggplot(data = d_plot_all %>% 
         filter(lc != "Shrubland" & lc != "Dryland" & scale != "global") %>%
         filter(!is.na(effect_significance)) %>%
         mutate(Label = factor(Label, levels = rev(fns_labels$Label))) %>%
         mutate(Label = ifelse(Group_function == "Function", as.character(Label), 
                               ifelse(Label_short == "AM fungi R", "Fungi (AM)",
                                      ifelse(Label_short == "EM fungi R", "Fungi (EM)",
                                             Organism)))) %>%
         mutate(lc_icon = ifelse(lc == "Cropland", "<img src='figures/icon_harvest.png' width='40'>",
                                    ifelse(lc == "Grassland", "<img src='figures/icon_grass.png' width='34'>",
                                           ifelse(lc == "Woodland", "<img src='figures/icon_forest.png' width='60'>", NA)))) %>%
         mutate(Label = factor(Label, rev(sort(unique(Label)))),
           Group_function = factor(Group_function, levels = rev(c("Dissimilarity", "Shannon", "Richness", "Function"))),
           lc_icon = factor(lc_icon, levels = c("<img src='figures/icon_harvest.png' width='40'>",
                                                           "<img src='figures/icon_grass.png' width='34'>",
                                                           "<img src='figures/icon_forest.png' width='60'>" ))),
       aes(y = effect_mean, x = Label))+
  
  geom_hline(aes(yintercept=0), linetype = "dashed")+  
  geom_pointinterval(aes(ymin = effect_ci_17, ymax = effect_ci_83, 
                         color = effect_significance),
                     linewidth = 12) +
  geom_pointinterval(aes(ymin = effect_ci_2.5, ymax = effect_ci_97.5,
                         fill = effect_significance, color = effect_significance), 
                     shape = 21, size = 10,
                     linewidth = 6) + 
  coord_flip()+
  ylab("Effect size")+
  facet_grid(Group_function ~ scale + lc_icon,
             scales = "free",
             axis.labels = "all_x",
             space = "free_y")+
  # facet_wrap2(vars(scale, lc_icon), nrow = 1, 
  #             #strip = ggh4x::strip_nested(clip = "off"),
  #             strip.position = "top")+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000",
                             "Dryland" = "#000000",
                             "not significant" = "grey50",
                             "positive" = "black",
                             "negative" = "black"))+
  scale_fill_manual(values = c("positive" = "#91bfdb", "negative" = "#fc8d59",
                               "not significant" = "white"),
                    name = "Effect direction")+
  theme_bw() + # use a white background
  guides(fill = guide_legend(override.aes = list(color = c("black", "black","grey50"))))+
  theme(legend.position = "right", axis.title.y =element_blank(),
        legend.title = element_text(size=18, color = "black"),
        legend.text = element_text(size=18, color = "black"),
        axis.text.y = element_text(size=18, color = "black"),  
        axis.text.x = element_text(size=10, color = "gray80"),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=18, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray80"),
        rect = element_blank(),  
        strip.background.x = element_rect(fill="white", color = NA), #chocolate4
        strip.background.y = element_rect(fill="white", color = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size=18, color = "black", hjust = 0),
        strip.text.x = ggtext::element_markdown(vjust = 0)) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_meanCI_allScales_fns.png"),
       plot = last_plot(),
       dpi = 300,
       width=16, height=10)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### Summarizing stats ####
# number significant effects
d_plot_all %>% filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% nrow() #38
# number ns
d_plot_all %>% filter(effect_significance=="not significant" & !is.na(effect_significance)) %>% nrow() #109

# number significant per lc
table(d_plot_all %>% filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% dplyr::select(scale, lc))
table(d_plot_all %>% filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% dplyr::select(effect_direction_c))
table(d_plot_all %>% filter(effect_significance!="not significant" & !is.na(effect_significance)) %>% dplyr::select(effect_direction_c, Group_function))
table(d_plot_all %>% filter(!is.na(effect_significance)) %>% dplyr::select(scale, lc))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
#### individual distributions ####

d_cont <- get(load(paste0(here::here(), "/results/d_1000_trails_continental.RData"))) #d_list
d_regi <- get(load(paste0(here::here(), "/results/d_1000_trails_regional.RData")))

d_cont <- do.call(rbind, d_cont)
d_cont$scale <- "continental"
d_regi <- do.call(rbind, d_regi)
d_regi$scale <- "regional"

d_all <- rbind(d_cont, d_regi)
str(d_all)

d_all <- d_all %>%
  right_join(fns_labels %>% dplyr::select(Label, Label_short, Function, Group_function, Organism), 
             by = c("fns" = "Function"))

ggplot(d_all[order(d_all$lower),] %>% filter(fns == "Soil_carbon_service") , aes(x = effect, y = run)) +
  geom_point(alpha = 1) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, alpha = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  facet_grid(vars(lc, scale))+
  coord_flip()+
  labs(x = "Effect size", y = "Trail")

p_dist <- 
  ggplot(d_all %>%
         mutate(Label = factor(Label, levels = rev(fns_labels$Label))) %>%
         mutate(Label = ifelse(Group_function == "Function", as.character(Label), 
                               ifelse(Label_short == "AM fungi R", "Fungi (AM)",
                                      ifelse(Label_short == "EM fungi R", "Fungi (EM)",
                                             Organism)))) %>%
          mutate(lc_icon = ifelse(lc == "Cropland", "<img src='figures/icon_harvest.png' width='40'>",
                               ifelse(lc == "Grassland", "<img src='figures/icon_grass.png' width='34'>",
                                      ifelse(lc == "Woodland", "<img src='figures/icon_forest.png' width='60'>", NA)))) %>%
         mutate(Label = factor(Label, rev(sort(unique(Label)))),
                Group_function = factor(Group_function, levels = rev(c("Dissimilarity", "Shannon", "Richness", "Function"))),
                lc_icon = factor(lc_icon, levels = c("<img src='figures/icon_harvest.png' width='40'>",
                                                     "<img src='figures/icon_grass.png' width='34'>",
                                                     "<img src='figures/icon_forest.png' width='60'>" ))) %>%
         filter(!is.na(lc_icon)), 
       aes(x = effect, y = Label)) +
  geom_point(alpha = 0.02, size = 4, shape = 15) +
  geom_linerange(aes(xmin = lower, xmax = upper), linewidth = 2, alpha = 0.01) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_bw() +
  #coord_flip()+
  facet_grid(Group_function ~ scale + lc_icon,
             scales = "free",
             axis.labels = "all_x",
             space = "free_y")+
  theme(legend.position = "right", axis.title.y =element_blank(),
        legend.title = element_text(size=18, color = "black"),
        legend.text = element_text(size=18, color = "black"),
        axis.text.y = element_text(size=18, color = "black"),  
        axis.text.x = element_text(size=10, color = "gray80"),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "gray80"),
        rect = element_blank(),  
        strip.background.x = element_rect(fill="white", color = NA), #chocolate4
        strip.background.y = element_rect(fill="white", color = NA),
        strip.placement = "outside",
        strip.text.y = element_text(size=18, color = "black", hjust = 0),
        strip.text.x = ggtext::element_markdown(vjust = 0)) #white
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_dist_allScales_fns.png"),
       plot = p_dist,
       dpi = 300,
       width=16, height=10)

#### individual distributions - stats ####

d_all %>% filter(sign(upper)==sign(lower)) %>%
  filter(Organism == "Bacteria") %>% count(lc, scale, Label)

d_all %>% filter(sign(upper)==sign(lower)) %>%
  filter(Organism == "Invertebrates") %>% count(lc, scale, Label)

d_all %>% filter(sign(upper)==sign(lower)) %>%
  filter(Organism == "Fungi") %>% count(lc, scale, Label)

d_all %>% filter(sign(upper)==sign(lower)) %>%
  filter(Organism == "Function") %>% count(lc, scale, Label)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SUPPLEMENTARY FIGURE 3 - Heatmap correlation Mahalanobis & difference ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

all_corr <- data.frame("LC" = "lc", "fns" = "fns", "correlation" = 1, "scale" = "scale", "p_value" = 1)[0,]

for(temp_scale in c("global", "continental", "regional")){

  if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names_all[lc_names_all != "Other" & lc_names_all != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
  }
  
  data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  
  # load pairs of PA and nonPA
  pa_pairs <- read_csv(file=sort(list.files(here::here("intermediates", temp_scale), pattern = "Pairs_paNonpa_1000trails_", full.names = TRUE), decreasing = TRUE)[1])
  
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
  
  if(temp_scale == "global"){ #necessary probably because of only one LC...
    correlation_matrix <- t(correlation_matrix)
    correlation_p_matrix <- t(correlation_p_matrix)
  } 

  
  # Convert the correlation matrix to a data frame for easier interpretation
  correlation_df <- as.data.frame(correlation_matrix)
  correlation_p_df <- as.data.frame(correlation_p_matrix)
  
  rownames(correlation_df) <- lc_names 
  colnames(correlation_df) <- fns 
  
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
         scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='35'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='35'>",
                                    ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='35'>", NA)))) %>%
  full_join(fns_labels %>% 
              filter(Function %in% fns),
            dplyr::select(Label_short, Function), 
            by = c("fns" = "Function")) %>%
  mutate(Label_short = factor(Label_short, levels = rev(fns_labels$Label_short))) 

# add d-values (effect size results)
all_corr_plot <- all_corr_plot %>%
  full_join(d_plot_all %>%
              filter(!is.na(effect_significance)) %>%
              rename(LC = lc) %>%
              dplyr::select(Label_short, LC, fns, Group_function, effect_significance, scale)) %>%
  filter(!is.na(correlation)) %>%
  bind_rows(data.frame(LC = NA, fns = "Soil_carbon_service", correlation = NA, scale = NA, p_value = NA, 
                       scale_icon = " ", Label = "Soil carbon service", Label_short="Soil carbon", Group_function ="Function", Organism = NA)) %>%
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_flag-Portugal.png' width='35'>",
                                                    "<img src='figures/icon_location-black.png' width='35'>",
                                                    " ",
                                                    "<img src='figures/icon_earth-globe-with-continents-maps.png' width='35'>")))

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
            dplyr::select(Group, Variable, Scale, Dryland, Cropland, Grassland, Woodland) %>%
            arrange(Group, Variable, Scale),
          paste0(here::here(), "/figures/Correlation_diff_mahal_allScales.csv"))

# plotting
ggplot(data = all_corr_plot %>%
         mutate(LC = factor(LC, levels = c("Cropland", "Grassland", "Woodland", "Dryland"))))+
  geom_tile(aes(x = LC, y = Label_short, fill = correlation))+
  
  geom_point(data = all_corr_plot %>% filter(p_value >= 0.05),
             aes(x = LC, y = Label_short), shape = 4)+
  
  geom_point(data = all_corr_plot %>% filter(effect_significance != "not significant"),
             aes(x = LC, y = Label_short), shape = 0, size = 8)+
  
  scale_shape_manual(values = c(NA, 4))+
  scale_fill_distiller(type = "div", na.value = "white", limits = c(-0.7, 0.7),
                       name = "Correlation")+
  scale_x_discrete(labels = c(
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>",
    "Dryland" = "<img src='figures/icon_land.png' width='33'>"
  ))+
  facet_grid(cols = vars(scale_icon), scales = "free_x", space = "free_x",
             drop = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 15),
        strip.text.x = ggtext::element_markdown(vjust = 0.5),
        strip.background = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text.x = ggtext::element_markdown(vjust = 0),
        panel.spacing.x = unit(1, "lines"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())
ggsave(paste0(here::here(), "/figures/Correlation_diff_mahal_allScales.png"),
       last_plot(),
       height = 10, width = 8)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Spatial distance and difference in attributes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for(temp_scale in c("global", "continental", "regional")){
  
  if(temp_scale == "global"){
    lc_names <- "Dryland" #lc_names_all[lc_names_all != "Other" & lc_names_all != "Cropland"]
  } 
  if(temp_scale == "continental"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
  }
  if(temp_scale == "regional"){
    lc_names <- lc_names_all[lc_names_all != "Other" & lc_names_all != "Shrubland" & lc_names_all != "Dryland"]
  }
  
  data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
  
  # load pairs of PA and nonPA
  pa_pairs <- read_csv(file=sort(list.files(here::here("intermediates", temp_scale), pattern = "Pairs_paNonpa_1000trails_", full.names = TRUE), decreasing = TRUE)[1])
  
  # calculate geographic distance
  pa_dist <- pa_pairs %>% 
    full_join(data_clean %>% dplyr::select(SampleID, Longitude, Latitude), by = c("nonPA" = "SampleID")) %>% 
    full_join(data_clean %>% dplyr::select(SampleID, Longitude, Latitude), by = c("ID" = "SampleID"), suffix = c(".pa", ".nonpa")) %>%
    mutate(distance = geosphere::distHaversine(
      cbind(Longitude.pa, Latitude.pa), # First site in the pair
      cbind(Longitude.nonpa, Latitude.nonpa)) /1000 ) %>% # Second site in the pair, transform from meters in kilometers
    filter(!is.na(nonPA) & !is.na(ID))
  
  # calculate difference for fns columns and mean for mahal.min column (= mahal.min because same for both) within each group
  pa_env <- pa_pairs %>% full_join(data_clean, by = c("nonPA" = "SampleID", "LC")) %>% mutate(data = "nonPA") %>%
    rbind(pa_pairs %>% full_join(data_clean, by = c("ID" = "SampleID", "LC")) %>% mutate(data = "PA")) %>%
    filter(!is.na(nonPA) & !is.na(ID)) %>%
    arrange(ID, nonPA, times) 
  
  pa_env <- pa_env %>% group_by(ID, nonPA, LC, times, n) %>% 
    summarize(across(all_of(fns), diff), #calculates nonPA-PA -> positive means PA larger/better
              across(mahal.min, mean)) %>% 
    dplyr::select(ID, nonPA, mahal.min, LC, all_of(fns))
  
  # merge with distances
  pa_env <- pa_env %>%
    full_join(pa_dist %>% dplyr::select(ID, nonPA, distance), by = c("nonPA", "ID")) %>% 
    unique() %>%
    ungroup() %>%
    pivot_longer(cols = all_of(fns), names_to = "Function")

  # Compute unique correlation values per Function
  cor_values <- pa_env %>%
    group_by(Function) %>%
    summarize(cor_value = cor.test(distance, abs(value), use = "complete.obs", method = "spearman")[["estimate"]],
              cor_p = cor.test(distance, abs(value), use = "complete.obs", method = "spearman")[["p.value"]],
              cor_t = cor.test(distance, abs(value), use = "complete.obs", method = "spearman")[["statistic"]])
  
  write_csv(cor_values, 
            paste0(here::here(), "/figures/Correlation_diff_distance_", temp_scale, ".csv"))
  
  # plot distance x difference
  ggplot(data = pa_env,
         aes(x = distance, y = abs(value)))+
    geom_point(cex = 0.5)+
    facet_wrap(vars(Function), scales = "free")+
    ggpubr::stat_cor(method = "spearman", color = "red", size = 2) + 
    #geom_abline(data = cor_values, aes(slope = cor_value, intercept = 0), linetype = "dashed")+ # Line from correlation
    xlab("Distance in km")+ 
    ylab("Absolute difference in attribute values")+
    theme_bw()
  ggsave(paste0(here::here(), "/figures/Correlation_diff_distance_", temp_scale, ".png"),
         last_plot(),
         height = 15, width = 15)
}
  
# merge cor_values table
corr_list <- list(
  read_csv(paste0(here::here(), "/figures/Correlation_diff_distance_", "global", ".csv")) %>% mutate(scale = "global"),
  read_csv(paste0(here::here(), "/figures/Correlation_diff_distance_", "continental", ".csv")) %>% mutate(scale = "continental"),
  read_csv(paste0(here::here(), "/figures/Correlation_diff_distance_", "regional", ".csv")) %>% mutate(scale = "regional")
)

corr_list <- do.call(rbind, corr_list) %>%
  dplyr::select(scale, Function, cor_value, cor_p, cor_t)

corr_list %>% filter(cor_p < 0.05) %>% count(scale) # number of significant cor
corr_list %>% filter(cor_p < 0.05) %>% count()
corr_list %>% filter(cor_value > 0 & cor_p < 0.05) %>% count(scale) #number positive cor
corr_list %>% filter(cor_value > 0 & cor_p < 0.05) %>% count()
corr_list %>% filter(cor_value < 0 & cor_p < 0.05) %>% count(scale) #number negative cor
corr_list %>% filter(cor_value < 0 & cor_p < 0.05) %>% count()

write_csv(corr_list, 
          paste0(here::here(), "/figures/Correlation_diff_distance_allScales.csv"))

corr_list <- read_csv(paste0(here::here(), "/figures/Correlation_diff_distance_allScales.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SUPPLEMENTARY FIGURE S4 & TABLE 2.3 - Pointrange plot grouped per estimate type ####
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
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_flag-Portugal.png' width='30'>", 
                                          "<img src='figures/icon_location-black.png' width='30'>", 
                             " ",
                             "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>"))) %>%
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
            filter(scale_icon != " ") %>%
            dplyr::select(-effect_ci2.5, -effect_ci92.5), 
          paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_grouped.csv"))

# plot
ggplot(data = d_df_grouped %>%
         bind_rows(data.frame(lc = NA, lower = NA, effect = NA, upper = NA, fns = NA, run = NA, Label = NA, Label_short = NA,
                   Group_function = "Function", Organism = NA, scale = NA, scale_icon = " ")),
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
                             "Other" = "#000000",
                             "Dryland" = "#000000"))+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000",
                             "Dryland" = "#000000"))+
  scale_x_continuous(breaks = c(-2, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 2))+
  scale_y_discrete(limits = c("<img src='figures/icon_flag-Portugal.png' width='30'>", 
                              "<img src='figures/icon_location-black.png' width='30'>", 
                              " ",
                              "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>"))+
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
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing.x = unit(40, "pt"))

#ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_groupedOrganisms.png"), #switch facet_wrap to Organism
ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianCI_allScales_grouped.png"),
       plot = last_plot(), 
       width = 2700, height = 2200,
       units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Power analysis for effect sizes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Cohen's D is difference between 2 means -> t.test function for power test

# sample size = number of protected areas
sample_size <- read_csv(paste0(here::here(), "/figures/Data_habitats_allScales.csv"))
sample_size <- sample_size %>%
  mutate(PA = ifelse(is.na(`Protection type`), 0, 1)) %>%
  group_by(`Habitat type`, PA) %>%
  summarize(across(c(regional:global),  \(x) sum(x, na.rm = TRUE)))

for(temp_scale in c("global", "continental", "regional")){
  print("################################################")
  print(temp_scale)
  d_sum <- read_csv(file=paste0(here::here(), "/figures/Results_d-value_summary_", temp_scale, ".csv"))
 
  d_sum <- d_sum %>%
    full_join(sample_size %>% filter(PA == 1) %>%
                pivot_longer(regional:global, names_to = "scale", values_to = "sample_size") %>%
                filter(scale == temp_scale) %>%
                dplyr::select(-PA, -scale),
              by = c("lc" = "Habitat type"))
  d_sum$power <- NA
  
  # power given sample size n and effect size d
  d_sum$power <- pwr::pwr.t.test(d = d_sum$effect_median, n = d_sum$sample_size, sig.level = 0.05, type = "two.sample", alternative = "two.sided")$power
  
  # number of effect sizes that meet power threshold of >=0.8
  print(sum(d_sum$power >= 0.8, na.rm=TRUE))
  print(sum(d_sum$power >= 0.6, na.rm=TRUE))
  
  write_csv(d_sum, file=paste0(here::here(), "/figures/Results_d-value_summary+power_", temp_scale, ".csv"))
}


d_power_glob <- read_csv(paste0(here::here(), "/figures/Results_d-value_summary+power_", "global", ".csv"))
d_power_cont <- read_csv(paste0(here::here(), "/figures/Results_d-value_summary+power_", "continental", ".csv"))
d_power_regi <- read_csv(paste0(here::here(), "/figures/Results_d-value_summary+power_", "regional", ".csv"))

d_power_all <- rbind(d_power_glob %>% mutate("scale" = "global"), 
                     d_power_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_power_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_power_glob, d_power_cont, d_power_regi)
d_power_all

#### Summary stats power analysis ####
# number significant effects
d_power_all %>% filter(sign(effect_ci_2.5)== sign(effect_ci_97.5)) %>% nrow() #31
# number ns
d_power_all %>% filter(sign(effect_ci_2.5)!= sign(effect_ci_97.5)) %>% nrow() #116

d_power_all <- d_power_all %>%
  mutate(effect_direction = sign(effect_median))

# number significant per lc
table(d_power_all %>% filter(sign(effect_ci_2.5)== sign(effect_ci_97.5)) %>% dplyr::select(scale, lc))
table(d_power_all %>% filter(sign(effect_ci_2.5)== sign(effect_ci_97.5) & power>=0.8) %>% dplyr::select(scale, lc))
table(d_power_all %>% filter(sign(effect_ci_2.5)== sign(effect_ci_97.5)) %>% dplyr::select(effect_direction))
table(d_power_all %>% filter(sign(effect_ci_2.5)== sign(effect_ci_97.5)) %>% dplyr::select(effect_direction , Group_function))
table(d_power_all %>% dplyr::select(scale, lc))

d_power_all %>% arrange(desc(power))

# # minimum detectable effect size
pwr::pwr.t.test(n = 14,
           power = 0.8, # typically desired power
           sig.level = 0.05, type = "two.sample", alternative = "two.sided")
pwr::pwr.t.test(n = 7,
                power = 0.8, # typically desired power
                sig.level = 0.05, type = "two.sample", alternative = "two.sided")


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Model comparison: random-slope & -intercept models ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
model_comparison_table <- data.frame()

for(temp_scale in c("global", "continental", "regional")){
  
  load(file=paste0(here::here(), "/intermediates/PAranks_ModelEval_", temp_scale, ".RData"))
  
  model_comparison <- lapply(model_comparison, function(x){
    
    # Extracting values
    x <- tibble(
      model = c("intercept", "slope"),
      elpd_loo = c(x$loo$intercept$elpd_loo, x$loo$slope$elpd_loo),
      p_loo = c(x$loo$intercept$p_loo, x$loo$slope$p_loo),
      looic = c(x$loo$intercept$looic, x$loo$slope$looic),  # Lower looic is better
      se_elpd_loo = c(x$loo$intercept$se_elpd_loo, x$loo$slope$se_elpd_loo), 
      se_p_loo = c(x$loo$intercept$se_p_loo, x$loo$slope$se_p_loo),
      se_looic = c(x$loo$intercept$se_looic, x$loo$slope$se_looic),
      elpd_loo_diff = c(x$loo$comparison[,1]), 
      se_loo_diff = c(x$loo$comparison[,2]), #If elpd_diff / se_diff > 2, the difference is statistically meaningful.
      elpd_waic = c(x$waic$intercept$elpd_waic, x$waic$slope$elpd_waic),
      p_waic = c(x$waic$intercept$p_waic, x$waic$slope$p_waic),
      waic = c(x$waic$intercept$waic, x$waic$slope$waic),
      se_elpd_waic = c(x$waic$intercept$se_elpd_waic, x$waic$slope$se_elpd_waic),
      se_p_waic = c(x$waic$intercept$se_p_waic, x$waic$slope$se_p_waic),
      se_waic = c(x$waic$intercept$se_waic, x$waic$slope$se_waic),
      r2 = c(x$r2$intercept[1], x$r2$slope[1]),
      r2_se = c(x$r2$intercept[2], x$r2$slope[2]),
      r2_q2.5 = c(x$r2$intercept[3], x$r2$slope[3]),
      r2_q97.5 = c(x$r2$intercept[4], x$r2$slope[4]))
    
    x <- x %>%
      mutate(loo_comp = ifelse(abs(x$elpd_loo_diff / x$se_loo_diff) > 2, 
                               "worse", "similar"), .before = 2)
    
    return(x)
  })
  
  
  model_comparison <- model_comparison %>% 
    imap_dfr(~ mutate(.x, fns = .y, .before = 1)) %>%
    mutate(scale = temp_scale, .before = 1)
  
  model_comparison_table <- rbind(model_comparison_table, model_comparison)
}

# look if there are "worse" models
model_comparison_table %>%
  filter(!is.na(loo_comp)) %>%
  arrange(desc(loo_comp))

# save output
write_csv(model_comparison_table, paste0(here::here(), "/results/PAranks_ModelEval_allScales.csv"))
model_comparison_table <- read_csv(paste0(here::here(), "/results/PAranks_ModelEval_allScales.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Random-slope model (PA ranks/ levels) ####
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

# # pointrange plot
# ggplot(pars_all %>%
#          full_join(fns_labels %>% dplyr::select(Function, Group_function, Label_short), by = c("fns" = "Function")) %>%
#          mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
#                                     ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
#                                            ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
#          mutate(Group_function = ifelse(Group_function=="Service", "Function", 
#                                         ifelse(Group_function=="Diversity", "Shannon", Group_function)),
#                 scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_flag-Portugal.png' width='70'>",
#                                                            "<img src='figures/icon_location-black.png' width='70'>",
#                                                            "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>"
#                                                            ))) %>%
#          mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
#          mutate(LC = factor(LC, levels = c("Dryland", "Woodland", "Shrubland", "Grassland", "Cropland","ns")),
#                 significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
#          mutate(significance = factor(significance, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
#          filter(!is.na(LC) & !is.na(scale)),
#        aes(x = PA_rank_rev.trend, y = Label_short,
#            xmin = lower.HPD, xmax = upper.HPD,
#            fill = significance, color = LC,
#            linetype = significance))+
#   
#   geom_vline(aes(xintercept=0), color="black")+
#   # geom_text(position = position_jitter(height = 0.3),
#   #           size = 2,
#   #           aes(label = fns),
#   #           hjust = 0.5) +
#   geom_pointrange(shape = 21,
#                   size = 0.5) +
# 
#   # geom_pointrange(data = pars_sum %>% filter(!is.na(LC) & !is.na(scale)) %>%
#   #                   mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
#   #                   mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='70'>",
#   #                                              ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='70'>",
#   #                                                     ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>", NA)))) %>%
#   #                   mutate(Group_function = ifelse(Group_function=="Service", "Function",
#   #                                                  ifelse(Group_function=="Diversity", "Shannon", Group_function)),
#   #                          scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_earth-globe-with-continents-maps.png' width='70'>",
#   #                                                                     "<img src='figures/icon_location-black.png' width='70'>",
#   #                                                                     "<img src='figures/icon_flag-Portugal.png' width='70'>" ))) %>%
#   #                   mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))) %>%
#   #                   mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","Dryland", "ns")),
#   #                          significance =  ifelse(sign(PA_rank_rev.trend_CI_lower)!= sign(PA_rank_rev.trend_CI_upper), "ns", as.character(as.factor(LC)))) %>%
#   #                   mutate(significance = factor(significance, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "ns"))),
#   # 
#   #                 aes(fill = significance, color = LC,
#   #                     y = Label_short, x = PA_rank_rev.trend_mean, 
#   #                     xmin = PA_rank_rev.trend_CI_lower, xmax = PA_rank_rev.trend_CI_upper),
#   # 
#   #                 position = position_dodge(width = 0.1),
#   #                 size = 1.5, shape = 21, linewidth = 1.5) +
# 
#   #coord_flip()+
#   coord_cartesian()+
#   ggh4x::facet_grid2(scale_icon + LC ~ Group_function, drop=TRUE, 
#              scales = "free", independent = "x", switch = "y",
#              space = "free",
#              shrink = FALSE)+
# 
#   #ylab("Effect size")+
#   scale_fill_manual(values=c("Cropland" = "#4A2040",
#                              "Grassland" = "#E69F00",
#                              "Shrubland" = "#0072B2", 
#                              "Woodland" = "#009E73", 
#                              "Other" = "#000000",
#                              "Dryland" = "#000000",
#                              "ns" = "white"),
#                     drop = FALSE)+
#   scale_linetype_manual(values=c("Cropland" = "solid",
#                              "Grassland" = "solid",
#                              "Shrubland" = "solid", 
#                              "Woodland" = "solid", 
#                              "Other" = "solid",
#                              "ns" = "longdash"),
#                     drop = FALSE)+
#   scale_color_manual(values=c("Cropland" = "#4A2040",
#                               "Grassland" = "#E69F00",
#                               "Shrubland" = "#0072B2", 
#                               "Woodland" = "#009E73", 
#                               "Other" = "#000000",
#                               "Dryland" = "#000000",
#                               "ns" = "grey"),
#                      drop = FALSE)+
#   theme_void()+
#   
#   guides(fill = guide_legend(reverse = F, override.aes = list(shape = 21, 
#                                                               color = c("Cropland" = "#4A2040",
#                                                                         "Grassland" = "#E69F00",
#                                                                         "Shrubland" = "#0072B2", 
#                                                                         "Woodland" = "#009E73",
#                                                                         "Dryland" = "#000000", 
#                                                                         "ns" = "grey"))), 
#          color = "none")+
#   theme(legend.position = "bottom", 
#         axis.title.y =element_blank(),
#         legend.text = element_text(size = 15),
#         legend.spacing.x = unit(0.5, "cm"),
#         legend.box.spacing = unit(c(1,0,0,0), "cm"),
#         legend.title = element_blank(),
#         plot.margin = unit(c(1, 2, 0.5, 0.5), "cm"),
#         axis.text.x = element_text(size = 15),
#         axis.ticks.x = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_line(color = "grey90"),
#         panel.spacing = unit(0.8, "cm"),
#         strip.background = element_rect(fill="white", color = "white"), #chocolate4
#         #strip.text.x.bottom = ggtext::element_markdown(vjust = 1),
#         strip.text.x = element_text(size = 30, hjust = 0, vjust = 1),
#         strip.text.y = ggtext::element_markdown(vjust = 0.5))
# 
# ggsave(filename=paste0(here::here(), "/figures/Results_slope_BayesianTrends_allScales_grouped.png"),
#        plot = last_plot(), 
#        width = 5000, height = 3000,
#        units = "px")

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
  mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","Dryland", "ns")),
         significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
  mutate(significance = factor(significance, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
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
      mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","Dryland", "ns")),
             significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
      mutate(significance = factor(significance, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "ns")))  %>% 
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

write_csv(pars_stats, paste0(here::here(), "/figures/Results_slope_BayesianTrends_stats_allScales.csv"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SUPPLEMENTARY FIGURE 3.2 -Heatmap Bayesian slopes ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## plot like heatmap
ggplot(data = pars_all %>%
         mutate(LC = factor(LC, levels = c("Woodland", "Shrubland", "Grassland", "Cropland","Dryland", "Other")),
                significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "not significant", as.character(as.factor(LC)))) %>%
         mutate(significance = factor(significance, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "Other", "not significant")))  %>%
         filter(!is.na(LC) & !is.na(scale)) %>%
         unique() %>%

         filter(LC %in% lc_names_all) %>%
         mutate(PA_trend_f = cut(abs(PA_rank_rev.trend),
                                    breaks=c(-Inf, -10, -1, 0, 1, 10, 100, Inf),
                                    labels=c("<-10", "<-1", "<0", "<1", "<10", "<100", ">100"))) %>%
         mutate(scale = factor(scale, levels = c("global","continental", "regional")),
                LC = factor(LC, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "Other", "Dryland"))) %>%
         full_join(expand.grid(scale = c("global", "continental", "regional"),
                               LC = c("Dryland", "Cropland", "Grassland", "Woodland"),
                               fns = fns_labels$Function)) %>%
         full_join(fns_labels %>% dplyr::select(Function, Group_function, Label_short, Label), by = c("fns" = "Function")) %>%
         mutate(effect_direction = as.factor(sign(PA_rank_rev.trend))) %>%
         mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
                                            ifelse(effect_direction=="1", "positive", "0"))) %>%
         mutate(Label = factor(Label, levels = labels_order)) %>%
         mutate(effect_significance = ifelse(significance == "not significant", "not significant", effect_direction_c)) %>%
         filter(LC != "Other" & !is.na(Label)) %>%
         mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "not significant"))) %>%
         mutate(effect_na = ifelse(is.na(PA_trend_f), "not available", NA)),
       aes(x = LC, y = scale))+

  geom_point(aes(size = PA_trend_f,
                 color = effect_direction_c,
                 fill= effect_significance,
                 shape = effect_na,
                 alpha = scale))+
  facet_wrap(vars(Label), ncol=6, drop=FALSE)+
  scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "not significant" = "white"),
                    name = "Direction of effect",
                    na.value = "black", drop = FALSE)+
  scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
                     name = "Direction of effect",
                     na.value = "black")+
  scale_size_manual(values =  c("<-10" = 10, "<-1" = 5, "<0" = 2, "<1" = 5, "<10" = 10, "<100" = 15, ">100" = 20), #c("marginal" = 2, "not significant" = 5, "small" = 5, "medium" = 10, "large" = 15),
                     name = "Slope estimate",
                     na.value = 1)+
  scale_shape_manual(values = c("not available" = 4),
                     name = "Missing data",
                     na.value = 21)+
  scale_alpha_manual(values = c("regional" = 1, "continental" = 1, "global" = 0.3))+
  scale_x_discrete(labels = c(
    "Dryland" = "<img src='figures/icon_land.png' width='33'>",
    "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
    "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
    "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
    "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
  ))+

  scale_y_discrete(labels = c(
    "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>",
    "continental" = "<img src='figures/icon_location-black.png' width='30'>",
    "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>"
  ))+
  xlab("")+ylab("")+
  theme_bw() + # use a white background

  guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black", "black"),
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
### Priors of Bayesian models ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# function to extract priors and merge in df
f_prior_merge <- function(list_priors, scale){
  list_priors <- lapply(list_priors, function(x){
    y <- x$intercept[, c("prior", "class", "coef", "group")]
    y$model <- "intercept"
    
    z <- x$slope[, c("prior", "class", "coef", "group")]
    z$model <- "slope"
    
    return(rbind(y, z))
  })
  list_priors <- lapply(names(list_priors), function(x){
    list_priors[[x]]$fns <- x
    list_priors[[x]]$scale <- scale
    list_priors[[x]] <- list_priors[[x]] %>%
      dplyr::select(model, scale, fns, class, prior)
    return(list_priors[[x]])
  })
  
  list_priors <- do.call(rbind, list_priors)
  
  return(list_priors)
}

# loop through spatial scales
df_prior <- data.frame()
for(temp_scale in c("continental", "regional", "global")){
  print(temp_scale)
  
  prior_list_noRF <- get(load(file=paste0(here::here(), "/intermediates/PAranks_Priors_", temp_scale, ".RData"))) #prior_list
  prior_list_noRF <- f_prior_merge(list_priors = prior_list_noRF, scale = temp_scale)
  prior_list_noRF$random_factors <- FALSE
  # remove flat priors
  prior_list_noRF <- prior_list_noRF %>% filter(class %in% c("Intercept", "sd", "sigma") & prior != "")
  
  if(temp_scale != "global"){
    prior_list_withRF <- get(load(file=paste0(here::here(), "/intermediates/PAranks_PriorsWithRandomFactor_", temp_scale, ".RData"))) #prior_list
    prior_list_withRF <- f_prior_merge(list_priors = prior_list_withRF, scale = temp_scale)
    prior_list_withRF$random_factors <- TRUE
    prior_list_withRF <- prior_list_withRF %>% filter(class %in% c("Intercept", "sd", "sigma") & prior != "")
  }
  
  df_prior <- rbind(df_prior, prior_list_noRF)
  df_prior <- rbind(df_prior, prior_list_withRF)
}
df_prior

df_prior %>% count(scale)

df_prior <- df_prior %>%
  group_by(scale, fns, prior, class) %>%
  mutate(
    rf = ifelse(random_factors, "withRF", "noRF"),
    group = paste0(model, " -", rf)) %>%
  summarise(model = paste(unique(group), collapse = ", ")) %>%
  arrange(scale, fns, class)

# check if sd = sigma everywhere (summary without model), if so, report together with note that sd only in withRF models

write_csv(prior_df, paste0(here::here(), "/results/PA_Bayesian_priors_all.csv"))
