library(tidyverse)
library(patchwork)
library(ggtext)

# load background map
world.inp <- map_data("world")
world.inp <- subset(world.inp, region %in% c("Albania", "Andorra", "Armenia", "Austria", "Azerbaijan",
                                               "Belarus", "Belgium", "Bosnia and Herzegovina", "Bulgaria",
                                               "Croatia", "Cyprus", "Czechia","Denmark","Estonia","Finland", 
                                               "France","Georgia", "Germany", "Greece","Hungary","Iceland", 
                                               "Ireland", "Italy","Kazakhstan", "Kosovo", "Latvia","Liechtenstein", 
                                               "Lithuania", "Luxembourg","Malta","Moldova","Monaco","Montenegro",
                                               "Macedonia", "Netherlands","Norway","Poland","Portugal","Romania",
                                               "Russia","San Marino","Serbia","Slovakia","Slovenia","Spain",
                                               "Sweden","Switzerland","Turkey","Ukraine","UK","Vatican"))


source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

data_locations_regi <- read_csv(file = paste0(here::here(), "/results/Locations_regional.csv"))
data_loc_bayes_regi <- read_csv(file = paste0(here::here(), "/results/Locations_PAranks_regional.csv"))

data_locations_cont <- read_csv(file = paste0(here::here(), "/results/Locations_continental.csv"))
data_loc_bayes_cont <- read_csv(file = paste0(here::here(), "/results/Locations_PAranks_continental.csv"))

# set limits for point maps
temp_limits_regi <- c(-9.02, -6.1, 40.66, 42.5)

# plot - FIGURE 1
p_regi <- ggplot()+
  geom_map(data = world.inp, map = world.inp,
           aes(map_id = region),  show.legend = FALSE,
           fill="white", color = "grey80", linewidth = 0.3) + #G:0.3, C: 0.6, R: 0.8

    geom_point(data=data_locations_regi, aes(x=Longitude, y=Latitude,
                                      shape = as.character(PA),
                                      color = LC,
                                      fill = as.character(PA),
                                      size = as.character(PA)),
             stroke = 2)+
  geom_point(data=data_loc_bayes_regi %>%
               filter(!(SampleID %in% unique(data_locations_regi$SampleID))),
             aes(x=Longitude, y=Latitude,
                 shape = as.character(PA),
                 color = LC,
                 fill = as.character(PA),
                 size = as.character(PA)),
             stroke = 2,  alpha = 0.15)+
  
  annotate(geom = "text",
           x = -8.85,
           y = 41,
           hjust = 0.5,
           color = "grey50",
           fontface = "italic",
           size = 5,
           label = "Atlantic\nOcean")+
  annotate(geom = "text",
           x = -7.85,
           y = 40.748,
           hjust = 0,
           color = "grey50",
           fontface = "italic",
           size = 5,
           label = "Portugal")+
  annotate(geom = "text",
           x = -6.5,
           y = 40.9,
           hjust = 0,
           color = "grey50",
           fontface = "italic",
           size = 5,
           label = "Spain")+
  
  annotate(geom = "rect",
           fill = "grey80",
           color = NA,
           xmin = temp_limits_regi[1],
           xmax = temp_limits_regi[2],
           ymin = temp_limits_regi[4]-(0.2 * (temp_limits_regi[4] - temp_limits_regi[3])),
           ymax = temp_limits_regi[4])+
  annotate(geom = "text",
           x = temp_limits_regi[1] + (0.22 * (temp_limits_regi[2] - temp_limits_regi[1])),
           y = temp_limits_regi[4]-((0.2 * (temp_limits_regi[4] - temp_limits_regi[3]))/2),
           hjust = 0,
           vjust = 0.5,
           color = "black",
           size = 10,
           label = "North of Portugal")+
  annotate(geom = "text",
           x = temp_limits_regi[2] - (0.1 * (temp_limits_regi[2] - temp_limits_regi[1])),
           y = temp_limits_regi[4]-((0.2 * (temp_limits_regi[4] - temp_limits_regi[3]))/2),
           hjust = 1,
           vjust = 0.5,
           color = "black",
           size = 5,
           label = paste0(table(data_locations_regi$PA)[2], "\n", table(data_locations_regi$PA)[1]))+
  annotate(geom = "point",
           x = temp_limits_regi[2] - (0.07 * (temp_limits_regi[2] - temp_limits_regi[1])),
           y = temp_limits_regi[4]-((0.153 * (temp_limits_regi[4] - temp_limits_regi[3]))/2),
           color = "black",
           fill = "black",
           size = 6,
           shape = 1,
           stroke = 2)+
  annotate(geom = "point",
           x = temp_limits_regi[2] - (0.07 * (temp_limits_regi[2] - temp_limits_regi[1])),
           y = temp_limits_regi[4]-((0.26 * (temp_limits_regi[4] - temp_limits_regi[3]))/2),
           color = "black",
           fill = "white",
           size = 2,
           shape = 21,
           stroke = 2)+
  ggtext::geom_richtext(aes(
    x = temp_limits_regi[1] + (0.05 * (temp_limits_regi[2] - temp_limits_regi[1])),
    y = temp_limits_regi[4]-((0.2 * (temp_limits_regi[4] - temp_limits_regi[3]))/2),
    hjust = 0,
    vjust = 0.5,
    label = "<img src='figures/icon_flag-Portugal.png' width='60'>"),
    fill = NA, label.color = NA)+ 
  
  #stroke = 1.4, color = "#000000")+ #increase circle line width; G: 2 (1.4), C+R:3
  scale_shape_manual(values = c("0" = 21, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" = 2, "1" = 6))+ #G: 1.4,4.5/0.3, 1, C+R:3,8/ 0.6,2
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000",
                                "Dryland" = "#000000"))+ 
  scale_fill_manual(values = c("0" = "white", "1" = "black"))+
  scale_x_continuous(expand = c(0,0), limits = c(temp_limits_regi[1], temp_limits_regi[2]))+
  scale_y_continuous(expand = c(0,0), limits = c(temp_limits_regi[3], temp_limits_regi[4]))+
  coord_map()+
  theme_bw()+
  theme(axis.title = element_blank(), 
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position ="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.background = element_rect(fill= "grey80"))
p_regi

# continental
temp_limits_cont <- c(-11, 34.8, 33, 71.38)
p_cont <- ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey80", linewidth = 0.3) + #G:0.3, C: 0.6, R: 0.8
  
  geom_point(data=data_locations_cont, aes(x=Longitude, y=Latitude, 
                                           shape = as.character(PA), 
                                           color = LC,
                                           fill = as.character(PA), 
                                           size = as.character(PA)),
             stroke = 1.5)+
  geom_point(data=data_loc_bayes_cont %>%
               filter(!(SampleID %in% unique(data_locations_cont$SampleID))), 
             aes(x=Longitude, y=Latitude,
                 shape = as.character(PA), 
                 color = LC,
                 fill = as.character(PA), 
                 size = as.character(PA)),
             stroke = 1.5,  alpha = 0.15)+
  
  annotate(geom = "text",
           x = -3,
           y = 63,
           hjust = 0.5,
           color = "grey50",
           fontface = "italic",
           size = 5,
           label = "Atlantic\nOcean")+
  annotate(geom = "text",
           x = 3,
           y = 56,
           hjust = 0.5,
           color = "grey50",
           fontface = "italic",
           size = 5, 
           label = "North\nSea")+
  annotate(geom = "text",
           x = 5,
           y = 35,
           hjust = 0, vjust = 1,
           color = "grey50",
           fontface = "italic",
           size = 5,
           label = "Mediterranean Sea")+
  
  annotate(geom = "rect",
           fill = "grey80",
           color = NA,
           xmin = temp_limits_cont[1],
           xmax = temp_limits_cont[2],
           ymin = temp_limits_cont[4]-(0.07 * (temp_limits_cont[4] - temp_limits_cont[3])),
           ymax = temp_limits_cont[4])+
  annotate(geom = "text",
           x = temp_limits_cont[1] + (0.2 * (temp_limits_cont[2] - temp_limits_cont[1])),
           y = temp_limits_cont[4]-((0.07 * (temp_limits_cont[4] - temp_limits_cont[3]))/2),
           hjust = 0,
           vjust = 0.5,
           color = "black",
           size = 10,
           label = "Europe")+
  annotate(geom = "text",
           x = temp_limits_cont[2] - (0.1 * (temp_limits_cont[2] - temp_limits_cont[1])),
           y = temp_limits_cont[4]-((0.07 * (temp_limits_cont[4] - temp_limits_cont[3]))/2),
           hjust = 1,
           vjust = 0.5,
           color = "black",
           size = 5,
           label = paste0("Protected sites:   ", table(data_locations_cont$PA)[2], "\n", 
                          "Unprotected sites: ", table(data_locations_cont$PA)[1]))+
  annotate(geom = "point",
           x = temp_limits_cont[2] - (0.07 * (temp_limits_cont[2] - temp_limits_cont[1])),
           y = temp_limits_cont[4]-((0.051 * (temp_limits_cont[4] - temp_limits_cont[3]))/2),
           color = "black",
           fill = "black",
           size = 6,
           shape = 1,
           stroke = 2)+
  annotate(geom = "point",
           x = temp_limits_cont[2] - (0.07 * (temp_limits_cont[2] - temp_limits_cont[1])),
           y = temp_limits_cont[4]-((0.091 * (temp_limits_cont[4] - temp_limits_cont[3]))/2),
           color = "black",
           fill = "white",
           size = 2,
           shape = 21,
           stroke = 2)+
  ggtext::geom_richtext(aes(
    x = temp_limits_cont[1] + (0.05 * (temp_limits_cont[2] - temp_limits_cont[1])),
    y = temp_limits_cont[4]-((0.07 * (temp_limits_cont[4] - temp_limits_cont[3]))/2),
    hjust = 0,
    vjust = 0.5,
    label = "<img src='figures/icon_location-black.png' width='60'>"),
    fill = NA, label.color = NA)+ 
  
  #stroke = 1.4, color = "#000000")+ #increase circle line width; G: 2 (1.4), C+R:3
  scale_shape_manual(values = c("0" = 21, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" = 1.7, "1" = 5))+ #G: 1.4,4.5/0.3, 1, C+R:3,8/ 0.6,2
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000",
                                "Dryland" = "#000000"))+ 
  scale_fill_manual(values = c("0" = "white", "1" = "black"))+
  scale_x_continuous(expand = c(0,0), limits = c(temp_limits_cont[1], temp_limits_cont[2]))+
  scale_y_continuous(expand = c(0,0), limits = c(temp_limits_cont[3], temp_limits_cont[4]))+
  coord_map()+
  theme_bw()+
  theme(axis.title = element_blank(), 
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position ="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black"),
        panel.background = element_rect(fill= "grey80"))
p_cont


## Bar plot 
data_locations_all <- data_locations_regi %>%
  mutate(scale = "regional") %>%
  rbind(data_locations_cont %>%
          mutate(scale = "continental"))

data_locations_all <- data_locations_all %>%
  mutate(LC_f = as.factor(LC)) %>%
  mutate(LC_f = factor(LC_f, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
  mutate(scale = factor(scale, levels = c("continental", "regional"))) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='30'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='30'>", NA))) %>%
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_location-black.png' width='30'>",
                                                    "<img src='figures/icon_flag-Portugal.png' width='30'>" )))

data_loc_bayes_all <- data_loc_bayes_regi %>%
  dplyr::select(SampleID, LC, Latitude, Longitude, PA) %>%
  mutate(scale = "regional") %>%
  rbind(data_loc_bayes_cont %>%
          dplyr::select(SampleID, LC, Latitude, Longitude, PA) %>%
          mutate(scale = "continental"))

data_loc_bayes_all <- data_loc_bayes_all %>%
  mutate(LC_f = as.factor(LC)) %>%
  mutate(LC_f = factor(LC_f, levels = c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
  mutate(scale = factor(scale, levels = c("continental", "regional"))) %>%
  mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='30'>",
                             ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='30'>", NA))) %>%
  mutate(scale_icon = factor(scale_icon, levels = c("<img src='figures/icon_location-black.png' width='30'>",
                                                    "<img src='figures/icon_flag-Portugal.png' width='30'>" )))

# plot
p_bar <- ggplot()+
  xlab("")+
  ylab("Number of sampling sites")+
  # geom_bar(data = data_loc_bayes_all, 
  #          aes(x = LC_f, color = LC_f), fill = "white", #fill = LC_f, 
  #          position = "dodge",
  #          alpha = 0.2,
  #          linewidth = 0.6,
  #          just = 0,
  #          linetype = "longdash",
  #          width = 0.5,
  #          na.rm=TRUE)+
  # geom_bar(data = data_loc_bayes_all %>%
  #            filter(PA == 0) %>%
  #            mutate(scale = factor(scale, levels = c("regional", "continental"))), 
  #          aes(x = LC_f, fill = LC_f, color = LC_f), #fill = "white", 
  #          position = "dodge",
  #          linewidth = 0.6,
  #          linetype = "longdash",
  #          just = 0,
  #          alpha = 0.2,
  #          width = 0.4,
  #          na.rm=TRUE)+
  geom_bar(data = data_locations_all, 
         aes(x = LC_f, color = LC_f),fill = "white", 
         position = "dodge",
         linewidth = 1,
         width = 0.8,
         na.rm=TRUE)+
  geom_bar(data = data_locations_all, 
           aes(x = LC_f, fill = LC_f, color = LC_f,
               alpha = "Protected"), #fill = "white", 
           position = "dodge",
           linewidth = 1,
           width = 0.8,
           na.rm=TRUE)+
  geom_bar(data = data_locations_all %>%
             filter(PA == 0) %>%
             mutate(scale = factor(scale, levels = c("regional", "continental"))), 
           aes(x = LC_f, color = LC_f, fill = LC_f,
               alpha = "Unprotected"), #fill = "white", 
           position = "dodge",
           linewidth = 1,
           width = 0.8,
           na.rm=TRUE)+
  
  scale_fill_manual(values = c("Cropland" = "#4A2040",
                               "Grassland" = "#E69F00",
                               "Shrubland" = "#0072B2", 
                               "Woodland" = "#009E73", 
                               "Other" = "#000000",
                               "Dryland" = "#000000"))+
  scale_color_manual(values = c("Cropland" = "#4A2040",
                               "Grassland" = "#E69F00",
                               "Shrubland" = "#0072B2", 
                               "Woodland" = "#009E73", 
                               "Other" = "#000000",
                               "Dryland" = "#000000"))+
  
  geom_text(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_location-black.png' width='30'>"),
           x = 1, y = 40,
           hjust = 0.5,
           label = c("Cropland"),
           angle = 90,
           color = "white",
           size = 6,
           inherit.aes = FALSE)+
  geom_text(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_location-black.png' width='30'>"),
            x = 2, y = 40,
            hjust = 0.5,
            label = c("Grassland"),
            angle = 90,
            color = "white",
            size = 6,
            inherit.aes = FALSE)+
  geom_text(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_location-black.png' width='30'>"),
            x = 3, y = 40,
            hjust = 0.5,
            label = c("Woodland"),
            angle = 90,
            color = "white",
            size = 6,
            inherit.aes = FALSE)+
  
  geom_point(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_location-black.png' width='30'>"),
             x = 3.35, y = 96.5,
             shape = 21,
             stroke = 2,
             color = "black",
             fill = "white",
             size = 3,
             inherit.aes = FALSE)+
  geom_point(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_location-black.png' width='30'>"),
             x = 3.35, y = 114.5,
             shape = 1,
             stroke = 2,
             color = "black",
             size = 6,
             inherit.aes = FALSE)+
  geom_point(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_flag-Portugal.png' width='30'>"),
             x = 3.35, y = 49,
             shape = 21,
             stroke = 2,
             color = "black",
             fill = "white",
             size = 3,
             inherit.aes = FALSE)+
  geom_point(data = subset(data_locations_all, scale_icon == "<img src='figures/icon_flag-Portugal.png' width='30'>"),
             x = 3.35, y = 61.5,
             shape = 1,
             stroke = 2,
             color = "black",
             size = 6,
             inherit.aes = FALSE)+
  
  scale_alpha_manual(values = c("Protected" = 0.1, "Unprotected" = 1))+
  guides(alpha = guide_legend(override.aes = list(fill = "black",
                                                  alpha = c(0.1, 1))))+
  scale_x_discrete(expand = c(0,0.7),
                   drop = TRUE,
                   labels = c(
                         #"Dryland" = "<img src='figures/icon_land.png' width='30'>",
                         "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
                         "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
                         #"Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
                         "Woodland" = "<img src='figures/icon_forest.png' width='30'>"
                       ))+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 333))+
  facet_grid(cols = vars(scale_icon), drop = TRUE, 
             scales = "free_y", 
             space = "free_y")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_line(color = "grey80"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey60"),
        strip.background = element_blank(),
        strip.text = ggtext::element_markdown(vjust = 0, hjust = 0.15),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(1, "cm"),
        axis.text.y = element_text(size = 15, vjust = 0.5),
        axis.title.y = element_text(size = 15),
        axis.text.x = ggtext::element_markdown(vjust = 0))
p_bar

# merge
p_design <- "
AAAABBBB
AAAABBBB
AAAABBBB
AAAACCCC
AAAACCCC"

p_fig1 <- (p_cont + p_regi + free(p_bar) & theme(legend.position = "none")) +
  plot_layout(design = p_design)+
  plot_annotation(tag_levels = "A") & 
  theme(
    plot.tag = element_text(face = "bold", size = 10)  # Bold and bigger tags
  )
ggsave(filename = paste0(here::here(), "/figures/FIGURE_1.png"),
       plot = p_fig1,
       dpi = 300,
       width=13, height=10)
