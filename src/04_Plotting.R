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

#temp_scale <- "global"
#temp_scale <- "continental"
temp_scale <- "regional"

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
if(temp_scale == "global") temp_date <- "2023-12-01"
if(temp_scale == "continental") temp_date <- "2023-12-14"
if(temp_scale == "regional") temp_date <- "2023-12-14"

if(temp_scale == "global") lc_names <- lc_names[lc_names != "Other"]
if(temp_scale == "continental") lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland"]
  min_size <- 7
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

#pars_sample <- pars_sample %>% group_by(fns, lc) %>% slice_sample(n=10000)

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
data_locations #G: nrow=235, C: 243, R: 152
nrow(data_locations %>% filter(PA==1)) #G: 56 PAs, C: 59, R: 38
nrow(data_locations %>% filter(PA==0)) #G: 179 PAs, C: 184, R: 124

# set limits for point maps
if(temp_scale == "global") temp_limits <- c(-180, 180, -180, 180)
if(temp_scale == "continental") temp_limits <- c(-10, 35, 35, 70)
if(temp_scale == "regional") temp_limits <- c(-9, -6, 40.5, 42.5)

ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey90", linewidth = 0.15) + #fill = "grey80", color="grey75"
  #xlim(-180, 180) +  ylim(-180, 180) + #global
  #xlim(-10, 35) +  ylim(35, 70) + #continental
  #xlim(-9.5, -6) +  ylim(40.5, 42.5) + #regional
  xlim(temp_limits[1], temp_limits[2])+
  ylim(temp_limits[3], temp_limits[4])+
  
  geom_point(data=data_locations, aes(x=Longitude, y=Latitude, 
                                      shape = as.character(PA), color=LC, 
                                      size = as.character(PA)),
             stroke = 3)+ #increase circle line width
  scale_shape_manual(values = c(19,1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c(3,8))+
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
## Pairing ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Save some numbers
sink(paste0(here::here(), "/results/Numbers_", temp_scale, ".txt"))

print(temp_scale)
# look what non-protected sites have (not) been paired to any PA
cat("Number of non-protected sites have never been paired to any PA:")
hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
length(setdiff(data_clean[data_clean$PA==0,]$SampleID, pa_pairs$nonPA))  
# nonPA sites never used: G 131

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

# cat("Mean number of nonPA sites per PA.")
# mean(table(pa_pairs$nonPA))
# 
# cat("SD number of nonPA sites per PA.")
# sd(table(pa_pairs$nonPA))
# 
# cat("Min/max number of nonPA sites per PA.")
# min(table(pa_pairs$nonPA))
# max(table(pa_pairs$nonPA))
# 
# cat("Mean/SD/min/max number of times individual PA were used.")
# mean(table(pa_pairs$ID))
# sd(table(pa_pairs$ID))
# min(table(pa_pairs$ID))
# max(table(pa_pairs$ID))

table(pa_pairs$ID)

cat("Number of sites (PA/ nonPA) per LC.")
pa_pairs %>% group_by(LC) %>% dplyr::select(ID) %>% unique() %>% count()
pa_pairs %>% group_by(LC) %>% dplyr::select(nonPA) %>% unique() %>% count()

sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplot mahalanobis distance ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

ggplot(pa_pairs, aes(x = LC, y = mahal.min, fill = LC))+
  geom_boxplot()+
  geom_violin(alpha = 0.3, adjust = 0.3)+
  theme_bw() +
  labs(x="Land-cover type",y="Mahalanobis distance") +
  theme(axis.text.x=element_text(size=15),text = element_text(size=20),  
        legend.position = "none", axis.text.y = element_text(size=15), legend.title = element_blank())+
  scale_fill_manual(values=c("Cropland" = "#4A2040",
                             "Grassland" = "#E69F00",
                             "Shrubland" = "#0072B2", 
                             "Woodland" = "#009E73", 
                             "Other" = "#000000")) #"gold3", "limegreen", "forestgreen"
ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_global.png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Effect size per LC type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
d_df <- do.call(rbind, d_list)
str(d_df)

d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## Violin plot effect size per LC type ####
# 
# ggplot(data=d_df, aes(x=round(effect,2), y=Label))+
#   geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
#   geom_vline(aes(xintercept=-0.8, linetype = "-0.8"), color="grey60")+ #large effect
#   geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
#   geom_vline(aes(xintercept=-0.5, linetype = "-0.5"), color="grey60")+ #medium effect
#   geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+ #small effect
#   geom_vline(aes(xintercept=-0.2, linetype = "-0.2"), color="grey60")+
#   #geom_vline(aes(xintercept=0, linetype = "0"))+
#   scale_linetype_manual(values = c("-0.8" = "dotted",
#                                    "-0.5" = "dashed", "-0.2" = "solid",
#                                    "0.2" = "solid","0.5" = "dashed",
#                                    "0.8" = "dotted"))+
#   
#   geom_violin(aes(fill = Organism), scale="width", width=0.5)+ #scale by width of violin bounding box
#   stat_summary(fun = "mean", geom = "pointrange", color = "black")+
#   facet_wrap(vars(lc))+
#   scale_x_continuous(labels = function(x) format(x, nsmall = 0))+#round c-value axis label
#   scale_fill_viridis_d()+
#   xlab("Effect size")+
#   theme_bw() + # use a white background
#   theme(legend.position = "right", axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),  
#         axis.text.x = element_text(size=7),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.spacing = unit(1, "lines"))
# ggsave(filename=paste0(here::here(), "/figures/Data_violin_d-value_global.png"),
#        plot = last_plot())

# save data for plot
write.csv(d_df, file=paste0(here::here(), "/figures/Data_pointrange_d-value_", temp_scale, ".csv"), row.names = FALSE)

d_df <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_", temp_scale, ".csv"))

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
write.csv(d_summary, file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_", temp_scale, ".csv"))

# mean per lc type
d_summary <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_", temp_scale, ".csv"))

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
  #n=nrow(d_summary %>% filter(abs(effect_mean) >= 0.2))
  ) #for print() command
sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pointrange effect size per LC type ####
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
ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_", temp_scale, ".png"),
       plot = last_plot(),
       width=5, height=4)

## different visualisation in grid
# ggplot(data = d_summary %>% 
#          filter(lc %in% lc_names), 
#        aes(y = effect_ci_17, x = effect_ci_83, 
#            #color=as.factor(sign(effect_median))
#            ))+
# 
#   # geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
#   # geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
#   # geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+
#   # scale_linetype_manual(values = c("0.8" = "dotted",
#   #                                  "0.5" = "dashed",
#   #                                  "0.2" = "solid"),
#   #                       name="Strength of effect")+
#   #annotate("rect", ymin = 0, ymax = 0.2, xmin=0, xmax=0.2, fill = "grey", alpha=0.1)+ #chocolate4
#   annotate("rect", ymin = -0.2, ymax = 0.2, xmin=-0.2, xmax=0.2, fill = "grey", alpha=0.3)+
#   annotate("rect", ymin = -0.5, ymax = 0.5, xmin=-0.5, xmax=0.5, fill = "grey", alpha=0.3)+
#   annotate("rect", ymin = -0.8, ymax = 0.8, xmin=-0.8, xmax=0.8, fill = "grey", alpha=0.3)+
#   
#   geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "grey")+
#   
#   # geom_pointrange(aes(xmin = abs(effect_ci_2.5), xmax = abs(effect_ci_97.5),
#   #                     shape = Group_function))+
#   # geom_linerange(aes(xmin = ifelse(sign(effect_ci_17) != sign(effect_ci_83), 
#   #                                  abs(effect_ci_17)), xmax = abs(effect_ci_83),
#   #                     #shape = Group_function
#   #                    ), linewidth=1, alpha = 0.2)+
#   geom_point(aes(color = Group_function), size = 5)+
#   ggrepel::geom_text_repel(aes(label = Label, color = Group_function), size = 3)+
#   #geom_text(aes(label = Label), size = 1.5, nudge_y = 0.005)+
#   
#   xlab("Upper CI (83%)")+ ylab("Lower CI (17%)")+
#   scale_y_continuous(breaks = c())+
#   scale_x_continuous(breaks = c())+
#   # scale_y_continuous(limits = c(0, 0.7),
#   #                    expand = c(0,0), breaks = c(0, 0.25, 0.5))+
#   #scale_color_manual(c(""))+
#   
#   facet_wrap(vars(lc), scales = "free")+
#   theme_bw() + # use a white background
#   theme(legend.position = "bottom", 
#         legend.direction = "vertical",
#         #axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),  
#         axis.text.x = element_text(size=10),
#         #panel.grid.major.y = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"), #chocolate4
#         strip.text = element_text(color="black")) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_medianSD_", temp_scale, ".png"),
#        plot = last_plot())

## one-sided
ggplot(data = d_summary %>%
         filter(lc %in% lc_names) %>%
         mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
                                         abs(effect_ci_17), 
                                         abs(effect_ci_83))),
       aes(y = abs(effect_ci_min66), x = abs(effect_median),
           color=as.factor(sign(effect_median))
           ))+

  # geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
  # geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
  # geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+
  # scale_linetype_manual(values = c("0.8" = "dotted",
  #                                  "0.5" = "dashed",
  #                                  "0.2" = "solid"),
  #                       name="Strength of effect")+
  #annotate("rect", ymin = 0, ymax = 0.2, xmin=0, xmax=0.2, fill = "grey", alpha=0.1)+ #chocolate4
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.2, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.5, fill = "grey", alpha=0.3)+
  annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.8, fill = "grey", alpha=0.3)+

  # geom_pointrange(aes(xmin = abs(effect_ci_2.5), xmax = abs(effect_ci_97.5),
  #                     shape = Group_function))+
  # geom_linerange(aes(xmin = ifelse(sign(effect_ci_17) != sign(effect_ci_83),
  #                                  abs(effect_ci_17)), xmax = abs(effect_ci_83),
  #                     #shape = Group_function
  #                    ), linewidth=1, alpha = 0.2)+
  geom_point(aes(shape = Group_function), size = 5)+
  ggrepel::geom_text_repel(aes(label = Label, shape = Group_function), size = 3)+
  #geom_text(aes(label = Label), size = 1.5, nudge_y = 0.005)+

  xlab("Median effect")+ ylab("Lower CI (17%)")+
  scale_y_continuous(breaks = c(0.2, 0.5, 0.8))+
  scale_x_continuous(breaks = c(0.2, 0.5, 0.8))+
  # scale_y_continuous(limits = c(0, 0.7),
  #                    expand = c(0,0), breaks = c(0, 0.25, 0.5))+
  #scale_color_manual(c(""))+

  facet_wrap(vars(lc), scales = "free")+
  theme_bw() + # use a white background
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        #axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        #panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white"), #chocolate4
        strip.text = element_text(color="black")) #white
ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_medianSD_", temp_scale, ".png"),
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
ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_medianSD_", temp_scale, ".png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Heatmap all 3 scales ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FIGURE 2 - Heatmap ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mean per lc type and all 3 scales
d_sum_glob <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_global.csv"))
d_sum_cont <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_continental.csv"))
d_sum_regi <- read.csv(file=paste0(here::here(), "/figures/Results_pointrange_d-value_summary_regional.csv"))

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
d_sum_all %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_median), 
                   list(mean=function(x) mean(abs(x)),
                        sd=function(x) sd(abs(x)))))

# heatmap
# ggplot(data = d_sum_all %>%
#          filter(lc %in% lc_names_all) %>%
#          mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
#                                          abs(effect_ci_17), 
#                                          abs(effect_ci_83))) %>%
#          mutate(effect_ci_min66f = cut(effect_ci_min66,
#                                        breaks=c(0, 0.2, 0.5, 0.8, Inf),
#                                        labels=c("ns", "small", "medium", "large"))) %>%
#          filter(!is.na(effect_ci_min66)) %>%
#          full_join(expand.grid(scale = c("global", "continental", "regional"), 
#                                lc = lc_names_all, 
#                                Label = fns_labels$Label)) %>%
#          mutate(scale = factor(scale, levels = c("global", "continental", "regional")),
#                 Label = factor(Label, levels = rev(fns_labels %>% arrange(Group_function, Label) %>% pull(Label))),
#                 lc = factor(lc, levels = c("Cropland", "Grassland", "Woodland", "Shrubland"))) %>%
#          filter(lc != "Other" & !is.na(Label)),
#        
#        aes(x = scale, y = Label, alpha=effect_ci_min66f, 
#            fill=as.factor(sign(effect_median))))+
#   
#   geom_tile()+
#   facet_grid(cols=vars(lc), drop=FALSE)+
#   scale_fill_manual(values = c("-1" = "#fc8d59", "0" = "#ffffbf", "1" = "#91bfdb"),
#                     name = "Direction of effect",
#                     na.value = "grey80")+
#   scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
#                      name = "Minimum effect size (66% CI)")+
#   theme_bw() + # use a white background
#   theme(legend.position = "right",
#         legend.direction = "vertical",
#         #axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"), #chocolate4
#         strip.text = element_text(color="black")) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_minCI66_allScales.png"),
#        plot = last_plot())


## switch lc and scales
d_plot_all <- d_sum_all %>%
  filter(lc %in% lc_names_all) %>%
  dplyr::select(-Label) %>%
  right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
             by = c("fns" = "Function")) %>%
  # mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
  #                                 abs(effect_ci_17), 
  #                                 abs(effect_ci_83))) %>%
  # mutate(effect_ci_min66 = ifelse(sign(effect_ci_17)!= sign(effect_ci_83), 0.01, effect_ci_min66)) %>%
  # mutate(effect_ci_min66f = cut(effect_ci_min66,
  #                               breaks=c(0, 0.2, 0.5, 0.8, Inf),
  #                               labels=c("marginal", "small", "medium", "large"))) %>%
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
  #mutate(effect_direction_c = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c)) %>%
  mutate(Label = factor(Label, levels = labels_order)) %>%
  mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c),
         effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
  mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "ns")))

ggplot(data = d_plot_all,
       aes(x = lc, y = scale))+
  
  # geom_tile(aes(alpha=effect_ci_min66f, 
  #                fill=as.factor(sign(effect_median))))+ 
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
  # scale_color_manual(values = c("negative" = "#fc8d59", "ns" = "black", "positive" = "#91bfdb"),
  #                   name = "Direction of effect",
  #                   na.value = "grey60")+
  scale_size_manual(values = c("marginal" = 2, "ns" = 5, "small" = 5, "medium" = 10, "large" = 15),
                     name = "Effect size",
                    na.value = 5)+
  # scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
  #                    name = "Effect size")+
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
  theme(#legend.position = c(0.96, -0.01),
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
ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_meanCI_allScales_fns.png"),
       plot = last_plot(), 
       width = 4400, height = 3800, units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summarizing stats ####
# number significant effects
d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% nrow()
# number ns
d_plot_all %>% filter(effect_significance=="ns" & !is.na(effect_significance)) %>% nrow()

# number significant per lc
table(d_plot_all %>% filter(effect_significance!="ns" & !is.na(effect_significance)) %>% dplyr::select(scale, lc))
table(d_plot_all %>% filter(!is.na(effect_significance)) %>% dplyr::select(scale, lc))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pointrange plot grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

d_df_glob <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_global.csv"))
d_df_cont <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_continental.csv"))
d_df_regi <- read_csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_regional.csv"))

d_df_all <- rbind(d_df_glob %>% mutate("scale" = "global"), 
                   d_df_cont %>% mutate("scale" = "continental")) %>%
  rbind(d_df_regi %>% mutate("scale" = "regional"))  %>%
  as_tibble()
rm(d_df_glob, d_df_cont, d_df_regi)


# d_plot_group <- d_df_all %>%
#   group_by(lc, scale, Organism) %>%
#   summarize(across(effect, .fns=list("mean" = function(x) mean(x, na.rm=T),
#                                      "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
#                                      "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
#                                      "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
#                                      "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE))))

ggplot(data = d_df_all %>% filter(!is.na(lc)) %>%
         mutate(scale = factor(scale, levels = c("regional", "continental", "global"))) %>%
         mutate(scale_icon = ifelse(scale == "regional", "<img src='figures/icon_flag-Portugal.png' width='30'>",
                                    ifelse(scale == "continental", "<img src='figures/icon_location-black.png' width='30'>",
                                           ifelse(scale == "global", "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>", NA)))) %>%
         mutate(Group_function = ifelse(Group_function=="Service", "Function", 
                                        ifelse(Group_function=="Diversity", "Shannon", Group_function))) %>%
        mutate(Group_function = factor(Group_function, levels = c("Function", "Richness", "Shannon", "Dissimilarity"))),
       aes(fill = lc, color = lc, 
           y = scale_icon, x = effect))+
  
  geom_vline(aes(xintercept=0), color="black")+
  ggdist::stat_pointinterval(fatten_point=1.2, shape=21,
                             position=position_dodgejust(width=0.5)) +
  coord_flip()+
  facet_wrap(vars(Group_function), drop=FALSE)+
  #ylab("Effect size")+
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
  # scale_x_discrete(labels = c(
  #   "global" = "<img src='figures/icon_earth-globe-with-continents-maps.png' width='30'>",
  #   "continental" = "<img src='figures/icon_location-black.png' width='30'>",
  #   "regional" = "<img src='figures/icon_flag-Portugal.png' width='30'>"
  # ))+
  scale_x_continuous(breaks = c(-2, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 2))+
  theme_void()+
  theme(legend.position = "right", axis.title.y =element_blank(),
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        plot.margin = unit(c(0, 0, 2, 0.5), "cm"),
        axis.text.y = element_text(size = 20, hjust = 1),  
        axis.text.x = ggtext::element_markdown(vjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        strip.background = element_rect(fill="white", color = "white"), #chocolate4
        strip.text = element_text(size = 40, hjust = 0))

ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_medianCI_allScales_grouped.png"),
       plot = last_plot(), 
       #width = 2000, height = 1500,
       units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplots values ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
  #scale_y_continuous(limits=c(-3,13)) + 
  
  xlab("")+ ylab("")+
  facet_wrap(vars(LC), ncol=1, nrow=4)+
  theme_bw() +
  theme(axis.text.x=element_text(size=5, angle=45, hjust=1),
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", legend.text = element_text(size=15), axis.text.y = element_text(size=15))
ggsave(filename=paste0(here::here(), "/figures/Results_boxplot_estimates_", temp_scale, ".png"),
       plot = last_plot())


# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## Combine all 3 scales into one ####
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# d_df_glob <- read.csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_global.csv"))
# d_df_cont <- read.csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_continental.csv"))
# d_df_regi <- read.csv(file=paste0(here::here(), "/figures/Data_pointrange_d-value_regional.csv"))
# 
# # heatmap
# 
# ggplot(data = d_df %>% 
#          filter(lc %in% lc_names) %>%
#          mutate(Label=factor(Label, levels = rev(fns_labels$Label))), 
#        aes(y = effect, x = Label))+
#   
#   geom_hline(aes(yintercept=0.8, linetype = "-0.8 / 0.8"), color="grey60")+
#   geom_hline(aes(yintercept=-0.8, linetype = "-0.8 / 0.8"), color="grey60")+ #large effect
#   geom_hline(aes(yintercept=0.5, linetype = "-0.5 / 0.5"), color="grey60")+
#   geom_hline(aes(yintercept=-0.5, linetype = "-0.5 / 0.5"), color="grey60")+ #medium effect
#   geom_hline(aes(yintercept=0.2, linetype = "-0.2 / 0.2"), color="grey90")+ #small effect
#   geom_hline(aes(yintercept=-0.2, linetype = "-0.2 / 0.2"), color="grey90")+
#   #geom_hline(aes(yintercept=0, linetype = "0"))+
#   scale_linetype_manual(values = c("-0.8 / 0.8" = "dotted",
#                                    "-0.5 / 0.5" = "dashed", 
#                                    "-0.2 / 0.2" = "solid"), name="Strength of effect")+
#   annotate("rect", ymin = -0.2, ymax = 0.2, xmin=-Inf, xmax=Inf, fill = "grey95")+
#   
#   ggdist::stat_pointinterval(fatten_point=1.2, shape=21) +
#   coord_flip()+
#   ylab("Effect size")+
#   facet_wrap(vars(lc))+
#   theme_bw() + # use a white background
#   theme(legend.position = "bottom", axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),  
#         axis.text.x = element_text(size=10),
#         panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(),
#         strip.text = element_text(color="white"))
# ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_d-value_allScales.png"),
#        plot = last_plot(),
#        width=5, height=4)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Bayesian results (PA types) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX S3 - Bayesian PA types ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
write_csv(pars_summary, file=paste0(here::here(), "/figures/Results_pointrange_parsBayesian_summary_", temp_scale, ".csv"))

# extract sample size
n_table <- data_clean %>% filter(LC!="Other") %>%
  group_by(LC, PA, PA_type, PA_rank) %>% count() %>%
  pivot_wider(id_cols = c("PA", "PA_type", "PA_rank"), 
              names_from = LC, values_from = n) %>% 
  arrange(PA_rank) %>% ungroup() %>%
  dplyr::select(-PA_rank)
n_table
write_csv(n_table, file=paste0(here::here(), "/figures/Results_pointrange_parsBayesian_nTable_", temp_scale, ".csv"))

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
                     #filter(Label_fns != "Water regulation") %>%%>%
                     filter(!is.na(Label_pa)) %>%
                     mutate(Label_fns = factor(Label_fns, levels = c(labels_order, "Number of sites"))),
       aes(y=Label_pa, x=value, color=lc))+
  
  ## adapt for scale
  annotate("rect", ymin = -Inf, ymax = 3+0.5, xmin=-Inf, xmax=Inf, fill = "grey90", alpha=0.5)+ #global: 4, C: 3, R: 2
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
ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_parsBayesian_", temp_scale, ".png"),
       plot = last_plot(),
       width=12, height=10)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Bayesian results (PA ranks) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

pred_list <- rbind(get(load(paste0(here::here(), "/results/PAranks_Bayesian_global_sample10k.RData"))), 
                   get(load(paste0(here::here(), "/results/PAranks_Bayesian_continental_sample10k.RData")))) %>% 
  rbind(get(load(paste0(here::here(), "/results/PAranks_Bayesian_regional_sample10k.RData"))))

for(temp_scale in c("global", "continental", "regional")){
  ggplot(data = pred_list %>% filter(scale == temp_scale) %>%
           right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
                      by = c("fns" = "Function")) %>%
           mutate(Label = factor(Label, levels = labels_order)), 
         
         aes(x = PA_rank, y = .epred, color = ordered(LC))) +
    stat_lineribbon() +
    #geom_point(data = pred_list) +
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
  ggsave(filename=paste0(here::here(), "/figures/Results_regressions_parsBayesian_", temp_scale,".png"),
         plot = last_plot(),
         width=15, height=10)
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Bayesian pointrange grouped per estimate type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
pars_glob <- read_csv(file=paste0(here::here(), "/results/PAranks_Bayesian_global_emtrends.csv"))
pars_cont <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_continental_emtrends.csv"))
pars_regi <- read_csv(paste0(here::here(), "/results/PAranks_Bayesian_regional_emtrends.csv"))

pars_all <- rbind(pars_glob, pars_cont) %>%
  rbind(pars_regi)  %>%
  as_tibble()
rm(pars_glob, pars_cont, pars_regi)
# 
# d_plot_all <- d_sum_all %>%
#   filter(lc %in% lc_names_all) %>%
#   dplyr::select(-Label) %>%
#   right_join(fns_labels %>% dplyr::select(Label, Label_short, Function), 
#              by = c("fns" = "Function")) %>%
#   # mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
#   #                                 abs(effect_ci_17), 
#   #                                 abs(effect_ci_83))) %>%
#   # mutate(effect_ci_min66 = ifelse(sign(effect_ci_17)!= sign(effect_ci_83), 0.01, effect_ci_min66)) %>%
#   # mutate(effect_ci_min66f = cut(effect_ci_min66,
#   #                               breaks=c(0, 0.2, 0.5, 0.8, Inf),
#   #                               labels=c("marginal", "small", "medium", "large"))) %>%
#   mutate(effect_mean_f = cut(abs(effect_mean),
#                              breaks=c(0, 0.2, 0.5, 0.8, Inf),
#                              labels=c("marginal", "small", "medium", "large"))) %>%
#   filter(!is.na(effect_mean)) %>% #!is.na(effect_ci_min66) & 
#   full_join(expand.grid(scale = c("global", "continental", "regional"), 
#                         lc = lc_names_all, 
#                         Label = fns_labels$Label)) %>%
#   mutate(scale = factor(scale, levels = rev(c("global", "continental", "regional"))),
#          Label = factor(Label, levels = rev(fns_labels %>% arrange(Group_function, Label) %>% pull(Label))),
#          lc = factor(lc, levels = c("Cropland", "Grassland", "Shrubland", "Woodland"))) %>%
#   filter(lc != "Other" & !is.na(Label)) %>% 
#   mutate(effect_direction = as.factor(sign(effect_mean)))%>%
#   mutate(effect_direction_c = ifelse(effect_direction=="-1", "negative",
#                                      ifelse(effect_direction=="1", "positive", "0"))) %>%
#   #mutate(effect_direction_c = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c)) %>%
#   mutate(Label = factor(Label, levels = labels_order)) %>%
#   mutate(effect_significance = ifelse(sign(effect_ci_2.5)!= sign(effect_ci_97.5), "ns", effect_direction_c),
#          effect_na = ifelse(is.na(effect_mean), "not available", NA)) %>%
#   mutate(effect_significance = factor(effect_significance, levels = c("negative", "positive", "ns")))
# 
# # plot
# ggplot(data = pars_sum,
#        aes(x = LC, y = scale))+
#   
#   # geom_tile(aes(alpha=effect_ci_min66f, 
#   #                fill=as.factor(sign(effect_median))))+ 
#   geom_point(aes(size = PA_rank.trend,
#                  color = effect_direction_c, 
#                  fill= effect_significance,
#                  shape = effect_na))+
#   facet_wrap(vars(Label), ncol=6, drop=FALSE)+
#   scale_fill_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb", "ns" = "white"),
#                     name = "Direction of effect",
#                     na.value = "black", drop = FALSE)+
#   scale_color_manual(values = c("negative" = "#fc8d59", "positive" = "#91bfdb"),
#                      name = "Direction of effect",
#                      na.value = "black")+
#   # scale_color_manual(values = c("negative" = "#fc8d59", "ns" = "black", "positive" = "#91bfdb"),
#   #                   name = "Direction of effect",
#   #                   na.value = "grey60")+
#   scale_size_manual(values = c("marginal" = 2, "ns" = 5, "small" = 5, "medium" = 10, "large" = 15),
#                     name = "Effect size",
#                     na.value = 5)+
#   # scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
#   #                    name = "Effect size")+
#   scale_shape_manual(values = c("not available" = 4),
#                      name = "Missing data",
#                      na.value = 21)+
#   scale_x_discrete(labels = c(
#     "Cropland" = "<img src='figures/icon_harvest.png' width='20'>",
#     "Grassland" = "<img src='figures/icon_grass.png' width='17'>",
#     "Shrubland" = "<img src='figures/icon_shrub-crop.png' width='35'>",
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
#   guides(fill = guide_legend(override.aes = list(color = c("#fc8d59", "#91bfdb", "black"), 
#                                                  shape = 21, size = 5)), #tell legend to use different point shape
#          color = "none", #don't show legend
#          shape = guide_legend(override.aes = list(size = 5)))+
#   theme(#legend.position = c(0.96, -0.01),
#     legend.position = "bottom", #c(0.96, -0.05),
#     legend.justification = c(1, 0),
#     legend.box = "horizontal",
#     legend.direction = "vertical",
#     #legend.key.size = unit(20, "pt"),
#     legend.title = element_text(size = 15),
#     legend.text = element_text(size = 15),
#     #axis.title.y =element_blank(),
#     axis.text.y = ggtext::element_markdown(hjust = 0),
#     axis.ticks = element_blank(),
#     axis.text.x = ggtext::element_markdown(vjust = 0),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     strip.background = element_rect(fill="white", color = "white"), #chocolate4
#     strip.text = element_text(color="black", size = 15, hjust = 0)) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_BayesianTrends_allScales_fns.png"),
#        plot = last_plot(), 
#        width = 4400, height = 3800, units = "px")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### FIGURE 3 - Bayesian summary ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# combined diversity and function estimates (x axis)
# scale on y axis
# habitat types as colors
# pointrange plot

# summarize into groups
pars_sum <- pars_all %>%
  full_join(fns_labels %>% dplyr::select(Function, Group_function), by = c("fns" = "Function")) %>%
  dplyr::select(-fns) %>%
  group_by(Group_function, LC, scale) %>%
  summarize(across(everything(), function(x) mean(x, na.rm=T)))

# plot
ggplot(data = pars_sum %>% filter(!is.na(LC)) %>%
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
               significance =  ifelse(sign(lower.HPD)!= sign(upper.HPD), "ns", as.character(as.factor(LC)))) %>%
         mutate(significance = factor(significance, levels = c("Cropland", "Grassland", "Shrubland", "Woodland", "ns"))),
       aes(fill = significance, color = LC, 
           y = LC, x = PA_rank.trend))+
  
  geom_vline(aes(xintercept=0), color="black")+
  geom_pointrange(aes(xmin = lower.HPD, xmax = upper.HPD),
                  position = position_dodge(width = 0.3),
                  size = 3, shape = 21) +
  #coord_flip()+
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

ggsave(filename=paste0(here::here(), "/figures/Results_pointrange_BayesianTrends_allScales_grouped.png"),
       plot = last_plot(), 
       #width = 2000, height = 1500,
       units = "px")
