#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Plot data and results             #
#            per spatial scale              #
#          author: Romy Zeiss               #
#            date: 2025-01-02               #
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

if(temp_scale == "global"){
  lc_names <- "Dryland" #lc_names[lc_names != "Other" & lc_names != "Cropland"]
} 
if(temp_scale == "continental"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
}
if(temp_scale == "regional"){
  lc_names <- lc_names[lc_names != "Other" & lc_names != "Shrubland" & lc_names != "Dryland"]
}
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
## Load soil biodiversity data ####
data_clean <- read_csv(paste0(here::here(), "/intermediates/Data_clean_", temp_scale, ".csv"))
data_clean

# exclude LC types if needed ####
data_clean <- data_clean %>% filter(LC %in% lc_names)

if(temp_scale == "global") data_clean <- data_clean %>% mutate(LC = "Dryland")

# define min_size
min_size <- min(table(data_clean$LC, 
                      data_clean$PA)[table(data_clean$LC, 
                                           data_clean$PA)
                                     >0])

# load pairs of PA and nonPA
pa_pairs <- read_csv(file=sort(list.files(here::here("intermediates", temp_scale), pattern = "Pairs_paNonpa_1000trails_", full.names = TRUE), decreasing = TRUE)[1])
head(pa_pairs)

# load effect sizes
load(file=paste0(here::here(), "/results/d_1000_trails_", temp_scale, ".RData")) #d_list

# Note: other output data are used in the 04b_Plotting_allScales.R script


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
data_locations #G: nrow=131, C: 316, R: 161
write_csv(data_locations, file = paste0(here::here(), "/results/Locations_", temp_scale, ".csv"))
nrow(data_locations %>% filter(PA==1)) #G: 28 PAs, G-together: 39; C: 48, R: 36
nrow(data_locations %>% filter(PA==0)) #G: 93 PAs, G-together: 92; C: 269, R: 125

# set limits for point maps
if(temp_scale == "global") temp_limits <- c(-180, 180, -180, 180)
if(temp_scale == "continental") temp_limits <- c(-10, 35, 35, 70)
if(temp_scale == "regional") temp_limits <- c(-9, -6, 40.5, 42.5)

# plot - FIGURE 1
ggplot()+
  geom_map(data = world.inp, map = world.inp, 
           aes(map_id = region),  show.legend = FALSE, 
           fill="white", color = "grey80", linewidth = 0.25) + #G:0.15, C+R:
  xlim(temp_limits[1], temp_limits[2])+
  ylim(temp_limits[3], temp_limits[4])+
  
  geom_point(data=data_locations, aes(x=Longitude, y=Latitude, 
                                      shape = as.character(PA), color=LC, 
                                      size = as.character(PA)),
             stroke = 2)+
             #stroke = 1.4, color = "#000000")+ #increase circle line width; G: 2 (1.4), C+R:3
  scale_shape_manual(values = c("0" = 19, "1" = 1))+ #label = c("Protected", "Unprotected")
  scale_size_manual(values = c("0" = 2, "1" = 6))+ #G: 1.4,4.5/0.3, 1, C+R:3,8/ 0.6,2
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000",
                                "Dryland" = "#000000"))+ 
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
       plot = last_plot(),
       width = 10, height = 8)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.1 - Maps for random-slope model ####
# extract list of sampling locations actually used in comparison
data_locations <- data_clean %>%
  filter(LC %in% lc_names)
data_locations #G: nrow=248; C: 745, R: 270
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
  scale_size_manual(values = c("0" =2, "1" = 6))+ #G:+C: 2,4; ,R:3,8 
  scale_color_manual(values = c("Cropland" = "#4A2040",
                                "Grassland" = "#E69F00",
                                "Shrubland" = "#0072B2", 
                                "Woodland" = "#009E73", 
                                "Other" = "#000000",
                                "Dryland" = "#000000"))+ 
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
       plot = last_plot(),
       width = 10, height = 8)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pairing ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Save some numbers
sink(paste0(here::here(), "/results/Numbers_d-value_", temp_scale, ".txt"))

print(temp_scale)
# look what non-protected sites have (not) been paired to any PA
cat("Number of non-protected sites have never been paired to any PA:")
#hist(table(pa_pairs$nonPA))  # frequency distribution of the use of sites from all runs
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

# ## one-sided
# ggplot(data = d_summary %>%
#          filter(lc %in% lc_names) %>%
#          mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
#                                          abs(effect_ci_17), 
#                                          abs(effect_ci_83))),
#        aes(y = abs(effect_ci_min66), x = abs(effect_median),
#            color=as.factor(sign(effect_median))
#            ))+
# 
#   annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.2, fill = "grey", alpha=0.3)+
#   annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.5, fill = "grey", alpha=0.3)+
#   annotate("rect", ymin = -Inf, ymax = Inf, xmin=-Inf, xmax=0.8, fill = "grey", alpha=0.3)+
# 
#   geom_point(aes(shape = Group_function), size = 5)+
#   ggrepel::geom_text_repel(aes(label = Label, shape = Group_function), size = 3)+
#   
#   xlab("Median effect")+ ylab("Lower CI (17%)")+
#   scale_y_continuous(breaks = c(0.2, 0.5, 0.8))+
#   scale_x_continuous(breaks = c(0.2, 0.5, 0.8))+
#   
#   facet_wrap(vars(lc), scales = "free")+
#   theme_bw() + # use a white background
#   theme(legend.position = "bottom",
#         legend.direction = "vertical",
#         axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"), #chocolate4
#         strip.text = element_text(color="black")) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianSD_", temp_scale, ".png"),
#        plot = last_plot())
# 
# ## heatmap
# ggplot(data = d_summary %>%
#          filter(lc %in% lc_names) %>%
#          mutate(effect_ci_min66 = ifelse(abs(effect_ci_17)<abs(effect_ci_83), 
#                                          abs(effect_ci_17), 
#                                          abs(effect_ci_83))) %>%
#          mutate(effect_ci_min66f = cut(effect_ci_min66,
#                                        breaks=c(0, 0.2, 0.5, 0.8, Inf),
#                                        labels=c("ns", "small", "medium", "large"))) %>%
#          filter(!is.na(effect_ci_min66)),
#        aes(x = lc, y = Label, alpha=effect_ci_min66f, 
#            fill=as.factor(sign(effect_median))))+
# 
#   geom_tile()+
#   scale_fill_manual(values = c("-1" = "#fc8d59", "0" = "#ffffbf", "1" = "#91bfdb"),
#                     name = "Direction of effect")+
#   scale_alpha_manual(values = c("ns" = 0.05, "small" = 0.3, "medium" = 0.65, "large" = 1),
#                      name = "Minimum effect size (66% CI)")+
#   theme_bw() + # use a white background
#   theme(legend.position = "bottom",
#         legend.direction = "vertical",
#         #axis.title.y =element_blank(),
#         axis.text.y = element_text(size=10),
#         axis.text.x = element_text(size=10),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill="white"), #chocolate4
#         strip.text = element_text(color="black")) #white
# ggsave(filename=paste0(here::here(), "/figures/Results_d-value_medianSD_", temp_scale, ".png"),
#        plot = last_plot())


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Random-intercept models (PA types) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.3 & 3.5 & 3.7 - Bayesian PA types ####

pars_intercept <- read_csv(file=paste0(here::here(), "/results/PAtypes_Bayesian_", temp_scale, "_emmeans.csv"))

pars_intercept <- pars_intercept %>%
  full_join(protect_type %>% 
                mutate(PA_rank = as.character(PA_rank)),
            by=c("PA_rank"))  %>%
    mutate("Label_pa" = ifelse(PA_rank=="Unprotected", "Unprotected", PA_type)) %>%
    mutate("Label_pa" = factor(Label_pa, levels = rev(c(protect_type$PA_type, "Unprotected")))) %>%
    full_join(fns_labels, by=c("fns"="Function")) %>%
    mutate("Label_fns" = factor(Label, fns_labels$Label))

if(temp_scale == "global") pars_intercept <- pars_intercept %>% dplyr::select(-LC)
  
## OLD: Model output only, we want emmeans from model
# # transform parameters to long format and assign labels
# pars_long <- pars_sample %>% 
#   dplyr::select(-sigma, -sigma_a, -mu_a) %>% 
#   pivot_longer(cols = colnames(.)[!(colnames(.) %in% c("lc", "fns"))])%>%
#   full_join(protect_type %>% 
#               mutate(PA_rank = as.character(PA_rank)), 
#             by=c("name" = "PA_rank"))  %>%
#   mutate("Label_pa" = ifelse(name=="11", "Unprotected", PA_type)) %>%
#   mutate("Label_pa" = factor(Label_pa, levels = rev(c(protect_type$PA_type, "Unprotected")))) %>%
#   full_join(fns_labels, by=c("fns"="Function")) %>%
#   mutate("Label_fns" = factor(Label, fns_labels$Label))
# head(pars_long)
# 
# # save summary (i.e., data from figure)
# pars_summary <- pars_long %>% group_by(lc, fns) %>% 
#   summarize(across(value, list("mean"  = function(x) mean(x, na.rm=TRUE), 
#                                "median" = function(x) median(x, na.rm=TRUE), 
#                                "ci_2.5" = function(x) quantile(x, 0.05, na.rm=TRUE), 
#                                "ci_17" = function(x) quantile(x, 0.17, na.rm=TRUE), 
#                                "ci_83" = function(x) quantile(x, 0.83, na.rm=TRUE), 
#                                "ci_97.5" = function(x) quantile(x, 0.975, na.rm=TRUE)))) %>%
#   arrange(lc, fns)
# pars_summary
# write_csv(pars_summary, file=paste0(here::here(), "/figures/Results_intercept_BayesianTrends_summary_", temp_scale, ".csv"))

# extract sample size
n_table <- data_clean %>% filter(LC!="Other") %>%
  group_by(LC, PA, PA_type, PA_rank) %>% count() %>%
  pivot_wider(id_cols = c("PA", "PA_type", "PA_rank"), 
              names_from = LC, values_from = n) %>% 
  arrange(PA_rank) %>% ungroup() %>%
  dplyr::select(-PA_rank)
n_table
write_csv(n_table, file=paste0(here::here(), "/figures/Results_intercept_BayesianTrends_nTable_", temp_scale, ".csv"))

ggplot(pars_intercept %>% filter(!is.na(Label)) %>% #filter(!is.na(PA_type)) %>%
                     # add number of sizes to plot
                     rbind(data_clean %>% filter(LC!="Other") %>%
                             group_by(PA_type) %>% count() %>%
                             rename(emmean=n) %>%
                             mutate(fns="Number of sites",
                                    PA_protected=NA,
                                    PA_rank = NA,
                                    scale = temp_scale,
                                    Label_pa=ifelse(is.na(PA_type), "Unprotected", PA_type),
                                    Label="Number of sites",
                                    Label_fns = Label,
                                    Label_short = Label,
                                    Group_function=NA,
                                    Organism=NA, 
                                    lower.HPD = NA,
                                    upper.HPD = NA) %>%
                             dplyr::select(colnames(pars_intercept))) %>%
                     filter(!is.na(Label_pa)) %>%
                     mutate(Label_fns = factor(Label_fns, levels = c(labels_order, "Number of sites")),
                            Label_n = ifelse(Label_fns == "Number of sites", "Number of sites", "Other fns")),
       aes(y=Label_pa, x=emmean, #color=LC,
           xmin = lower.HPD, xmax = upper.HPD, shape = Label_n))+
  
  ## adapt for scale
  annotate("rect", ymin = -Inf, ymax = 4+0.5, xmin=-Inf, xmax=Inf, fill = "grey90", alpha=0.5)+ #global: 4, C: 3, R: 2
  annotate("rect", ymin = -Inf, ymax = 1+0.5, xmin=-Inf, xmax=Inf, fill = "grey85", alpha=0.5)+

  geom_pointrange(position=position_dodgejust(width=0.5))+  
  # ggdist::stat_pointinterval(fatten_point=1, shape=3, 
  #                            position=position_dodgejust(width=0.5))+ 
  facet_wrap(vars(Label_fns), scales = "free_x", ncol=6, drop=FALSE)+
  
  scale_shape_manual(values=c("Number of sites" = 3, "Other fns" = 19))+
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000",
                              "Dryland" = "#000000"), name="Habitat type")+
  ylab("")+ xlab("")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 13))
ggsave(filename=paste0(here::here(), "/figures/Results_intercept_BayesianTrends_", temp_scale, ".png"),
       plot = last_plot(),
       width=12, height=10)


# ggdist to show whole distribution and not just emmeans
pred_list <- get(load(paste0(here::here(), "/results/PAtypes_Bayesian_", temp_scale, "_sample10k.RData")))

pred_list <- pred_list %>%
  ungroup() %>%
  rowwise() %>%
  mutate(value = sum(c_across(c(Soil_carbon_service, OM_decomposition_service:Invertebrate_JaccDist_av)), na.rm=TRUE)) %>%
  dplyr::select(LC:fns, value, -Soil_carbon_service) %>%
  full_join(protect_type %>% 
              mutate(PA_rank = as.character(PA_rank)),
            by=c("PA_rank"))  %>%
  mutate("Label_pa" = ifelse(PA_rank=="Unprotected", "Unprotected", PA_type)) %>%
  mutate("Label_pa" = factor(Label_pa, levels = rev(c(protect_type$PA_type, "Unprotected")))) %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label_fns" = factor(Label, fns_labels$Label))

ggplot(data = pred_list %>%
         mutate(Label = factor(Label, levels = labels_order)), 
       
       aes(x = Label_pa, y = .epred, color = ordered(LC))) +
  stat_slab(normalize = "panels", position = "dodge", scale = 1.5) +
  facet_wrap(vars(Label), scales = "free", ncol=6)+
  coord_flip()+
  #scale_fill_brewer(palette = "Greys") +
  scale_color_manual(values=c("Cropland" = "#4A2040",
                              "Grassland" = "#E69F00",
                              "Shrubland" = "#0072B2", 
                              "Woodland" = "#009E73", 
                              "Other" = "#000000",
                              "Dryland" = "#000000"), name="Habitat type")+
  theme_void()+
  theme(axis.text = element_text(),
        panel.grid.minor.x =  element_line(color = "grey"),
        strip.text = element_text(size = 15, hjust=0),
        legend.position = c(0.8, 0.1),
        legend.box = "horizontal")
ggsave(filename=paste0(here::here(), "/figures/Results_intercept_parsBayesian_", temp_scale,".png"),
       plot = last_plot(),
       width=15, height=10)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Random-slope model (PA ranks/ levels) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### APPENDIX FIGURE 3.3 & 3.6 & 3.8 (PA ranks/ levels) ####

pred_list <- get(load(paste0(here::here(), "/results/PAranks_Bayesian_", temp_scale, "_sample10k.RData")))

ggplot(data = pred_list %>%
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
                              "Other" = "#000000",
                              "Dryland" = "#000000"), name="Habitat type")+
  scale_x_continuous(limits = c(1, 11), breaks = c(2, 10), minor_breaks = c(2,4,6,8, 10))+
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
ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_estimates_", temp_scale, ".png"),
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
                              "Other" = "#000000",
                              "Dryland" = "#000000"), name="Habitat type")+
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

