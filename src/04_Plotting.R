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
library(ggridges) #to plot distributions (Bayesian)

# load background map
world.inp <- map_data("world")

source("src/00_Parameters_functions.R")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

# load pairs of PA and nonPA
pa_pairs <- read.csv(file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_2023-02-22.csv"))
head(pa_pairs)

# load p values and effect sizes
load(file=paste0(here::here(), "/results/p_1000_trails.RData")) #p_list
load(file=paste0(here::here(), "/results/d_1000_trails.RData")) #d_list

# load Bayesian results
load(file=paste0(here::here(), "/results/delta_1000_trails_Bayesian_sample100.RData")) #delta_sample

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Sampling locations ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# extract list of sampling locations actually used in comparison
data_locations <- data_glob %>% 
  filter(Order_ID %in% unique(pa_pairs$Order_ID) | 
           Order_ID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Longitude_c,Latitude_c,Order_ID, PA, LC)
data_locations #nrow=319
nrow(data_locations %>% filter(PA==1)) #120 PAs

ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80", show.legend = FALSE) +
  xlim(-180, 180) +
  ylim(-180, 180) +
  
  geom_point(data=data_locations, aes(x=Longitude_c, y=Latitude_c, fill=as.factor(PA), shape=LC), alpha=0.5)+
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
ggsave(filename=paste0(here::here(), "/figures/Data_locations_global.png"),
       plot = last_plot())


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
  scale_fill_manual(values=c("gold3", "limegreen", "forestgreen"))
ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_mahal.distance_global.png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## P-value per land use type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

p_df <- do.call(rbind, p_list)
str(p_df)

p_df <- p_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Violin plot p-value per land use type ####

ggplot(data=p_df, aes(x=round(p_value,2), y=Label))+
  geom_violin(aes(fill = Organism), scale="width", width=0.5)+ #scale by width of violin bounding box
  geom_vline(aes(xintercept=0.05, linetype = "p = 0.05"))+
  geom_vline(aes(xintercept=0, linetype = "p = 0"))+
  scale_color_manual(values = c("p = 0.05" = "dashed", "p = 0" = "solid"))+
  stat_summary(fun = "mean", geom = "point", color = "black")+
  facet_wrap(vars(lc))+
  scale_x_continuous(labels = function(x) format(x, nsmall = 2))+#round p-value axis label
  scale_fill_viridis_d()+
  xlab("p-value")+
  theme_bw() + # use a white background
  theme(legend.position = "right", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=7),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(filename=paste0(here::here(), "/figures/Data_violin_p-value_global.png"),
       plot = last_plot())

# save data for plot
write.csv(p_df, file=paste0(here::here(), "/figures/Data_violin_p-value_global.csv"))

p_summary <- p_df %>% 
  dplyr::select(-run) %>%
  #pivot_longer(cols = c(p_value, ci_lower, ci_upper, t_stats),
  #             names_to = "metric") %>%
  group_by(lc, fns, Label, Group_function, Organism) %>%
  summarize(across(everything(), .fns = list("mean"=mean, "SD"=sd)))
p_summary
write.csv(p_summary, file=paste0(here::here(), "/figures/Data_violin_p-value_summary_global.csv"))

# check for significant mean p-values
#p_summary %>% arrange(p_value_mean)

sink(file=paste0(here::here(), "/results/P-values_global.txt"))
cat("#################################################", sep="\n")
cat("##  Significant p-values from global analysis  ##", sep="\n")
cat(paste0("##  Sys.Date() ", Sys.Date(), "  ##"), sep="\n")
cat("#################################################", sep="\n")
print(p_summary %>% ungroup() %>%
        filter(p_value_mean <= 0.05) %>%
        dplyr::select(lc, fns, starts_with("p_value"), starts_with("ci"),
                      starts_with("t_stats")),
      n=nrow(p_summary %>% filter(p_value_mean <= 0.05)))
sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Effect size per LC type ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
d_df <- do.call(rbind, d_list)
str(d_df)

d_df <- d_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Violin plot effect size per LC type ####

ggplot(data=d_df, aes(x=round(effect,2), y=Label))+
  geom_vline(aes(xintercept=0.8, linetype = "0.8"), color="grey60")+
  geom_vline(aes(xintercept=-0.8, linetype = "-0.8"), color="grey60")+ #large effect
  geom_vline(aes(xintercept=0.5, linetype = "0.5"), color="grey60")+
  geom_vline(aes(xintercept=-0.5, linetype = "-0.5"), color="grey60")+ #medium effect
  geom_vline(aes(xintercept=0.2, linetype = "0.2"), color="grey60")+ #small effect
  geom_vline(aes(xintercept=-0.2, linetype = "-0.2"), color="grey60")+
  #geom_vline(aes(xintercept=0, linetype = "0"))+
  scale_linetype_manual(values = c("-0.8" = "dotted",
                                   "-0.5" = "dashed", "-0.2" = "solid",
                                   "0.2" = "solid","0.5" = "dashed",
                                   "0.8" = "dotted"))+
  
  geom_violin(aes(fill = Organism), scale="width", width=0.5)+ #scale by width of violin bounding box
  stat_summary(fun = "mean", geom = "point", color = "black")+
  facet_wrap(vars(lc))+
  scale_x_continuous(labels = function(x) format(x, nsmall = 0))+#round c-value axis label
  scale_fill_viridis_d()+
  xlab("Effect size")+
  theme_bw() + # use a white background
  theme(legend.position = "right", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=7),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(filename=paste0(here::here(), "/figures/Data_violin_d-value_global.png"),
       plot = last_plot())

# save data for plot
write.csv(d_df, file=paste0(here::here(), "/figures/Data_violin_d-value_global.csv"))

d_summary <- d_df %>% 
  dplyr::select(-run) %>%
  #pivot_longer(cols = c(p_value, ci_lower, ci_upper, t_stats),
  #             names_to = "metric") %>%
  group_by(lc, fns, Label, Group_function, Organism) %>%
  summarize(across(everything(), .fns = list("mean"=mean, "SD"=sd)))
d_summary
write.csv(d_summary, file=paste0(here::here(), "/figures/Data_violin_d-value_summary_global.csv"))

# mean per lc type
d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_mean, lower_mean, upper_mean), 
                   function(x) mean(x)))
d_summary %>% ungroup() %>% group_by(lc) %>% 
  summarize(across(c(effect_mean, lower_mean, upper_mean), 
                   function(x) mean(abs(x))))

# mean per Group_function
d_summary %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_mean, lower_mean, upper_mean), 
                   function(x) mean(abs(x))))
d_summary %>% ungroup() %>% group_by(Group_function) %>% 
  summarize(across(c(effect_mean, lower_mean, upper_mean), 
                   function(x) mean(x)))

# check for significant mean p-values
#d_summary %>% arrange(p_value_mean)

sink(file=paste0(here::here(), "/results/D-values_global.txt"))
cat("#################################################", sep="\n")
cat("##  Significant d-values from global analysis  ##", sep="\n")
cat(paste0("##  Sys.Date() ", Sys.Date(), "  ##"), sep="\n")
cat("#################################################", sep="\n")
print(d_summary %>% ungroup() %>%
  filter(abs(effect_mean) >= 0.2) %>%
  dplyr::select(lc, fns, starts_with("effect"), starts_with("lower"),
                starts_with("upper")),
  n=nrow(d_summary %>% filter(abs(effect_mean) >= 0.2)))
sink()

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Pointrange effect size per LC type ####
# with confidence intervals 

ggplot(data = d_summary, 
       aes(y = effect_mean, x = Label, ymin = lower_mean, ymax = upper_mean))+
  
  geom_hline(aes(yintercept=0.8, linetype = " 0.8"), color="grey60")+
  geom_hline(aes(yintercept=-0.8, linetype = "-0.8"), color="grey60")+ #large effect
  geom_hline(aes(yintercept=0.5, linetype = " 0.5"), color="grey60")+
  geom_hline(aes(yintercept=-0.5, linetype = "-0.5"), color="grey60")+ #medium effect
  geom_hline(aes(yintercept=0.2, linetype = " 0.2"), color="grey60")+ #small effect
  geom_hline(aes(yintercept=-0.2, linetype = "-0.2"), color="grey60")+
  #geom_hline(aes(yintercept=0, linetype = "0"))+
  scale_linetype_manual(values = c("-0.8" = "dotted",
                                   "-0.5" = "dashed", "-0.2" = "solid",
                                   " 0.2" = "solid"," 0.5" = "dashed",
                                   " 0.8" = "dotted"))+
  
  geom_linerange(size=1, color="black") +
  geom_point(aes(color=Group_function), shape=19, size=3)+
  #scale_color_viridis_d()+
  coord_flip()+
  ylab("Effect size")+
  facet_wrap(vars(lc))+
  theme_bw() + # use a white background
  theme(legend.position = "right", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.spacing = unit(1, "lines"))
ggsave(filename=paste0(here::here(), "/figures/Data_pointrange_d-value_global.png"),
       plot = last_plot())



#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplots values ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# extract list of sampling locations actually used in comparison
data_values <- data_glob %>% 
  filter(Order_ID %in% unique(pa_pairs$Order_ID) | 
         Order_ID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Order_ID, PA, LC, fns) %>%
  mutate(across(fns, .fns = function(x) { (x - mean(x)) / sd(x)})) %>%
  
  pivot_longer(cols = c(fns), names_to = "fns") %>%
  full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))
data_values

ggplot(data_values, aes(x = Label, y = value)) +
  geom_hline(yintercept = 0, lty="dashed", color="grey60") +
  geom_boxplot(aes(fill=as.factor(PA)))+#,fatten = NULL) + # activate to remove the median lines
  scale_fill_manual(values=c("orange","darkgreen"),name = NULL, labels = c("Non-Protected","Protected")) +
  #scale_y_continuous(limits=c(-3,13)) + 
  
  xlab("")+ ylab("")+
  facet_wrap(vars(LC), ncol=1, nrow=3)+
  theme_bw() +
  theme(axis.text.x=element_text(size=5, angle=45, hjust=1),
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", legend.text = element_text(size=15), axis.text.y = element_text(size=15))
ggsave(filename=paste0(here::here(), "/figures/Data_boxplot_values_global.png"),
       plot = last_plot())

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Bayesian results (distributions) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

delta_df <- delta_sample %>% 
  bind_rows(.id="LC") %>% 
  pivot_longer(cols=Soil_carbon_service:Diss_invert_std, names_to = "fns") %>%
  mutate("value"=value*(-1)) #convert all values to get difference PA-nonPA
head(delta_df)

# add categories of functions
delta_df <- delta_df %>% full_join(fns_labels, by=c("fns"="Function")) %>%
  mutate("Label" = factor(Label, levels = rev(fns_labels$Label)),
         "Organism" = factor(Organism, levels = unique(fns_labels$Organism)))

ggplot(delta_df, aes(x = value, y = Label, color=Group_function, fill=Group_function)) +
  geom_vline(xintercept=0, color="grey60")+
  ggridges::geom_density_ridges() +
  theme_ridges() + 
  facet_wrap(vars(LC), ncol=3)+
  #xlim(-10, 10)+
  xlab ("")+ ylab("")+
  #theme_bw()+
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())
ggsave(filename=paste0(here::here(), "/figures/Data_distr_delta1000_global.png"),
       plot = last_plot())


# save CI for delta_summary
delta_summary <- delta_df %>% group_by(fns, LC) %>% 
  summarize(across(value, list("mean"= mean,
                               "CI_05"=function(x) quantile(x,0.05),
                               "CI_95"=function(x) quantile(x,0.95)
                               )))
delta_summary
write.csv(delta_summary, file=paste0(here::here(), "/figures/Data_distr_delta1000_global.csv"))
