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
library(ggridges) #to plot distributions (Bayesian)
library(ggdist) #to plot distributions (Bayesian)
library(ggpubr) #plot multiple ggplots

# load background map
world.inp <- map_data("world")

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load soil biodiversity data ####
data_glob <- read_csv(paste0(here::here(), "/intermediates/Data_global.csv"))
data_glob

# load pairs of PA and nonPA
pa_pairs <- read.csv(file=paste0(here::here(), "/intermediates/Pairs_paNonpa_1000trails_global.csv"))
head(pa_pairs)

# load effect sizes
load(file=paste0(here::here(), "/results/d_1000_trails_global.RData")) #d_list

# load Bayesian results from PA_type comparison
load(file=paste0(here::here(), "/results/pars_PAtypes_Bayesian_df_global.RData")) #pars_sample

pars_sample <- pars_sample %>% group_by(fns, lc) %>% slice_sample(n=1000)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Sampling locations ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# extract list of sampling locations actually used in comparison
data_locations <- data_glob %>% 
  filter(Order_ID %in% unique(pa_pairs$ID) | 
         Order_ID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Longitude_c,Latitude_c,Order_ID, PA, LC)
data_locations #nrow=262
nrow(data_locations %>% filter(PA==1)) #68 PAs

ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region),  show.legend = FALSE, fill = "grey80", color="grey75") +
  xlim(-180, 180) +
  ylim(-180, 180) +
  
  geom_point(data=data_locations, aes(x=Longitude_c, y=Latitude_c, shape=as.character(PA), color=LC), size=2)+
  scale_shape_manual(values = c(1,4))+ #label = c("Protected", "Unprotected")
  scale_color_manual(values = c("gold3", "limegreen", "forestgreen"))+
  
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

ggplot(data = d_df, 
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

  ggdist::stat_pointinterval(fatten_point=1.2, shape=21) +
  coord_flip()+
  ylab("Effect size")+
  facet_wrap(vars(lc))+
  theme_bw() + # use a white background
  theme(legend.position = "bottom", axis.title.y =element_blank(),
        axis.text.y = element_text(size=10),  
        axis.text.x = element_text(size=10),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(filename=paste0(here::here(), "/figures/Data_pointrange_d-value_global.png"),
       plot = last_plot())


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Boxplots values ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

# extract list of sampling locations actually used in comparison
data_values <- data_glob %>% 
  filter(Order_ID %in% unique(pa_pairs$ID) | 
         Order_ID %in% unique(pa_pairs$nonPA)) %>%
  dplyr::select(Order_ID, PA, LC, all_of(fns)) %>%
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
## Bayesian results (PA types) ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -

pars_df <- pars_sample %>% dplyr::select(-sigma) %>% pivot_longer(cols = c("a[1]":"a[11]"))
a <- ggplot(pars_df,
       aes(y=name, x=value))+
  ggdist::stat_halfeye(normalize="xy")+ 
  facet_wrap(vars(lc))+
  theme_bw()

a %+% subset(pars_df %>% filter(fns==fns[6]))
