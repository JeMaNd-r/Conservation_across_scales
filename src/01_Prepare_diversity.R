#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Prepare diversity data           #
#          author: Romy Zeiss               #
#            date: 2023-10-26               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

data_wd <- "T:/_students/Romy/Sequences/Rarified_data_2024-06-28/"

library(tidyverse)

temp_scale <- "global"
#temp_scale <- "continental"
#temp_scale <- "regional"


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Load raw data ####

if(temp_scale == "global") {
  raw_data <- list(
    "div_bac" = read_delim(paste0(data_wd, "Drylands/Dry16S_rarefyd_DivMetrics.csv")),
    "div_euk" = read_delim(paste0(data_wd, "Drylands/Dry18S_rarefyd_DivMetrics.csv")),
    
    "dis_bac" = read_delim(paste0(data_wd, "Drylands/Dry16SBac_rarefyd_JacDist.csv")),
    "dis_euk" = read_delim(paste0(data_wd, "Drylands/Dry18sEukaryotes_rarefyd_JacDist.csv")),
    "dis_fun" = read_delim(paste0(data_wd, "Drylands/Dry18sFungi_rarefyd_JacDist.csv")),
    "dis_inv" = read_delim(paste0(data_wd, "Drylands/Dry18sInvertebrate_rarefyd_JacDist.csv")),
    "dis_pro" = read_delim(paste0(data_wd, "Drylands/Dry18sProtists_rarefyd_JacDist.csv"))
  )
}

## for alternative dataset
# if(temp_scale == "global") {
#   raw_data <- list(
#     "div_bac" = read_delim(here::here("Global/Glob16S_rarefyd_DivMetrics.csv")),
#     "div_euk" = read_delim(here::here("Global/Glob18s_rarefyd_DivMetrics.csv")),
#     
#     "dis_bac" = read_delim(here::here("Global/Glob16SBac_rarefyd_JacDist.csv")),
#     "dis_euk" = read_delim(here::here("Global/Glob18sEukaryotes_rarefyd_JacDist.csv")),
#     "dis_fun" = read_delim(here::here("Global/Glob18sFungi_rarefyd_JacDist.csv")),
#     "dis_inv" = read_delim(here::here("Global/Glob18sInvertebrate_rarefyd_JacDist.csv")),
#     "dis_pro" = read_delim(here::here("Global/Glob18sProtists_rarefyd_JacDist.csv"))
#   )
# }

if(temp_scale == "continental") {
  raw_data <- list(
    "div_bac" = read_delim(paste0(data_wd, "Europe/Europ16S_rarefyd_DivMetrics.csv")),
    "div_euk" = read_delim(paste0(data_wd, "Europe/Europ18S_rarefyd_DivMetrics.csv")),
    
    "dis_bac" = read_delim(paste0(data_wd, "Europe/Europ16SBac_rarefyd_JacDist.csv")),
    "dis_euk" = read_delim(paste0(data_wd, "Europe/Europ18sEukaryotes_rarefyd_JacDist.csv")),
    "dis_fun" = read_delim(paste0(data_wd, "Europe/Europ18sFungi_rarefyd_JacDist.csv")),
    "dis_inv" = read_delim(paste0(data_wd, "Europe/Europ18sInvertebrate_JacDist.csv")),
    "dis_pro" = read_delim(paste0(data_wd, "Europe/Europ18sProtists_JacDist.csv"))
  )
}

if(temp_scale == "regional"){
  raw_data <- list(
    "div_bac" = read_delim(paste0(data_wd, "Portugal/Port16S_rarefyd_DivMetrics.csv")),
    "div_euk" = read_delim(paste0(data_wd, "Portugal/Port18S_rarefyd_DivMetrics.csv")),
    
    "dis_bac" = read_delim(paste0(data_wd, "Portugal/Port16SBac_rarefyd_JacDist.csv")),
    "dis_euk" = read_delim(paste0(data_wd, "Portugal/Port18sEukaryotes_rarefyd_JacDist.csv")),
    "dis_fun" = read_delim(paste0(data_wd, "Portugal/Port18sFungi_rarefyd_JacDist.csv")),
    "dis_inv" = read_delim(paste0(data_wd, "Portugal/Port18sInvertebrate_rarefyd_JacDist.csv")),
    "dis_pro" = read_delim(paste0(data_wd, "Portugal/Port18sProtists_rarefyd_JacDist.csv"))
  )
}

raw_data

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Calculate mean Jaccard distance for each sample ID ####

# select dissimilarity data
dis_data <- raw_data[names(raw_data)[str_detect(names(raw_data), "dis_")]]

# remove Control sites in European dataset
if(temp_scale == "continental") {
  dis_data <- lapply(dis_data, function(x) {
    x <- x %>% 
      rename_with(toupper, starts_with("Lucas")) %>%
      dplyr::select(colnames(.)[!(str_detect(colnames(.), "Contr"))]) %>% #exclude Control plots
      mutate("...1" = toupper(`...1`)) %>%
      mutate("...1" = str_replace(`...1`, "LUCAS", "")) %>%
      mutate("...1" = as.numeric(`...1`))  %>%
      filter(!is.na(`...1`))
  })
}

# calculate mean and sd distance for each sample
dis_data <- lapply(dis_data, function(x) {
  temp_mean <- apply(x[,-1], 1, function(y) mean(y, na.rm=T))
  temp_sd <- apply(x[,-1], 1, function(y) sd(y, na.rm=T))
  tibble("SampleID" = pull(x[,1]), 
         temp_mean, temp_sd)
})

# add column with name (organism)
dis_data <- lapply(setNames(names(dis_data), names(dis_data)), function(nameindex) {
  if(str_detect(nameindex, "bac")) temp_name <- "Bac"
  if(str_detect(nameindex, "euk")) temp_name <- "Euk"
  if(str_detect(nameindex, "fun")) temp_name <- "Fungi"
  if(str_detect(nameindex, "inv")) temp_name <- "Invertebrate"
  if(str_detect(nameindex, "pro")) temp_name <- "Protist"
  
  dis_data[[nameindex]]$metric <- paste0(temp_name, "_JaccDist")
  dis_data[[nameindex]]
  })

dis_data <- do.call(rbind, dis_data)

# merge results into 1 table
dis_data <- dis_data %>% 
  rename("av" = temp_mean, 
         "sd" = temp_sd) %>%
  pivot_wider(id_cols = "SampleID", names_from = "metric", 
              values_from = c("av", "sd"), names_glue = "{metric}_{.value}")
dis_data

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Merge diversity and dissimilarity measures ####

# prepare diversity measures
div_data <- raw_data[names(raw_data)[str_detect(names(raw_data), "div_")]]
div_data

# remove Control sites in European dataset
if(temp_scale == "continental") {
  div_data <- lapply(div_data, function(x) {
    x <- x %>% 
      rename_with(toupper, starts_with("Lucas")) %>%
      dplyr::select(colnames(.)[!(str_detect(colnames(.), "Contr"))]) %>% #exclude Control plots
      mutate("sampleID" = toupper(sampleID)) %>%
      mutate("sampleID" = str_replace(sampleID, "LUCAS", "")) %>%
      mutate("Sample_ID" = as.numeric(sampleID))  %>%
      dplyr::select(Sample_ID, everything(), -sampleID) %>%
      filter(!is.na(Sample_ID))
  })
}

div_data <- div_data[[1]] %>% full_join(div_data[[2]], by="Sample_ID")
div_data

# merge dissimilarity & div. measures
data <- dis_data %>% full_join(div_data, by=c("SampleID" = "Sample_ID"))
data

# rename sites in Portuguese dataset
if(temp_scale == "regional") {
  data <- data %>% mutate("SampleID" = str_replace(SampleID, "X", ""))
}

rm(div_data, dis_data)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Save ####
write_csv(data, file=paste0(data_wd, "Diversity_", temp_scale, ".csv"))


