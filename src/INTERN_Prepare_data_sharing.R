#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#     What data to share & what not?        #
#          author: Romy Zeiss               #
#            date: 2025-03-25               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# This script aims to compare the published versions of the global datasets
# to see, what rows and columns can be published/shared with the paper and
# which data is still under embargo.

library(here)
library(tidyverse)
library(readxl)

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## GLOBAL ----------------------------------------- ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Load published datasets ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# download file from Maestre et al. 2015
download.file(url = "https://figshare.com/ndownloader/files/2365595",
              destfile = here::here("data_raw"), 
              mode = "wb")

# download file from Maestre et al. 2022 (BIODESERT)
download.file(url = "https://figshare.com/ndownloader/files/32115179",
              destfile = here::here("data_raw"), 
              mode = "wb")

# download zip file from Maestre et al. 2022b 
download.file(url = "https://figshare.com/ndownloader/files/37544623",
              destfile = here::here("data_raw"), 
              mode = "wb")
# manually: unzip!

# load all
maestre2015 <- readxl::read_xlsx(here::here("data_raw", "Data_Maestre_14_10_2015.xlsx"))
maestre2022a <- readxl::read_xlsx(here::here("data_raw", "Global_Biodesert_15122021_coordinates.xlsx"))
maestre2022b <- read_delim(here::here("data_raw", "Main_Data_code/Data_Source_Maestre_Science.txt"))

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Start checking newer data ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
coords1 <- maestre2022b %>% dplyr::select(LAT, LONG, ID) %>% distinct()
coords2 <- maestre2022a %>% dplyr::select(Lat_decimal, Long_decimal, ID) %>% 
  rename(LAT = Lat_decimal, LONG = Long_decimal) %>% distinct()

coords1 #326
coords2 #326

inner_join(coords1, coords2, by = c("LAT", "LONG")) #308
anti_join(coords1, coords2, by = c("LAT", "LONG")) #18
anti_join(coords2, coords1, by = c("LAT", "LONG")) #18 (exactly the same 18 sites)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## SUMMARY: BIODESERT data and Science data share the same sites, 
## but probably just some minor changes in digits.
## BIODESERT has the coordinates only, Science paper also soil functions and
## summarized diversity data (richness of selected groups).

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Check 2015 data ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
coords3 <- maestre2015 %>% dplyr::select(LAT, LON, NUM) %>% 
  rename(LONG = LON) %>% distinct()

coords3 #80

inner_join(coords1, coords3, by = c("LAT", "LONG")) #0
anti_join(coords1, coords3, by = c("LAT", "LONG")) #326
anti_join(coords3, coords1, by = c("LAT", "LONG")) #80

coords_both <- full_join(coords1, coords3)

# check with map
ggplot()+
  geom_point(data = coords1, 
             aes(x = LONG, y = LAT, fill = "BIODESERT"), alpha = 0.5, color = "grey")+
  geom_point(data = coords3, 
             aes(x = LONG, y = LAT, fill = "2015"), shape = 1, size = 5, alpha = 0.5)


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## SUMMARY: 2022 and 2015 data are complementary and do not share sites.
## Same for map and manual check using View(coords_both): not the same sites,
## but ~11+ km away from each other (0.1 digits and more).

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Compare with updated dataset shared internally ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
data_new <- readxl::read_xlsx(paste0(here::here(), "/data_raw/Global_Atlas_drylands_V2.xlsx"),
                                sheet = "Database")
data_new #551 

# remove sites without barcoding data
data_new <- data_new %>% filter(!is.na(ID_sequencing_16S))

coords0 <- data_new %>% dplyr::select(Latitude_c, Longitude_c , Plot_ID) %>% 
  rename(LAT = Latitude_c, LONG = Longitude_c) %>% distinct()
coords0 #338

inner_join(coords_both, coords0, by = c("LAT", "LONG")) #319
anti_join(coords_both, coords0, by = c("LAT", "LONG")) #87 present in coords_both but not coords0
data_unpublished <- anti_join(coords0, coords_both, by = c("LAT", "LONG"))
data_unpublished #19 present in coords0 but not coords_both

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Prepare raw data for sharing ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# load clean data
data_clean <- read_csv(here::here("intermediates", "Data_clean_global.csv"))
data_clean

## 1. remove rows that have embargo, i.e. sites are not yet shared elsewhere
# sampleIDs to remove
data_unpublished <- data_unpublished %>%
  inner_join(data_clean %>% dplyr::select(Latitude, Longitude, SampleID),
            by = c("LAT" = "Latitude", "LONG" = "Longitude")) 
data_unpublished #18 only... 1 (Plot_ID==4) for unknown reason missing

data_sharing <- data_clean %>% 
  filter(!(SampleID %in% unique(data_unpublished$SampleID)))
data_sharing #230 = 248-18

# Note: BIOCOM samples are not all included because in only 80 of their sites, 
#       barcoding was done. Therefore Plot_ID < 1000 missing

## 2. remove columns that have embargo, i.e. variables are not yet shared elsewhere
# check which variables were shared
colnames(maestre2015) #pH, C, Veg, diversity of microbes
colnames(maestre2022a) #just ID, country, site name & ID, grazing intensity and coordinates
colnames(maestre2022b) #sand, invertebrate richness, diverse soil properties
colnames(data_clean)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## SUMMARY: Data were shared before, even though we did not find matching values.
## Probably due to fixed or different calculations? Checked for DON only.

## 3. save for sharing
write_csv(data_sharing, here::here("data_raw", "Data_clean_global_SHARED.csv"))

# manually added METADATA file to data_sharing on 2025-03-26

#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## CONTINENTAL ------------------------------------ ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Continental data are available through request at ESDAC. Coordinates cannot
# be shared but are available on request to ESDAC or to authors.

# load clean data
data_clean <- read_csv(here::here("intermediates", "Data_clean_continental.csv"))
data_clean

# remove coordinates
data_sharing <- data_clean %>%
  mutate(Latitude = "available on request",
         Longitude = "available on request")

## save for sharing
write_csv(data_sharing, here::here("data_raw", "Data_clean_continental_SHARED.csv"))

# manually added METADATA file to data_sharing on 2025-03-26


#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
## REGIONAL --------------------------------------- ####
#- - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Regional data are not yet shared but just available on request as described
# in Duarte et al. 2024. Data used here is under embargo but can be shared 
# without coordinates, which are available on request.

# load clean data
data_clean <- read_csv(here::here("intermediates", "Data_clean_regional.csv"))
data_clean

# remove coordinates
data_sharing <- data_clean %>%
  mutate(Latitude = "available on request",
         Longitude = "available on request")

## save for sharing
write_csv(data_sharing, here::here("data_raw", "Data_clean_regional_SHARED.csv"))

# manually added METADATA file to data_sharing on 2025-03-26

