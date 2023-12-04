#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#     Prepare protected area network        #
#          author: Romy Zeiss               #
#            date: 2022-11-04               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(here)
library(tidyverse)
library(terra)

# change temporary directory
terraOptions(tempdir =  "D:/EIE_Macroecology/00_datasets/Trash")

source(paste0(here::here(), "/src/00_Parameters.R"))
source(paste0(here::here(), "/src/00_Functions.R"))

#- - - - - - - - - - - - - - - - - - - - - -
## load all files provided by protectedplanet.net  ####
shp0 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_0/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp1 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_1/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp2 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_2/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")

#- - - - - - - - - - - - - - - - - - - - - -
## merge all files
shp_all <- rbind(shp0, shp1)
rm(shp0, shp1)

shp_all <- rbind(shp_all, shp2)
rm(shp2)

shp_all
#plot(shp_all)

#terra::writeVector(shp_all, "D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged.shp") #too big
save(shp_all, file="D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged_new.RData")
gc()

load("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged_new.RData") #shp_all
#shp_all <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged.shp") #= merged_new.RData but not all saved due to size

#- - - - - - - - - - - - - - - - - - - - - -
## Extract protection information based on point coordinates ####
#- - - - - - - - - - - - - - - - - - - - - -

### global ####
# load points
data_glob <- readxl::read_xlsx(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v2.xlsx"),
                               sheet = "Data")
data_glob

data_pa <- f_extract_pa(data = data_glob, 
             col_lon = "Longitude_c", col_lat = "Latitude_c",
             col_id = "Order_ID",
             shp = shp_all, col_pa = "WDPAID")

table(data_pa$PA)

write_csv(data_pa, paste0(here::here(), "/intermediates/PA_assignment_global.csv"))

#- - - - - - - - - - - - - - - - - - - - - -
### continental ####
# load points
data_cont <- read_csv(paste0(here::here(), "/data_raw/LUCAS_2018_iDiv_20221018.csv"))
data_cont

data_pa <- f_extract_pa(data = data_cont,
                        col_lon = "Longitude", col_lat = "Latitude",
                        col_id = "BARCODE_ID",
                        shp = shp_all, col_pa = "WDPAID")

table(data_pa$PA)

data_pa <- data_pa %>% rename(SampleID = BARCODE_ID)

write_csv(data_pa, paste0(here::here(), "/intermediates/PA_assignment_continental.csv"))

#- - - - - - - - - - - - - - - - - - - - - -
### regional ####
# load points
data_regi <- read_csv(paste0(here::here(), "/data_raw/SoilReCon_Data_211123.csv"))
data_regi

# transform coordinates
data_regi

data_pa <- f_extract_pa(data = data_regi, 
                        col_lon = "Longitude", col_lat = "Latitude",
                        col_id = "SampleID",
                        shp = shp_all, col_pa = "WDPAID")

table(data_pa$PA)

write_csv(data_pa, paste0(here::here(), "/intermediates/PA_assignment_regional.csv"))
