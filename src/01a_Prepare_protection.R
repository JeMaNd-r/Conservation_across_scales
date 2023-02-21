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
terraOptions(tempdir =  "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - - -
## load all files provided by protectedplanet.net
shp0 <- terra::vect("D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_0/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp1 <- terra::vect("D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_1/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp2 <- terra::vect("D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_2/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")

#- - - - - - - - - - - - - - - - - - - - - -
## merge all files
shp_all <- rbind(shp0, shp1)
rm(shp0, shp1)

shp_all <- rbind(shp_all, shp2)
rm(shp2)

shp_all
#plot(shp_all)

#terra::writeVector(shp_all, "D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged.shp")
save(shp_all, file="D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged.RData")
gc()

load("D:/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_merged.RData") #shp_all

#- - - - - - - - - - - - - - - - - - - - - -
## Extract protection information based on point coordinates
#- - - - - - - - - - - - - - - - - - - - - -
# load points
data_glob <- read_csv(paste0(here::here(), "/data_raw/GlobalAtlasv2_conservation_heterogeneity_papers_v1.csv"))
data_glob

data_vect <- terra::vect(data_glob[,c("Longitude_c", "Latitude_c")], geom=c("Longitude_c", "Latitude_c"))
#plot(data_vect)

# extract protection status
data_pa <- terra::extract(shp_all, data_vect, xy=TRUE)
data_pa

data_pa$PA <- 0
data_pa[!is.na(data_pa$WDPAID), "PA"] <- 1

data_pa
table(data_pa$PA)

write_csv(data_pa, paste0(here::here(), "/intermediates/PA_assignment_global.csv"))
