#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Parameters & functions             #
#          author: Romy Zeiss               #
#            date: 2023-02-22               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## re-define sample function (if only 1 value in sample(), it considers it as vector)
resample <- function(x, ...) x[sample.int(length(x), ...)]


## define functions to be compared
fns <- c("Soil_carbon_service", "OM_decomposition_service", "Water_regulation_service", 
        "Soil_stability_service", "Nutrient_service", "Pathogen_control", 
        "Richness_bacteria", "Fungi_18s_richness", "Invertebrate_18s_richness", 
        "Protist_18s_richness", "Ectomycorrhizal_18s_richness", 
        "Arbuscular_mycorrhizal_18s_richness", "Decomposers_18s_richness", 
        "Diss_Bacteria_std", "Diss_Fungi_std", "Diss_Protists_std", "Diss_invert_std")

fns_labels <- read.csv(paste0(here::here(), "/data_raw/LABELS_functions.csv"),
                       fileEncoding = "UTF-8-BOM")

## define variables to be compared between PA and nonPA, & threshold
mahal_vars <- c("Latitude_c", "Longitude_c", "Elevation", "AnnualPrec", "AnnualTemp", 
                #"MonthlyPrecipSum","MonthlyMeanTemp", 
                "Soil_pH", "Soil_salinity", "Clay_silt_c")

mahal_thres <- qchisq(.975, df=length(mahal_vars)) #21.92005

mahal_vars_z <- paste0(mahal_vars, ".z")

# define each land cover type
lc_names <- c("Grassland", "Shrubland", "Woodland")

# number of samples/ sites that should be paired per LC type
min_size <- 20