#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Parameter settings               #
#          author: Romy Zeiss               #
#            date: 2023-03-22               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

set.seed(1759) #& sample PA ID in pairing [last option used]
#set.seed(5738375)

## define functions to be compared
fns <- c("Soil_carbon_service", "OM_decomposition_service", #"Water_regulation_service", #water only for global
        "Soil_stability_service", "Nutrient_service", "Pathogen_control", 
        "Bac_richness", "Fungi_richness", "Invertebrate_richness", 
        "Protist_richness", "Nematode_richness", 
        "Bac_shannonDiv", "Fungi_shannonDiv", "Invertebrate_shannonDiv", "Protist_shannonDiv",
        "Ectomycorrhizal_richness", "Arbuscularmycorrhizal_richness", "Decomposer_richness",
        "Bac_JaccDist_av", "Fungi_JaccDist_av", "Protist_JaccDist_av",
        "Invertebrate_JaccDist_av")

fns_labels <- read.csv(paste0(here::here(), "/data_raw/LABELS_functions.csv"),
                       fileEncoding = "UTF-8-BOM")

## define variables to be compared between PA and nonPA, & threshold
mahal_vars <- c("Latitude", "Longitude", "Elevation", "AnnualPrec", "AnnualTemp", 
                #"MonthlyPrecipSum","MonthlyMeanTemp", 
                "Soil_pH", "Soil_salinity", "Soil_texture")

mahal_thres <- qchisq(.975, df=length(mahal_vars)) #21.92005

mahal_vars_z <- paste0(mahal_vars, ".z")

# define each land cover type
lc_names <- c("Dryland", "Cropland", "Grassland", "Shrubland", "Woodland", "Other")

# radius for continental site pairing
radius_thres <- 500000 #in m; 500000 for continental

# minimum number of nonPA sites per PA that can be paired
min_nonPA <- 10

# number of randomization = number of pairings 
number_times <- 1000

# define protected area ranking
protect_type <- data.frame("PA_type" = c("Ia", "Ib", "II", "III", "IV", "V", "VI",
                                         "Not Reported", "Not Assigned", "Not Applicable"),
                           "PA_rank" = 1:10,
                           "PA_protected" = c(rep(1,7), rep(0,3)))
protect_type

