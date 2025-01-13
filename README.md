# **Comparison of protected and unprotected soils across scales**
This repository contains the scripts, input and output for the comparison of soil diversity and functioning estimates in protected and unprotected sampling sites across three spatial scales (global, continental: Europe, regional: North-Portugal). 

This README.md file was written by Romy Zeiss and last updated on 2025-01-13. 
The structure of the README file is based on [this one](https://github.com/JWicquart/fish_growth/blob/master/README.md)

**The respective study is presented in the following article**:

[Zeiss et al (in Review) ... .](DOI)


## Table of Contents
- [1. How to download this project?](#1-how-to-download-this-project)
- [2. Description of the project](#2-description-of-the-project)
  - [2.1 Project organization](#21-project-organization)
  - [2.2 Datasets description](#22-datasets-description)
    - [2.2.1 Raw data](#221-raw-data)
    - [2.2.2 Intermediate results](#222-intermediate-results)
    - [2.2.3 Results](#223-results)
  - [2.3 Code description](#23-code-description)
- [3. How to report issues?](#3-how-to-report-issues)
- [4. Reproducibility parameters](#4-reproducibility-parameters)
- [5. Contributing](#5-contributing)
- [6. License](#6-license)
- [7. Use of ChatGPT](#7-use-of-chatgpt-in-project-development)


## 1. How to download this project?

On the project main page on GitHub, click on the green button `Code` and then click on `Download ZIP`.
Alternatively, you can `clone` the repository.


## 2. Description of the project

The aim of the analysis was the comparison of soil biodiversity and functioning attributes in protected and unprotected sites at three spatial scales and different habitat types.

Spatial scales are global, continental (Europe) and regional (north Portugal) and included as `temp_scale` or `scale` and with 4 letters (glob, cont, regi) or whole scale names in the scripts and output.
Habitat types are Dryland at global scale, and Cropland, Grassland, and Woodland at continental and regional scale. They are included as `LC` or `lc` (land cover types) individually or `lc_names` and `lc_names_all` as whole in the analysis and names by their full names.

At some point in the analysis, `temp_date` is used to track different versions of the intermediate datasets (mainly for randomized pairing and later comparisons based on these set of pairs).

### 2.1 Project organization

This project is divided in 6 folders:

* :open_file_folder:	`data_raw` folder contains input 11 datasets, 3 metadata files and one reference table (see part _2.2 Datasets description_).
* :open_file_folder:	`intermediates` folder serves as temporary storage place for intermediate results (see part _2.2 Datasets description_).
* :open_file_folder:	`results` folder contains output files (csv, RData, txt). It also contains intermediate results, figures and final results for the sensitivity analysis (see part _2.2 Datasets description_).
* :open_file_folder:	`src` folder contains X _src_ codes and a `functions` folder (see part _2.3 Code description_).
* :open_file_folder: `figures` folder contains figures and their underlying data (see part _2.3 Code description_).
* :open_file_folder: `docs` folder contains files to create Appendix as Quarto markdown report (*Appendix.qmd* and *Appendix.html*), including for example description of response variables. 


### 2.2 Datasets description

#### 2.2.1 Raw data

The spatial rasters **CHELSA_bio1_1981-2010_V.2.1.tif** and **CHELSA_bio12_1981-2010_V.2.1.tif** were downloaded from [CHELSA](https://chelsa-climate.org/downloads/) on 2022-04-10. 
Reference: Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017): Climatologies at high resolution for the Earth land surface areas. *Scientific Data*. 4 170122. https://doi.org/10.1038/sdata.2017.122

The datasets **Diversity_[temp_scale].csv** with `temp_scale` being one of three (global, continental, regional) correspond to the pre-processed sequencing data from the three original datasets described in the associated publication. The dataset **Diversity_global_atlas.csv** corresponds to the alternative global dataset that was analysed in a previous version of the study. This dataset covered very diverse habitats worldwide that are not necessarily comparable, which is why we decided to focus on the results from the dryland dataset only.   

- `SampleID` Unique code identifying each sample
- The remaining columns contain the richness (`richness`), Shannon index (`shannon`) and mean and standard deviation of the Jaccard Dissimilarity (`JaccDist_av` and `JaccDist_sd`) of the investigated soil taxa, namely bacteria (`Bac`), eukaryotes (`Euk`), fungi (`Fungi`), invertebrates (`Invertebrate`), protists (`Protist`), nematodes (`Nematode`), ecto- and arbuscular mycorrhizal fungi (`Ectomycorrhizal` and `Arbuscularmycorrhizal`), fungal plant pathogens (`Plant_pathogen`) and decomposers (`Decomposer`). The additional columns containing the richness of other functional groups were not considered in this analysis.

The spatial raster **Elevation_Worldclim_Global.tif** was downloaded from [WorldClim](https://www.worldclim.org/data/worldclim21.html) on 2022-12-21. The associated file **Elevation_Worldclim_Global_METADATA.txt** contains basic information to the raster file.
Reference: Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: new 1‐km spatial resolution climate surfaces for global land areas. *International journal of climatology*, 37(12), 4302-4315. https://doi.org/10.1002/joc.5086"

The dataset **Global_Atlas_drylands_V2.xlsx** corresponds to the global drylands dataset used for analysis.
References:  
- Maestre FT, Delgado-Baquerizo M, Jeffries TC, Eldridge DJ, Ochoa V, Gozalo B, et al. (2015) Increasing aridity reduces soil microbial diversity and abundance in global drylands. Proceedings of the National Academy of Sciences 112: 15684–15689.
- Maestre FT, Le Bagousse-Pinguet Y, Delgado-Baquerizo M, Eldridge DJ, Saiz H, Berdugo M, et al. (2022) Grazing and ecosystem service delivery in global drylands. Science 378: 915–920.

The dataset **GlobalAtlasv2_conservation_heterogeneity_papers_v2.xlsx** corresponds to the alternative global dataset used in a previous version of the study. 
References: 
- Delgado-Baquerizo M, Bardgett RD, Vitousek PM, Maestre FT, Williams MA, Eldridge DJ, et al. (2019) Changes in belowground biodiversity during ecosystem development. Proceedings of the National Academy of Sciences 116: 6891–6896.
- Delgado-Baquerizo M, Reich PB, Bardgett RD, Eldridge DJ, Lambers H, Wardle DA, et al. (2020a) The influence of soil age on ecosystem structure and function across biomes. Nature Communications 11: 4721.
- Guerra CA, Berdugo M, Eldridge DJ, Eisenhauer N, Singh BK, Cui H, et al. (2022) Global hotspots for soil nature conservation. Nature: 1–6.
- Maestre FT, Le Bagousse-Pinguet Y, Delgado-Baquerizo M, Eldridge DJ, Saiz H, Berdugo M, et al. (2022) Grazing and ecosystem service delivery in global drylands. Science 378: 915–920.

The table **LABELS_functions.csv** serves as reference table to link raw columns names of response variables (`Function`) to shortnames and labels used in figures.

- `Function` Column names of response variables in raw data
- `Label` Corresponding labels for figures 
- `Label_short` Short labels for figures
- `Group_function` Group of response variable defined by the authors of the study. Response variables were grouped into Function, Richness, Shannon (index) and Dissimilarity measures.
- `Organism` Group of response variables based on taxonomic classification defined by the authors of the study. Response variables were grouped into Function, Bacteria, Fungi, Invertebrates, Protists, Nematodes, and Decomposers.

The dataset **LUCAS_2018_iDiv_20221018.csv** corresponds to the continental dataset used for analysis. This data was collected in the European Land Cover/Area frame Survey (LUCAS) soil module in 2018. The associated file **LUCAS_2018_iDiv_metadata.txt** contains information about the columns and sampling/analysis methods.
References: 
- European Commission. Joint Research Centre. (2022) LUCAS 2018 soil module: presentation of dataset and results. LU: Publications Office.
- Hiederer R (2018) Data Evaluation of LUCAS Soil Component Laboratory Data for Soil Organic Carbon. JRC Technical Reports: 73.
- Orgiazzi A, Ballabio C, Panagos P, Jones A, Fernández‐Ugalde O (2018) LUCAS Soil, the largest expandable soil dataset for Europe: a review. European Journal of Soil Science 69: 140–153.
- Orgiazzi A, Panagos P, Fernández-Ugalde O, Wojda P, Labouyrie M, Ballabio C, et al. (2022) LUCAS Soil Biodiversity and LUCAS Soil Pesticides, new tools for research and policy development. European Journal of Soil Science 73: e13299.

The dataset **SoilReCon_Data_211123.csv** corresponds to the regional dataset used for analysis. This data was collected in the project "Soil Ecosystems in the XXI Century: pressures, conservation and future scenarios" in the North of Portugal. The associated file **SoilReCon_MetaData_211123.csv** contains information about the columns and sampling/analysis methods. 
Project number: PTDC/BIA-CBI/2340/2020 (FCT-Portugal)


#### 2.2.2 Intermediate results

For each `temp_scale`, there is one folder containing the pairs of protected and unprotected sites (`Pairs_paNonpa...csv`) and protected sites that could not been paired with environmentally similar unprotected sites(`Unpaired_protected_sites...csv`). Latter ones should be empty, as sites that cannot be paired are excluded now beforehand. Creating and saving this output is an artefact of previous code versions. The `Pairs_paNonpa...csv` instead is used to compare protected and unprotected sites by calculating effect sizes (d-values). 
In addition, the scale-specific folders contain `Geographically_close_sites.RData` for at least continental scale, which is an R list object of unprotected sites with distances smaller or equal to the defined distance threshold for pairing.

The `intermediate` folder also contains the clean datasets (`Data_clean...csv`), the sites actually used in randomized pairing for calculating effect sizes (`Data_paired...csv`, subset of `Data_clean....csv` but with scaled variables), the assignment of protection status to the sampling locations (`PA_assignment...csv`), and intermediate results from random-slope (`PAranks_Bayesian...RData`) and random-intercept models (`pars_PAtypes_Bayesian...RData` from rstan models, and `PAtypes_Bayesian...RData` from brms models).


#### 2.2.3 Results

The `results` folder contains output files of the analyses (csv, RData, txt). It also contains intermediate results, figures and final results for the sensitivity analysis, namely the analysis of the alternative dataset (`sensitivity_globalAlternativeDataset` [outdated version]) and the analysis of all global bacterial sites (`sensitivity_globalBacteria`). The naming conventions within the `sensitivity_` folders follow those from the `intermediate`, `results` and `figures` folders.

Naming conventions in the `results` folder:
- `corMat` Correlation matrix calculated as *Pearson* or *Spearman* rank correlation coefficient of the environmental variables considered for environmental similarity.
- `d_1000_trails` Effect sizes from the 1000 randomizations.
- `D-values` Significant effect sizes summarized across the 1000 randomited pairings and sorted by absolute median effect size (`effect_median`).
- `Locations` Sampling locations for effect size analysis (no additional naming) and random-intercept and -slope models (`PAranks`). Sampling locations for effect size analysis are a subset of the ones from random-intercept and -slope models.
- `Numbers` Number of sampling sites (not) used in pairing and per habitat type.
- `PAranks_Bayesian` Results of random-slope models, either as estimated marginal means (`emtrends`) or as subset of 10,000 randomly samples slope estimates (`sample10k`).
- `pars_PAtypes_Bayesian_df` Results of the random-intercept models
- `VIF_envVars` Variable inflation factor of the environmental variables calculated using the *usdm* package.


### 2.3 Code description

#### 00_Functions.R
Script to define functions for preparing and analysing given data.

##### Prepare data
1. **f_extract_pa:**
   This function extracts information about the protection status of locations from a given dataset using a shapefile containing polygons of protected areas. The input parameters include the dataset (`data`), column names for longitude (`col_lon`), latitude (`col_lat`), and sample ID (`col_id`). It also requires a shapefile (`shp`) with protection area information and the column name (`col_pa`) in the shapefile indicating the protection status. The function returns a modified dataset with added information about whether each location is within a protected area.

2. **f_extract_env:**
   This function extracts environmental data such as elevation, annual temperature, and annual precipitation for given locations from raster datasets. The input parameters are similar to the previous function, including the dataset (`data`), and column names for longitude (`col_lon`) and latitude (`col_lat`). The function uses raster data for elevation and climate variables, extracts the corresponding values for each location, and adds them to the dataset.

3. **f_add_protect:**
   This function adds information about the protection level of locations based on a provided dataset and information about protection status (`data_pa`). It includes double-entry removal for cases with multiple protected area types for a single protected area, selecting the highest protection level, and adding this information to the original dataset (`data`). The string `col_id` is used to acces the ID column in both `data` and `data_pa`.

4. **f_colinearity:**
   The purpose of this function is to check for collinearity between environmental variables in a given dataset (`data`). It calculates the Variable Inflation Factor (VIF) and correlation matrices (Spearman and Pearson) to assess the correlation and potential multicollinearity between variables. The columns to calculate correlation are given by the input parameters `col_lon` (name of longitude column), `col_lat` (name of latitude column), and `vars_env` (column names of environmental variables).

5. **f_scale_vars:**
   This function scales specified variables within a dataset. It takes the dataset (`data`) and a vector of column names for variables to be scaled (`vars`). It checks if the specified variables are present in the dataset and scales each variable individually. Scaling is done using the base R function scale().

##### Pairing
6. **f_check_pairs:**
   This function checks the pairing of protected and non-protected sites for a given dataset (`data`) based on Mahalanobis distance. It assesses the Mahalanobis distance values for each non-protected site compared to protected sites and identifies unpaired protected sites. It also provides information about the number of non-protected sites per protected area. The string `col_id` is used to acces the ID column in `data`. The parameters `col_lc` and `vars_z` refer to the column names of the habitat type (or land cover, lc) and scaled variables selected to calculate the Mahalanobis distance.

7. **f_resample:**
   This is a utility function used in the pairing process to resample values from a vector. It is employed within the `f_pairing` function.

8. **f_pairing:**
   This function is designed to pair protected and non-protected sites based on Mahalanobis distance. `f_pairing` takes the same input parameters as `f_check_pairs`. It iteratively pairs sites for each land cover type, considering environmental variables. The process involves randomly ordering protected area (PA) sites using the function `f_resample`, calculating Mahalanobis distances to non-protected sites, and selecting pairs based on distance thresholds. If not enough pairs are found in a run, the function repeats the process to ensure an adequate number of paired sites. It returns information about unpaired protected sites, the number of non-protected sites per protected area, and the Mahalanobis distance values for analysis. For `temp_scale` == "continental", unprotected sites with a geographical distance > `radius_thres` value are removed. Geographical distances have to be given in a list called `list_dist`.

##### Comparison
9. **f_compare_pa_nonpa:**
   This function compares the effect sizes (Cohen's d) between protected (PA) and non-protected sites for each land cover type and response variable of a given dataset (`data`). It performs individual tests for each combination of land cover type and response variable, considering multiple runs of pairings. The function utilizes the `psych::cohen.d` function to calculate Cohen's d, and the results are stored in a list for further analysis. This function is useful for understanding the impact of protection on different environmental variables across various land cover types. The input parameter `data_pairs` is a list of a length equal to the number of times with pairings (`number_times`). The input parameters `col_id` and `col_fns` provide the name of the ID and response variable columns, respectively.

10. **f_compare_pa_types:**
    This function conducts Bayesian analysis to compare different protection levels, including non-protected sites, within each land cover type for a given dataset (`data`). It fits a Bayesian model for each response variable and land cover type combination using the `rstan` package. The function utilizes Stan models to estimate parameters such as group means and variances, providing a Bayesian perspective on the differences between protection levels. The results are stored in a list for further examination. The input parameter `protection_levels` is a dataframe containing the translation of protection type levels into ranks. `col_fns` provides the names of the response variable columns.

11. **f_compare_pa_dummy:**
    This function conducts Bayesian analysis to compare protected and unprotected sites, within each land cover type for a given dataset (`data`). It is based on the function `f_compare_pa_types` but only considered 2 groups, namely "protected" and "unprotected" instead of the protected area type categories. 

12. **f_combine_pars_list:**
    This function combines the output of Bayesian analyses (`data`) into a single dataframe for easier interpretation and analysis. It takes the list of Bayesian parameter samples for each land cover type and response variable, organizing them into a tidy dataframe. This consolidated dataframe can then be used for further statistical analysis or visualization of the results obtained from the Bayesian comparisons. The input parameter `col_fns` provides the names of the response variable columns.

#### 00_Parameters.R
Script to define parameters used throughout the analysis, such as the response variables (`fns`) and its labels for figures (`fns_labels`), environmental variables for similarity (`mahal_vars` and scaled `mahal_vars_z`) and the Mahalobis distance threshold (`mahal_thres`), individually considered land cover types (`lc_names`), and the protected area type ranking (`protect_type`). It also sets user-defined values, namely `radius_thres` as maximum geographical distance for continental site pairing, `min_nonPA` as minimum number of unprotected sites per protected site that can be paired, or `number_times` as number of pairings during randomization process. Note: some parameters are defined in each script, e.g. `temp_scale`.

#### 01_Prepare_diversity.R
Script to prepare soil diversity data: load pre-processed diversity data, calculate missing variables (e.g. mean Jaccard distance for each sample ID), merge dissimilarity and diversity measures.

#### 01a_Prepare_protection.R
Script to extract the protection status of the investigated sampling sites by merging shapefiles downloaded from protectedplanet.net and extracting protected area types for each coordinate for all three spatial scales.

#### 01b_Prepare_soil_data.R
Script to merge diversity data prepared in script `01_Prepare_diversity.R`, soil functioning data, and environmental data for each spatial scale separately. Additional variables are calculated or renamed if necessary, and environmental data is extracted using the function `f_extract_env`. This script also contains calculation of colinearity using the function `f_colinearity`. 

#### 02_Pairing.R
Script to build pairs of protected and environmentally similar unprotected sites for each spatial scale (`temp_scale`). At the beginning, `temp_date` has to be set to the date when this script is run to guarantee access of correct intermediate resuls file. For `temp_scale` == "continental", unprotected sites with a geographical distance > `radius_thres` value will be removed, which is why geographical distances are calculated and saved in `list_dist`.

#### 03a_Compare_PA_nonPA.R
Script to compare protected (`PA`) and paired unprotected (`nonPA`) sampling sites for each spatial scale (`temp_scale`) using the function `f_compare_pa_nonpa`. At the beginning, `temp_date` has to be set to the date when this script is run to guarantee access of correct intermediate resuls file. The output is saves in the `results` folder in `d_1000_trails_[temp_scale].RData`.

#### 03b_Compare_PA_types_Bayesian.R
Script to compare sampling sites in different protected area types for each spatial scale (`temp_scale`) using the function `f_compare_pa_types`. The first part corresponds to the random-intercept model using *rstan*, the second part to the random-slope model using the R packages *brms*, *modelr*, and *tidybayes*. The *emmeans* R package is used to calculate estimated marginal means of the random-slope model. At the beginning and in the second part of the script, `temp_date` has to be set to the date when this script is run to guarantee access of correct intermediate resuls file. The output for the random-intercept model is saved in the `results` folder in `pars_PAtypes_Bayesian_df_[temp_scale].RData` for each `temp_scale`. An intermediate result named `pars_PAtypes_Bayesian_[temp_scale].RData` is created for each `temp_scale`. The output for the random-slope model is saved in the `results` folder in `PAranks_Bayesian_[temp_scale]_sample10k.RData` and `PAranks_Bayesian_[temp_scale]_emtrends.csv` for each `temp_scale`. An intermediate result named `PAranks_Bayesian_[temp_scale].RData` and `PAranks_Bayesian_[temp_scale]_summary.RData` is created for each `temp_scale`.

#### 04_Plotting.R
Script to create figures saved in the folder `figures` and the underlying data behind the figures. The script also creates summarizing statistics that are saved in the `figures` folder, such as the whole data for effect sizes (`Data_d-value_[temp_scale]_.csv`) and a summarized table that is used to create **Figure 2** and **Appendix Figure 2.3** of the publication (i.e. mean and confidence intervals (CI), table `Results_d-value_summary_[temp_scale].csv` corresponding to **Appendix Table 2.3**). In addition, TXT files are saved in the `results` folder summarizing numbers and significant results, namely `D-values_[temp_scale].txt` and `Numbers_d-values_[temp_scale].txt`. Note: Some figures are created for each `temp_scale` separately, while others are created in a `for` loop across all three spatial scales overwriting the R object `temp_scale`.

Naming conventions in the `figures` folder:
- `icon_[...]` Icons downloaded from Flaticon.com.
- `Data_[...]` Raw data and supporting data.
- `Results_[...]` Figures and underlying data of the results.
- `[...]_allScales` Combination of results for all three spatial scales.
- `[...]_grouped` Results for grouped response variables, namely Function, Richness, Shannon (index) and Dissimilarity measures.

Specific file naming conventions in the `figures` folder:
- `Correlation_` Correlation between Mahalanobis distances and differences in protected and unprotected sampling sites. If not indicated differently, correlation was calculated as Spearman's rank correlation coefficient.
- `Data_boxplot_mahal.distance_` Mahalanobis distance of pairs and of all sites visualized as boxplots.
- `Data_clean_PAranks_fns` Raw data visualized as scatter plots with regression lines.
- `Data_d-value` Effect sizes of all randomizations with 95 % confidence intervals.
- `Data_locations` Maps of the sampling locations for effect size analysis and random-intercept and -slope models (`PAranks`).
- `icon_` Icons downloaded from Flaticon.com: “Grass”, “Globe” and “Flag” created by Freepik, “Europe” (`icon_location-black.png`) modified after icon by amoghdesign, “Forest” by Vitaly Gorbachev, “Bush” (`icon_shrub-crop.png`) by Buandesign, and “Agriculture” (`icon_harvest.png`) by Vectors Tank.
- `Locations_connection` Network of sampling sites with lines connecting the pairs of protected and unprotected sites. At continental scale, different thresholds were tested (500, 600, 700 km, no threshold).
- `Results_boxplot_estimates` Raw data visualized in boxplots comparing protected and unprotected sites.
- `Results_d-values` Effect sizes visualized as heatmap (means and medians, `meanCI` and `medianSD`), and pointrange plots for individual and grouped response variables (`medianCI` and `medianCI_grouped`).
- `Results_intercept_parsBayesian` Results of random-intercept models with number of sampling sites documented in `_nTable` CSV files.
- `Results_slope_BayesianTrends` Results of random-slope models.


### 2.4 How to reproduce the final datasets?

To reproduce the results and figures, open the **Conservation_across_scales.Rproj** file, and/or open and run successively the R codes _01_Prepare_diversity.R_ to _04_Plotting.R_. The working directory should be in a folder that contains all 5 subfolders (i.e. `data_raw`, `figures`, `intermediates`, `results`, `src`). Make sure that all required packages were previously installed. 

For each `temp_scale`, that is global, continental and regional, parts of the code have to be run repeatately by changing or out-commenting the assignment of the `temp_scale` R object (usually at the beginning of the scripts).

Due to the elevated number of iterations from the Bayesian model, the codes necessitate a certain amount of time to run. Please note that a full reproduction of the results is not possible, as Bayesian estimations are different at each run.


## 3. How to report issues?

Please report any bugs or issues [HERE](https://github.com/JeMaNd-r/Conservation_across_scales/issues).


## 4. Reproducibility parameters

```R
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows Server 2019 x64 (build 17763)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

time zone: Europe/Berlin
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] emmeans_1.8.9       rstan_2.32.3        StanHeaders_2.26.28 magick_2.8.1        corrplot_0.92       ggh4x_0.2.8         ggtext_0.1.2        ggrepel_0.9.3      
 [9] gridExtra_2.3       ggdist_3.3.0        terra_1.7-29        lubridate_1.9.2     forcats_1.0.0       stringr_1.5.0       dplyr_1.1.2         purrr_1.0.1        
[17] readr_2.1.4         tidyr_1.3.0         tibble_3.2.1        ggplot2_3.4.2       tidyverse_2.0.0     here_1.0.1         

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3   tensorA_0.36.2.1     rstudioapi_0.15.0    jsonlite_1.8.5       magrittr_2.0.3       estimability_1.4.1   farver_2.1.1         ragg_1.2.5          
  [9] vctrs_0.6.3          base64enc_0.1-3      htmltools_0.5.5      distributional_0.3.2 broom_1.0.5          tidybayes_3.0.6      raster_3.6-20        cellranger_1.1.0    
 [17] htmlwidgets_1.6.2    plyr_1.8.8           zoo_1.8-12           commonmark_1.9.0     igraph_1.5.0         mime_0.12            lifecycle_1.0.4      pkgconfig_2.0.3     
 [25] colourpicker_1.3.0   Matrix_1.5-4         R6_2.5.1             fastmap_1.1.1        shiny_1.7.4          digest_0.6.31        colorspace_2.1-0     ps_1.7.5            
 [33] rprojroot_2.0.3      brms_2.20.4          textshaping_0.3.6    crosstalk_1.2.1      labeling_0.4.3       fansi_1.0.4          timechange_0.2.0     mgcv_1.8-42         
 [41] abind_1.4-5          compiler_4.3.0       bit64_4.0.5          withr_2.5.2          backports_1.4.1      inline_0.3.19        shinystan_2.6.0      psych_2.3.6         
 [49] QuickJSR_1.0.7       pkgbuild_1.4.1       maps_3.4.1           gtools_3.9.4         loo_2.6.0            tools_4.3.0          httpuv_1.6.11        threejs_0.3.3       
 [57] glue_1.6.2           callr_3.7.3          nlme_3.1-162         promises_1.2.0.1     gridtext_0.1.5       grid_4.3.0           checkmate_2.3.0      reshape2_1.4.4      
 [65] generics_0.1.3       gtable_0.3.4         tzdb_0.4.0           hms_1.1.3            sp_1.6-1             xml2_1.3.4           utf8_1.2.3           pillar_1.9.0        
 [73] markdown_1.11        vroom_1.6.3          posterior_1.5.0      later_1.3.1          splines_4.3.0        lattice_0.21-8       bit_4.0.5            tidyselect_1.2.0    
 [81] miniUI_0.1.1.1       arrayhelpers_1.1-0   xfun_0.41            stats4_4.3.0         bridgesampling_1.1-2 matrixStats_1.0.0    DT_0.31              stringi_1.7.12      
 [89] codetools_0.2-19     usdm_1.1-18          cli_3.6.1            RcppParallel_5.1.7   shinythemes_1.2.0    xtable_1.8-4         systemfonts_1.0.4    munsell_0.5.0       
 [97] processx_3.8.1       modelr_0.1.11        Rcpp_1.0.10          readxl_1.4.2         mapproj_1.2.11       png_0.1-8            coda_0.19-4          svUnit_1.0.6        
[105] parallel_4.3.0       rstantools_2.3.1.1   ellipsis_0.3.2       prettyunits_1.1.1    dygraphs_1.1.1.6     bayesplot_1.10.0     Brobdingnag_1.2-9    mvtnorm_1.2-2       
[113] scales_1.2.1         xts_0.13.2           crayon_1.5.2         rlang_1.1.1          mnormt_2.1.1         shinyjs_2.1.0
```

## 5. Contributing

Contributors names and contact info:

Romy Zeiss  
[@RomyZeiss](https://twitter.com/romyzeiss)

Calderón-Sanou
Concha Cano-Díaz
Manuel Delgado-Baquerizo
Paulo Fernandes
Fernando T. Maestre
Susana Mendes
Bala Singavarapu
Brajesh K. Singh
Carlos A. Guerra


## 6. License

This project is licensed under the MIT License - see the LICENSE.md file for details.


## 7. Use of ChatGPT in Project Development

This project was partially developed with the assistance of OpenAI's ChatGPT, a large language model designed to support developers in coding, debugging, and explaining complex concepts. ChatGPT was used to help with tasks such as generating code snippets, providing algorithmic insights, refining documentation, and ensuring proper project organization. While the model provided suggestions, all code has been thoroughly reviewed, tested, and modified by the project maintainers to ensure accuracy, efficiency, and relevance to the project’s objectives.

