# Conservation_across_scales
 Compare soil biodiversity in- and outside of protected areas across scales (global, Europe, country)

Bayesian script: 250 runs take >18h

### Functions.R

1. **f_extract_pa:**
   This function extracts information about the protection status of locations from a given dataset using a shapefile containing polygons of protected areas. The input parameters include the dataset (`data`), column names for longitude (`col_lon`), latitude (`col_lat`), and sample ID (`col_id`). It also requires a shapefile (`shp`) with protection area information and the column name (`col_pa`) in the shapefile indicating the protection status. The function returns a modified dataset with added information about whether each location is within a protected area.

2. **f_extract_env:**
   This function extracts environmental data such as elevation, annual temperature, and annual precipitation for given locations from raster datasets. The input parameters are similar to the previous function, including the dataset (`data`), and column names for longitude (`col_lon`) and latitude (`col_lat`). The function uses raster data for elevation and climate variables, extracts the corresponding values for each location, and adds them to the dataset.

3. **f_add_protect:**
   This function adds information about the protection level of locations based on a provided dataset and information about protection status (`data_pa`). It includes double-entry removal for cases with multiple protected area types for a single protected area, selecting the highest protection level, and adding this information to the original dataset.

4. **f_colinearity:**
   The purpose of this function is to check for collinearity between environmental variables in a given dataset. It calculates the Variable Inflation Factor (VIF) and correlation matrices (Spearman and Pearson) to assess the correlation and potential multicollinearity between variables.

5. **f_scale_vars:**
   This function scales specified variables within a dataset. It takes the dataset (`data`) and a vector of column names for variables to be scaled (`vars`). It checks if the specified variables are present in the dataset and scales each variable individually.

6. **f_check_pairs:**
   This function checks the pairing of protected and non-protected sites based on Mahalanobis distance. It assesses the Mahalanobis distance values for each non-protected site compared to protected sites and identifies unpaired protected sites. It also provides information about the number of non-protected sites per protected area.

7. **f_resample:**
   This is a utility function used in the pairing process to resample values from a vector. It is employed within the `f_pairing` function.

8. **f_pairing:**
   This function is designed to pair protected and non-protected sites based on Mahalanobis distance. It iteratively pairs sites for each land cover type, considering environmental variables. The process involves randomly ordering protected area (PA) sites, calculating Mahalanobis distances to non-protected sites, and selecting pairs based on distance thresholds. If not enough pairs are found in a run, the function repeats the process to ensure an adequate number of paired sites. It returns information about unpaired protected sites, the number of non-protected sites per protected area, and the Mahalanobis distance values for analysis.

9. **f_compare_pa_nonpa:**
   This function compares the effect sizes (Cohen's d) between protected (PA) and non-protected sites for each land cover type and response variable. It performs individual tests for each combination of land cover type and response variable, considering multiple runs of pairings. The function utilizes the `psych::cohen.d` function to calculate Cohen's d, and the results are stored in a list for further analysis. This function is useful for understanding the impact of protection on different environmental variables across various land cover types.

10. **f_compare_pa_types:**
    This function conducts Bayesian analysis to compare different protection levels, including non-protected sites, within each land cover type. It fits a Bayesian model for each response variable and land cover type combination using the `rstan` package. The function utilizes Stan models to estimate parameters such as group means and variances, providing a Bayesian perspective on the differences between protection levels. The results are stored in a list for further examination.

11. **f_combine_pars_list:**
    This function combines the output of Bayesian analyses into a single dataframe for easier interpretation and analysis. It takes the list of Bayesian parameter samples for each land cover type and response variable, organizing them into a tidy dataframe. This consolidated dataframe can then be used for further statistical analysis or visualization of the results obtained from the Bayesian comparisons.

