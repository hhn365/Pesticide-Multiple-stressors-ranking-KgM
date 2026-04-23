## Step 3. Analyze the effects of pesticide mixture versus multiple stressors on ecological responses
#  - Import and check raw data quality of each stressor
#  - Filter all stressors restricted to date prior or on the same date of macroinvertebrate sampling
#  - Test models with selected four pesticide metrics from Step 2 script
#  - Use correlation (> 0.7) and VIF (> 2) to reduce redundant predictors in multiple stressor models
#  - Check linearity assumption
#  - Use stepAIC backward selection to define optimum model from the initial full model with all stressors
#  - Run model validation 

# -------------------------

# -------------------------
## Import packages
library(readr)    
library(openxlsx) 
library(here)      
library(readxl)
library(writexl)
library(tidyverse)
library(dplyr)
library(car)       
library(MASS)      
library(relaimpo)
library(tibble)
library(broom)
library(purrr)
library(lubridate) 
library(ggplot2)
library(corrplot)
library(GGally)
library(patchwork)

# -------------------------

# -------------------------
## Import abiotic and biotic metrics published in the main manuscript
# Included variables: study site IDs, macroinvertebrate metrics, pesticide metrics, and other available abiotic variables in the main manuscript Liess et al. (2021)
Metrics_in_the_manuscript <- read_excel("Input_data/1-s2.0-S0043135421004607-mmc2.xlsx", sheet = "Site Parameters")
names(Metrics_in_the_manuscript)

# -------------------------
## All available abiotic parameters:
Abiotic_metrics <- Metrics_in_the_manuscript |>  
  dplyr::select(1, 2, 4,5, 9:21) |>                # Including WWTP, pH, and Morphology. Other abiotic variables are specified in the rename())
  dplyr::rename(ID = "Site ID",
         year = Year,
         AgriLand_percent = "Agricultural land cover [%]",
         NO2_origin = "NO2 [mg/L]",
         NO3_origin = "NO3 [mg/L]",
         NH4_origin = "NH4 [mg/L]",
         PO4_origin = "PO4 [mg/L]",
         TN_origin = "TN [mg/L]",
         TP_origin = "TP [mg/L]",
         TUmetal_origin = "TUmetal [log]",
         TMP_origin = "Temperature [°C]",
         Flow_origin = "Flow velocity [m/s]",
         DO_origin = "Dissolved O2 [mg/L]",
         pH_origin = pH,
         Habitat = "Bed habitat structure")  # Did not add '*_origin' for variables that had single value per site across years 
                                             # and no data aggregation and comparison to data aggregation choice in the original manuscript required

names(Abiotic_metrics)

# Add stream width and depth provided by authors of the original manuscript
multiple <- read.delim("Input_data/multipleRegression_data.txt", header = TRUE, sep = "\t")

# Joined abiotic dataset
Abiotic_metrics <- Abiotic_metrics |> 
  left_join(multiple |> 
              dplyr::select(siteID, year, Stream_width  = stream_width_m, Stream_depth  = water_depth_cm),
    by = c("ID" = "siteID", "year" = "year"))

names(Abiotic_metrics)

# -------------------------
## Macroinvertebrate metrics: Extract three biotic metrics selected from the main manuscript Liess et al. 2021 (i.e., those with highest explained variance with multiple stressors)
Invertebrate_metrics <- Metrics_in_the_manuscript |>  
  dplyr::select(1, 2, 22:24) |>                 
  rename(ID = "Site ID", 
         year = Year, 
         speartaxa = SPEARpesticides,
         eptProz= "EPT%",
         sapr = "Saprobic index")

# -------------------------

# -------------------------
## Study sites: Extract 101 IDs to filter relevant sites for macroinvertebrate and multiple stressor analyses
unique_IDs <- unique(Invertebrate_metrics$ID)
unique_IDs

# -------------------------
## Macroinvertebrate sampling dates: Extraction from raw invertebrate dataset 
Invertebrates <- read.table("Input_data/5_1_Invertebrate_Abundances.txt", header = TRUE, sep = "\t")
Invertebrates_date <- Invertebrates |>
  mutate(date_Invertebrates = as.Date(date)) |>  
  dplyr::select(ID, date_Invertebrates) |>
  distinct()

# Extract invertebrate sampling dates within study areas of 101 IDs
Invertebrates_date_subset <- Invertebrates_date[Invertebrates_date$ID %in% unique_IDs, ]
unique(Invertebrates_date_subset$ID) # Note: Within 101 sites, 4 sites collected macroinvertebrate samples in both years (S74, S76, S78, S116) 

# -------------------------

# -------------------------
## Import and analyze each stressor from raw data 
# Filter all stressor data prior to invertebrate sampling

# -------------------------
## Pesticide metrics: Four metrics from step 2.4 
Pesticide_metrics <- read_excel("Recalculated_metrics/2.TUsum_all_joined.xlsx")
Pesticide_metrics <- Pesticide_metrics |> 
  dplyr::select(ID, year, 
                TUsum_full, TUsum_extreme_removed, 
                TUsum_before_inver, TUsum_before_inver_extreme_removed)

# Check correlations of 4 metrics (all above 0.8)
cor(Pesticide_metrics$TUsum_full,  Pesticide_metrics$TUsum_extreme_removed,  use = "pairwise.complete.obs")
cor(Pesticide_metrics$TUsum_full,  Pesticide_metrics$TUsum_before_inver,  use = "pairwise.complete.obs")
cor(Pesticide_metrics$TUsum_full,  Pesticide_metrics$TUsum_before_inver_extreme_removed,  use = "pairwise.complete.obs")
cor(Pesticide_metrics$TUsum_extreme_removed,  Pesticide_metrics$TUsum_before_inver_extreme_removed,  use = "pairwise.complete.obs")
cor(Pesticide_metrics$TUsum_before_inver,  Pesticide_metrics$TUsum_before_inver_extreme_removed,  use = "pairwise.complete.obs")

# -------------------------
## pH:
# Data collection in the manuscript: "PH was measured with every grab samples"
# Data aggregation in the manuscript: "Mean of values measured per site"
# Approach of recalculation: Filter for data at 101 study site IDs, 
#                            select only pH data before macroinvertebrate sampling dates, 
#                            and use the mean aggregation approach (following the method in the main manuscript)

# Import the raw data
pH_Field <- read.table("Input_data/4_3_pH_Field.txt", header = TRUE, sep = "\t")
unique(pH_Field$ID) # 124 IDs

# Format date and filter pH sampling to 101 study site IDs
pH_Field_date <- pH_Field |> mutate(date = as.Date(date))
pH_Field_subset <- pH_Field_date[pH_Field_date$ID %in% unique_IDs, ]

# Filter pH data to keep only samples prior or on the day of invertebrate sampling
pH_Field_valid <- pH_Field |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date <= date_Invertebrates) 

unique(pH_Field_valid$ID)  
# Note: number of samples reduced from 852 to 551, but data covers all 101/101 IDs

# Calculate means per sample (year x ID)
pH_metric <- pH_Field_valid |> 
  dplyr::select(ID, date, pH) |>  
  mutate(date = as.Date(date),
         year = as.integer(format(date, "%Y"))) |> 
  group_by(ID, year) |> 
  summarise(pH = mean(pH, na.rm = TRUE), .groups = "drop")

# Join with the original pH table for comparison
pH_metric_comparison <- pH_metric |> 
  left_join(Abiotic_metrics, by = c("ID", "year")) |> 
  dplyr::select(ID, year, pH, pH_origin)

R2 <- cor(pH_metric_comparison$pH, pH_metric_comparison$pH_origin,  
          use = "pairwise.complete.obs")

ggplot(pH_metric_comparison, aes(x = pH_origin, y = pH)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  annotate("text", 
           x = min(pH_metric_comparison$pH_origin, na.rm = TRUE),
           y = max(pH_metric_comparison$pH, na.rm = TRUE),
           hjust = 0, vjust = 1, label = paste0("R² = ", round(R2, 2)))

# Compare the matching pH values for ID x year combinations: 
sum(round(pH_metric_comparison$pH, 1) == round(pH_metric_comparison$pH_origin, 1), na.rm = TRUE)

# Note: Filtering sampling dates (i.e., reduced number of samples across IDs) and recalculating the mean value
#       resulted in only 31/116 samples matching values, but correlation is high 0.89

# -------------------------
## Nutrients: 
# Manuscript Liess et al. 2021: "Ortho-phosphate, nitrate, nitrite and ammonium concentrations were determined in all grab and EDS samples"
# Supplementary Information of Liess et al. 2021: "Data aggregation: Maximum value measured per site (Grab and EDS samples)"
# "Total nitrogen and total phosphorous concentrations may be lower than specific nutrient (NO2, NO3, NH4, PO4) concentrations in some cases due to different analytical methods (see Methods)"
# Selected approach: 
#   - Raw data extract prior to ecological sampling (analogue to data processing for Pesticide data) 
#   - Using mean value per site instead of maximum, as the study aims to associate with macroinvertebrate responses (i.e., indicator of annual pattern of stream systems rather than extreme events), 
#     , which differ between nutrients and pesticide toxicity 

# Import Nitrate, nitrite, ammonium, ortho-phosphate data
Nutrients_Field <- read.table("Input_data/4_2_Nutrients_Field.txt", header = TRUE, sep = "\t")
unique(Nutrients_Field$ID)  # 124/124 full IDs

# Import total nitrogen and total phosphorus data
Nutrients_Lab <- read.table("Input_data/4_2_Nutrients_Metals_Lab.txt", header = TRUE, sep = "\t")
unique(Nutrients_Lab$ID)  # 121/124 IDs 
# Note: 3 sites S80, S86, and S112 are missing TN and TP values compared to other nutrient forms
# Consultation with coauthors of Liess et al. 2021: Data at these sites were collected in 2019, with location switching downstream and renamed IDs from the initial locations S79, S85, and S111 
# The missing 2019 data of these sites are erroneously aggregated in the dataset under site IDs S79, S85, and S111
# A similar missing data at the above three sites for trace elements

# Correct data values for three sites 
Nutrients_Lab <- Nutrients_Lab |> 
  mutate(ID = case_when(ID == "S79"  & year(date) == 2019 ~ "S80",
                        ID == "S85"  & year(date) == 2019 ~ "S86",
                        ID == "S111" & year(date) == 2019 ~ "S112",
                        TRUE ~ ID))

unique(Nutrients_Lab$ID)  # 124/124 IDs 

# Format date and filter nutrients sampling to 101 study site IDs
Nutrients_Field_date <- Nutrients_Field |> mutate(date_Nutrients_Field = as.Date(date))
Nutrients_Field_subset <- Nutrients_Field_date[Nutrients_Field_date$ID %in% unique_IDs, ]

Nutrients_Lab_date <- Nutrients_Lab |> mutate(date_Nutrients_Lab = as.Date(date))
Nutrients_Lab_subset <- Nutrients_Lab_date[Nutrients_Lab_date$ID %in% unique_IDs, ]

# Filter nutrients to keep only samples prior or on the day of invertebrate sampling
Nutrients_Field_valid <- Nutrients_Field_subset |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date_Nutrients_Field <= date_Invertebrates) 

Nutrients_Lab_valid <- Nutrients_Lab_subset |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date_Nutrients_Lab <= date_Invertebrates) 

unique(Nutrients_Field_valid$ID)  # 101/101 IDs
unique(Nutrients_Lab_valid$ID)    # 101/101 IDs
setdiff(unique_IDs, unique(Nutrients_Lab_valid$ID))

# Correlation among nutrient species
cor(Nutrients_Field_valid$NO2,  Nutrients_Field_valid$NO2_N,  use = "pairwise.complete.obs")
cor(Nutrients_Field_valid$NO3,  Nutrients_Field_valid$NO3_N,  use = "pairwise.complete.obs")
cor(Nutrients_Field_valid$PO4,  Nutrients_Field_valid$PO4_P,  use = "pairwise.complete.obs")
cor(Nutrients_Field_valid$NH4,  Nutrients_Field_valid$NH4_N,  use = "pairwise.complete.obs")
# Note: all R2 ~ 1.00, thus no need to use both  
#       Select NO2_N, NO3_N, PO4_P, NH4_N variables to better link to ecological impacts than their molecular ions NO2, NO3, NH4, and PO4

# Calculate metric means per sample (year x ID)
Nutrients_Field_metrics <- Nutrients_Field_valid |> 
  dplyr::select(ID, date_Nutrients_Field, NO2_N, NO3_N, NH4_N, PO4_P) |>  
  rename(date = date_Nutrients_Field) |>                         
  mutate(date = as.Date(date),               
    year = format(date, "%Y")) |> 
  group_by(ID, year) |> 
  summarise(across(c(NO2_N, NO3_N, NH4_N, PO4_P), ~ mean(.x, na.rm = TRUE)), 
            .groups = "drop")

Nutrients_Lab_metrics <- Nutrients_Lab_valid |> 
  dplyr::select(ID, date_Nutrients_Lab, TNb_mgL, TP_mgL) |>  
  rename(date = date_Nutrients_Lab,
         TP = "TP_mgL",
         TN = "TNb_mgL") |>                         
  mutate(date = as.Date(date),               
    year = format(date, "%Y")) |> 
  group_by(ID, year) |> 
  summarise(across(c(TN, TP), ~ mean(.x, na.rm = TRUE)), 
            .groups = "drop")

# Join all nutrient variables
Nutrients_metrics <- Nutrients_Field_metrics |> 
  left_join(Nutrients_Lab_metrics, by = c("ID", "year"))  

## Plot the changes in correlation of recalculated nutrient metrics (using mean method) with the metrics provided in the main manuscript (using maximum values)
# Joint table for plotting
Nutrients_metrics_plotting <- Nutrients_metrics |> 
  dplyr::mutate(year = as.integer(year)) |> 
  inner_join(Abiotic_metrics, by = c("ID", "year")) |> 
  dplyr::select(ID, year,
         NO2_N, NO3_N, NH4_N, PO4_P, TN, TP,
         NO2_origin, NO3_origin, NH4_origin, PO4_origin, TN_origin, TP_origin) 

# Variable pairs for plotting
nutrient_vars <- c("NO2_N", "NO3_N", "NH4_N", "PO4_P", "TN", "TP")
origin_vars   <- c("NO2_origin", "NO3_origin", "NH4_origin",
                   "PO4_origin", "TN_origin", "TP_origin")
# Correlation
pairs_vars <- c(nutrient_vars, origin_vars)
data_corr  <- Nutrients_metrics_plotting[, pairs_vars]
ggpairs(data_corr,
        upper = list(continuous = wrap("cor", size = 3)),
        lower = list(continuous = wrap("points", alpha = 0.5, size = 1)))

# Compare the matching nutrient values for ID x year combinations: 
sum(round(Nutrients_metrics_plotting$NO2_N, 1) == round(Nutrients_metrics_plotting$NO2_origin, 1), na.rm = TRUE) # match = 19/112, R2 = 0.671
sum(round(Nutrients_metrics_plotting$NO3_N, 1) == round(Nutrients_metrics_plotting$NO3_origin, 1), na.rm = TRUE) # match = 3/112, R2 = 0.860
sum(round(Nutrients_metrics_plotting$NH4_N, 1) == round(Nutrients_metrics_plotting$NH4_origin, 1), na.rm = TRUE) # match = 30/112, R2 = 0.762
sum(round(Nutrients_metrics_plotting$PO4_P, 1) == round(Nutrients_metrics_plotting$PO4_origin, 1), na.rm = TRUE) # match = 4/112, R2 = 0.657
sum(round(Nutrients_metrics_plotting$TN, 1) == round(Nutrients_metrics_plotting$TN_origin, 1), na.rm = TRUE) # match = 1/112, R2 = 0.850
sum(round(Nutrients_metrics_plotting$TP, 1) == round(Nutrients_metrics_plotting$TP_origin, 1), na.rm = TRUE) # match = 37/112, R2 = 0.801

# Note: All pair metrics (recalculated mean with reduced sampling vs max original) 
#       do not numerically match (i.e. exact values), due to different species selection (with *_N and *_P), different aggregation method (mean instead of max), and filtering date prior to macroinvertebrate sampling
#       but all recalculated metrics have correlation > 0.65

# -------------------------
##### Trace elements:
# Manuscript Liess et al. 2021: "The toxicity of trace elements was calculated using literature LC50 values. The local maximum of summed TUs (TUsum) including all trace elements per sample is considered in the multiple linear regression."
# Supplementary Information of Liess et al. 2021: ""Toxic Unit aggregated from trace elements to quantify their toxicity (lead, mercury, copper, arsenic, cadmium, zinc"
# Data aggregation: "Maximum of summed Toxic Units per site"
# Selected approach: Use analogue approach of metric calculation and max aggregation of TUsum of trace elements as for TUsum of pesticides in the main manuscript, 
#                    but first correcting values at three sites and filtering the valid sampling prior to macroinvertebrate sampling 

# Extract trace element data provided in the nutrient_lab file 
names(Nutrients_Lab_subset)

Trace_elements <- Nutrients_Lab_subset |> 
  dplyr::select(ID, date, As_ugL, Pb_ugL, Cd_ugL, Cu_ugL, Hg_ugL, Zn_ugL) |> 
  mutate(date = as.Date(date))

# Filter trace element data to keep only samples prior or on the day of macroinvertebrate sampling
Trace_elements_valid <- Trace_elements |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date <= date_Invertebrates) 

unique(Trace_elements_valid$ID)        # Filtering reduced sample size from 715 to 525 observations, but data available at 101/101 IDs
setdiff(unique_IDs, unique(Trace_elements_valid$ID))

# Reformat trace element data for TU calculation
# LC50 information obtained from the Table SI3, page 19/36 of the Supplementary Information of the main manuscript Liess et al. 2021
Trace_elements_filtered <- Trace_elements_valid |> 
  mutate(year = as.integer(format(date, "%Y"))) |> 
  dplyr::select(-date_Invertebrates) |> 
  pivot_longer(cols = c(As_ugL, Pb_ugL, Cd_ugL, Cu_ugL, Hg_ugL, Zn_ugL),
    names_to = "substance",
    values_to = "concentration_µgL") |> 
  dplyr::mutate(LC50_Daphnia_µgL = case_when(    
      substance == "As_ugL" ~ 3800,
      substance == "Pb_ugL" ~ 398,
      substance == "Cd_ugL" ~ 55,
      substance == "Cu_ugL" ~ 26.3,
      substance == "Hg_ugL" ~ 56,
      substance == "Zn_ugL" ~ 617,
      TRUE ~ NA_real_))

# Trace element toxic unit metric
Trace_elements_metric <- Trace_elements_filtered |> 
  mutate(TUsubstance_unlog = concentration_µgL / LC50_Daphnia_µgL) |>
  group_by(ID, year, date) |>
  summarise(TUsum_sample_unlog = sum(TUsubstance_unlog, na.rm = TRUE), 
            .groups = "drop") |>
  mutate(TUsum_sample = log10(TUsum_sample_unlog)) |>
  group_by(ID, year) |>
  summarise(TUsum_Trace_elements = max(TUsum_sample, na.rm = TRUE), 
            .groups = "drop")

# Join with the table having original trace element calculation to compare recalculated vs original metric values
Trace_element_metrics_comparison <- Trace_elements_metric |> 
  left_join(Abiotic_metrics, by = c("ID", "year")) |> 
  dplyr::select(ID, year, TUsum_Trace_elements, TUmetal_origin)

R2_metal <- cor(Trace_element_metrics_comparison$TUsum_Trace_elements, Trace_element_metrics_comparison$TUmetal_origin,  
          use = "pairwise.complete.obs")

ggplot(Trace_element_metrics_comparison, aes(x = TUmetal_origin, y = TUsum_Trace_elements)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  annotate("text", 
           x = min(Trace_element_metrics_comparison$TUmetal_origin, na.rm = TRUE),
           y = max(Trace_element_metrics_comparison$TUsum_Trace_elements, na.rm = TRUE),
           hjust = 0, vjust = 1, label = paste0("R² = ", round(R2_metal, 2)))

# Compare the matching trace element values for ID x year combinations: 
sum(round(Trace_element_metrics_comparison$TUsum_Trace_elements, 1) == round(Trace_element_metrics_comparison$TUmetal_origin, 1), na.rm = TRUE)

# Note: Filtering sampling dates (i.e., reduced number of samples across IDs) + recalculating from raw data
#       resulted in only 12/116 samples matching values, but correlation with the metric in the original manuscript is high (0.82) 

# -------------------------
##### Temperature, DO, and flow:
# Manuscript Liess et al. 2021: "Oxygen, temperature, water level was continuously measured throughout the sampling period from April to June"
# Supplementary Information of Liess et al. 2021: 
#   - Dissolved oxygen: "Data aggregation: 25th quantile of all values measured per site, missing values were replaced by average of similar sites"
#   - Temperature "Data aggregation: 75th quantile of all measuring points in data series, missing values were replaced by average of similar sites"
#   - Flow velocity: "Data aggregation: Geometric mean of all values measured per site"
# Selected approach: 
# - For flow velocity: check data completeness, filter to 101 study sites and sampling prior to macroinvertebrate sampling, and aggregate the mean values for consistencies with other abiotic metrics
# - For temperature and DO: data summary suggests significant missing data and unclear approach to fill NA from nearby sites. 
#                           Check the data gap and consider to remove these two metrics.

# Import flow velocity data
Flow_data <- read.table("Input_data/3_5_Sites_Flow_Velocity.txt", header = TRUE, sep = "\t")
unique(Flow_data$ID)  # 124/124 full IDs
summary(Flow_data)    # No NAs, data range from 0 to 1.3 

# Import temperature and DO
Temp_DO_data <- read.table("Input_data/4_1_Temp_ElectrConduct_Oxygen_Pressure.txt", header = TRUE, sep = "\t")
unique(Temp_DO_data$ID)  # 112/124 IDs
summary(Temp_DO_data)    # Note: Temperature: 52470/5526741 ~ 1% NAs, OK 
                         #       DO:          5004015/5526741 ~ 91% NAs, remove DO from further analyses
Temp_data <- Temp_DO_data |> 
  dplyr::select(ID, date, temperature)

# Format date and filter flow and temperature sampling to 101 study site IDs
Flow_date <- Flow_data |> mutate(date_Flow = as.Date(date))
Flow_subset <- Flow_date[Flow_date$ID %in% unique_IDs, ]

Temp_date <- Temp_data |> mutate(date_Temp = as.Date(date))
Temp_subset <- Temp_date[Temp_date$ID %in% unique_IDs, ]

# Filter flow and temperature to keep only samples prior or on the day of invertebrate sampling
Flow_valid <- Flow_subset |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date_Flow <= date_Invertebrates) 

Temp_valid <- Temp_subset |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date_Temp <= date_Invertebrates) 

unique(Flow_valid$ID)  # 100/101 IDs, missing one site S52 (because raw flow data was in late June and July, which are after macroinvertebrate sampling)
setdiff(unique_IDs, unique(Flow_valid$ID))

unique(Temp_valid$ID)  # 66/101 ~ 65% IDs => remove Temperature from further analyses as the focus is on effects of pesticides at 101 sites in relation to other stressors with available data

# Calculate mean flow per sample (year x ID)
Flow_metric <- Flow_valid |> 
  dplyr::select(ID, date, flow_velocity_ms) |>  
  mutate(date = as.Date(date),
         year = as.integer(format(date, "%Y"))) |> 
  group_by(ID, year) |> 
  summarise(Flow = mean(flow_velocity_ms, na.rm = TRUE), .groups = "drop")

# Join with the original flow table and calculate correlation
Flow_metric_comparison <- Flow_metric |> 
  left_join(Abiotic_metrics, by = c("ID", "year")) |> 
  dplyr::select(ID, year, Flow, Flow_origin)

R2 <- cor(Flow_metric_comparison$Flow, Flow_metric_comparison$Flow_origin,  
          use = "pairwise.complete.obs")

ggplot(Flow_metric_comparison, aes(x = Flow_origin, y = Flow)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue")+
  annotate("text", 
           x = min(Flow_metric_comparison$Flow_origin, na.rm = TRUE),
           y = max(Flow_metric_comparison$Flow, na.rm = TRUE),
           hjust = 0, vjust = 1, label = paste0("R² = ", round(R2, 2)))

# Compare the matching pH values for ID x year combinations: 
sum(round(Flow_metric_comparison$Flow, 1) == round(Flow_metric_comparison$Flow_origin, 1), na.rm = TRUE)

# Note: Filtering sampling dates (i.e., reduced number of samples across IDs) + recalculating from raw data the mean (unchanged aggregation method as in main manuscript) 
#       resulted in 61/115 samples matching values, but the correlation is high (0.82) and no extreme pattern

# -------------------------
## Join tables of recalculated metrics (with filtered sampling dates prior to macroinvertebrate sampling) and save the file for further analyses

# Extract metrics used in the main manuscript
Abiotic_metrics_selection <- Abiotic_metrics |>
  dplyr::select(ID, year, AgriLand_percent, Morphology, Habitat, Stream_width, Stream_depth) |>
  distinct(ID, year, .keep_all = TRUE)

# Check that other tables have distinct ID and year for joining
Pesticide_metrics_selection <- Pesticide_metrics |> 
  mutate(year = as.integer(year)) |> 
  distinct(ID, year, .keep_all = TRUE)
pH_metric_selection <- pH_metric |> 
  mutate(year = as.integer(year)) |> 
  distinct(ID, year, .keep_all = TRUE)
Nutrients_metrics_selection <- Nutrients_metrics |> 
  dplyr::select(ID, year, TN, TP) |> 
  mutate(year = as.integer(year)) |> 
  distinct(ID, year, .keep_all = TRUE)
Trace_elements_metric_selection <- Trace_elements_metric |> 
  mutate(year = as.integer(year)) |> 
  distinct(ID, year, .keep_all = TRUE)
Flow_metric_selection <- Flow_metric |> 
  mutate(year = as.integer(year)) |> 
  distinct(ID, year, .keep_all = TRUE)

# Join all five abiotic variable tables to the Invertebrate_metrics table
All_metrics <- Invertebrate_metrics |>
  left_join(Abiotic_metrics_selection,        by = c("ID", "year")) |>
  left_join(Pesticide_metrics_selection,      by = c("ID", "year")) |>
  left_join(pH_metric_selection,              by = c("ID", "year")) |>
  left_join(Nutrients_metrics_selection,      by = c("ID", "year")) |>
  left_join(Trace_elements_metric_selection,  by = c("ID", "year")) |>
  left_join(Flow_metric_selection,            by = c("ID", "year"))

# Save recalculated metrics
write_xlsx(All_metrics, "Recalculated_metrics/3.All_biotic_abiotic_metrics_mean_aggregation.xlsx")

# -------------------------

# -------------------------
# Multiple linear regression 
All_metrics <- read_excel("Recalculated_metrics/3.All_biotic_abiotic_metrics_mean_aggregation.xlsx")

## Aggregate the mean of year 2018 and 2019 for all metrcis to use in the multiple regression models
All_metrics_mean <- All_metrics |>
  mutate(Flow = if_else(is.na(Flow), mean(Flow, na.rm = TRUE), Flow),  # impute one NA in Flow values
    TN = log10(TN),                                                    # log transform nutrients for linearity model assumption
    TP = log10(TP)) |>
  dplyr::select(-year) |>
  group_by(ID) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  arrange(as.numeric(sub("S", "", ID)))

# ---------------------

# ---------------------
# Check correlation among predictor & response variables 
dat <- All_metrics_mean

# List of variables included in the correlation
cor_vars <- c("speartaxa", "sapr", "eptProz",
              "AgriLand_percent", "Morphology", "Habitat",
              "Stream_width", "Stream_depth",
              "TUsum_full", "TUsum_extreme_removed", "TUsum_before_inver", "TUsum_before_inver_extreme_removed", 
              "TUsum_Trace_elements",
              "pH", "TN", "TP", "Flow")

# Correlation data table
cor_data <- dat |>
  dplyr::select(dplyr::all_of(cor_vars))

# Correlation matrix
cor_mat <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")

# Function for the matrix with p-values of pairwise correlations
cor_test_matrix <- function(mat, ...) {
  mat    <- as.matrix(mat)
  n      <- ncol(mat)
  p_mat  <- matrix(NA_real_, n, n)
  diag(p_mat) <- 0

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value}}

  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  p_mat}

# Matrix of p-values
p_mat <- cor_test_matrix(cor_data)

# Plot correlation matrix
png("correlation_plot.png", width = 2300, height = 1500, res = 300)

corrplot(cor_mat, order = "original", method = "number",
         tl.col = "black", number.cex = 0.55, tl.srt = 45,
         p.mat = p_mat, sig.level = 0.05,
         insig = "blank")

dev.off()

# ---------------------

# ---------------------
# Check linearity: predictor–response scatterplots
# Response variables 
response_vars <- c("speartaxa", "sapr", "eptProz")

# Predictor variables 
predictor_vars <- c("TUsum_full", 
                    "TUsum_extreme_removed", 
                    "TUsum_before_inver", 
                    "TUsum_before_inver_extreme_removed", 
                    "TUsum_Trace_elements",
                    "AgriLand_percent",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow")

# Scatter plot for linearity check
for (resp in response_vars) {
  dat_long <- All_metrics_mean |>
    dplyr::select(dplyr::all_of(c(resp, predictor_vars))) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(predictor_vars),
      names_to  = "predictor",
      values_to = "x")

  p <- ggplot(dat_long, aes(x = x, y = .data[[resp]])) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, colour = "red") +
    facet_wrap(~ predictor, scales = "free_x") +
    labs(title = paste("Linearity check for response:", resp),
         x = "Predictor value", y = resp) +
    theme_minimal()

  print(p)
}

# Note: Linearity assumption valid for all three response variables
#---------------------

# ---------------------
# Standardize predictors for comparability of estimates across continous predictors
All_metrics_scaled <- All_metrics_mean

All_metrics_scaled[, 5:ncol(All_metrics_mean)] <- scale(All_metrics_mean[, 5:ncol(All_metrics_mean)])

summary(All_metrics_scaled)

#---------------------
#---------------------
# Model for each pesticide toxicity metric
# Include: - "TUsum_full"
#          - "TUsum_extreme_removed", 
#          - "TUsum_before_inver", 
#          - "TUsum_before_inver_extreme_removed"
# ---------------------

# ---------------------
## Model 1 "TUsum_full"
# Model for SPEARpesticides (speartaxa)
resp <- "speartaxa"

predictor_vars <- c("TUsum_full", 
                    "AgriLand_percent",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow",
                    "TUsum_Trace_elements")

dat1.1 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model1.1 <- lm(as.formula(paste(resp, "~ .")), data = dat1.1)
summary(full_model1.1)

vif1.1 <- car::vif(full_model1.1)
vif1.1

keep_pred1.1 <- names(vif1.1)[vif1.1 <= 2]
keep_pred1.1

dat1.1_red <- dat1.1 |>
  dplyr::select(dplyr::all_of(c(resp, keep_pred1.1)))

full_model1.1_red <- lm(as.formula(paste(resp, "~ .")), data = dat1.1_red)
summary(full_model1.1_red)

final_model1.1 <- MASS::stepAIC(full_model1.1_red, direction = "backward", trace = FALSE)
summary(final_model1.1)

# Model residual check (validation)
par(mfrow = c(2, 2))
plot(final_model1.1)
par(mfrow = c(1, 1))
dev.off()

# Relative importance
rel1.1  <- relaimpo::calc.relimp(final_model1.1, type = "lmg", rela = FALSE)
coef1.1 <- broom::tidy(final_model1.1)

res1.1 <- data.frame(indicator = resp,
                   stressor  = names(rel1.1$lmg),
                   r2        = as.numeric(rel1.1$lmg),
                   estimate  = coef1.1$estimate[match(names(rel1.1$lmg), coef1.1$term)],
                   signifi   = coef1.1$p.value[match(names(rel1.1$lmg), coef1.1$term)])

res1.1 <- rbind(res1.1, data.frame(indicator = resp,
                               stressor  = "Full model",
                               r2        = rel1.1$R2,
                               estimate  = NA,
                               signifi   = NA))

# ---------------------
# Model for %EPT
resp <- "eptProz"

dat1.2 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model1.2 <- lm(as.formula(paste(resp, "~ .")), data = dat1.2)
summary(full_model1.2)

vif1.2 <- car::vif(full_model1.2)
vif1.2 

keep_pred1.2 <- names(vif1.2)[vif1.2 <= 2]   
keep_pred1.2

dat1.2_red <- dat1.2 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred1.2)))

full_model1.2_red <- lm(as.formula(paste(resp, "~ .")), data = dat1.2_red)
summary(full_model1.2_red)

final_model1.2    <- stepAIC(full_model1.2_red, direction = "backward", trace = FALSE)
summary(final_model1.2)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model1.2)
par(mfrow = c(1,1))
dev.off()

rel1.2  <- calc.relimp(final_model1.2, type = "lmg", rela = FALSE)
coef1.2 <- broom::tidy(final_model1.2)

res1.2 <- data.frame(indicator = resp,
                   stressor  = names(rel1.2$lmg),
                   r2        = as.numeric(rel1.2$lmg),
                   estimate  = coef1.2$estimate[match(names(rel1.2$lmg), coef1.2$term)],
                   signifi   = coef1.2$p.value[match(names(rel1.2$lmg), coef1.2$term)])

res1.2 <- rbind(res1.2, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel1.2$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Model for Saprobic index (same steps as for SPEARpesticide)
resp <- "sapr"

dat1.3 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model1.3 <- lm(as.formula(paste(resp, "~ .")), data = dat1.3)
summary(full_model1.3)

vif1.3 <- car::vif(full_model1.3)
vif1.3

keep_pred1.3 <- names(vif1.3)[vif1.3 <= 2] 
keep_pred1.3

dat1.3_red <- dat1.3 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred1.3)))

full_model1.3_red <- lm(as.formula(paste(resp, "~ .")), data = dat1.3_red)
summary(full_model1.3_red)

final_model1.3 <- stepAIC(full_model1.3_red, direction = "backward", trace = FALSE)
summary(final_model1.3)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model1.3)
par(mfrow = c(1,1))
dev.off()

rel1.3  <- calc.relimp(final_model1.3, type = "lmg", rela = FALSE)

coef1.3 <- broom::tidy(final_model1.3)

res1.3 <- data.frame(indicator = resp, 
                   stressor  = names(rel1.3$lmg),
                   r2        = as.numeric(rel1.3$lmg),
                   estimate  = coef1.3$estimate[match(names(rel1.3$lmg), coef1.3$term)],
                   signifi   = coef1.3$p.value[match(names(rel1.3$lmg), coef1.3$term)])

res1.3 <- rbind(res1.3, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel1.3$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Combine model results for all three response varialbes
groesult_simple1 <- dplyr::bind_rows(res1.1, res1.2, res1.3)|>
  dplyr::select(indicator, stressor, signifi, r2, estimate)

# Save results
write_xlsx(groesult_simple1, "Recalculated_metrics/3.1.Final_model_results_TUsum_full.xlsx")

#----------------------
# Plot
# List of stressors in the plot
stress_levels_simple1 <- c("Best-fit model",
                          "Pesticide\ntoxicity",          # TUsum_Pesticide (changed)
                          "Trace\nelements",              # TUsum_Trace_elements
                          "%Agricultural\nland uses",     # AgriLand_percent (changed)
                          "Hydromorphological\ndegradation", # Morphology (changed)
                          "Bed habitat\nstructure",       # Habitat (changed)
                          "pH",                           # pH
                          "TP",                           # TP
                          "Flow\nvelocity",               # Flow
                          "Stream depth")                 # Stream depth


# Labels for y-axis
y_labels1 <- c("Best-fit model"              = " Best-fit model",
              "Pesticide\ntoxicity"          = "Pesticide\ntoxicity",
              "Trace\nelements"              = "Trace\nelements",
              "%Agricultural\nland uses"     = "%Agricultural\nland uses",
              "Hydromorphological\ndegradation" = "Hydromorphological\ndegradation",
              "Bed habitat\nstructure"       = "Bed habitat\nstructure",
              "pH"                           = "pH",
              "TP"                           = "TP",
              "Flow\nvelocity"               = "Flow\nvelocity",
              "Stream_depth"                 = "Stream depth")

# Table for plotting

groesult_new1 <- groesult_simple1 |> 
  # Significance labels
  dplyr::mutate(signifi = dplyr::case_when(signifi <= 0.001 ~ "***",
                                           signifi <= 0.01  ~ "**",
                                           signifi <= 0.05  ~ "*",
                                           TRUE             ~ NA_character_),

    # Recode stressors for plot labels 
    stressor = dplyr::recode(stressor,
      "Full model"           = "Best-fit model",
      "TUsum_full"           = "Pesticide\ntoxicity",
      "TUsum_Trace_elements" = "Trace\nelements",
      "AgriLand_percent"     = "%Agricultural\nland uses",
      "Morphology"           = "Hydromorphological\ndegradation",
      "Habitat"              = "Bed habitat\nstructure",
      "pH"                   = "pH",
      "TP"                   = "TP",
      "Flow"                 = "Flow\nvelocity",
      "Stream_depth"         = "Stream depth"),

    # Rename & reorder response variables
    indicator = dplyr::recode(indicator,
      "speartaxa" = "SPEARpesticides",
      "sapr"      = "Saprobic\nindex",
      "eptProz"   = "%EPT"),

    # Colours: black for Full model, red for negative estimate, blue for positive
    Col = dplyr::case_when(stressor == "Best-fit model" ~ "black",
      estimate <= 0                      ~ "firebrick3",
      TRUE                               ~ "dodgerblue3"),
    
    # Effects
    direction = dplyr::case_when(
      stressor == "Best-fit model" ~ "Best-fit model",
      estimate < 0 ~ "Negative",
      estimate > 0 ~ "Positive",
      TRUE ~ NA_character_),

    r2     = round(as.numeric(r2), 2),
    labels = dplyr::if_else(!is.na(signifi), paste(r2, signifi), as.character(r2)),
    labels = dplyr::if_else(labels %in% c("NA", "NA NA"), NA_character_, labels)) 

# Stress levels: only keep those that are actually present in the final model
used_levels <- stress_levels_simple1[stress_levels_simple1 %in% groesult_new1$stressor]

groesult_new1$stressor <- factor(groesult_new1$stressor,
                                levels = rev(used_levels))  
groesult_new1$indicator <- factor(groesult_new1$indicator,
                                 levels = c("SPEARpesticides", "%EPT", "Saprobic\nindex"))

n_stress <- nlevels(groesult_new1$stressor)
full_row <- which(levels(groesult_new1$stressor) == "Best-fit model")

# ---------------------------------------------------------------------------
# Plotting 
# ---------------------------------------------------------------------------

(Fig_multiple_stressors1 <- ggplot(groesult_new1, aes(x = indicator, y = stressor, size = r2, label = labels, colour = Col)) +
  
  # Light horizontal grid lines for each stressor row
  geom_hline(yintercept = seq(1.5, n_stress - 0.5, 1), colour = "grey80", linewidth = 0.25) +
  
  # Black lines separating sections 
  geom_hline(yintercept = c(9.5, 8.5, 0.5), linewidth = 0.5, colour = "#000000") +
  
  # R² bubbles with adjusted transparency
  geom_point(alpha = 0.75, position = position_nudge(y = +0.16)) +
  
  # Text labels (R² + significance stars within the table) 
  geom_text(stat = "identity", 
            size = 4.5,  
            colour = "#000000", 
            family = "Arial",  
            vjust = 1, 
            alpha = 1, 
            position = position_nudge(y = -0.15)) +
  
  # X axis – indicators
  scale_x_discrete(position = "top",
                   labels = c("SPEARpesticides" = expression(SPEAR[pesticides]),
                              "%EPT" = "%EPT",
                              "Saprobic\nindex" = "Saprobic\nindex")) +
  
  # Y axis – stressors
  scale_y_discrete(drop = FALSE, labels = y_labels1) +
  
  # Color and size scales
  scale_colour_identity(guide = "none") +
  scale_size(range = c(3, 14)) +  
  
  # Theme 
  theme_classic(base_size = 12, base_family = "Arial") +  
  theme(
    # Background
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Axes
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Axis text 
    axis.text.x = element_text(size = 14, 
                               vjust = 0.5, 
                               hjust = 0.5,
                               colour = "#000000",
                               face = "plain"),
    axis.text.y = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 1,
                               colour = "#000000"),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 11, face = "plain"),  
    legend.text = element_text(size = 11),  
    legend.key.size = unit(0.6, "cm"),  
    legend.background = element_blank(),
    
    # Margins - tighter for publication
    plot.margin = unit(c(5, 2, 2, 2), "mm")) +
  
  coord_cartesian(xlim = c(0.5, 3), clip = "off"))+

  # Effect sign legend 
  annotate("point", x = 3.45, y = 0.7, size = 4.5, colour = "firebrick3") +
  annotate("text",  x = 3.55, y = 0.7, label = "Negative", hjust = 0, size = 3.6, family = "Arial")  +

  annotate("point", x = 3.45, y = 1.2, size = 4.5, colour = "dodgerblue3") +
  annotate("text",  x = 3.55, y = 1.2, label = "Positive", hjust = 0, size = 3.6, family = "Arial")


# ----------------------------------

# ----------------------------------
## Model 2 "TUsum_extreme_removed"
# Model for SPEARpesticides (speartaxa)
resp <- "speartaxa"

predictor_vars <- c("TUsum_extreme_removed", 
                    "AgriLand_percent",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow",
                    "TUsum_Trace_elements")

dat2.1 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model2.1 <- lm(as.formula(paste(resp, "~ .")), data = dat2.1)
summary(full_model2.1)

vif2.1 <- car::vif(full_model2.1)
vif2.1

keep_pred2.1 <- names(vif2.1)[vif2.1 <= 2]
keep_pred2.1

dat2.1_red <- dat2.1 |>
  dplyr::select(dplyr::all_of(c(resp, keep_pred2.1)))

full_model2.1_red <- lm(as.formula(paste(resp, "~ .")), data = dat2.1_red)
summary(full_model2.1_red)

final_model2.1 <- MASS::stepAIC(full_model2.1_red, direction = "backward", trace = FALSE)
summary(final_model2.1)

# Model residual check (validation)
par(mfrow = c(2, 2))
plot(final_model2.1)
par(mfrow = c(1, 1))
dev.off()

# Relative importance
rel2.1  <- relaimpo::calc.relimp(final_model2.1, type = "lmg", rela = FALSE)
coef2.1 <- broom::tidy(final_model2.1)

res2.1 <- data.frame(indicator = resp,
                   stressor  = names(rel2.1$lmg),
                   r2        = as.numeric(rel2.1$lmg),
                   estimate  = coef2.1$estimate[match(names(rel2.1$lmg), coef2.1$term)],
                   signifi   = coef2.1$p.value[match(names(rel2.1$lmg), coef2.1$term)])

res2.1 <- rbind(res2.1, data.frame(indicator = resp,
                               stressor  = "Full model",
                               r2        = rel2.1$R2,
                               estimate  = NA,
                               signifi   = NA))

# ---------------------
# Model for %EPT
resp <- "eptProz"

dat2.2 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model2.2 <- lm(as.formula(paste(resp, "~ .")), data = dat2.2)
summary(full_model2.2)

vif2.2 <- car::vif(full_model2.2)
vif2.2 

keep_pred2.2 <- names(vif2.2)[vif2.2 <= 2]   # identical to first model, TP is excluded
keep_pred2.2

dat2.2_red <- dat2.2 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred2.2)))

full_model2.2_red <- lm(as.formula(paste(resp, "~ .")), data = dat2.2_red)
summary(full_model2.2_red)

final_model2.2    <- stepAIC(full_model2.2_red, direction = "backward", trace = FALSE)
summary(final_model2.2)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model2.2)
par(mfrow = c(1,1))
dev.off()

rel2.2  <- calc.relimp(final_model2.2, type = "lmg", rela = FALSE)
coef2.2 <- broom::tidy(final_model2.2)

res2.2 <- data.frame(indicator = resp,
                   stressor  = names(rel2.2$lmg),
                   r2        = as.numeric(rel2.2$lmg),
                   estimate  = coef2.2$estimate[match(names(rel2.2$lmg), coef2.2$term)],
                   signifi   = coef2.2$p.value[match(names(rel2.2$lmg), coef2.2$term)])

res2.2 <- rbind(res2.2, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel2.2$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Model for Saprobic index (same steps as for SPEARpesticide)
resp <- "sapr"

dat2.3 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model2.3 <- lm(as.formula(paste(resp, "~ .")), data = dat2.3)
summary(full_model2.3)

vif2.3 <- car::vif(full_model2.3)
vif2.3

keep_pred2.3 <- names(vif2.3)[vif2.3 <= 2] 
keep_pred2.3

dat2.3_red <- dat2.3 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred2.3)))

full_model2.3_red <- lm(as.formula(paste(resp, "~ .")), data = dat2.3_red)
summary(full_model2.3_red)

final_model2.3 <- stepAIC(full_model2.3_red, direction = "backward", trace = FALSE)
summary(final_model2.3)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model2.3)
par(mfrow = c(1,1))
dev.off()

rel2.3  <- calc.relimp(final_model2.3, type = "lmg", rela = FALSE)

coef2.3 <- broom::tidy(final_model2.3)

res2.3 <- data.frame(indicator = resp, 
                   stressor  = names(rel2.3$lmg),
                   r2        = as.numeric(rel2.3$lmg),
                   estimate  = coef2.3$estimate[match(names(rel2.3$lmg), coef2.3$term)],
                   signifi   = coef2.3$p.value[match(names(rel2.3$lmg), coef2.3$term)])

res2.3 <- rbind(res2.3, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel2.3$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Combine model results for all three response varialbes
groesult_simple2 <- dplyr::bind_rows(res2.1, res2.2, res2.3)|>
  dplyr::select(indicator, stressor, signifi, r2, estimate)

# Save results
write_xlsx(groesult_simple2, "Recalculated_metrics/3.2. Final_model_results_TUsum_extreme_removed.xlsx")

#----------------------
#----------------------
# Plot
# List of stressors in the plot
stress_levels_simple2 <- c("Best-fit model",
                          "Pesticide\ntoxicity",          # TUsum_Pesticide 
                          "Trace\nelements",              # TUsum_Trace_elements
                          "%Agricultural\nland uses",     # AgriLand_percent 
                          "Hydromorphological\ndegradation", # Morphology 
                          "Bed habitat\nstructure",       # Habitat 
                          "pH",                           # pH
                          "TP",                           # TP
                          "Flow\nvelocity",               # Flow
                          "Stream depth")                 # Stream depth


# Labels for y-axis
y_labels2 <- c("Best-fit model"              = " Best-fit model",
              "Pesticide\ntoxicity"          = "Pesticide\ntoxicity",
              "Trace\nelements"              = "Trace\nelements",
              "%Agricultural\nland uses"     = "%Agricultural\nland uses",
              "Hydromorphological\ndegradation" = "Hydromorphological\ndegradation",
              "Bed habitat\nstructure"       = "Bed habitat\nstructure",
              "pH"                           = "pH",
              "TP"                           = "TP",
              "Flow\nvelocity"               = "Flow\nvelocity",
              "Stream_depth"                 = "Stream depth")

# Table for plotting

groesult_new2 <- groesult_simple2 |> 
  # Significance labels
  dplyr::mutate(signifi = dplyr::case_when(signifi <= 0.001 ~ "***",
                                           signifi <= 0.01  ~ "**",
                                           signifi <= 0.05  ~ "*",
                                           TRUE             ~ NA_character_),

    # Recode stressors for plot labels (save space)
    stressor = dplyr::recode(stressor,
      "Full model"           = "Best-fit model",
      "TUsum_extreme_removed"  = "Pesticide\ntoxicity",
      "TUsum_Trace_elements" = "Trace\nelements",
      "AgriLand_percent"     = "%Agricultural\nland uses",
      "Morphology"           = "Hydromorphological\ndegradation",
      "Habitat"              = "Bed habitat\nstructure",
      "pH"                   = "pH",
      "TP"                   = "TP",
      "Flow"                 = "Flow\nvelocity",
      "Stream_depth"         = "Stream depth"),

    # Rename & reorder response variables
    indicator = dplyr::recode(indicator,
      "speartaxa" = "SPEARpesticides",
      "sapr"      = "Saprobic\nindex",
      "eptProz"   = "%EPT"),

    # Colours: black for Full model, red for negative estimate, blue for positive
    Col = dplyr::case_when(stressor == "Best-fit model" ~ "black",
      estimate <= 0                      ~ "firebrick3",
      TRUE                               ~ "dodgerblue3"),

    r2     = round(as.numeric(r2), 2),
    labels = dplyr::if_else(!is.na(signifi), paste(r2, signifi), as.character(r2)),
    labels = dplyr::if_else(labels %in% c("NA", "NA NA"), NA_character_, labels)) 

# Stress levels: only keep those that are actually present in the final model
used_levels <- stress_levels_simple2[stress_levels_simple2 %in% groesult_new2$stressor]

groesult_new2$stressor <- factor(groesult_new2$stressor,
                                levels = rev(used_levels))  
groesult_new2$indicator <- factor(groesult_new2$indicator,
                                 levels = c("SPEARpesticides", "%EPT", "Saprobic\nindex"))

n_stress <- nlevels(groesult_new2$stressor)
full_row <- which(levels(groesult_new2$stressor) == "Best-fit model")

# ---------------------------------------------------------------------------
# Plotting 
# ---------------------------------------------------------------------------

(Fig_multiple_stressors2 <- ggplot(groesult_new2, aes(x = indicator, y = stressor, size = r2, label = labels, colour = Col)) +
  
  # Light horizontal grid lines for each stressor row
  geom_hline(yintercept = seq(1.5, n_stress - 0.5, 1), colour = "grey80", linewidth = 0.25) +
  
  # Black lines separating sections 
  geom_hline(yintercept = c(9.5, 8.5, 0.5), linewidth = 0.5, colour = "#000000") +
  
  # R² bubbles with adjusted transparency
  geom_point(alpha = 0.75, position = position_nudge(y = +0.16)) +
  
  # Text labels (R² + significance stars within the table) 
  geom_text(stat = "identity", 
            size = 4.5,  
            colour = "#000000", 
            family = "Arial",  
            vjust = 1, 
            alpha = 1, 
            position = position_nudge(y = -0.15)) +
  
  # X axis – indicators
  scale_x_discrete(position = "top",
                   labels = c("SPEARpesticides" = expression(SPEAR[pesticides]),
                              "%EPT" = "%EPT",
                              "Saprobic\nindex" = "Saprobic\nindex")) +
  
  # Y axis – stressors
  scale_y_discrete(drop = FALSE, labels = y_labels2) +
  
  # Color and size scales
  scale_colour_identity(guide = "none") +
  scale_size(range = c(3, 14)) +  
  
  # Theme 
  theme_classic(base_size = 12, base_family = "Arial") + 
  theme(
    # Background
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Axes
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Axis text 
    axis.text.x = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 0.5,
                               colour = "#000000",
                               face = "plain"),
    axis.text.y = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 1,
                               colour = "#000000"),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 11, face = "plain"),  
    legend.text = element_text(size = 11), 
    legend.key.size = unit(0.6, "cm"),  
    legend.background = element_blank(),
    
    # Margins - tighter for publication
    plot.margin = unit(c(5, 2, 2, 2), "mm")) +
  
  coord_cartesian(xlim = c(0.5, 3), clip = "off"))

# ---------------------
## Model 3 "TUsum_before_inver"
# Model for SPEARpesticides (speartaxa)
resp <- "speartaxa"

predictor_vars <- c("TUsum_before_inver", 
                    "AgriLand_percent",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow",
                    "TUsum_Trace_elements")

dat3.1 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model3.1 <- lm(as.formula(paste(resp, "~ .")), data = dat3.1)
summary(full_model3.1)

vif3.1 <- car::vif(full_model3.1)
vif3.1

keep_pred3.1 <- names(vif3.1)[vif3.1 <= 2]
keep_pred3.1

dat3.1_red <- dat3.1 |>
  dplyr::select(dplyr::all_of(c(resp, keep_pred3.1)))

full_model3.1_red <- lm(as.formula(paste(resp, "~ .")), data = dat3.1_red)
summary(full_model3.1_red)

final_model3.1 <- MASS::stepAIC(full_model3.1_red, direction = "backward", trace = FALSE)
summary(final_model3.1)

# Model residual check (validation)
par(mfrow = c(2, 2))
plot(final_model3.1)
par(mfrow = c(1, 1))
dev.off()

# Relative importance
rel3.1  <- relaimpo::calc.relimp(final_model3.1, type = "lmg", rela = FALSE)
coef3.1 <- broom::tidy(final_model3.1)

res3.1 <- data.frame(indicator = resp,
                   stressor  = names(rel3.1$lmg),
                   r2        = as.numeric(rel3.1$lmg),
                   estimate  = coef3.1$estimate[match(names(rel3.1$lmg), coef3.1$term)],
                   signifi   = coef3.1$p.value[match(names(rel3.1$lmg), coef3.1$term)])

res3.1 <- rbind(res3.1, data.frame(indicator = resp,
                               stressor  = "Full model",
                               r2        = rel3.1$R2,
                               estimate  = NA,
                               signifi   = NA))

# ---------------------
# Model for %EPT
resp <- "eptProz"

dat3.2 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model3.2 <- lm(as.formula(paste(resp, "~ .")), data = dat3.2)
summary(full_model3.2)

vif3.2 <- car::vif(full_model3.2)
vif3.2 

keep_pred3.2 <- names(vif3.2)[vif3.2 <= 2]  
keep_pred3.2

dat3.2_red <- dat3.2 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred3.2)))

full_model3.2_red <- lm(as.formula(paste(resp, "~ .")), data = dat3.2_red)
summary(full_model3.2_red)

final_model3.2    <- stepAIC(full_model3.2_red, direction = "backward", trace = FALSE)
summary(final_model3.2)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model3.2)
par(mfrow = c(1,1))
dev.off()

rel3.2  <- calc.relimp(final_model3.2, type = "lmg", rela = FALSE)
coef3.2 <- broom::tidy(final_model3.2)

res3.2 <- data.frame(indicator = resp,
                   stressor  = names(rel3.2$lmg),
                   r2        = as.numeric(rel3.2$lmg),
                   estimate  = coef3.2$estimate[match(names(rel3.2$lmg), coef3.2$term)],
                   signifi   = coef3.2$p.value[match(names(rel3.2$lmg), coef3.2$term)])

res3.2 <- rbind(res3.2, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel3.2$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Model for Saprobic index (same steps as for SPEARpesticide)
resp <- "sapr"

dat3.3 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model3.3 <- lm(as.formula(paste(resp, "~ .")), data = dat3.3)
summary(full_model3.3)

vif3.3 <- car::vif(full_model3.3)
vif3.3

keep_pred3.3 <- names(vif3.3)[vif3.3 <= 2] 
keep_pred3.3

dat3.3_red <- dat3.3 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred3.3)))

full_model3.3_red <- lm(as.formula(paste(resp, "~ .")), data = dat3.3_red)
summary(full_model3.3_red)

final_model3.3 <- stepAIC(full_model3.3_red, direction = "backward", trace = FALSE)
summary(final_model3.3)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model3.3)
par(mfrow = c(1,1))
dev.off()

rel3.3  <- calc.relimp(final_model3.3, type = "lmg", rela = FALSE)

coef3.3 <- broom::tidy(final_model3.3)

res3.3 <- data.frame(indicator = resp, 
                   stressor  = names(rel3.3$lmg),
                   r2        = as.numeric(rel3.3$lmg),
                   estimate  = coef3.3$estimate[match(names(rel3.3$lmg), coef3.3$term)],
                   signifi   = coef3.3$p.value[match(names(rel3.3$lmg), coef3.3$term)])

res3.3 <- rbind(res3.3, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel3.3$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Combine model results for all three response varialbes
groesult_simple3 <- dplyr::bind_rows(res3.1, res3.2, res3.3)|>
  dplyr::select(indicator, stressor, signifi, r2, estimate)

# Save results
write_xlsx(groesult_simple3, "Recalculated_metrics/3.3.Final_model_results_TU_before_invertebrates.xlsx")

#----------------------
# Plot
# List of stressors in the plot
stress_levels_simple3 <- c("Best-fit model",
                          "Pesticide\ntoxicity",          # TUsum_Pesticide 
                          "Trace\nelements",              # TUsum_Trace_elements
                          "%Agricultural\nland uses",     # AgriLand_percent 
                          "Hydromorphological\ndegradation", # Morphology 
                          "Bed habitat\nstructure",       # Habitat 
                          "pH",                           # pH
                          "TP",                           # TP
                          "Flow\nvelocity",               # Flow
                          "Stream depth")                 # Stream depth


# Labels for y-axis
y_labels3 <- c("Best-fit model"              = " Best-fit model",
              "Pesticide\ntoxicity"          = "Pesticide\ntoxicity",
              "Trace\nelements"              = "Trace\nelements",
              "%Agricultural\nland uses"     = "%Agricultural\nland uses",
              "Hydromorphological\ndegradation" = "Hydromorphological\ndegradation",
              "Bed habitat\nstructure"       = "Bed habitat\nstructure",
              "pH"                           = "pH",
              "TP"                           = "TP",
              "Flow\nvelocity"               = "Flow\nvelocity",
              "Stream_depth"                 = "Stream depth")

# Table for plotting
groesult_new3 <- groesult_simple3 |> 
  # Significance labels
  dplyr::mutate(signifi = dplyr::case_when(signifi <= 0.001 ~ "***",
                                           signifi <= 0.01  ~ "**",
                                           signifi <= 0.05  ~ "*",
                                           TRUE             ~ NA_character_),

    # Recode stressors for plot labels (save space)
    stressor = dplyr::recode(stressor,
      "Full model"           = "Best-fit model",
      "TUsum_before_inver"   = "Pesticide\ntoxicity",
      "TUsum_Trace_elements" = "Trace\nelements",
      "AgriLand_percent"     = "%Agricultural\nland uses",
      "Morphology"           = "Hydromorphological\ndegradation",
      "Habitat"              = "Bed habitat\nstructure",
      "pH"                   = "pH",
      "TP"                   = "TP",
      "Flow"                 = "Flow\nvelocity",
      "Stream_depth"         = "Stream depth"),

    # Rename & reorder response variables
    indicator = dplyr::recode(indicator,
      "speartaxa" = "SPEARpesticides",
      "sapr"      = "Saprobic\nindex",
      "eptProz"   = "%EPT"),

    # Colours: black for Full model, red for negative estimate, blue for positive
    Col = dplyr::case_when(stressor == "Best-fit model" ~ "black",
      estimate <= 0                      ~ "firebrick3",
      TRUE                               ~ "dodgerblue3"),

    r2     = round(as.numeric(r2), 2),
    labels = dplyr::if_else(!is.na(signifi), paste(r2, signifi), as.character(r2)),
    labels = dplyr::if_else(labels %in% c("NA", "NA NA"), NA_character_, labels)) 

# Stress levels: only keep those that are actually present in the final model
used_levels <- stress_levels_simple3[stress_levels_simple3 %in% groesult_new3$stressor]

groesult_new3$stressor <- factor(groesult_new3$stressor,
                                levels = rev(used_levels))  
groesult_new3$indicator <- factor(groesult_new3$indicator,
                                 levels = c("SPEARpesticides", "%EPT", "Saprobic\nindex"))

n_stress <- nlevels(groesult_new3$stressor)
full_row <- which(levels(groesult_new3$stressor) == "Best-fit model")

# ---------------------------------------------------------------------------
# Plotting 
# ---------------------------------------------------------------------------

(Fig_multiple_stressors3 <- ggplot(groesult_new3, aes(x = indicator, y = stressor, size = r2, label = labels, colour = Col)) +
  
  # Light horizontal grid lines for each stressor row
  geom_hline(yintercept = seq(1.5, n_stress - 0.5, 1), colour = "grey80", linewidth = 0.25) +
  
  # Black lines separating sections 
  geom_hline(yintercept = c(9.5, 8.5, 0.5), linewidth = 0.5, colour = "#000000") +
  
  # R² bubbles with adjusted transparency
  geom_point(alpha = 0.75, position = position_nudge(y = +0.16)) +
  
  # Text labels (R² + significance stars within the table) 
  geom_text(stat = "identity", 
            size = 4.5,  # Increased from 3.0
            colour = "#000000", 
            family = "Arial",  
            vjust = 1, 
            alpha = 1, 
            position = position_nudge(y = -0.15)) +
  
  # X axis – indicators
  scale_x_discrete(position = "top",
                   labels = c("SPEARpesticides" = expression(SPEAR[pesticides]),
                              "%EPT" = "%EPT",
                              "Saprobic\nindex" = "Saprobic\nindex")) +
  
  # Y axis – stressors
  scale_y_discrete(drop = FALSE, labels = y_labels2) +
  
  # Color and size scales
  scale_colour_identity(guide = "none") +
  scale_size(range = c(3, 14)) + 
  
  # Theme 
  theme_classic(base_size = 12, base_family = "Arial") +  
  theme(
    # Background
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Axes
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Axis text 
    axis.text.x = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 0.5,
                               colour = "#000000",
                               face = "plain"),
    axis.text.y = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 1,
                               colour = "#000000"),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 11, face = "plain"),  
    legend.text = element_text(size = 11), 
    legend.key.size = unit(0.6, "cm"),  
    legend.background = element_blank(),
    
    # Margins - tighter for publication
    plot.margin = unit(c(5, 2, 2, 2), "mm")) +
  
  coord_cartesian(xlim = c(0.5, 3), clip = "off"))

# ----------------------------------

# ---------------------
## Model "TUsum_before_inver_extreme_removed"
# Model for SPEARpesticides (speartaxa)
resp <- "speartaxa"

predictor_vars <- c("TUsum_before_inver_extreme_removed", 
                    "AgriLand_percent",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow",
                    "TUsum_Trace_elements")

dat4.1 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model4.1 <- lm(as.formula(paste(resp, "~ .")), data = dat4.1)
summary(full_model4.1)

vif4.1 <- car::vif(full_model4.1)
vif4.1

keep_pred4.1 <- names(vif4.1)[vif4.1 <= 2]
keep_pred4.1

dat4.1_red <- dat4.1 |>
  dplyr::select(dplyr::all_of(c(resp, keep_pred4.1)))

full_model4.1_red <- lm(as.formula(paste(resp, "~ .")), data = dat4.1_red)
summary(full_model4.1_red)

final_model4.1 <- MASS::stepAIC(full_model4.1_red, direction = "backward", trace = FALSE)
summary(final_model4.1)

# Model residual check (validation)
par(mfrow = c(2, 2))
plot(final_model4.1)
par(mfrow = c(1, 1))
dev.off()

# Relative importance
rel4.1  <- relaimpo::calc.relimp(final_model4.1, type = "lmg", rela = FALSE)
coef4.1 <- broom::tidy(final_model4.1)

res4.1 <- data.frame(indicator = resp,
                   stressor  = names(rel4.1$lmg),
                   r2        = as.numeric(rel4.1$lmg),
                   estimate  = coef4.1$estimate[match(names(rel4.1$lmg), coef4.1$term)],
                   signifi   = coef4.1$p.value[match(names(rel4.1$lmg), coef4.1$term)])

res4.1 <- rbind(res4.1, data.frame(indicator = resp,
                               stressor  = "Full model",
                               r2        = rel4.1$R2,
                               estimate  = NA,
                               signifi   = NA))

# ---------------------
# Model for %EPT
resp <- "eptProz"

dat4.2 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model4.2 <- lm(as.formula(paste(resp, "~ .")), data = dat4.2)
summary(full_model4.2)

vif4.2 <- car::vif(full_model4.2)
vif4.2 

keep_pred4.2 <- names(vif4.2)[vif4.2 <= 2]  
keep_pred4.2

dat4.2_red <- dat4.2 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred4.2)))

full_model4.2_red <- lm(as.formula(paste(resp, "~ .")), data = dat4.2_red)
summary(full_model4.2_red)

final_model4.2    <- stepAIC(full_model4.2_red, direction = "backward", trace = FALSE)
summary(final_model4.2)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model4.2)
par(mfrow = c(1,1))
dev.off()

rel4.2  <- calc.relimp(final_model4.2, type = "lmg", rela = FALSE)
coef4.2 <- broom::tidy(final_model4.2)

res4.2 <- data.frame(indicator = resp,
                   stressor  = names(rel4.2$lmg),
                   r2        = as.numeric(rel4.2$lmg),
                   estimate  = coef4.2$estimate[match(names(rel4.2$lmg), coef4.2$term)],
                   signifi   = coef4.2$p.value[match(names(rel4.2$lmg), coef4.2$term)])

res4.2 <- rbind(res4.2, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel4.2$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Model for Saprobic index (same steps as for SPEARpesticide)
resp <- "sapr"

dat4.3 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp, predictor_vars)))

full_model4.3 <- lm(as.formula(paste(resp, "~ .")), data = dat4.3)
summary(full_model4.3)

vif4.3 <- car::vif(full_model4.3)
vif4.3

keep_pred4.3 <- names(vif4.3)[vif4.3 <= 2] 
keep_pred4.3

dat4.3_red <- dat4.3 |> 
  dplyr::select(dplyr::all_of(c(resp, keep_pred4.3)))

full_model4.3_red <- lm(as.formula(paste(resp, "~ .")), data = dat4.3_red)
summary(full_model4.3_red)

final_model4.3 <- stepAIC(full_model4.3_red, direction = "backward", trace = FALSE)
summary(final_model4.3)

# Model residual check (validation)
par(mfrow = c(2,2))
plot(final_model4.3)
par(mfrow = c(1,1))
dev.off()

rel4.3  <- calc.relimp(final_model4.3, type = "lmg", rela = FALSE)

coef4.3 <- broom::tidy(final_model4.3)

res4.3 <- data.frame(indicator = resp, 
                   stressor  = names(rel4.3$lmg),
                   r2        = as.numeric(rel4.3$lmg),
                   estimate  = coef4.3$estimate[match(names(rel4.3$lmg), coef4.3$term)],
                   signifi   = coef4.3$p.value[match(names(rel4.3$lmg), coef4.3$term)])

res4.3 <- rbind(res4.3, data.frame(indicator = resp, 
                               stressor = "Full model",
                               r2 = rel4.3$R2, 
                               estimate = NA, 
                               signifi = NA))

# ---------------------
# Combine model results for all three response varialbes
groesult_simple4 <- dplyr::bind_rows(res4.1, res4.2, res4.3)|>
  dplyr::select(indicator, stressor, signifi, r2, estimate)

# Save results
write_xlsx(groesult_simple4, "Recalculated_metrics/3.4.Final_model_results_TU_before_invertebrates_extreme_removed.xlsx")

#----------------------
# Plot
# List of stressors in the plot
stress_levels_simple4 <- c("Best-fit model",
                          "Pesticide\ntoxicity",          # TUsum_Pesticide 
                          "Trace\nelements",              # TUsum_Trace_elements
                          "%Agricultural\nland uses",     # AgriLand_percent 
                          "Hydromorphological\ndegradation", # Morphology 
                          "Bed habitat\nstructure",       # Habitat 
                          "pH",                           # pH
                          "TP",                           # TP
                          "Flow\nvelocity",               # Flow
                          "Stream depth")                 # Stream depth


# Labels for y-axis
y_labels4 <- c("Best-fit model"              = " Best-fit model",
              "Pesticide\ntoxicity"          = "Pesticide\ntoxicity",
              "Trace\nelements"              = "Trace\nelements",
              "%Agricultural\nland uses"     = "%Agricultural\nland uses",
              "Hydromorphological\ndegradation" = "Hydromorphological\ndegradation",
              "Bed habitat\nstructure"       = "Bed habitat\nstructure",
              "pH"                           = "pH",
              "TP"                           = "TP",
              "Flow\nvelocity"               = "Flow\nvelocity",
              "Stream_depth"                 = "Stream depth")

# Table for plotting
groesult_new4 <- groesult_simple4 |> 
  # Significance labels
  dplyr::mutate(signifi = dplyr::case_when(signifi <= 0.001 ~ "***",
                                           signifi <= 0.01  ~ "**",
                                           signifi <= 0.05  ~ "*",
                                           TRUE             ~ NA_character_),

    # Recode stressors for plot labels 
    stressor = dplyr::recode(stressor,
      "Full model"           = "Best-fit model",
      "TUsum_before_inver_extreme_removed"  = "Pesticide\ntoxicity",
      "TUsum_Trace_elements" = "Trace\nelements",
      "AgriLand_percent"     = "%Agricultural\nland uses",
      "Morphology"           = "Hydromorphological\ndegradation",
      "Habitat"              = "Bed habitat\nstructure",
      "pH"                   = "pH",
      "TP"                   = "TP",
      "Flow"                 = "Flow\nvelocity",
      "Stream_depth"         = "Stream depth"),

    # Rename & reorder response variables
    indicator = dplyr::recode(indicator,
      "speartaxa" = "SPEARpesticides",
      "sapr"      = "Saprobic\nindex",
      "eptProz"   = "%EPT"),

    # Colours: black for Full model, red for negative estimate, blue for positive
    Col = dplyr::case_when(stressor == "Best-fit model" ~ "black",
      estimate <= 0                      ~ "firebrick3",
      TRUE                               ~ "dodgerblue3"),

    r2     = round(as.numeric(r2), 2),
    labels = dplyr::if_else(!is.na(signifi), paste(r2, signifi), as.character(r2)),
    labels = dplyr::if_else(labels %in% c("NA", "NA NA"), NA_character_, labels)) 

# Stress levels: only keep those that are actually present in the final model
used_levels <- stress_levels_simple4[stress_levels_simple4 %in% groesult_new4$stressor]

groesult_new4$stressor <- factor(groesult_new4$stressor,
                                levels = rev(used_levels))  
groesult_new4$indicator <- factor(groesult_new4$indicator,
                                 levels = c("SPEARpesticides", "%EPT", "Saprobic\nindex"))

n_stress <- nlevels(groesult_new4$stressor)
full_row <- which(levels(groesult_new4$stressor) == "Best-fit model")

# ---------------------------------------------------------------------------
# Plotting 
# ---------------------------------------------------------------------------

(Fig_multiple_stressors4 <- ggplot(groesult_new4, aes(x = indicator, y = stressor, size = r2, label = labels, colour = Col)) +
  
  # Light horizontal grid lines for each stressor row
  geom_hline(yintercept = seq(1.5, n_stress - 0.5, 1), colour = "grey80", linewidth = 0.25) +
  
  # Black lines separating sections 
  geom_hline(yintercept = c(9.5, 8.5, 0.5), linewidth = 0.5, colour = "#000000") +
  
  # R² bubbles with adjusted transparency
  geom_point(alpha = 0.75, position = position_nudge(y = +0.16)) +
  
  # Text labels (R² + significance stars within the table) 
  geom_text(stat = "identity", 
            size = 4.5, 
            colour = "#000000", 
            family = "Arial",  
            vjust = 1, 
            alpha = 1, 
            position = position_nudge(y = -0.15)) +
  
  # X axis – indicators
  scale_x_discrete(position = "top",
                   labels = c("SPEARpesticides" = expression(SPEAR[pesticides]),
                              "%EPT" = "%EPT",
                              "Saprobic\nindex" = "Saprobic\nindex")) +
  
  # Y axis – stressors
  scale_y_discrete(drop = FALSE, labels = y_labels3) +
  
  # Color and size scales
  scale_colour_identity(guide = "none") +
  scale_size(range = c(3, 14)) + 
  
  # Theme 
  theme_classic(base_size = 12, base_family = "Arial") +  
  theme(
    # Background
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    
    # Axes
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    # Axis text 
    axis.text.x = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 0.5,
                               colour = "#000000",
                               face = "plain"),
    axis.text.y = element_text(size = 14,  
                               vjust = 0.5, 
                               hjust = 1,
                               colour = "#000000"),
    
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 11, face = "plain"),  
    legend.text = element_text(size = 11), 
    legend.key.size = unit(0.6, "cm"),  
    legend.background = element_blank(),
    
    # Margins - tighter for publication
    plot.margin = unit(c(5, 2, 2, 2), "mm")) +
  
  coord_cartesian(xlim = c(0.5, 3), clip = "off"))

#--------------------------------------------
# Combined plots (Figure 3)
# 1-TUsum_full and 2-TUsum_extreme_removed
combined_plot1_2 <- Fig_multiple_stressors1 + Fig_multiple_stressors2
combined_plot1_2
ggsave("combined_plot1_2.svg", combined_plot1_2, width = 18, height = 12)

# 3-TUsum_before_inver and 4-TUsum_before_inver_extreme_removed
combined_plot3_4 <- Fig_multiple_stressors3 + Fig_multiple_stressors4
combined_plot3_4
ggsave("combined_plot3_4.svg", combined_plot3_4, width = 18, height = 12)

# 1 + 2 + 3 + 4
combined_plot_all <- (Fig_multiple_stressors1 + Fig_multiple_stressors2)/
                     (Fig_multiple_stressors3 + Fig_multiple_stressors4)
combined_plot_all

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("combined_plot_all.svg", combined_plot_all, width = 320, height = 340, units = "mm")
