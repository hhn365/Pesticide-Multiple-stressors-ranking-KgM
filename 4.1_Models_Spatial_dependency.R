## Step 4. Spatial dependency check

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
library(sf)
library(spdep)  
library(gstat)

# -------------------------

# -------------------------
## Import data
# Site Longitude and Latitude
sites <- read.table("Input_data/1_1_Sites_MetaData.txt", header = TRUE, sep = "\t")

# Biotic and abiotic variables
All_metrics <- read_excel("Recalculated_metrics/3.All_biotic_abiotic_metrics_mean_aggregation.xlsx")
names(All_metrics)

# -------------------------

# -------------------------
## Join tables and preprocess data for spatial dependency analyses
# IDs of study sites
unique_IDs <- unique(All_metrics$ID)
unique_IDs

# Filter relevant sites
sites <- sites |> 
  filter(ID %in% unique_IDs) |> 
  distinct(ID, Latitude, Longitude) |> 
  select(ID, Latitude, Longitude)

All_metrics_update <- All_metrics |> 
  left_join(sites, by = "ID", relationship = "many-to-one") |>            
  relocate(Latitude, Longitude, .after = ID)

# Preprocess data for models
All_metrics_mean <- All_metrics_update |>
  mutate(Flow = if_else(is.na(Flow), mean(Flow, na.rm = TRUE), Flow),  # impute one NA in Flow values
    TN = log10(TN),                                                    # log transform nutrients for linearity relationship
    TP = log10(TP)) |>
  dplyr::select(-year) |>                             # exclude year after averaging, and TUmax that is highly correlated to TUsum
  group_by(ID) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  rename(Pesticide = TUsum_before_inver_extreme_removed,               # Use one of four pesticide mixture estimates for the test
         Trace_elements = TUsum_Trace_elements,
         Agriculture = AgriLand_percent) |>
  arrange(as.numeric(sub("S", "", ID)))

# -------------------------

# -------------------------
# Check spatial dependency - Moran I and variogram residuals
# Create  a spatial object 
dat_sf <- All_metrics_mean |> 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)  # Geographic coordinates EPSG:4326

# Convert lat/long to UTM Zone 32N(in meters)
dat_utm <- st_transform(dat_sf, 25832)

# Extract numeric X-Y coordinates as a matrix
xy <- st_coordinates(dat_utm)

# Calculate Euclidean distance between all site pairs in a distance matrix
dmat <- as.matrix(dist(xy))

# Identify the max range of closest neighbour distances of across sites
max_min_dist <- max(apply(dmat, 1, function(x) min(x[x > 0])))

# Define spatial neighbour and convert to spatial weights
nb <- dnearneigh(xy, d1 = 0, d2 = max_min_dist * 1.05)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
                               
# Moran I on a response variable
moran.test(dat_utm$speartaxa, lw, zero.policy = TRUE)  # p-value > 0.05, not significant
moran.test(dat_utm$eptProz, lw, zero.policy = TRUE)    # p-value > 0.05, not significant
moran.test(dat_utm$sapr, lw, zero.policy = TRUE)       # p-value > 0.05, not significant

# Permutation-based Moran I - similar insignificant p-values as above
moran.mc(dat_utm$speartaxa, lw, nsim = 999, zero.policy = TRUE)   # p-value > 0.05, not significant
moran.mc(dat_utm$eptProz, lw, nsim = 999, zero.policy = TRUE)     # p-value > 0.05, not significant
moran.mc(dat_utm$sapr, lw, nsim = 999, zero.policy = TRUE)        # p-value > 0.05, not significant

# -------------------------

# -------------------------
# Check spatial dependency on residuals of final models (test with the final model in script 3.1) using Moran I and Variogram

## SPEARpesticide model => mo spatial pattern
# Residual Moran I
m0_spear <- lm(speartaxa ~ Pesticide + Agriculture + Morphology + Habitat + TP, data = dat_utm)
moran.test(residuals(m0_spear), lw, zero.policy = TRUE)                # p = 0.14, not significant

# Residual variogram
dat_utm$resid_m0_spear <- residuals(m0_spear)
v_res_spear <- variogram(resid_m0_spear ~ 1, data = dat_utm)
plot(v_res_spear)    # no distance-dependent data structure pattern 

# -------------------------
## %EPT model => no spatial pattern
# Residual Moran I
m0_EPT <- lm(eptProz ~ Pesticide + Trace_elements + Agriculture + Morphology + Habitat + pH + TP + Flow, data = dat_utm)
moran.test(residuals(m0_EPT), lw, zero.policy = TRUE)                # p = 0.18, not significant

# Residual variogram
dat_utm$resid_m0_EPT <- residuals(m0_EPT)
v_res_EPT <- variogram(resid_m0_EPT ~ 1, data = dat_utm)
plot(v_res_EPT)  

# -------------------------
## Saprobic index model - mo spatial pattern
# Residual Moran I
m0_sapr <- lm(sapr ~ Pesticide + Agriculture + TP + Flow, data = dat_utm)
moran.test(residuals(m0_sapr), lw, zero.policy = TRUE)                # p = 0.75, not significant

# Residual variogram
dat_utm$resid_m0_sapr <- residuals(m0_sapr)
v_res_sapr <- variogram(resid_m0_sapr ~ 1, data = dat_utm)
plot(v_res_sapr)  

# Note: Results suggest no spatial dependency in data structure and no need for use of a spatial model

