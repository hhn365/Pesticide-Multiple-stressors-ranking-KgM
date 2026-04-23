## Step 1.Reproduce pesticide metric (TUmax) from raw data
## Include:
# - Import raw data from PANEGEAN dataset Liess et al. 2021b (https://doi.org/10.1594/PANGAEA.931673)
# - Calculate TUmax metric following steps described in the main manuscript Liess et al. 2021a (https://doi.org/10.1016/j.watres.2021.117262)
# - Compare values of reproduced and original (presented in the main manuscript) TUmax metrics for reproducibility 
# - Contact main authors of the KgM manuscript to verify potential mismatch steps in metric calculation approaches
# - Additional: Plot Figure 3A (including the 90% prediction interval, line a95% and b5% benchmarks) to associate reproduced TUmax with SPEARpesticide

# -------------------------

# -------------------------

# Load libraries
library(dplyr)
library(readr) 
library(writexl)
library(readxl)
library(cowplot)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(stabs)

#-------------------------

# -------------------------
## Import pesticide data
d <- read_tsv("Input_data/4_4_Pesticide_TargetAnalytics.txt", locale = locale(encoding = "latin1"),  
                          quote = "\"", trim_ws = TRUE, show_col_types = FALSE)

# Check the dataset
dplyr::glimpse(d)

# -------------------------

# -------------------------
# Pre-process pesticide data table: rename and reformat some variables for individual TU per substance calculation
# Format date and year
d <- d |> 
  rename(siteID = ID,
         conc = concentration_µgL) |> 
  mutate(sampleDate = date,            
         year = format(date, "%Y")) |> 
  mutate(CODE = paste(siteID, year, sep = "_"),              # 'CODE' variable is combination of siteID and year
         labelShort = paste(siteID,                          # 'labelShort' variable is combination of siteID, sampleDate, and method variables
                           format(sampleDate, "%Y%m%d"), 
                           method, sep = "_"))|> 
  select(siteID, labelShort, CODE, sampleDate, year, substance,
    method, conc,                                            # 'conc' is substance concentration above limit of detection (LOD)
    limit_of_quantification_µgL,                        
    concentration_greater_limit_of_quantification_µgL)       # LOQ = limit of quantification LOQ, this variable is provided but not used in the manuscript (i.e., consider conc >= LOD, but accept cases conc < LOQ)

# Check the list of substance names for consistency before joining with LC50 for Toxic unit (TU) calculations:
# Remove the name spinosad_D and only keep spinosad_A as spinosad across samples 
d <- d |> filter(is.na(substance) | substance != "spinosad_D") |>            
          mutate(substance = if_else(substance == "spinosad_A",              
                             "spinosad", substance))

# -------------------------

# -------------------------
# Join pesticide raw data table (d) with the pesticide table having LC50 values 
# Import LC50 values for individual substances
Pesticide_parameters <- read_tsv("Input_data/6_1_Pesticides_data.txt", locale = locale(encoding = "latin1"),  
                                quote = "\"", trim_ws = TRUE, show_col_types = FALSE)

# Verify that 'substance' names match between the two tables before joining the two tables
typeof(d$substance) == typeof(Pesticide_parameters$substance)
param_names_upd <- sort(unique(Pesticide_parameters$substance))
data2_names_upd <- sort(unique(d$substance))
setdiff(param_names_upd, data2_names_upd)

# Join pesticide data (d) and pesticide parameter (Pesticide_parameters) tables
d2 <- d |> 
  left_join(Pesticide_parameters |> 
      select(substance, LC50_invertebrates_µgL),
             by = c("substance"))|>
  select(siteID, labelShort, CODE, year, sampleDate, substance, 
    method, conc, 
    limit_of_quantification_µgL, 
    concentration_greater_limit_of_quantification_µgL, 
    LC50_invertebrates_µgL) |> 
  # Calculate TU for each substance having LC50 values
  mutate(tu = ifelse(is.na(LC50_invertebrates_µgL), 
                     NA_real_, 
                     conc / LC50_invertebrates_µgL)) 

#---------------

#---------------
# Aggregate TU per sample (distinct grab and event samples)
# Define the joining variable labelShort(incorporating siteID, sampleDate, and method) as factor
d2 <- d2 |> 
  mutate(labelShort = as.factor(labelShort))

# Define TUmax per sample at each study site and year
d3 <- d2 |> 
  group_by(labelShort) |> 
  summarise(CODE = first(CODE),          
            sampleDate = first(sampleDate),    
            year = first(year),          
            siteID = first(siteID), 
            method = first(method),
            tuMax = max(tu, na.rm = TRUE),
            tuMaxSub = {i <- which.max(tu)               
                        if (all(is.na(tu))) NA_real_ 
                        else substance[i]},
            Nsub = sum(conc > 0, na.rm = TRUE),
            .groups = "drop")

#---------------

#---------------
# Calculate TUmax per site and year to compare with original TUmax in the KgM manuscript
# Including the filtering of 'outliers'* as defined in the main manuscript
# * Outlier: is defined by a toxicity higher by a factor of a magnitude of 2 (which equals to a factor of 100 in unlog scale)
outlier_threshold = 2 

# TUmax per site and year (grouping variable = CODE)
tu <- d3 |> 
  group_by(CODE) |> 
  arrange(desc(tuMax)) |> 
  dplyr::summarise(siteID = siteID[1],
            year = year[1],
            TUmax_full = log10(max(tuMax, na.rm = T)),
            substance_TUmax_full = as.character(tuMaxSub[1]),
            second_TUmax = log10(tuMax[2]),
            second_TUmaxSub = as.character(tuMaxSub[2]),
            n_samples = n(),
            n_eds_samples = sum(method=="eds"),
            mean_next = mean(log10(tuMax[2:6]), na.rm = T),                                    # mean log TUmax of up to five subsequent samples with lower TUmax values
            max_is_outlier = ifelse((TUmax_full - mean_next)>= outlier_threshold, 1, 0)) |>    # check if difference exceeds the outlier threshold of 2
  as.data.frame()

# If TUmax is outlier, take the second highest TUmax value
tu$TUmax_extreme_removed = ifelse(tu$max_is_outlier==0 | is.na(tu$max_is_outlier), tu$TUmax_full, tu$second_TUmax)
tu$substance_TUmax_extreme_removed = ifelse(tu$max_is_outlier==0 | is.na(tu$max_is_outlier), tu$substance_TUmax_full, tu$second_TUmaxSub)

# Round values below -5 to -5 (assuming values below -5 is likely to not contribute to overall toxicity)
tu$TUmax_full[tu$TUmax_full<(-5)] <- (-5) 
tu$TUmax_extreme_removed[tu$TUmax_extreme_removed<(-5)] <- (-5)  

#---------------

#---------------
# Compare TUmax metric calculations
# Import the list of 101 IDs (including 11 sites monitored in both years) in the main manuscript to filter relevant sites for a comparison
metrics_in_manuscript <- readxl::read_excel("Input_data/1-s2.0-S0043135421004607-mmc2.xlsx", sheet = "Site Parameters") |> 
  select(`Site ID`, Year, `TUmax [log]`, SPEARpesticides) |> 
  distinct(`Site ID`, Year, .keep_all = TRUE)

# Extract unique 101 study site IDs in the manuscript for comparison of the metric reproducibility
ids_101 <- unique(metrics_in_manuscript$'Site ID')  

# Extract TUmax recalculation table and filter for 101 site IDs
tu_exclude_extreme <- tu
tu_exclude_extreme_101sites <- tu_exclude_extreme[tu_exclude_extreme$siteID %in% ids_101, ]  
unique(tu_exclude_extreme_101sites$siteID) 

# Join TUmax recalculated and original for comparison
tu_exclude_extreme_101sites_comparison <- metrics_in_manuscript |> 
  mutate(Year = as.integer(Year)) |> 
  left_join(tu_exclude_extreme_101sites |> mutate(year = as.integer(year)),
            by = c("Site ID" = "siteID", "Year" = "year")) |> 
  select("Site ID", Year, CODE,
         substance_TUmax_full,TUmax_full, 
         second_TUmaxSub, second_TUmax,
         n_samples, n_eds_samples,mean_next, max_is_outlier,
         substance_TUmax_extreme_removed, TUmax_extreme_removed,
         "TUmax [log]")

# Check number of cases with matching values of the TUmax (reproduced) and TUmax [log] (in the manuscript, file "SI_2_Site_Parameters.xlsx")
tu_exclude_extreme_101sites_comparison$Match <- 
  round(tu_exclude_extreme_101sites_comparison$TUmax_extreme_removed, 2) ==
  round(tu_exclude_extreme_101sites_comparison$`TUmax [log]`, 2)

tu_exclude_extreme_101sites_comparison$Match
# Conclusion: 112/112 matching results 

# Export the TUmax
tu_exclude_extreme <- tu_exclude_extreme_101sites_comparison 
write_xlsx(list(tu_max_full = tu, tu_exclude_extreme = tu_exclude_extreme), 
           "Recalculated_metrics/1.1.Reproduced_TUmax_using-raw-data.xlsx")

# -------------------------

# -------------------------
# Part 2. Plot relationship of SPEARpesticide and reproduced TUmax

## Import SPEAR_pesticide data - File '1-s2.0-S0043135421004607-mmc2' (Supporting Information - Liess et al. 2021a)
input_metrics <- metrics_in_manuscript 
names(input_metrics)

## Import reproduced TUmax
TUmax <- tu_exclude_extreme
names(TUmax)

## Extract variables to reproduce linear correlation as in Figure 3A of the main manuscript Liess et al. (2021)
Fig3_metrics <- input_metrics |>
  left_join(TUmax |> 
              select(`Site ID`, Year, TUmax_extreme_removed),
              by = c("Site ID", "Year")) |>
  mutate(TUmax_extreme_removed = round(TUmax_extreme_removed, 2))|>
  select(Site_ID = `Site ID`,
         Year,
         TUmax_log = TUmax_extreme_removed,
         SPEAR_pesticides = SPEARpesticides)

## Aggregate mean values of metrics per Site_ID for 101 sites monitored from 2018 to 2019
Fig3_metrics_agg <- Fig3_metrics |> 
  group_by(Site_ID) |> 
  summarise(mean_TUmax = mean(TUmax_log, na.rm = TRUE),
            mean_SPEAR = mean(SPEAR_pesticides, na.rm = TRUE),
            .groups = "drop")

# -------------------------
# Linear regression model used for Figure 3A
lm_model <- lm(mean_SPEAR ~ mean_TUmax, data = Fig3_metrics_agg)
# Model statistics
(summary_lm <- summary(lm_model))        # Note: Significant negative relationships
(r_squared <- summary_lm$r.squared)      # R2 = 0.4247 ~ 0.43, matches the reported value in Fig 3A
(p_value <- summary_lm$coefficients[2, 4])
(intercept <- coef(lm_model)[1])
(slope <- coef(lm_model)[2])
(sigma <- summary_lm$sigma)

# -------------------------
## Figure 3A plotting 
# a95% benchmark - as specified in Fig. 3. caption
good_moderate_benchmark <- 0.6  
(a95_benchmark <- good_moderate_benchmark - 1.645 * sigma)    # a95% = 0.345

# b5%
(b5_threshold <- (good_moderate_benchmark - intercept)/slope) # b5% = -3.27

# 90% prediction interval
pred_data <- data.frame(mean_TUmax = seq(min(Fig3_metrics_agg$mean_TUmax), max(Fig3_metrics_agg$mean_TUmax), length.out = 100))
pred <- predict(lm_model, newdata = pred_data, interval = "prediction", level = 0.90)
pred_data <- cbind(pred_data, fit = pred[,1], lo = pred[,2], up = pred[,3])

# Plot Figure 3A: SPEAR_pesticide vs TUmax(log)
plot_FigA1 <- ggplot() +
  geom_point(data = Fig3_metrics_agg, aes(x = mean_TUmax, y = mean_SPEAR), size = 2, color = "black") +
  geom_ribbon(data = pred_data, aes(x = mean_TUmax, ymin = lo, ymax = up), fill = "#87CEEB", alpha = 0.2) +
  
  geom_line(data = pred_data, aes(x = mean_TUmax, y = fit), color = "#4A90E2", linewidth = 1) +
  geom_vline(xintercept = -3.27, linetype = "dotted", color = "red") +
  geom_hline(yintercept = a95_benchmark, linetype = "dotted", color = "red") +
  geom_hline(yintercept = c(0.2, 0.4, 0.6, 0.8), color = "lightgrey", linewidth = 0.3) +

  annotate("text", x = max(Fig3_metrics_agg$mean_TUmax) + 0.2, y = 1.05, 
           label = sprintf("R² = %.2f\np < 0.001", ceiling(r_squared * 100) / 100), hjust = 1, vjust = 1, size = 4) +
  annotate("text", x = max(Fig3_metrics_agg$mean_TUmax) + 0.2, y = a95_benchmark + 0.05, 
           label = "a95%", size = 4, color = "red") +
  annotate("text", x = -3.27 + 0.2, y = 0.01, 
           label = "b5%", size = 4, color = "red") +
  
  labs(x = "Estimated pesticide toxicity", y = "SPEARpesticides", 
       title = "Field-based adaptive cause-effect relationship for pesticide") +
  
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(breaks = seq(-5, 1, by = 1)) +
  coord_cartesian(xlim = c(-5, 1), ylim = c(0, NA)) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),  
        panel.grid.major.y = element_blank(),  
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(0.2, "cm"),  
        plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"))

plot_FigA1

# Export the reanalysis plot comparing TUmax and SPEARpesticide
# Export as *.svg
ggsave("Appendix_Figure_A1.svg", plot = plot_FigA1, device = "svg", width = 8, height = 6)

#-----------------------------------------

#-----------------------------------------
# Conclusion: 
# Figure 3A is reproducible using the reproduced TUmax 
# and associating with the SPEARpesticide from the Supplementary data '1-s2.0-S0043135421004607-mmc2' (Liess et al. 2021) 
