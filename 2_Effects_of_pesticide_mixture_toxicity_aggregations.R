## Step 2. Different data aggregations of pesticide toxicity estimates
# Aim: Analyze the changes in the relationship of pesticide toxicity and ecological responses (using the SPEARpesticides metric) 
# Include comparisons of:
# 2.1. Reproduced TUsum (log) vs TUmax (log) based on full raw data 
# 2.2. Reproduced TUsum (log) separated for grab and event sampling types versus both sampling types (i.e. full data)
# 2.3. Reproduced TUsum (log) excluding high toxicity peaks and with full dataset
# 2.4. Reproduced TUsum (log) for the full time series vs only before macroinvertebrate sampling, and additionally with full dataset vs excluding high peaks (from 2.3.)
# Note: Reuse the SPEAR_pesticides as provided in Liess et al. 2021 for comparisons with reproduced pesticide toxicity metrics

# -------------------------

# -------------------------
## Import packages
library(here)
library(readr)     
library(openxlsx)
library(readxl)
library(writexl)
library(tidyverse)
library(ggplot2)
library(lubridate) 
library(patchwork) 
library(tibble)
library(grid)

# -------------------------

# -------------------------
# Import metrics in the main manuscript to extract biotic metrics and site IDs for a comparison
metrics_in_manuscript <- readxl::read_excel("Input_data/1-s2.0-S0043135421004607-mmc2.xlsx", sheet = "Site Parameters") |> 
  select(1, 2, 22) |>                     # Include: Site.ID, Year, and SPEAR_pesticide
  distinct(`Site ID`, Year, .keep_all = TRUE)|>            
  mutate(Year = as.character(Year)) |> 
  rename(ID = `Site ID`,
    year = 'Year',
    SPEAR_pesticides = SPEARpesticides)

# Extract 101 unique study site IDs 
unique_IDs <- unique(metrics_in_manuscript$ID)  
unique_IDs

# -------------------------

# -------------------------
# Join pesticide data and parameters for calculating pesticide toxicity unit (TU) per substance (repeat procedure as in Step 1)
# Import pesticide data
Pesticide_data <- read_tsv("Input_data/4_4_Pesticide_TargetAnalytics.txt", locale = locale(encoding = "latin1"),  
                          quote = "\"", trim_ws = TRUE, show_col_types = FALSE)

names(Pesticide_data) 

# Subset data for 101 sampling sites IDs reported in the main manuscript
Pesticide_data_subset <- Pesticide_data[Pesticide_data$ID %in% unique_IDs, ]

# Format date and year
Pesticide_data2 <- Pesticide_data_subset |> 
  rename(conc = concentration_µgL) |> 
  mutate(sampleDate = as.Date(date),
         year = format(sampleDate, "%Y")) |> 
  mutate(CODE = paste(ID, year, sep = "_"),         # 'CODE' variable is the combination of siteID and year variables
         labelShort = paste(ID,                     # 'labelShort' variable is the combination of siteID, sampleDate, and method variables
                           format(sampleDate, "%Y%m%d"), 
                           method, sep = "_"))|> 
  select(ID, labelShort, CODE, sampleDate, year, substance,
    method, conc,
    limit_of_quantification_µgL,                        # 'conc' is substance concentration above limit of detection (LOD)
    concentration_greater_limit_of_quantification_µgL)  # import limit of quantification LOQ, but do not use in the manuscript 
                                                        # (i.e., method adapted in consultation with coauthors of the main manuscript, which considers the criteria conc >= Limit Of Detection (conc is not NA), but preserve cases conc < LOQ)
  
# Check substance name list:
# Remove substance spinosad_D and rename spinosad_A as spinosad across samples (method adapted in communication with coauthors of the main manuscript for reproducibility) 
Pesticide_data2 <- Pesticide_data2 |> 
  filter(is.na(substance) | substance != "spinosad_D") |>            
          mutate(substance = if_else(substance == "spinosad_A",              
                             "spinosad", substance))

# -------------------------
# Join pesticide raw data table with the table with LC50 values 
# Import LC50 values for individual substances
Pesticide_parameters <- read_tsv("Input_data/6_1_Pesticides_data.txt", locale = locale(encoding = "latin1"),  
                                quote = "\"", trim_ws = TRUE, show_col_types = FALSE)

# Verify that 'substance' names match between the two tables
typeof(Pesticide_data2$substance) == typeof(Pesticide_parameters$substance)
param_names_upd <- sort(unique(Pesticide_parameters$substance))
data2_names_upd <- sort(unique(Pesticide_data2$substance))
setdiff(param_names_upd, data2_names_upd)

# Join pesticide data (Pesticide_data2) and pesticide parameter (Pesticide_parameters) tables
# And calculate TU values for individual substances
Pesticide_joined <- Pesticide_data2 |> 
  left_join(Pesticide_parameters |> 
      select(substance, LC50_invertebrates_µgL),
    by = c("substance"))|>
  select(ID, labelShort, CODE, year, sampleDate, substance, 
    method, conc, 
    limit_of_quantification_µgL, concentration_greater_limit_of_quantification_µgL, 
    LC50_invertebrates_µgL) |> 
  mutate(labelShort = as.factor(labelShort))  |> 
  # Calculate TU for each substance
  mutate(tu = ifelse(is.na(LC50_invertebrates_µgL), NA_real_,conc / LC50_invertebrates_µgL))      

# -------------------------

# -------------------------
## 2.1. Compare TUsum versus TUmax (full dataset)

## TUsum
# First calculate sum of TU substance per sampling date, then max per year for each site ID 
# this follows the method description of TUsum of trace elements in the main manuscript Liess et al. (2021)

# Sum of TUsum per sampling date
TU_sum_calc <- Pesticide_joined |> 
  group_by(labelShort) |>
  summarise(ID = first(ID),
            CODE = first(CODE),
            year = first(year),
            sampleDate = first(sampleDate),
            method = first(method),
            TUsum_sample  = sum(tu, na.rm = TRUE), 
            .groups = "drop") 
  
# Max of TUsum of each year and site ID 
TU_sum_calc2 <- TU_sum_calc |>
  group_by(CODE) |>
  summarise(ID = first(ID),
            year = first(year),
            TUsum_full  = log10(max(TUsum_sample,  na.rm = TRUE)), .groups = "drop")

# Round values below -5 to -5 (assuming values below -5 is likely do not contributing to the overall toxicity)
TU_sum_calc2$TUsum_full[TU_sum_calc2$TUsum_full<(-5)] <- (-5) 

#---------------
## TUmax
# Aggregate TUmax per sampling date (using 'labelShort' == year x date x method), 
# keep separated between grab and event sampling types
TU_max_calc <- Pesticide_joined |> 
  group_by(labelShort) |> 
  summarise(CODE = first(CODE),          
            sampleDate = first(sampleDate),    
            year = first(year),          
            ID = first(ID), 
            method = first(method),
            tuMax = max(tu, na.rm = TRUE),
            tuMaxSub = {i <- which.max(tu)              
                        if (all(is.na(tu))) NA_real_ 
                        else substance[i]},
            Nsub = sum(conc > 0, na.rm = TRUE),
            .groups = "drop")

# Calculate TUmax per site and year 
TU_max_calc2 <- TU_max_calc |> 
  group_by(CODE) |> 
  arrange(desc(tuMax)) |> 
  dplyr::summarise(ID = ID[1],
            year = year[1],
            TUmax_full = log10(max(tuMax, na.rm = T)),
            substance_TUmax_full = as.character(tuMaxSub[1])) |>  
  as.data.frame()

# Round values below -5 to -5 (assuming values below -5 is likely do not contributing to the overall toxicity)
TU_max_calc2$TUmax_full[TU_max_calc2$TUmax_full<(-5)] <- (-5) 

# Join tables of TUsum and TUmax with SPEARpesticide table, and save for use in later calculations
TUmax_sum_join <- TU_max_calc2 |>
  mutate(year = as.integer(year)) |>
  left_join(TU_sum_calc2 |> 
              transmute(ID, year = as.integer(year), TUsum_full),
    by = c("ID", "year"))

TUmax_sum_join <- metrics_in_manuscript |>
  left_join(select(TU_sum_calc2, - CODE), by = c("ID", "year")) %>%
  left_join(select(TU_max_calc2, - CODE), by = c("ID", "year"))

write_xlsx(TUmax_sum_join, here("Recalculated_metrics", "2.1.Pesticide_TU_max_sum_full_comparisons.xlsx"))

# -------------------------
# Aggregate mean values to get metrics for 101 sites
TU_max_sum_comparison_years_mean <- TUmax_sum_join |> 
  select(-substance_TUmax_full) |> 
  group_by(ID) |> 
  summarise(mean_SPEAR      = round(mean(SPEAR_pesticides, na.rm = TRUE), 2),
            mean_TUsum_full = round(mean(TUsum_full, na.rm = TRUE), 2),
            mean_TUmax_full = round(mean(TUmax_full, na.rm = TRUE), 2),
            .groups = "drop")

summary(TU_max_sum_comparison_years_mean)

# -------------------------
# Plot TUsum and TUmax values (Figure 2a)
dist_TU_max_sum_full_table <- TU_max_sum_comparison_years_mean |>
  select(ID, mean_TUsum_full, mean_TUmax_full) |>
  pivot_longer(cols = c(mean_TUsum_full, mean_TUmax_full),
               names_to = "TU_type", values_to = "TU_value") |> 
  mutate(TU_type = factor(TU_type, levels = c("mean_TUmax_full", "mean_TUsum_full")))

TU_violin_box_plot <- ggplot(dist_TU_max_sum_full_table,
  aes(x = TU_value, y = TU_type, fill = TU_type)) +
  geom_violin(trim = TRUE, width = 0.85, color = NA, alpha = 0.9, na.rm = TRUE) +
  geom_boxplot(width = 0.16, fill = "white", color = "black", linewidth = 0.5, outlier.shape = NA, na.rm = TRUE) +
  scale_fill_manual(values = c(mean_TUmax_full = "#4F8FC6", mean_TUsum_full = "#BDBDBD")) +
  scale_y_discrete(labels = c(mean_TUmax_full = "TUmax", mean_TUsum_full = "TUsum")) +
  scale_x_continuous(breaks = seq(-5, 1, by = 1)) +
  labs(x = "Estimated pesticide toxicity", y = NULL) +
  theme_classic(base_size = 15) +
  theme(text = element_text(),
    axis.text = element_text(size = rel(1.1), color = "black"),
    axis.title = element_text(size = 16),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"))

TU_violin_box_plot

#-----------------------------

#-----------------------------
## 2.2.Compare reproduced TUsum for grab, event, and both sample types

# TUsum full data
TU_sum_method_mixed <- TUmax_sum_join |>
  mutate(year = as.character(year)) |>
  select(ID, year, TUsum_full)

# Calculate TUsum for event (eds) and grab methods separately 
TU_sum_method_separated <- TU_sum_calc |>
  filter(!is.na(TUsum_sample), 
         TUsum_sample > 0,
         method %in% c("grab","eds")) |>
  group_by(ID, year, method) |>
  summarise(TUsum = log10(max(TUsum_sample, na.rm = TRUE)), 
            .groups = "drop") |>
  pivot_wider(names_from  = method, 
              values_from = TUsum, 
              names_glue  = "TUsum_{method}")|> 
  mutate(TUsum_grab = ifelse(TUsum_grab < -5, -5, TUsum_grab),
    TUsum_eds  = ifelse(TUsum_eds  < -5, -5, TUsum_eds))

# Join three TUsum metrics
TU_sum_joined <- full_join(TU_sum_method_mixed, TU_sum_method_separated, by = c("ID","year")) |> 
                            select(ID, year, TUsum_full, TUsum_grab, TUsum_eds)

# Add SPEARpesticide to the joined TUsum table
Pesticide_metrics_recalc_sum_grab_eds_separated <- list(metrics_in_manuscript, TU_sum_joined) |> 
                                    reduce(full_join, by = c("ID", "year"))|> 
                                    select(ID, year, SPEAR_pesticides, 
                                           TUsum_full, TUsum_grab, TUsum_eds) |> 
                                    filter(!(is.na(SPEAR_pesticides))) 

summary(Pesticide_metrics_recalc_sum_grab_eds_separated)

# Export the table of main recalculated metrics
write_xlsx(Pesticide_metrics_recalc_sum_grab_eds_separated, here("Recalculated_metrics", "2.2.Pesticide_TUsum_grab_EDS_comparisons.xlsx"))

# Aggregate mean for 101 site IDs
Pesticide_metrics_recalc_sum_grab_eds <- Pesticide_metrics_recalc_sum_grab_eds_separated |>
  group_by(ID) |>
  summarise(mean_SPEAR      = round(mean(SPEAR_pesticides, na.rm = TRUE), 2),
            mean_TUsum_full = round(mean(TUsum_full,       na.rm = TRUE), 2),
            mean_TUsum_grab = round(mean(TUsum_grab,       na.rm = TRUE), 2),
            mean_TUsum_eds  = round(mean(TUsum_eds,        na.rm = TRUE), 2), 
            .groups = "drop")

summary(Pesticide_metrics_recalc_sum_grab_eds)

## Plot data ranges of TUsum based on both and grab and event (eds) sampling separately (Figure 2b)
TUsum_different_methods_plot <- Pesticide_metrics_recalc_sum_grab_eds |>
  select(ID, mean_TUsum_full, mean_TUsum_grab, mean_TUsum_eds) |>
  pivot_longer(cols = c(mean_TUsum_full, mean_TUsum_grab, mean_TUsum_eds),
               names_to = "metric", values_to = "value") |>
  mutate(metric = factor(metric,
                         levels = c("mean_TUsum_grab","mean_TUsum_eds","mean_TUsum_full"),
                         labels = c("Only grab","Only event","Both"))) |>
  filter(is.finite(value))

TUsum_violin_box_plot2 <- ggplot(TUsum_different_methods_plot,
                                 aes(x = value, y = metric, fill = metric)) +
  geom_violin(trim = TRUE, width = 0.85, color = NA, alpha = 0.9, na.rm = TRUE) +
  geom_boxplot(width = 0.16, fill = "white", color = "black", linewidth = 0.5, outlier.shape = NA, na.rm = TRUE) +
  # colors 
  scale_fill_manual(values = c("Only grab"  = "#A8D5A5", "Only event" = "#31a354", "Both" = "#BDBDBD")) +
  scale_x_continuous(breaks = seq(-5, 0, by = 1), limits = c(-5, 0)) +
  labs(x = "Estimated pesticide toxicity", y = NULL) +
  theme_classic(base_size = 15) +
  theme(text = element_text(), 
        axis.text = element_text(size = rel(1.1), color = "black"),
        axis.title = element_text(size = 16),
        legend.position = "none",
        axis.line = element_line(linewidth = 0.8),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"))

TUsum_violin_box_plot2

# -----------------------------------------
# Count cases with TUsum values higher, equal, and smaller between grab sampling (101 values) and event (eds) sampling types
TUsum_per_method_counts <- Pesticide_metrics_recalc_sum_grab_eds |>
  mutate(comparison = case_when(mean_TUsum_grab >  mean_TUsum_eds ~ "Grab > EDS",
                                mean_TUsum_grab <  mean_TUsum_eds ~ "Grab < EDS",
                                mean_TUsum_grab == mean_TUsum_eds ~ "Grab = EDS")) |>
  count(comparison)

TUsum_per_method_counts  
# Note: Similar results to TUmax, with 29 cases grab > eds and 63 cases eds > grab 
#       9 NAs are the 9/101 NA's of mean_TUsum_eds

#------------------------------------

#------------------------------------
## 2.3. TUsum (log) from raw data excluding high peaks versus full dataset
# Two cases: full dataset and filtering of 'outliers'*: 
# * An 'outlier' is defined by a toxicity higher by a factor of a magnitude of 2 (which equals to a factor of 100 in unlog scale). 
outlier_threshold = 2 

TU_sum_calc2 <- TU_sum_calc |> 
  group_by(CODE) |> 
  arrange(desc(TUsum_sample)) |> 
  dplyr::summarise(ID = ID[1],
            year = year[1],
            TUsum_full = log10(max(TUsum_sample, na.rm = TRUE)),
            second_TUsum = log10(TUsum_sample[2]),
            n_samples = n(),
            n_eds_samples = sum(method=="eds"),
            mean_next = mean(log10(TUsum_sample[2:6]), na.rm = TRUE),                                 # mean log TUsum of up to five subsequent samples with lower values
            max_is_outlier = ifelse((TUsum_full - mean_next)>=outlier_threshold, 1, 0),
            .groups = "drop") |>                                                                      # check if difference exceeds the outlier threshold of 2
  as.data.frame()

# If TUsum is outlier, take the second highest TUsum value within the site x year
TU_sum_calc2$TUsum_extreme_removed = ifelse(
  TU_sum_calc2$max_is_outlier==0 | is.na(TU_sum_calc2$max_is_outlier), 
  TU_sum_calc2$TUsum_full, 
  TU_sum_calc2$second_TUsum)

# Round values below -5 to -5 
TU_sum_calc2$TUsum_full[TU_sum_calc2$TUsum_full<(-5)] <- (-5) 
TU_sum_calc2$TUsum_extreme_removed[TU_sum_calc2$TUsum_extreme_removed<(-5)] <- (-5) 

# -------------------------------------
# Join with SPEARpesticide for linear regression modelling
TUsum_high_peak_removed <- metrics_in_manuscript |>
  left_join(TU_sum_calc2 |> select(ID, year, CODE, TUsum_full, TUsum_extreme_removed),
    by = c("ID", "year"))

# Export
write_xlsx(TUsum_high_peak_removed, here("Recalculated_metrics", "2.3.Pesticide_TUsum_full_vs_High_peak_removed.xlsx"))

# Aggregate mean of two years
TUsum_high_peak_removed_years_mean <- TUsum_high_peak_removed |>
  group_by(ID) |>
  summarise(mean_SPEAR = round(mean(SPEAR_pesticides, na.rm = TRUE), 2),
    mean_TUsum_full = round(mean(TUsum_full, na.rm = TRUE), 2),
    mean_TUsum_extreme_removed = round(mean(TUsum_extreme_removed, na.rm = TRUE), 2),
    .groups = "drop")

# Note: The plot of TUsum metrics in 2.3. is joined with TUsum metrics in the subsequent section 2.4.

#------------------------------------

#------------------------------------
## 2.4. Reproduced TUsum (log) for the full time series vs only before macroinvertebrate sampling, 
#  and additionally with full dataset vs excluding high peaks (from 2.3.)
# Import pesticide and macroinvertebrate sampling dataset from Liess et al. 2021 - PANGEAN
# Pesticides data
Pesticides_date <- Pesticide_data|> 
  mutate(date_Pesticides = as.Date(date)) |>        
  mutate(Year = year(date_Pesticides),
         Day = yday(date_Pesticides)) 

# Macroinvertebrates data
Invertebrates <- read.table("Input_data/5_1_Invertebrate_Abundances.txt", 
                            header = TRUE, sep = "\t")

Invertebrates_date <- Invertebrates |>
  mutate(date_Invertebrates = as.Date(date)) |>  
  mutate(Year = year(date_Invertebrates),
         Day = yday(date_Invertebrates)) |> 
  select(ID, date_Invertebrates, Year, Day) |>
  distinct()

# Check sampling IDs of two sampling types
unique(Pesticides_date$ID)    # 124 sites
unique(Invertebrates_date$ID) # 122 sites - missing two sites S122 and S124 compared to pesticide sampling,
setdiff(unique(Pesticides_date$ID), unique(Invertebrates_date$ID)) # both sites S122 and S124 do not belong to the list of 112 sites in the main manuscript

## Reduce both Pesticide and macroinvertebrate datasets to 101 sites included in the main manuscript
# Subset the pesticide and macroinvertebrate data to 101 study sites
Pesticides_date_subset <- Pesticides_date[Pesticides_date$ID %in% unique_IDs, ]
Invertebrates_date_subset <- Invertebrates_date[Invertebrates_date$ID %in% unique_IDs, ]

(Pest_IDs_compare <- unique(Pesticides_date_subset$ID))    # 101 sites
(Inver_IDs_compare <- unique(Invertebrates_date_subset$ID)) # 101 sites

setdiff(Pest_IDs_compare, Inver_IDs_compare)
# Note: Both pesticides and macroinvertebrates were sampled at all 101 sites of focus in the main manuscript

# -------------------------
# Plot date to compare the distribution of pesticide and macroinvertebrate sampling dates
# Create a table with Sample_date (1:365) and count IDs of pesticide and macroinvertebrate sampling over a year (not distinct between 2018 and 2019)
Pesticides_counts <- Pesticides_date_subset |> 
  group_by(Day) |> summarise(n_IDs_Pesticides = n_distinct(ID)) 

Invertebrates_counts <- Invertebrates_date_subset |> 
  group_by(Day) |> summarise(n_IDs_Invertebrates = n_distinct(ID)) 

Pest_Inver_counts <- full_join(Pesticides_counts, Invertebrates_counts, by = "Day")

full_days <- tibble(Sample_date = 1: 365)

Pest_Inver_Sampling <- full_days |> 
  left_join(Pest_Inver_counts, by = c("Sample_date" = "Day")) |> 
  replace_na(list(n_IDs_Pesticides = 0, n_IDs_Invertebrates = 0)) |> 
  mutate(n_IDs_Invertebrates_neg = - n_IDs_Invertebrates)

# Plot the sampling dates
Sampling_date_plot <- ggplot(Pest_Inver_Sampling, aes(x = Sample_date)) +
  # Pesticide bars above 
  geom_col(aes(y = n_IDs_Pesticides), fill = "#E63946", alpha = 0.85, width = 1) +
  # Invertebrate bars below 
  geom_col(aes(y = n_IDs_Invertebrates_neg), fill = "#457B9D", alpha = 0.85, width = 1) +
  # Baseline with thinner line
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  # Updated labels
  labs(x = "Sampling days (1-365)",
    y = "Number of sampling sites") +
  # Theme
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14, hjust = 0),
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 15, 10, 10)) +
  # Annotation (for clarity)
  annotate("text", x = 30, y = max(Pest_Inver_Sampling$n_IDs_Pesticides) * 0.85,
           label = "Pesticides", color = "#E63946", fontface = "bold", hjust = 0) +
  annotate("text", x = 30, y = min(Pest_Inver_Sampling$n_IDs_Invertebrates_neg) * 0.85,
           label = "Macroinvertebrates", color = "#457B9D", fontface = "bold", hjust = 0) +
  # Add breaks for x and y axis
  scale_x_continuous(breaks = seq(0, 365, by = 50), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

Sampling_date_plot

# ggsave("Sampling_date_plot.svg", Sampling_date_plot, width = 8.8, height = 4)

# -------------------------
# Filter pesticides to only include samples that have sampling dates on or before macroinvertebrate sampling dates
Pesticides_valid <- Pesticides_date_subset |> 
  inner_join(Invertebrates_date_subset, by = "ID", relationship = "many-to-many") |> 
  filter(date_Pesticides <= date_Invertebrates) 

unique(Pesticides_valid$ID)  # Note: valid rows of sampling observations decreased from 91,103 to 75,735,
                             #       but all 101 sites (S1 to S86, S104 to S118) remain valid for the study

write_xlsx(Pesticides_valid, here("Recalculated_metrics","2.4.Pesticides_Ecological_date_match.xlsx"))

# -------------------------
## Repeat steps in 2.1 to calculate toxic unit (TU) per sample, with data from the variable 'Pesticides_valid'
Pesticide_data_filtered <- Pesticides_valid |> 
  mutate(sampleDate = as.Date(date),
         year = format(date, "%Y")) |> 
  rename(conc = concentration_µgL) |> 
  mutate(CODE = paste(ID, year, sep = "_"),         # 'CODE' variable is combined from siteID and year variables
         labelShort = paste(ID,                     # 'labelShort' variable is combined from siteID, sampleDate, and method variables
                           format(sampleDate, "%Y%m%d"), 
                           method, sep = "_"))|> 
  select(ID, labelShort, CODE, sampleDate, year, substance,
    method, conc,
    limit_of_quantification_µgL,                        # 'conc' is substance concentration above limit of detection (LOD)
    concentration_greater_limit_of_quantification_µgL)  # import limit of quantification LOQ but do not use in the manuscript (i.e., method adapted in consultation with authors of the main manuscript, which considers the criteria conc >= Limit Of Detection (conc is not NA), but preserve cases conc < LOQ)

# Remove substance spinosad_D and only keep spinosad_A as spinosad across samples (method adapted in consultation with authors of the main manuscript) 
Pesticide_data3 <- Pesticide_data_filtered |> 
  filter(is.na(substance) | substance != "spinosad_D") |>                    # remove spinosad_D
          mutate(substance = if_else(substance == "spinosad_A",              # keep spinosad_A and change the name to spinosad
                             "spinosad", substance))

# Verify that 'substance' names match between the two tables
typeof(Pesticide_data3$substance) == typeof(Pesticide_parameters$substance)
param3_names_upd <- sort(unique(Pesticide_parameters$substance))
data3_names_upd <- sort(unique(Pesticide_data3$substance))
setdiff(param3_names_upd, data3_names_upd)

# Join pesticide data (Pesticide_data3) and pesticide parameter (Pesticide_parameters) tables
# And calculate TU values for individual substances
Pesticide_joined3 <- Pesticide_data3 |> 
  left_join(Pesticide_parameters |> 
      select(substance, LC50_invertebrates_µgL),
    by = c("substance"))|>
  select(ID, labelShort, CODE, year, sampleDate, substance, 
    method, conc, 
    limit_of_quantification_µgL, concentration_greater_limit_of_quantification_µgL, 
    LC50_invertebrates_µgL) |> 
  mutate(labelShort = as.factor(labelShort)) |>                                                    # Define the joining variable labelShort as factor
  mutate(tu = ifelse(is.na(LC50_invertebrates_µgL), NA_real_,conc / LC50_invertebrates_µgL))       # Calculate TU for each substance

# -------------------------
# Repeat steps from 2.3 to calculate TUsum for the case of filtered pesticide sampling dates prior to macroinvertebrate sampling to compare with TUsum at step 2.3
# Sum of TU per sampling date
TU_sum_calc3 <- Pesticide_joined3|> 
  group_by(labelShort) |>
  summarise(ID = first(ID),
            CODE = first(CODE),
            year = first(year),
            sampleDate = first(sampleDate),
            method = first(method),
            TUsum_sample  = sum(tu,  na.rm = TRUE), 
            .groups = "drop") 

# Two cases: before macroinvertebrate sampling, and additionally filtering of 'outliers'*: 
outlier_threshold = 2 

TU_sum_calc4 <- TU_sum_calc3 |> 
  group_by(CODE) |> 
  arrange(desc(TUsum_sample)) |> 
  dplyr::summarise(ID = ID[1],
            year = year[1],
            TUsum_before_inver = log10(max(TUsum_sample, na.rm = TRUE)),
            second_TUsum_before_inver = log10(TUsum_sample[2]),
            n_samples = n(),
            n_eds_samples = sum(method=="eds"),
            mean_next = mean(log10(TUsum_sample[2:6]), na.rm = TRUE),                                 # mean log TUsum of up to five subsequent samples with lower values
            max_is_outlier = ifelse((TUsum_before_inver - mean_next)>=outlier_threshold, 1, 0),
            .groups = "drop") |>  # check if difference exceeds the outlier threshold of 2
  as.data.frame()

# If TUsum is outlier, take the second highest TUsum value within the site x year
TU_sum_calc4$TUsum_before_inver_extreme_removed = ifelse(
  TU_sum_calc4$max_is_outlier==0 | is.na(TU_sum_calc4$max_is_outlier), 
  TU_sum_calc4$TUsum_before_inver, 
  TU_sum_calc4$second_TUsum_before_inver)

# Round values below -5 to -5 
TU_sum_calc4$TUsum_before_inver[TU_sum_calc4$TUsum_before_inver<(-5)] <- (-5) 
TU_sum_calc4$TUsum_before_inver_extreme_removed[TU_sum_calc4$TUsum_before_inver_extreme_removed <(-5)] <- (-5) 

# -------------------------------------
# Join with SPEARpesticide for metric comparisons
TUsum_high_peak_removed_before_inver <- metrics_in_manuscript |>
  left_join(TU_sum_calc4 |> select(ID, year, CODE, TUsum_before_inver, TUsum_before_inver_extreme_removed),
    by = c("ID", "year"))

# Join with the table of two TUsum at 2.3 excluding high peaks and full dataset
TUsum_high_peak_removed_before_inver2 <- TUsum_high_peak_removed |>
  left_join(TUsum_high_peak_removed_before_inver  |> select(ID, year, TUsum_before_inver, TUsum_before_inver_extreme_removed),
    by = c("ID", "year"))

# Export
write_xlsx(TUsum_high_peak_removed_before_inver2, here("Recalculated_metrics", "2.4.Pesticide_TUsum_Full_HighPeakRemoved_BeforInver.xlsx"))

# -------------------------------------
# Aggregate mean of two years
TUsum_high_peak_removed_before_inver2_years_mean <- TUsum_high_peak_removed_before_inver2 |>
  group_by(ID) |>
  summarise(mean_SPEAR = round(mean(SPEAR_pesticides, na.rm = TRUE), 2),
    mean_TUsum_full = round(mean(TUsum_full, na.rm = TRUE), 2),
    mean_TUsum_extreme_removed = round(mean(TUsum_extreme_removed, na.rm = TRUE), 2),
    mean_TUsum_before_inver = round(mean(TUsum_before_inver, na.rm = TRUE), 2),
    mean_TUsum_before_inver_extreme_removed = round(mean(TUsum_before_inver_extreme_removed, na.rm = TRUE), 2),
    .groups = "drop")

# ----------------------------
# Plotting four TUsum in steps 2.3 and 2.4. (Figure 2c)
dist_TUsum_4types <- TUsum_high_peak_removed_before_inver2_years_mean |>
  select(ID, mean_TUsum_full,
         mean_TUsum_extreme_removed,
         mean_TUsum_before_inver,
         mean_TUsum_before_inver_extreme_removed) |>
  pivot_longer(cols = c(mean_TUsum_full,
             mean_TUsum_extreme_removed,
             mean_TUsum_before_inver,
             mean_TUsum_before_inver_extreme_removed),
    names_to = "TU_type",
    values_to = "TU_value") |>
  mutate(TU_type = factor(TU_type,
      levels = c("mean_TUsum_before_inver_extreme_removed",
        "mean_TUsum_before_inver",
        "mean_TUsum_extreme_removed",
        "mean_TUsum_full"))) |>
  filter(is.finite(TU_value))

TUsum_violin_box_plot3 <- ggplot(dist_TUsum_4types,
  aes(x = TU_value, y = TU_type, fill = TU_type)) +
  geom_violin(trim = TRUE, width = 0.85, color = NA, alpha = 0.9, na.rm = TRUE) +
  geom_boxplot(width = 0.16, fill = "white",color = "black", linewidth = 0.5, outlier.shape = NA, na.rm = TRUE) +
  scale_fill_manual(values = c(
    mean_TUsum_full = "#BDBDBD",                          
    mean_TUsum_extreme_removed = "#d6604d",               
    mean_TUsum_before_inver = "#dfc27d",                  
    mean_TUsum_before_inver_extreme_removed = "#a6611a")) +  
  
  # labels
  scale_y_discrete(labels = c(
    mean_TUsum_full = "All observations",
    mean_TUsum_extreme_removed = "Exclude high peaks",
    mean_TUsum_before_inver = "Before invertebrates",
    mean_TUsum_before_inver_extreme_removed = "Before invertebrates\n& Exclude high peaks"
  )) +
  scale_x_continuous(breaks = seq(-6, 1, by = 1)) +
  labs(x = "Estimated pesticide toxicity", y = NULL) +
  theme_classic(base_size = 15) +
  theme(text = element_text(),
    axis.text = element_text(size = rel(1.1), color = "black"),
    axis.title = element_text(size = 16),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"))

TUsum_violin_box_plot3

# Combined TUsum plots (Figure 2a, b, c)
TU_combined_row <- TU_violin_box_plot +
  TUsum_violin_box_plot2 +
  TUsum_violin_box_plot3 +
  plot_layout(nrow = 1)

TU_combined_row
# ---------------------------------------

#----------------------------------------
## Join all pesticide toxicity estimates
TUsum_all_joined <- TUmax_sum_join |>
  left_join(Pesticide_metrics_recalc_sum_grab_eds_separated |>
      select(ID, year, TUsum_grab, TUsum_eds),by = c("ID", "year")) |>
  left_join(TUsum_high_peak_removed |>
      select(ID, year, TUsum_extreme_removed), by = c("ID", "year")) |>
  left_join(TUsum_high_peak_removed_before_inver2 |>
      select(ID, year,
             TUsum_before_inver,
             TUsum_before_inver_extreme_removed),by = c("ID", "year")) |>
  select(ID, year, SPEAR_pesticides,
    TUsum_full, TUmax_full,
    TUsum_grab, TUsum_eds,
    TUsum_extreme_removed,
    TUsum_before_inver,
    TUsum_before_inver_extreme_removed)

# Export
write_xlsx(TUsum_all_joined, here("Recalculated_metrics", "2.TUsum_all_joined.xlsx"))

#-----------------------
# Read the table with all seven toxicity metrics
# TUsum_all_joined <- read_excel("Recalculated_metrics/2.TUsum_all_joined.xlsx")

# Calculate mean over two years for models and plotting
TUsum_all_years_mean <- TUsum_all_joined |>
  group_by(ID) |>
  summarise(mean_SPEAR = round(mean(SPEAR_pesticides, na.rm = TRUE), 2),
    mean_TUsum_full = round(mean(TUsum_full, na.rm = TRUE), 2),
    mean_TUmax_full = round(mean(TUmax_full, na.rm = TRUE), 2), 
    mean_TUsum_grab = round(mean(TUsum_grab, na.rm = TRUE), 2),
    mean_TUsum_eds  = round(mean(TUsum_eds,  na.rm = TRUE), 2),
    mean_TUsum_extreme_removed = round(mean(TUsum_extreme_removed, na.rm = TRUE), 2),
    mean_TUsum_before_inver = round(mean(TUsum_before_inver, na.rm = TRUE), 2),
    mean_TUsum_before_inver_extreme_removed = round(mean(TUsum_before_inver_extreme_removed, na.rm = TRUE), 2),
    .groups = "drop")

summary(TUsum_all_years_mean)

# ---------------------------------------
# Run models with all pesticide toxicity metrics
# Model 1: TUsum 
model1_TUsum_full <- lm(mean_SPEAR ~ mean_TUsum_full, data = TUsum_all_years_mean)
(summary_model1 <- summary(model1_TUsum_full))
(r_squared1 <- summary_model1$r.squared)             # R2 = 0.3582 
(slope1 <- coef(model1_TUsum_full)[2])               # Slope = -0.1370

# Model 2: TUmax 
model2_TUmax_full <- lm(mean_SPEAR ~ mean_TUmax_full, data = TUsum_all_years_mean)
(summary_model2 <- summary(model2_TUmax_full))
(r_squared2 <- summary_model2$r.squared)             # R2 = 0.3340
(slope2 <- coef(model2_TUmax_full)[2])               # Slope = -0.1276

# Model 3: TUsum (grab)
model3_TUsum_grab <- lm(mean_SPEAR ~ mean_TUsum_grab, data = TUsum_all_years_mean)
(summary_model3 <- summary(model3_TUsum_grab))
(r_squared3 <- summary_model3$r.squared)             # R2 = 0.3106
(slope3 <- coef(model3_TUsum_grab)[2])               # Slope = -0.1315

# Model 4: TUsum (event)
model4_TUsum_eds <- lm(mean_SPEAR ~ mean_TUsum_eds, data = TUsum_all_years_mean)
(summary_model4 <- summary(model4_TUsum_eds))
(r_squared4 <- summary_model4$r.squared)             # R2 = 0.2318
(slope4 <- coef(model4_TUsum_eds)[2])                # Slope = -0.0957

# Model 5: TUsum (excluded high peaks)
model5_TUsum_extreme_removed <- lm(mean_SPEAR ~ mean_TUsum_extreme_removed, data = TUsum_all_years_mean)
(summary_model5 <- summary(model5_TUsum_extreme_removed))
(r_squared5 <- summary_model5$r.squared)             # R2 = 0.4593
(slope5 <- coef(model5_TUsum_extreme_removed)[2])    # Slope = -0.1426

# Model 6: TUsum (before macroinvertebrates)
model6_TUsum_before_inver <- lm(mean_SPEAR ~ mean_TUsum_before_inver, data = TUsum_all_years_mean)
(summary_model6 <- summary(model6_TUsum_before_inver))
(r_squared6 <- summary_model6$r.squared)             # R2 = 0.2731
(slope6 <- coef(model6_TUsum_before_inver)[2])       # Slope = -0.1090

# Model 7: TUsum (before macroinvertebrates + removed peak removed)
model7_TUsum_before_inver_extreme_removed <- lm(mean_SPEAR ~ mean_TUsum_before_inver_extreme_removed,
  data = TUsum_all_years_mean)
(summary_model7 <- summary(model7_TUsum_before_inver_extreme_removed))
(r_squared7 <- summary_model7$r.squared)                       # R2 = 0.3814
(slope7 <- coef(model7_TUsum_before_inver_extreme_removed)[2]) # Slope = -0.1138

#--------------------------------
# Plot all sevel toxicity estimates (Figure 2d)
# Long-pivot table
plot_TUsum_7models_table <- TUsum_all_years_mean |>
  select(ID, mean_SPEAR,
    mean_TUsum_full,
    mean_TUmax_full,
    mean_TUsum_grab,
    mean_TUsum_eds,
    mean_TUsum_extreme_removed,
    mean_TUsum_before_inver,
    mean_TUsum_before_inver_extreme_removed) |>
  pivot_longer(cols = -c(ID, mean_SPEAR),
    names_to = "TU_type",
    values_to = "TU_value") |>
  filter(is.finite(TU_value), is.finite(mean_SPEAR)) |>
  mutate(TU_type = factor(TU_type,
      levels = c("mean_TUsum_full",
        "mean_TUmax_full",
        "mean_TUsum_grab",
        "mean_TUsum_eds",
        "mean_TUsum_extreme_removed",
        "mean_TUsum_before_inver",
        "mean_TUsum_before_inver_extreme_removed"),
      labels = c("TUsum",
        "TUmax",
        "TUsum_grab",
        "TUsum_eds",
        "TUsum_exclude_high_peaks",
        "TUsum_before_inv",
        "TUsum_before_inv_exclude_high_peaks")))

# R2 label table (place them stacked top-left)
x_label_pos <- quantile(plot_TUsum_7models_table$TU_value, 0.5, na.rm = TRUE)
R2_TU_df <- tibble(TU_type = factor(c("TUsum", "TUmax", "TUsum_grab", "TUsum_eds",
      "TUsum_exclude_high_peaks", "TUsum_before_inv", "TUsum_before_inv_exclude_high_peaks"),
    levels = levels(plot_TUsum_7models_table$TU_type)),
  label = c(sprintf("TUsum: R² = %.2f", r_squared1),
    sprintf("TUmax: R² = %.2f", r_squared2),
    sprintf("Only grab: R² = %.2f", r_squared3),
    sprintf("Only event: R² = %.2f", r_squared4),
    sprintf("Exclude high peaks: R² = %.2f", r_squared5),
    sprintf("Before invertebrates: R² = %.2f", r_squared6),
    sprintf("Before invertebrates\n& Exclude high peaks: R² = %.2f", r_squared7)),
  x = x_label_pos,
  y = c(0.97, 0.91, 0.85, 0.79, 0.73, 0.67, 0.61))

# Color schemes and line types
TU_colors <- c("TUsum" = "#BDBDBD",                         
               "TUmax" = "#4F8FC6",                          
               "TUsum_grab" = "#80cdc1",                          
               "TUsum_eds"  = "#018571",                         
               "TUsum_exclude_high_peaks" = "#d6604d",              
               "TUsum_before_inv" = "#dfc27d",                 
               "TUsum_before_inv_exclude_high_peaks" = "#a6611a") 

TU_linetypes <- c("TUsum" = "solid",
                  "TUmax" = "solid",
                  "TUsum_grab" = "solid",
                  "TUsum_eds"  = "solid",
                  "TUsum_exclude_high_peaks" = "solid",
                  "TUsum_before_inv" = "solid",
                  "TUsum_before_inv_exclude_high_peaks" = "solid")

# Plot
TU_7models_plot <- ggplot(plot_TUsum_7models_table,
  aes(x = TU_value, y = mean_SPEAR, color = TU_type, linetype = TU_type)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, fullrange = TRUE, na.rm = TRUE) +
  scale_color_manual(values = TU_colors) +
  scale_linetype_manual(values = TU_linetypes) +
  labs(x = "Estimated pesticide toxicity", y = "SPEARpesticides", title = "") +
  scale_x_continuous(breaks = seq(-6, 1, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_classic(base_size = 15) +
  theme(text = element_text(size = 14),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(size = 16),
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"),
    legend.position = "none") +
  geom_text(data = R2_TU_df,
    aes(x = x, y = y, label = label, color = TU_type),
    hjust = 0, vjust = 1,
    size = 4.5,
    show.legend = FALSE,
    inherit.aes = FALSE)

TU_7models_plot

#------------------------------------------
# Combined plot (2×2 layout)
combined_plot <- TU_combined_row /
  (plot_spacer() + TU_7models_plot + plot_spacer() +
     plot_layout(widths = c(1, 8, 1)))

combined_plot

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape/Affinity)
ggsave("combined_plot.svg", combined_plot, width = 310, height = 200, units = "mm",  device = "svg")

#------------------------------------

#------------------------------------
# Scatter plots (Supplementary material - Figure A3)
# R2 label
r2_lookup <- c(mean_TUsum_full = r_squared1,
  mean_TUmax_full = r_squared2,
  mean_TUsum_grab = r_squared3,
  mean_TUsum_eds = r_squared4,
  mean_TUsum_extreme_removed = r_squared5,
  mean_TUsum_before_inver = r_squared6,
  mean_TUsum_before_inver_extreme_removed = r_squared7)

# Plot labels
metric_labels <- c(mean_TUsum_full = "TUsum",
  mean_TUmax_full = "TUmax",
  mean_TUsum_grab = "Only grab",
  mean_TUsum_eds = "Only event",
  mean_TUsum_extreme_removed = "Exclude high peaks",
  mean_TUsum_before_inver = "Before invertebrates",
  mean_TUsum_before_inver_extreme_removed = "Before invertebrates\n& Exclude high peaks")

# Plot colors 
metric_colors <- c(mean_TUsum_full = "#000000",   
                   mean_TUmax_full = "#56A0C9",   
                   mean_TUsum_grab = "#4DB6AC",   
                   mean_TUsum_eds = "#148E7C",    
                   mean_TUsum_extreme_removed = "#d6604d",  
                   mean_TUsum_before_inver = "#F28E2B",     
                   mean_TUsum_before_inver_extreme_removed = "#a6611a")  

# Plot panels
panel_titles <- c(mean_TUmax_full = "(a)",
  mean_TUsum_eds = "(b)",
  mean_TUsum_grab = "(c)",
  mean_TUsum_extreme_removed = "(d)",
  mean_TUsum_before_inver = "(e)",
  mean_TUsum_before_inver_extreme_removed = "(f)")

# Use a function to generate coupled plots (TUsum versus one alternative metric, vs SPEARpesticides)
make_TUsum_pair_plot <- function(data, compare_var, panel_title,
                                 x_breaks = seq(-6, 1, by = 1),
                                 y_breaks = seq(0, 1.1, by = 0.2),
                                 y_limits = c(0, 1.1),
                                 x_label_pos = 1,
                                 y_label_pos = c(0.95, 0.88)) {
  
  #-----------------------------------------
  # Build pairwise long table
  pair_table <- data |>
    select(ID, mean_SPEAR, mean_TUsum_full, all_of(compare_var)) |>
    pivot_longer(cols = c(mean_TUsum_full, all_of(compare_var)),
                 names_to = "TU_type",
                 values_to = "TU_value") |>
    filter(is.finite(TU_value), is.finite(mean_SPEAR))
  
  #-----------------------------------------
  # R² label table
  R2_df <- tibble(
    TU_type = c("mean_TUsum_full", compare_var),
    label = c(sprintf("%s: R² = %.2f",
              metric_labels[["mean_TUsum_full"]],
              r2_lookup[["mean_TUsum_full"]]),
      sprintf("%s: R² = %.2f",
              metric_labels[[compare_var]],
              r2_lookup[[compare_var]])),
    x = x_label_pos,
    y = y_label_pos)
  
  #-----------------------------------------
  # Plot
  ggplot(pair_table,
         aes(x = TU_value, y = mean_SPEAR,
             color = TU_type, linetype = TU_type)) +
    
    # Background points: TUsum
    geom_point(data = subset(pair_table, TU_type == "mean_TUsum_full"),
      aes(x = TU_value, y = mean_SPEAR),
      color = metric_colors[["mean_TUsum_full"]],
      size = 2,
      inherit.aes = FALSE) +
    
    # Points of the compared alternative metric
    geom_point(data = subset(pair_table, TU_type == compare_var),
      aes(x = TU_value, y = mean_SPEAR),
      color = metric_colors[compare_var],
      size = 4,
      alpha = 0.5,
      inherit.aes = FALSE) +
    
    # Line of a linear regression model
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, fullrange = TRUE) +
    
    scale_color_manual(values = c(mean_TUsum_full = metric_colors[["mean_TUsum_full"]],
                                  setNames(metric_colors[[compare_var]], compare_var))) +
    
    scale_linetype_manual(values = c(mean_TUsum_full = "solid",
                                     setNames("dashed", compare_var))) +
    
    labs(x = "Estimated pesticide toxicity",
      y = "SPEARpesticides",
      title = panel_title) +
    
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    
    theme_minimal() +
    theme(text = element_text(size = 12),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      plot.margin = margin(t = 10, r = 10, b = 5, l = 5, unit = "pt"),
      legend.position = "none") +
    
    geom_text(data = R2_df,
      aes(x = x, y = y, label = label, color = TU_type),
      hjust = 1, vjust = 0.1,
      size = 4.5,
      show.legend = FALSE,
      inherit.aes = FALSE)
}

# Plot six scatter plots
plot_TUsum_vs_TUmax <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUmax_full",
  panel_title = panel_titles[["mean_TUmax_full"]])
plot_TUsum_vs_TUmax

plot_TUsum_vs_eds <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUsum_eds",
  panel_title = panel_titles[["mean_TUsum_eds"]])
plot_TUsum_vs_eds

plot_TUsum_vs_grab <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUsum_grab",
  panel_title = panel_titles[["mean_TUsum_grab"]])
plot_TUsum_vs_grab

plot_TUsum_vs_extreme_removed <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUsum_extreme_removed",
  panel_title = panel_titles[["mean_TUsum_extreme_removed"]])
plot_TUsum_vs_extreme_removed

plot_TUsum_vs_before_inver <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUsum_before_inver",
  panel_title = panel_titles[["mean_TUsum_before_inver"]])
plot_TUsum_vs_before_inver

plot_TUsum_vs_before_inver_extreme_removed <- make_TUsum_pair_plot(
  data = TUsum_all_years_mean,
  compare_var = "mean_TUsum_before_inver_extreme_removed",
  panel_title = panel_titles[["mean_TUsum_before_inver_extreme_removed"]])
plot_TUsum_vs_before_inver_extreme_removed

# Combine individual scatter plots into one joint plot
Appendix_pesticide_toxicity_plots <-
  (plot_TUsum_vs_TUmax | plot_TUsum_vs_eds | plot_TUsum_vs_grab) /
  (plot_TUsum_vs_extreme_removed | plot_TUsum_vs_before_inver | plot_TUsum_vs_before_inver_extreme_removed)

Appendix_pesticide_toxicity_plots

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("Appendix_pesticide_toxicity_plots.svg", Appendix_pesticide_toxicity_plots, width = 270, height = 180, units = "mm",  device = "svg")
