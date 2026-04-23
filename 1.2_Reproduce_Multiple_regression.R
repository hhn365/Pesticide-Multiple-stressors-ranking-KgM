## Multiple Linear Regression in Loop - using KgM dataset (Liess et al., 2021)

# This code is rewritten from the shared workflow script from co-authors of Liess et al. (2021) to reproduce results 
# in the Table 1 of the manuscript for the three response variables: SPEARpesticide, %EPT, and Saprobic Index
# -------------------------

# -------------------------
# Load packages
library(tidyr)
library(tibble)
library(data.table)
library(dplyr)
library(stringr)
library(broom)
library(MASS)
library(car)
library(relaimpo)
library(corrplot)
library(ggplot2)
library(openxlsx) 
library(readxl)
library(writexl)

# ---------------------

# ---------------------
# Import data (provided by co-authors of Liess et al. 2021)
multiple <- read.delim("Input_data/multipleRegression_data.txt", header = TRUE, sep = "\t")

# Check variable names
colnames(multiple)

# Calculate mean of two years across variables
multiple_mean <- multiple |>
  group_by(siteID) |>
  dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# ---------------------

# ---------------------
# Check correlation among predictors & responses

# List of variables included in the correlation
cor_vars <- c("speartaxa", "sapr", "eptProz",
  "TUmax", "TUmetalpeak", "fliessgeoMean", 
  "Agriculture", "area_km2", "stream_width_m", "water_depth_cm",
  "Gesamt", "grabby", "pHquant9", 
  "NH4", "NO2", "NO3", "PO4", "TN_mgLlog", "TP_mgLlog", 
  "meanO2p25", "meantempP75_nt")

# Correlation data table
cor_data <- multiple_mean |> 
  dplyr::select(all_of(cor_vars))

# Correlation matrix
cor_mat <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")

# Function for the pairwise correlation matrix (including p-values)
cor_test_matrix <- function(mat, ...) {
  mat    <- as.matrix(mat)
  n      <- ncol(mat)
  p_mat  <- matrix(NA_real_, n, n)
  diag(p_mat) <- 0

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- stats::cor.test(mat[, i], mat[, j], ...)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
      }
  }

  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  p_mat
}

# Matrix results
p_mat <- cor_test_matrix(cor_data)

# Plot correlation matrix (21 variables)
png("correlation_plot.png", width = 2300, height = 1500, res = 300)

corrplot(cor_mat, order="original", method="number",
         tl.col="black", number.cex=0.55,tl.srt=45,
         p.mat = p_mat, sig.level = 0.05, 
         insig = "blank")

dev.off()

# ---------------------------------

# ---------------------------------
# Subset the data frame and reshape the table for multiple regression modelling
# Response variables (3)
response_vars <- c("speartaxa", "sapr", "eptProz")

# Predictor variables (including 13 out of 18 variables)
# Note: Excluded predictors are "Agriculture", "area_km2", "NO3", "PO4", "TN_mgLlog". 
predictor_vars <- c("TUmax", "TUmetalpeak", "fliessgeoMean",
  "Gesamt", "grabby", "pHquant9",
  "NH4", "NO2", "TP_mgLlog",
  "meanO2p25", "stream_width_m",
  "water_depth_cm", "meantempP75_nt")

# Reshape the subset data (wide to long)
argrup <- multiple_mean |>
  dplyr::select(siteID, all_of(response_vars), all_of(predictor_vars)) |>
  tidyr::pivot_longer(
    cols      = all_of(response_vars),
    names_to  = "group",   # three biotic metric groups: speartaxa, sapr, eptProz
    values_to = "index") |>
  dplyr::relocate(siteID, group, index) |> 
  as.data.frame()

names(argrup) 

# ---------------------------------

# ---------------------------------
# Multiple linear regression models

# Drop ID and grouping variable
drop_cols <- c("siteID", "group")

# Add a vector: used to identify the single predictor (when only one predictor remains in the model)
namono <- c("index", "Residuals")

# Empty data frame to store results
gresult <- data.frame()

# Run the loop model for each response variable
#----------------------------------------------
# The model loop include four main steps:
# Step 1 – Data preparation & VIF filter of relevant predictors
#   - Run model for each biotic variable (GROUP) separately
#   - Subset argrup to that GROUP, then drop "siteID" and "group" variables
#   - Run a full main-effects model to remove highly collinear predictors using a VIF threshold (VIF ≤ 2)
#
# Step 2 – Stepwise model selection with interactions
#   - Run a "biggest" model with main effects and 2-way interactions
#   - Run a "null" model with no predictors
#   - Use stepAIC (forward) from a null model to the biggest model
#
# Step 3 – Keep only “important” main terms
#   - From the ANOVA table of the selected model, keep only terms with p < 0.05
#   - Apply the rule: keep interaction terms only if BOTH of their main effects are also significant; 
#     otherwise drop those interactions
#
# Step 4 – Final models and results
#   - If multiple predictors remain: refit a model with those terms,
#     compute relative importance (lmg) and merge with estimates and p-values, plus a “Full model” R² row
#   - If only one predictor remains: fit a simple linear model
#     and store its adj.R², estimate, p-value, and a “Full model” row
#   - Append all results for this GROUP to the final results table "gresult"
#----------------------------------------
for (GROUP in unique(argrup$group)) {
  #---------------------------------------------------------
  # Step 1: Prepare data for individual biotic GROUP and VIF for predictor filtering
  #---------------------------------------------------------
  # Extract data per biotic metric
  tempdf <- argrup[argrup$group == GROUP, ]
  
  # Filter out cases where the biotic metric is absent too often 
  # Criteria: when number of zeros <= 50
  if (length(tempdf$index[tempdf$index == 0]) <= 50) {
    rownames(tempdf) <- tempdf$siteID
    tempdf <- tempdf[, !names(tempdf) %in% drop_cols] # Drop ID and grouping columns, keep the index and predictors
    
    full.model <- lm(index ~ ., data = tempdf)        # Run full model: only to check VIF to identify high collinear predictors
    
    # VIF-based filtering of predictors: 
    # Note: run filter inside tryCatch to avoid car::vif() errors on singular fits
    tryCatch({
      v <- car::vif(full.model)
      vif <- data.frame(variable = names(v), 
                        vif = as.numeric(v), 
                        row.names = NULL)
      
      # Criteria: Keep predictors with VIF ≤ 2 (low collinearity)
      vif <- subset(vif, vif <= 2)       
      vif <- vif$variable
      
      # If at least one predictor passes VIF, reduce the tempdf to index and those predictor columns 
      if (length(vif) > 1) {
        tempdf <- tempdf[, colnames(tempdf) %in% c("index", vif)]
        names(tempdf) <- make.names(names(tempdf), unique = TRUE)
      }
    }, error = function(e) { 
      cat("ERROR :", conditionMessage(e), "\n") })
    
    # ------------------------------------------------
    # STEP 2: Full model with interactions
    # ------------------------------------------------
    # Biggest model with all main effects and two-way interactions
    biggest   <- lm(index ~ .^2, data = tempdf)
    # Min.model is the null model with only an intercept
    min.model <- lm(index ~ 1,   data = tempdf)
    # Forward stepwise selection from null to biggest model
    step.model3 <- stepAIC(min.model, 
                           direction = "forward", 
                           scope = formula(biggest), 
                           trace = FALSE)
    
    # ------------------------------------------------
    # STEP 3: Keep only significant terms and apply rules for interactions
    # ------------------------------------------------
    # Get ANOVA table from selected model and keep only rows with p.value and only significant predictors (p < 0.05) 
    colls  <- broom::tidy(aov(step.model3))          
    colls  <- colls[complete.cases(colls$p.value), ] 
    collsp <- colls[colls$p.value < 0.05, ]          
    
    # Interaction term: keep only if both main effects are significant:
    # - collsp: all significant terms (main effects + interactions)
    # - collsoint: the significant main effects
    # - collsmint: only the significant interaction term
    
    # Significant main effects 
    collsoint <- collsp$term[!grepl(":", collsp$term)]  
    
    # Interactions are considered when there are at least 2 (i.e., > 1) main effects
    if (length(collsoint) > 1) {
      # significant interaction term
      collsmint <- collsp$term[grepl(":", collsp$term)] 
      
      # Identify interaction components that are NOT significant as main effects
      # e.g. A:B is significant, but A (or B) is not in collsoint
      coco <- setdiff(unlist(stringr::str_split(collsmint, ":")), collsoint)
      
      # If such components exist, drop any term containing them from collsp
      if (length(coco) > 0 && any(coco %in% colls$term)) {
        pattern <- paste(coco, collapse = "|")
        collsp  <- collsp[!grepl(pattern, collsp$term), ]
      }
    } else {
      # If there are fewer than 2 significant main effects, drop interactions entirely
      collsp <- collsp[collsp$term %in% collsoint, ]
    }
    
    # Vector of significant terms (for subsequent steps) 
    collsp <- c("index", collsp$term)
    
    # Coefficients from the stepwise-selected model
    este <- broom::tidy(step.model3)
    este <- este[este$term %in% collsp, c("term","estimate")]
    
    # ------------------------------------------------
    # STEP 4: Calculate the relative importance and assemble results
    # Note: Two cases - Case 1 with multiple predictors vs Case 2 with only a single predictor
    # ------------------------------------------------
    # Case 1: more than one significant term and fewer than 20 terms in total
    if (length(collsp) > 2 & length(collsp) < 20) {
      
      # Reduce data to the selected terms and put index first
      tempdf2 <- tempdf[, colnames(tempdf) %in% collsp]
      tempdf2 <- dplyr::relocate(tempdf2, index)
      
      # Case 1.1: Singificant term include at least one interaction
      if (any(grepl(":", collsp))) {
        
        # Split into main terms and interaction terms
        main_terms <- collsp[!grepl(":", collsp)]
        main_terms <- main_terms[main_terms != "index"] # drop response
        inter_terms <- collsp[grepl(":", collsp)]
        
        # Build model formula: index ~ main1 + main2 + ...+ inter1 + inter2 + ...
        varias <- paste(c(
          paste(main_terms, collapse = "+"),
          paste(inter_terms, collapse = "+")), 
          collapse = "+")
        
        form     <- as.formula(paste("index ~", varias))
        modelint <- lm(form, data = tempdf2)
        
        # Relative importance for all terms (no grouping of interactions)
        calc <- calc.relimp(modelint, 
                            type = "lmg", 
                            groups = NULL, groupnames = NULL, always = NULL, 
                            rela = FALSE)
        
        # Assemble results: including adj.r2 (lmg) + estimates + p-values + group
        result <- data.frame(term = names(calc$lmg), 
                             adj.r2 = as.numeric(calc$lmg), 
                             row.names = NULL)
        result <- merge(result, este, by = "term", all.x = TRUE)
        result <- merge(result, colls[, c("term","p.value")], by = "term", all.x = TRUE)
        result$group <- GROUP
        
        # Add the full model row with R²
        fulltemp <- data.frame(term = "Full model", 
                               adj.r2 = calc$R2, 
                               estimate = NA, p.value = NA, 
                               group = GROUP)
        
        # Append to the final Results data frame
        gresult  <- rbind(gresult, rbind(result, fulltemp))
      
      # Case 1.2: no interactions among significant terms   
      } else {
        
        # Fit a main-effects-only model based on selected predictor variables 
        full.model.noint <- lm(index ~ ., data = tempdf2)
        
        # Relative importance for main effects
        calc <- calc.relimp(full.model.noint, 
                            type = "lmg", 
                            groups = NULL, groupnames = NULL, always = NULL, rela = FALSE)
        
        # Stepwise selection (both directions) to align estimates with selected terms
        step.model <- stepAIC(full.model.noint, 
                              direction = "both", trace = FALSE)
        
        este <- broom::tidy(step.model)
        este <- este[este$term %in% collsp, c("term","estimate")]
        
        # Assemble results table
        result <- data.frame(term = names(calc$lmg), 
                             adj.r2 = as.numeric(calc$lmg), 
                             row.names = NULL)
        result <- merge(result, este, by = "term", all.x = TRUE)
        result <- merge(result, colls[, c("term","p.value")], by = "term", all.x = TRUE)
        result$group <- GROUP
        
        fulltemp <- data.frame(term = "Full model", 
                               adj.r2 = calc$R2, 
                               estimate = NA, p.value = NA, 
                               group = GROUP)
        gresult  <- rbind(gresult, rbind(result, fulltemp))
      }
     
    # Case 2: exactly one explanatory term (i.e., index + 1 predictor)
    } else if (length(collsp) == 2) {  
      
      # Keep only the index and that single predictor
      tempdf3 <- tempdf[, colnames(tempdf) %in% collsp]
      
      # Name the single predictor
      varname <- setdiff(collsp, namono)  # the single predictor name
      
      # Simple linear model: index ~ predictor
      mod <- summary(lm(stats::reformulate(varname, response = "index"), data = tempdf3))
      
      # Assemble result table - one predictor and including full model row
      result2 <- data.frame(term     = varname,
        adj.r2   = mod$adj.r.squared,
        estimate = mod$coefficients[2, 1],
        p.value  = mod$coefficients[2, 4],
        group    = GROUP)
      
      result2 <- rbind(result2, 
                       data.frame(term = "Full model", 
                                  adj.r2 = mod$adj.r.squared, 
                                  estimate = NA, p.value = NA, 
                                  group = GROUP))
      
      gresult <- rbind(gresult, result2)
    }
  }
  # Show the model progress on which GROUP has been processed
  print(GROUP)
}

#----------------------------------------

#----------------------------------------
# Post-processing of the results

# Copy full results before filtering
groesult <- gresult

# Apply threshold adj.r2 to filter results, but keep Full model
groesult_filt <- groesult[groesult$adj.r2 >= 0.05 | groesult$term == "Full model", ]

# Define interaction term
groesult_filt$is_interaction <- grepl(":", groesult_filt$term)

# Define interaction term per GROUP in the final model
groesult_filt$has_interaction <- ave(
  groesult_filt$is_interaction,
  groesult_filt$group,
  FUN = any)

# Interaction term in the FINAL table
interaction_result <- subset(groesult_filt, is_interaction)

# Final output table
groesult_out <- groesult_filt[, c("group","term","p.value","adj.r2","estimate","has_interaction")]
colnames(groesult_out) <- c("indicator","stressor","signifi","r2","estimate","has_interaction")

write_xlsx(groesult_out, "Recalculated_metrics/1.2.Model_results_with_interactions.xlsx")

# Note: No interaction term in the final models
# -----------------------------------

# -----------------------------------
# Prepare the table for plotting
groesult <- groesult_out |>
  
  # Significance lables (***, **, *, or NA)
  mutate(signifi = case_when(signifi <= 0.001 ~ "***",
                             signifi <= 0.01  ~ "**",
                             signifi <= 0.05  ~ "*",
                             TRUE             ~ NA_character_)) |>
  
  # Rename stressors
  mutate(stressor = dplyr::recode(stressor,
      "Full model"      = "Full model\nStressors",
      "TUmax"           = "Pesticide\npressure",
      "meantempP75_nt"  = "Temperature",
      "meanO2p25"       = "O2 deficiency",
      "grabby"          = "Poor bed\nhabitat\nstructure",
      "fliessgeoMean"   = "Flow\nvelocity",
      "Gesamt"          = "Poor hydro-\nmorphology",
      "TUmetalpeak"     = "TUmetal",
      "pHquant9"        = "pH",
      "TN_mgL"          = "Total\nnitrate",
      "TP_mgLlog"       = "Total\nphosphorous",
      "stream_width_m"  = "Stream width"),
      
      # Reorder stressors for plotting
      stressor = factor(stressor, levels = rev(c("Full model\nStressors",
        "Pesticide\npressure", "O2 deficiency", "Poor hydro-\nmorphology",
        "Poor bed\nhabitat\nstructure", "pH", "NO2", "NH4", "Flow\nvelocity",
        "Stream width", "TUmetal", "Total\nphosphorous")))) |>
  
  # Rename and reorder ecological metrics
  mutate(indicator = dplyr::recode(indicator,
      "speartaxa" = "SPEARpesticides",
      "sapr"      = "Saprobic\nindex",
      "eptProz"   = "%EPT"),
    indicator = factor(indicator, levels = c("SPEARpesticides", "%EPT", "Saprobic\nindex"))) |>
  
  # Define colours for bars and points
  mutate(Col = case_when(
      stressor == "Full model\nStressors" ~ "black",
      estimate <= 0                    ~ "firebrick3",
      TRUE                             ~ "dodgerblue3")) |>
  
  # Add R2 and label text
  mutate(r2 = round(as.numeric(r2), 2),
         labels = if_else(!is.na(signifi), paste(r2, signifi), as.character(r2)),
         labels = if_else(labels == "NA NA", NA_character_, labels))

# -------------------------------------

# -------------------------------------
# Plotting
(Fig <- ggplot(groesult, aes(x = indicator, y = stressor, size  = r2, label = labels, colour = Col)) +
  # Horizontal grid lines
  geom_hline(yintercept = seq(1.5, 16.5, 1), colour = "black", alpha = 0.16) +
  geom_hline(yintercept = c(7.5, 6.5, 0.5)) +
  # Points (R² bubbles)
  geom_point(alpha = 0.7, position = position_nudge(y = +0.16)) +
  # Text labels (R² plus significance)
  geom_text(stat  = "identity", size = 4, family = "Calibri Light",
    colour = "black", vjust = 1, alpha = 1, position = position_nudge(y = -0.15)) +  
  # X axis - biotic metrics
  scale_x_discrete(position = "top", 
                   labels = c("SPEARpesticides" = expression(SPEAR[pesticides]),
                              "%EPT" = "%EPT",
                              "Saprobic\nindex" = "Saprobic\nindex")) +
  # Y axis - Stressors with named labels as in the manuscript
  scale_y_discrete(labels = c(
      "Full model\nStressors"       = "Best-fit model",
      "Pesticide\npressure"         = "Pesticide\ntoxicity",
      "O2 deficiency"               = expression(O[2]),
      "Poor hydro-\nmorphology"     = "Hydromorphological\ndegradation",
      "Poor bed\nhabitat\nstructure"= "Bed\nhabitat\nstructure",
      "NO2"                         = expression(NO[2]),
      "NH4"                         = expression(NH[4]),
      "Flow\nvelocity"              = "Flow\nvelocity")) +
  # Colors
  scale_colour_identity(guide = "none") +
  scale_size(range = c(3, 16)) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 14, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 14, vjust = 1.2, hjust = 0.5),
    axis.title.y = element_blank(),
    legend.position  = "none",
    plot.margin = unit(c(1, 2, 0.5, 1), "cm")) +
  # Annotation
  coord_cartesian(xlim = c(0.5, 3), clip = "off"))

# Save plot as SVG for further refinement in vector graphics softwares (Affinity)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("MLR_matrix_Indices.svg", Fig, width = 8, height = 10)

#-----------------------------------------

#-----------------------------------------
# Main conclusions
# - The code and models could be reproduced partly as results in Table 1 (Liess et al., 2021) using data and code provided by main authors
# - The reused code produced exact estimate and explained variance values as in the main manuscript, but mismatch in estimate signs of some variables
# - The main issue with reproducibility was due to deviating dataset used in the study to the published dataset. 
# - Additionally, the model developed in the manuscript is statistically fragile:
#      + Too complex model structure (13 predictors + 2-way interactions), which in the end produced no interaction outputs
#      + For a small sample size (101 sites), this complicated model structure causes high risk of overfitting and unstable estimates
#      + Inconsistent use of stepwise selection (direction = "forward" for model with interaction vs direction = "both" for model without interaction)
#      + Keep only significant stressors (p < 0.05), which causes a risk of false positives
# - Suggested alternative workflow: Develop multiple linear regression without interaction term and reduce number of predictors for a more stable/reliable outcomes
