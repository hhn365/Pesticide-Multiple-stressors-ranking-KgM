# Step 4. Stability analysis for the main model AIC-based linear regression using boostrap and additional 4 models
# 1- Elastic Net - elastic net regression
# 2- Lasso - least absolute shrinkage and selection operator regression
# 3- MCP - minimax convex penalty
# 4- Triangulated combined model (Reference: https://github.com/cran/stabiliser/tree/master)

# -------------------------

# -------------------------
## Import packages
library(dplyr)
library(purrr)
library(MASS)        
library(stabiliser)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(svglite)
library(readr)    
library(openxlsx) 
library(readxl)
library(writexl)

# -------------------------

# -------------------------
# Import data for modelling
All_metrics <- read_excel("Recalculated_metrics/3.All_biotic_abiotic_metrics_mean_aggregation.xlsx")
names(All_metrics)

# Preprocess data for models
All_metrics_mean <- All_metrics |>
  mutate(Flow = if_else(is.na(Flow), mean(Flow, na.rm = TRUE), Flow),                                           # impute one NA in Flow values
    TN = log10(TN),                                                                                             # log transform nutrients for linearity relationship
    TP = log10(TP)) |>
  dplyr::select(-year, -TUsum_full, -TUsum_extreme_removed, -TUsum_before_inver) |>                             # exclude the year variable after averaging
  group_by(ID) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  rename(Pesticide = TUsum_before_inver_extreme_removed,                                                        # Model with one of the four pesticide metrics from step 3
         Trace_elements = TUsum_Trace_elements,
         Agriculture = AgriLand_percent) |>
  arrange(as.numeric(sub("S", "", ID)))

# Scale all predictors for comparability of model estimates 
All_metrics_scaled <- All_metrics_mean
All_metrics_scaled[, 5:ncol(All_metrics_mean)] <- scale(All_metrics_mean[, 5:ncol(All_metrics_mean)])

summary(All_metrics_scaled)

# List of response and predictor variables
# Response variables 
response_vars <- c("speartaxa", "sapr", "eptProz")

# Predictor variables 
predictor_vars <- c("Pesticide",
                    "Trace_elements",
                    "Agriculture",
                    "Morphology",
                    "Habitat",
                    "Stream_width",
                    "Stream_depth",
                    "pH",
                    "TN",
                    "TP",
                    "Flow")

#---------------------
## AIC-based linear regression bootstrap analysis to address whether variable selection is robust or not
# Compare with other shrinkage-based robustness using 4 bootstrap penalized regression stability: Lasso / Elastic Net / MCP / triangulate model 
# Include separate models for three biotic metrics: speartaxa, eptProz, and sapr
# ---------------------

# ---------------------
# Model for SPEARpesticides (speartaxa)

# Input data
# Response variable
resp1 <- "speartaxa"

# Data table with the specific response variable and all relevant predictors
dat1 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp1, predictor_vars)))

# Predictor names (excluding the response column)
predictors <- setdiff(names(dat1), resp1)

# Number of bootstrap resamples
B <- 1000   

# A list storing results of each bootstrap iteration
boot_results <- vector("list", B)

# Bootstrap and stepwise AIC loop 
set.seed(365)  # For reproducibility of of bootstrap resampling

for (b in seq_len(B)) {
  # Bootstrap sample: sample rows with replacement to create a resampled dataset
  df_boot <- dat1[sample(nrow(dat1), replace = TRUE), ]

  # Full model: formula response ~ all predictors
  form_full <- as.formula(paste(resp1, "~", paste(predictors, collapse = "+")))
  fit_full <- lm(form_full, data = df_boot)

  # Backward stepwise selection using AIC (start from the full model with all predictors and remove predictors to have reduced AIC)
  fit_step <- MASS::stepAIC(fit_full, direction = "backward", trace = FALSE)
  
  # Store coefficients of the final selected models
  boot_results[[b]] <- coef(fit_step)
}

# Check example outputs of the bootstrap model (e.g., b = 1 and 10)
boot_results[[1]]
boot_results[[10]]

# Create summary table of bootstrap results, including stability index
stability_AIC1 <- tibble(variable = predictors,        # include all predictors in the final models, each predictor per row
                        mean_coefficient = NA_real_,   # mean coefficient (among times the predictor was selected)
                        ci_lower = NA_real_,           # 2.5% bootstrap quantile
                        ci_upper = NA_real_,           # 97.5% bootstrap quantile
                        bootstrap_p = NA_real_,        # Bootstrap_p instability measure (closer to 0 indicates lower instability or higher stability)
                        stability = NA_real_,          # % of bootstrap samples where predictor was selected
                        stable = NA_character_)        # Flag when the predictor selection is stable (i.e., above threshold)

# For each predictor (v), compute the above-listed metrics
for (v in predictors) {
  # If predictor v was selected in the bootstrap model, extract its coefficient estimate, otherwise store NA
  coef_values <- map_dbl(boot_results,
    ~ ifelse(v %in% names(.x), .x[[v]], NA_real_))

  selected <- !is.na(coef_values)                                 # A logical vector, TRUE when v was selected
  stab <- mean(selected) * 100                                    # selection frequency in percent
  stability_AIC1$stability[stability_AIC1$variable == v] <- stab  # store stability value stab in the corresponding row of the summary table
  
  # Compute coefficient summaries for all cases when predictor v was selected
  if (any(selected)) {
    vals <- coef_values[selected]                                 # coefficient estimates
    m <- mean(vals)                                               # mean coefficient across selected runs

    stability_AIC1$mean_coefficient[stability_AIC1$variable == v] <- m
    stability_AIC1$ci_lower[stability_AIC1$variable == v] <- quantile(vals, .025)
    stability_AIC1$ci_upper[stability_AIC1$variable == v] <- quantile(vals, .975)
    stability_AIC1$bootstrap_p[stability_AIC1$variable == v] <- mean(sign(vals) != sign(m))  
  }
}

# Stability threshold
perm_thresh_AIC <- 70

# Identify stable stressors (i.e., selected in >= threshold of bootstrap samples)
stability_AIC1$stable <- ifelse(stability_AIC1$stability >= perm_thresh_AIC, "*", NA)

#------------------------
# Other 4 stability analyses
set.seed(365)

stable_all1 <- stabilise(data = dat1,
                        outcome = resp1,
                        models = c("enet", "lasso", "mcp"),
                        boot_reps  = 1000)

triangulated1 <- triangulate(stable_all1)

#------------------------
# Joint table of all five model results
# Function for joining tables
extract_model_table <- function(stab_obj, model_name) {
  stab_obj$stability |>
    dplyr::transmute(model = model_name, 
                     variable,
                     mean_coefficient,
                     ci_lower,
                     ci_upper,
                     bootstrap_p,
                     stability,
                     stable)}

stab_AIC1 <- stability_AIC1 |> dplyr::mutate(model = "stepwise AIC", .before = 1)
stab_enet1  <- extract_model_table(stable_all1$enet,  "Enet")
stab_lasso1 <- extract_model_table(stable_all1$lasso, "Lasso")
stab_mcp1   <- extract_model_table(stable_all1$mcp,   "MCP")
stab_comb1  <- extract_model_table(triangulated1$combi, "Combined")

all_models_table1 <- dplyr::bind_rows(stab_AIC1,
                                     stab_enet1,
                                     stab_lasso1,
                                     stab_mcp1,
                                     stab_comb1)

all_models_table1

#------------------------
# Plot stability 
plot_stability <- function(stab_tbl, perm_thresh, title) {
  
  # Identify 5 most stable stressors
  top5 <- stab_tbl |> dplyr::arrange(desc(stability)) |>
    dplyr::slice(1:5) |>
    dplyr::pull(variable)
  
  stab_tbl |> dplyr::filter(!is.na(bootstrap_p)) |>
    dplyr::mutate(side = dplyr::if_else(stability < perm_thresh, "Unstable", "Stable")) |>
    
    ggplot(aes(x = stability, y = bootstrap_p, colour = side)) +
    # Points
    geom_jitter(height = 0.01, width = 1, size = 2.8, alpha = 0.9) +
    geom_vline(xintercept = perm_thresh, linetype = "dashed") +
    
    # Labels only for the top 5 stable stressors
    ggrepel::geom_text_repel(data = ~ dplyr::filter(.x, variable %in% top5),
      aes(label = variable),
      size = 4,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.size = 0.3,
      show.legend = FALSE) +
    
    # Axes
    scale_y_reverse(expand = expansion(mult = c(0.01, 0.05))) +
    coord_cartesian(ylim = c(0.4, 0)) +
    scale_x_continuous(limits = c(20, 100)) +
    
    # Colours
    scale_colour_manual(values = c("Stable" = "dodgerblue3", "Unstable" = "firebrick3")) +
    
    # Labels
    labs(x = "Stability (%)",
      y = "Bootstrap-p",
      title = title,
      colour = "Variable stability") +
    
    theme_minimal(base_size = 12, base_family = "Helvetica") +
    theme(
      plot.title   = element_text(size = 13, face = "bold", hjust = 0),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 14, colour = "black"),
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      plot.margin  = margin(5, 5, 5, 5))}

p_AIC1   <- plot_stability(stability_AIC1, perm_thresh_AIC,           "AIC-based linear regression")
p_enet1  <- plot_stability(stable_all1$enet$stability,  stable_all1$enet$perm_thresh,  "Enet")
p_lasso1 <- plot_stability(stable_all1$lasso$stability, stable_all1$lasso$perm_thresh, "Lasso")
p_mcp1   <- plot_stability(stable_all1$mcp$stability,   stable_all1$mcp$perm_thresh,   "MCP")
p_comb1  <- plot_stability(triangulated1$combi$stability, triangulated1$combi$perm_thresh, "Triangulated\nElastic Net + Lasso + MCP")

final_plot1.1 <- (p_lasso1| p_enet1) / ( p_mcp1| p_comb1) 

final_plot1.2 <- (p_AIC1 | final_plot1.1 ) +
  plot_annotation(tag_levels = "A", title = "Stability analysis across five model-selection methods")

final_plot1.2  

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("Stability_plot_SPEAR.svg", final_plot1.2, width = 540, height = 240, units = "mm")

# Note on the plot results: Stability plot identifies thresholds that are repeatedly selected and directionally robust across model resampling
# x-axis: Stability (%): proportion of bootstrap samples in which a stressor was retained in the model and contributed to variation of the response variable
# y-axis: Bootstrap-p: the proportion of bootstrap estimates whose sign differs from the mean bootstrap coefficient, i.e., quantifies the directional uncertainty of each stressor effect ( 0 means high consistency, 1 indicate unstable/inconsistent effects)
# dash line: stability threshold derived from permutation-based null models (penalised approach) or AIC criterion, i.e. predictor above this threshold are stable

#----------------------

# ---------------------
# Model for %EPT (analogue method to SPEAR model)
resp2 <- "eptProz"

dat2 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp2, predictor_vars)))

predictors <- setdiff(names(dat2), resp2)

# Bootstrap samples
B <- 1000   

boot_results <- vector("list", B)

# Fit a model with bootstrap
set.seed(72)

for (b in seq_len(B)) {
  # Bootstrap resample
  df_boot <- dat2[sample(nrow(dat2), replace = TRUE), ]

  # Full model
  form_full <- as.formula(
    paste(resp2, "~", paste(predictors, collapse = "+"))
  )
  fit_full <- lm(form_full, data = df_boot)

  # Stepwise AIC
  fit_step <- MASS::stepAIC(fit_full, direction = "backward", trace = FALSE)

  boot_results[[b]] <- coef(fit_step)
}

# Stability table
stability_AIC2 <- tibble(variable = predictors,
                        mean_coefficient = NA_real_,
                        ci_lower = NA_real_,
                        ci_upper = NA_real_,
                        bootstrap_p = NA_real_,
                        stability = NA_real_,
                        stable = NA_character_)

for (v in predictors) {
  coef_values <- map_dbl(
    boot_results,
    ~ ifelse(v %in% names(.x), .x[[v]], NA_real_))

  selected <- !is.na(coef_values)
  stab <- mean(selected) * 100
  stability_AIC2$stability[stability_AIC2$variable == v] <- stab

  if (any(selected)) {
    vals <- coef_values[selected]
    m <- mean(vals)

    stability_AIC2$mean_coefficient[stability_AIC2$variable == v] <- m
    stability_AIC2$ci_lower[stability_AIC2$variable == v] <- quantile(vals, .025)
    stability_AIC2$ci_upper[stability_AIC2$variable == v] <- quantile(vals, .975)
    stability_AIC2$bootstrap_p[stability_AIC2$variable == v] <- mean(sign(vals) != sign(m))
  }
}

# Stability threshold
perm_thresh_AIC <- 70

# Identify stable stressors from stepAIC method
stability_AIC2$stable <- ifelse(stability_AIC2$stability >= perm_thresh_AIC, "*", NA)

#------------------------
# Other 4 stability approaches
set.seed(72)

stable_all2 <- stabilise(data = dat2,
                        outcome = resp2,
                        models = c("enet", "lasso", "mcp"),
                        boot_reps  = 1000)

triangulated2 <- triangulate(stable_all2)

#------------------------
# Joint table of all five model results
stab_AIC2 <- stability_AIC2 |> dplyr::mutate(model = "AIC", .before = 1)
stab_enet2  <- extract_model_table(stable_all2$enet,  "Elastic Net")
stab_lasso2 <- extract_model_table(stable_all2$lasso, "Lasso")
stab_mcp2   <- extract_model_table(stable_all2$mcp,   "MCP")
stab_comb2  <- extract_model_table(triangulated2$combi, "Combined")

all_models_table2 <- dplyr::bind_rows(stab_AIC2,
                                     stab_enet2,
                                     stab_lasso2,
                                     stab_mcp2,
                                     stab_comb2)

all_models_table2

#------------------------
# Plot stability 
p_AIC2   <- plot_stability(stability_AIC2, perm_thresh_AIC,           "AIC-based linear regression")
p_enet2  <- plot_stability(stable_all2$enet$stability,  stable_all2$enet$perm_thresh,  "Enet")
p_lasso2 <- plot_stability(stable_all2$lasso$stability, stable_all2$lasso$perm_thresh, "Lasso")
p_mcp2   <- plot_stability(stable_all2$mcp$stability,   stable_all2$mcp$perm_thresh,   "MCP")
p_comb2  <- plot_stability(triangulated2$combi$stability, triangulated2$combi$perm_thresh, "Triangulated\nElastic Net + Lasso + MCP")

final_plot2.1 <- (p_lasso2| p_enet2) / ( p_mcp2| p_comb2) 

final_plot2.2 <- (p_AIC2 | final_plot2.1 ) +
  plot_annotation(tag_levels = "A", title = "Stability analysis across five model-selection methods")

final_plot2.2

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("Stability_plot_EPT.svg", final_plot2.2, width = 540, height = 240, units = "mm")


#----------------------

# ---------------------
# Model for  Saprobic index (analogue method to SPEAR model)
resp3 <- "sapr"

dat3 <- All_metrics_scaled |>
  dplyr::select(dplyr::all_of(c(resp3, predictor_vars)))

predictors <- setdiff(names(dat3), resp3)

# Bootstrap samples
B <- 1000   

boot_results <- vector("list", B)

# Fit a model with bootstrap
set.seed(58)

for (b in seq_len(B)) {
  # Bootstrap resample
  df_boot <- dat3[sample(nrow(dat3), replace = TRUE), ]

  # Full model
  form_full <- as.formula(
    paste(resp3, "~", paste(predictors, collapse = "+"))
  )
  fit_full <- lm(form_full, data = df_boot)

  # Stepwise AIC
  fit_step <- MASS::stepAIC(fit_full, direction = "backward", trace = FALSE)

  boot_results[[b]] <- coef(fit_step)
}

# Stability table
stability_AIC3 <- tibble(variable = predictors,
                        mean_coefficient = NA_real_,
                        ci_lower = NA_real_,
                        ci_upper = NA_real_,
                        bootstrap_p = NA_real_,
                        stability = NA_real_,
                        stable = NA_character_)

for (v in predictors) {
  coef_values <- map_dbl(
    boot_results,
    ~ ifelse(v %in% names(.x), .x[[v]], NA_real_))

  selected <- !is.na(coef_values)
  stab <- mean(selected) * 100
  stability_AIC3$stability[stability_AIC3$variable == v] <- stab

  if (any(selected)) {
    vals <- coef_values[selected]
    m <- mean(vals)

    stability_AIC3$mean_coefficient[stability_AIC3$variable == v] <- m
    stability_AIC3$ci_lower[stability_AIC3$variable == v] <- quantile(vals, .025)
    stability_AIC3$ci_upper[stability_AIC3$variable == v] <- quantile(vals, .975)
    stability_AIC3$bootstrap_p[stability_AIC3$variable == v] <- mean(sign(vals) != sign(m))
  }
}

# Stability threshold
perm_thresh_AIC <- 70

# Identify stable stressors from stepAIC method
stability_AIC3$stable <- ifelse(stability_AIC3$stability >= perm_thresh_AIC, "*", NA)

#------------------------
# Other 4 stability approaches
set.seed(58)

stable_all3 <- stabilise(data = dat3,
                        outcome = resp3,
                        models = c("enet", "lasso", "mcp"),
                        boot_reps  = 1000)

triangulated3 <- triangulate(stable_all3)

#------------------------
# Joint table of all five model results
stab_AIC3 <- stability_AIC3 |> dplyr::mutate(model = "AIC", .before = 1)
stab_enet3  <- extract_model_table(stable_all3$enet,  "Enet")
stab_lasso3 <- extract_model_table(stable_all3$lasso, "Lasso")
stab_mcp3   <- extract_model_table(stable_all3$mcp,   "MCP")
stab_comb3  <- extract_model_table(triangulated3$combi, "Combined")

all_models_table3 <- dplyr::bind_rows(stab_AIC3,
                                     stab_enet3,
                                     stab_lasso3,
                                     stab_mcp3,
                                     stab_comb3)

all_models_table3

#------------------------
# Plot stability 
p_AIC3   <- plot_stability(stability_AIC3, perm_thresh_AIC,           "AIC-based linear regression")
p_enet3  <- plot_stability(stable_all3$enet$stability,  stable_all3$enet$perm_thresh,  "Enet")
p_lasso3 <- plot_stability(stable_all3$lasso$stability, stable_all3$lasso$perm_thresh, "Lasso")
p_mcp3   <- plot_stability(stable_all3$mcp$stability,   stable_all3$mcp$perm_thresh,   "MCP")
p_comb3  <- plot_stability(triangulated3$combi$stability, triangulated3$combi$perm_thresh, "Triangulated\nElastic Net + Lasso + MCP")

final_plot3.1 <- (p_lasso3| p_enet3) / ( p_mcp3| p_comb3) 

final_plot3.2 <- (p_AIC3 | final_plot3.1 ) +
  plot_annotation(tag_levels = "A", title = "Stability analysis across five model-selection methods")

final_plot3.2

# Save plot as SVG for further refinement in vector graphics softwares (Inkscape)
# Note: Final figure in the manuscript differs slightly from these outputs
# due to manual adjustments of labels, annotations, and layout for publication
ggsave("Stability_plot_SaprobicIndex.svg", final_plot3.2, width = 540, height = 240, units = "mm")

#-------------------------
# Save Stability tables of three biotic indicators
write_xlsx(list(SPEAR = all_models_table1,
                percentEPT = all_models_table2,
                SaprobicIndex = all_models_table3),"Recalculated_metrics/4.Stability_outputs_combined.xlsx")
