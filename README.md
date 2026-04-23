# Pesticide-Multiple-stressors-ranking-KgM

## R code and data for the publication: 
**How do choices in data analysis influence the ranking of pesticide and other environmental stressors as primary drivers of stream macroinvertebrates?**

*Prepared by Hanh H. Nguyen, revised by Ralf B. Schäfer*

*Contact: honghanh.nguyen@uni-due.de*

## Overview
This repository contains the data and R scripts to reproduce all analyses presented in the manuscript. The study evaluates how data aggregation and modelling choices influence conclusions about pesticide toxicity and multiple stressor rankings in 101 small agricultural streams across Germany, using the Kleingewässermonitoring (KgM) dataset (Liess et al., 2021).

## Repository structure

| File/Folder | Description |
|---|---|
| `Input_data/` | Raw data for analyses and modelling |
| `Recalculated_metrics/` | Recalculated pesticide toxicity metrics and model outputs |
| `1.1_Reproduce_TUmax.R` | Reproduction of pesticide toxicity metric (TUmax)|
| `1.2_Reproduce_Multiple_regression.R` | Reproduction of multiple linear regressions of pesticide and multiple stressor effects on stream macroinvertebrates |
| `2_Effects_of_pesticide_mixture_toxicity_aggregations.R` | Pesticide toxicity metrics |
| `3_Models_Multiple_linear_regression_TUsum.R` | Multiple stressor models with AIC-based variable selection |
| `4.1_Models_Spatial_dependency.R` | Spatial dependency check |
| `4.2_Models_Stability_check.R` | Bootstrap stability analysis across five statistical models |

## Requirements
- Required R packages are listed at the top of each script

## Data source
The original KgM monitoring data are from Liess et al. (2021):  
Liess, M. et al. (2021). The lowland stream monitoring dataset (KgM). PANGAEA.  https://doi.org/10.1594/PANGAEA.931673
