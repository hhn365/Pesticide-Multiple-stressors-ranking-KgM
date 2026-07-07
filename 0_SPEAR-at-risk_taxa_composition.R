# 0.Check the percentage of EPT and Spear-at-risk taxa overlap
# Aim: to define the percentage of taxa overlaps contributing to high correlation between SPEARpesticides and %EPT metrics
# -------------------------

# -------------------------
## Import packages
library(dplyr)
library(readr)

# -------------------------

# -------------------------
## Import data
df <- read.table("Input_data/5_1_Invertebrate_Abundances.txt", 
                 quote = '"', na = c("", "NA"),
                 header = TRUE, sep = "\t")

# Inspect dataset
names(df)
glimpse(df)

# Extract data of two variables of interest
df1 <- df |> 
  mutate(SPEciesAtRisk = as.numeric(SPEciesAtRisk),
         EPT = group %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) |> 
  filter(!is.na(SPEciesAtRisk))   # drop records lacking SPEAR trait classification 
                                  # (e.g., taxa identified at coarser level such as Dicranota sp. or Nematoda Gen. sp.)

# -------------------------

# -------------------------
## Taxon percentage overlap analysis
tx <- df1 |>
  distinct(taxa, group, EPT, SPEciesAtRisk)

tab_taxon <- table(EPT = tx$EPT, AtRisk = tx$SPEciesAtRisk)

print(tab_taxon)                                          # raw counts
print(round(prop.table(tab_taxon, margin = 1) * 100, 1))  # count % within each row (i.e., within EPT and non-EPT separately) 

# Note: Results suggest that:
#
#   Raw counts:            non-EPT:  249 not-at-risk, 17 at-risk
#                          EPT:      64 not-at-risk, 155 at-risk
#   Row % (within EPT/non-EPT):      70.8% of EPT taxa ARE at risk
#                                    6.4% of non-EPT taxa ARE at risk
#                                    29.2% of EPT taxa are NOT at risk
#                                    93.6% of non-EPT taxa are NOT at risk          
