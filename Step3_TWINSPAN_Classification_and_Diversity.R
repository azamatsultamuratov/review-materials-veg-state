# ==============================================================================
# Script Name: Syntaxonomic Classification and Alpha/Beta Diversity Analysis
# Description: This script performs the final Modified TWINSPAN classification 
#              to establish the hierarchical syntaxonomic structure (6 upper and 
#              8 lower clusters). It identifies indicator/diagnostic species, 
#              calculates species constancy, and computes key data variability 
#              statistics including alpha, gamma, and partitioned beta diversity 
#              (turnover vs. nestedness).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(twinspanR)    
library(dplyr)        
library(tidyr)        
library(vegan)        
library(indicspecies) 
library(tibble)       
library(betapart)     

# ------------------------------------------------------------------------------
# 1. LOAD AND PREPARE INPUT DATA
# ------------------------------------------------------------------------------
# Load the frequency-filtered matrices (from Script 01)
# mtws_f4: sqrt-transformed matrix used for clustering to reduce dominance effect
# comm_raw_f4: raw matrix used for calculating true cover means and constancy
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
comm_raw_f4 <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)


# Safety check: Ensure strict alignment of plots and species between both matrices
common_plots <- intersect(rownames(mtws_f4), rownames(comm_raw_f4))
common_spp   <- intersect(colnames(mtws_f4), colnames(comm_raw_f4))

mtws_use <- mtws_f4[common_plots, common_spp, drop = FALSE]
raw_use  <- comm_raw_f4[common_plots, common_spp, drop = FALSE]

# ------------------------------------------------------------------------------
# 2. FINAL MODIFIED TWINSPAN CLASSIFICATION
# ------------------------------------------------------------------------------
# Standard pseudospecies cut levels used in traditional and modified TWINSPAN
cut_use <- c(0, 0.3, 0.7, 2, 4)


valid_plots <- rowSums(mtws_use) > 0
mtws_use <- mtws_use[valid_plots, , drop = FALSE]
raw_use  <- raw_use[valid_plots, , drop = FALSE]

# Generate Upper Hierarchical Level (6 Clusters)
tw_upper <- twinspan(
  com = mtws_use,
  modif = TRUE,
  clusters = 6,
  cut.levels = cut_use,
  min.group.size = 5,
  diss = "bray",
  quiet = TRUE,
  show.output.on.console = FALSE
)

# Generate Lower Hierarchical Level (8 Clusters)
tw_lower <- twinspan(
  com = mtws_use,
  modif = TRUE,
  clusters = 8,
  cut.levels = cut_use,
  min.group.size = 5,
  diss = "bray",
  quiet = TRUE,
  show.output.on.console = FALSE
)

# Extract cluster membership assignments
upper_class <- as.factor(cut(tw_upper))
lower_class <- as.factor(cut(tw_lower))

# Compile cluster memberships into a data frame
class_df <- data.frame(
  plot_id = rownames(mtws_use),
  upper_class = upper_class,
  lower_class = lower_class
)

# Export plot memberships to specific syntaxa
write.csv(class_df, "mtws_membership_upper6_lower8.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. INDICATOR (DIAGNOSTIC) SPECIES ANALYSIS
# ------------------------------------------------------------------------------
# Function to compute Indicator Value (IndVal) and extract significant species
get_indicator_table <- function(comm, groups, alpha = 0.05, nperm = 999) {
  groups <- as.factor(groups)
  
  ind <- multipatt(
    comm,
    groups,
    func = "IndVal.g",
    duleg = TRUE,
    control = how(nperm = nperm)
  )
  
  sign_tab <- ind$sign
  if(is.null(sign_tab) || nrow(sign_tab) == 0) return(NULL)
  
  # Format the output table for publication readiness
  sign_tab <- sign_tab %>%
    rownames_to_column("species") %>%
    mutate(
      p.value = ifelse(is.null(p.value), NA, p.value)
    ) %>%
    filter(p.value <= alpha)
  
  return(sign_tab)
}

# Run indicator analysis for both hierarchical levels
ind_upper <- get_indicator_table(mtws_use, upper_class, alpha = 0.05, nperm = 999)
ind_lower <- get_indicator_table(mtws_use, lower_class, alpha = 0.05, nperm = 999)

# Export significant diagnostic species
if(!is.null(ind_upper)) write.csv(ind_upper, "indicator_species_upper6.csv", row.names = FALSE)
if(!is.null(ind_lower)) write.csv(ind_lower, "indicator_species_lower8.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. PREPARE LONG-FORMAT DATA FOR SYNTAXONOMIC SUMMARIES
# ------------------------------------------------------------------------------
# Convert raw matrix to long format and merge with cluster assignments
raw_long <- raw_use %>%
  rownames_to_column("plot_id") %>%
  pivot_longer(-plot_id, names_to = "species", values_to = "cover") %>%
  left_join(class_df, by = "plot_id")

# ------------------------------------------------------------------------------
# 5. SPECIES SUMMARY BY UPPER CLASSES (k=6)
# ------------------------------------------------------------------------------
summary_upper <- raw_long %>%
  group_by(upper_class, species) %>%
  summarise(
    mean_cover = mean(cover, na.rm = TRUE),
    total_cover = sum(cover, na.rm = TRUE),
    constancy_pct = mean(cover > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(upper_class, desc(constancy_pct), desc(mean_cover))

write.csv(summary_upper, "species_summary_upper6_raw.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. SPECIES SUMMARY BY LOWER CLASSES (k=8)
# ------------------------------------------------------------------------------
summary_lower <- raw_long %>%
  group_by(lower_class, species) %>%
  summarise(
    mean_cover = mean(cover, na.rm = TRUE),
    total_cover = sum(cover, na.rm = TRUE),
    constancy_pct = mean(cover > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(lower_class, desc(constancy_pct), desc(mean_cover))

write.csv(summary_lower, "species_summary_lower8_raw.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. EXTRACT DOMINANT SPECIES PER SYNTAXON
# ------------------------------------------------------------------------------
# Isolate the top 10 dominant species (by mean cover) for each cluster
dominants_upper <- summary_upper %>%
  group_by(upper_class) %>%
  slice_max(order_by = mean_cover, n = 10, with_ties = FALSE) %>%
  ungroup()

dominants_lower <- summary_lower %>%
  group_by(lower_class) %>%
  slice_max(order_by = mean_cover, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(dominants_upper, "dominant_species_upper6_top10.csv", row.names = FALSE)
write.csv(dominants_lower, "dominant_species_lower8_top10.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. CLUSTER SIZE DISTRIBUTION (NUMBER OF PLOTS PER SYNTAXON)
# ------------------------------------------------------------------------------
upper_sizes <- data.frame(upper_class = names(table(upper_class)),
                          n_plots = as.integer(table(upper_class)))
lower_sizes <- data.frame(lower_class = names(table(lower_class)),
                          n_plots = as.integer(table(lower_class)))

write.csv(upper_sizes, "upper6_cluster_sizes.csv", row.names = FALSE)
write.csv(lower_sizes, "lower8_cluster_sizes.csv", row.names = FALSE)

cat("\n====================================\n")
cat("CLUSTER SIZES SUMMARY\n")
cat("====================================\n")
print(upper_sizes)
print(lower_sizes)

# ------------------------------------------------------------------------------
# 9. DIVERSITY STATISTICS (ALPHA, GAMMA, BETA PARTITIONING)
# ------------------------------------------------------------------------------

# 9.1 Convert community matrix to Presence/Absence for binary metrics
pa_use <- raw_use
pa_use[pa_use > 0] <- 1
pa_use <- as.data.frame(pa_use)

# 9.2 Compute Alpha, Gamma, and Whittaker's Beta Diversity
alpha_i <- rowSums(pa_use)
mean_alpha <- mean(alpha_i, na.rm = TRUE) # Mean species richness per plot

gamma_total <- sum(colSums(pa_use) > 0)   # Total species richness in the dataset

# Whittaker's original beta diversity formula: (Gamma / Mean Alpha) - 1
whittaker_beta <- gamma_total / mean_alpha - 1 

# 9.3 Compute Mean Pairwise Compositional Dissimilarities
# Mean Jaccard dissimilarity
mean_jaccard <- mean(vegdist(pa_use, method = "jaccard", binary = TRUE))

# Mean Sørensen dissimilarity (equivalent to Bray-Curtis on binary data)
mean_sorensen <- mean(vegdist(pa_use, method = "bray", binary = TRUE))

# 9.4 Partition Beta Diversity (Turnover vs. Nestedness)
# Utilizing the betapart package to understand the drivers of beta diversity
beta_core <- betapart.core(pa_use)
beta_pair <- beta.pair(beta_core, index.family = "sorensen")

# Simpson turnover component (species replacement)
mean_simpson <- mean(as.matrix(beta_pair$beta.sim)[upper.tri(as.matrix(beta_pair$beta.sim))], na.rm = TRUE)

# Nestedness-resultant component (species loss/gain)
mean_nestedness <- mean(as.matrix(beta_pair$beta.sne)[upper.tri(as.matrix(beta_pair$beta.sne))], na.rm = TRUE)

# Verification check: Overall Sørensen = Simpson + Nestedness
mean_sorensen_check <- mean(as.matrix(beta_pair$beta.sor)[upper.tri(as.matrix(beta_pair$beta.sor))], na.rm = TRUE)

# 9.5 Display Diversity Results in Console
cat("\n====================================\n")
cat("DATA VARIABILITY & DIVERSITY STATISTICS\n")
cat("====================================\n")
cat(sprintf("Mean alpha richness: %.2f\n", mean_alpha))
cat(sprintf("Gamma richness: %d\n", gamma_total))
cat(sprintf("Whittaker beta-diversity: %.2f\n\n", whittaker_beta))

cat("Pairwise compositional dissimilarity (Presence/Absence):\n")
cat(sprintf("Jaccard dissimilarity: %.2f\n", mean_jaccard))
cat(sprintf("Sørensen dissimilarity: %.2f\n", mean_sorensen))
cat(sprintf("Simpson turnover (Replacement): %.2f\n", mean_simpson))
cat(sprintf("Nestedness-resultant (Loss/Gain): %.2f\n", mean_nestedness))

# 9.6 Export Diversity Metrics to CSV (Removed "_Q1" from filename)
variability_out <- data.frame(
  metric = c(
    "Mean alpha richness",
    "Gamma richness",
    "Whittaker beta-diversity",
    "Jaccard dissimilarity",
    "Sorensen dissimilarity",
    "Simpson turnover",
    "Nestedness-resultant"
  ),
  value = c(
    mean_alpha,
    gamma_total,
    whittaker_beta,
    mean_jaccard,
    mean_sorensen,
    mean_simpson,
    mean_nestedness
  )
)

write.csv(variability_out, "data_variability_stats.csv", row.names = FALSE)
cat("\nProcess successfully completed. All output files saved.\n")

