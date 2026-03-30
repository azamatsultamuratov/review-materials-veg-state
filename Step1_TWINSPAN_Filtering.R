# ==============================================================================
# Script Name: Data Preprocessing and Filtering for Modified TWINSPAN
# Description: This script cleans raw geobotanical plot data, removes empty 
#              plots/species, and generates frequency-filtered matrices. 
#              While thresholds of 3, 4, and 5 are computed for sensitivity 
#              testing, freq >= 4 is established as the optimal threshold for 
#              the final analysis to balance noise reduction and rare syntaxa 
#              preservation.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(dplyr) # For efficient data manipulation

# ------------------------------------------------------------------------------
# 1. READ AND PREPARE RAW DATA
# ------------------------------------------------------------------------------
# Read the dataset. The first column must contain the unique plot IDs.
df <- read.csv("ready_data_AD_DM.csv", header = TRUE, check.names = FALSE)

# Extract plot IDs
plot_id <- df[[1]]

# Isolate the species matrix by dropping the ID column
comm <- df[, -1, drop = FALSE]

# Ensure all community data columns are numeric; replace NAs with 0
comm[] <- lapply(comm, function(x) as.numeric(as.character(x)))
comm[is.na(comm)] <- 0

# Assign plot IDs as row names for community matrix operations
rownames(comm) <- plot_id

# explicitly convert to data.frame
comm <- as.data.frame(comm)

# ------------------------------------------------------------------------------
# 2. BASIC DIAGNOSTICS OF RAW DATA
# ------------------------------------------------------------------------------
n_plots   <- nrow(comm)
n_species <- ncol(comm)

# Identify plots with no species and species occurring in no plots
empty_plots   <- rowSums(comm) == 0
empty_species <- colSums(comm) == 0

# Print raw data summary
cat("====================================\n")
cat("RAW DATA DIAGNOSTICS\n")
cat("====================================\n")
cat("Number of plots   :", n_plots, "\n")
cat("Number of species :", n_species, "\n")
cat("Empty plots       :", sum(empty_plots), "\n")
cat("Empty species     :", sum(empty_species), "\n")

# Calculate and print matrix sparsity (proportion of zeros)
sparsity <- sum(comm == 0) / (nrow(comm) * ncol(comm))
cat("Sparsity (%)      :", round(sparsity * 100, 2), "\n")

# Calculate species frequencies (number of plots each species occurs in)
sp_freq <- colSums(comm > 0)

cat("Species in <=1 plot :", sum(sp_freq <= 1), "\n")
cat("Species in <=2 plots:", sum(sp_freq <= 2), "\n")
cat("Species in <=3 plots:", sum(sp_freq <= 3), "\n")
cat("Species in <=4 plots:", sum(sp_freq <= 4), "\n")
cat("Species in <=5 plots:", sum(sp_freq <= 5), "\n")

# ------------------------------------------------------------------------------
# 3. REMOVE EMPTY PLOTS AND SPECIES
# ------------------------------------------------------------------------------
# Filter out plots with zero total cover and species with zero occurrences
comm_clean <- comm[!empty_plots, , drop = FALSE]
comm_clean <- comm_clean[, colSums(comm_clean) > 0, drop = FALSE]

cat("\n====================================\n")
cat("AFTER REMOVING EMPTY PLOTS/SPECIES\n")
cat("====================================\n")
cat("Plots   :", nrow(comm_clean), "\n")
cat("Species :", ncol(comm_clean), "\n")

# ------------------------------------------------------------------------------
# 4. FREQUENCY FILTERING FOR SENSITIVITY ANALYSIS
# ------------------------------------------------------------------------------
# Recalculate frequencies on the cleaned matrix
sp_freq_clean <- colSums(comm_clean > 0)

# Generate matrices with different frequency thresholds.
# Note: comm_f4 is selected for the primary analysis in the manuscript. 
# It effectively eliminates stochastic noise without losing diagnostic species 
# of rare plant communities, which is crucial for ICPN syntaxonomic classification.
comm_f3 <- comm_clean[, sp_freq_clean >= 3, drop = FALSE]
comm_f4 <- comm_clean[, sp_freq_clean >= 4, drop = FALSE]
comm_f5 <- comm_clean[, sp_freq_clean >= 5, drop = FALSE]

cat("\n====================================\n")
cat("AFTER FREQUENCY FILTERING\n")
cat("====================================\n")
cat("freq >= 3 :", ncol(comm_f3), "species\n")
cat("freq >= 4 :", ncol(comm_f4), "species\n")
cat("freq >= 5 :", ncol(comm_f5), "species\n")

# ------------------------------------------------------------------------------
# 5. DATA TRANSFORMATION
# ------------------------------------------------------------------------------
# Square-root transformation is applied to reduce the influence of highly 
# dominant species, a standard recommendation for TWINSPAN classification.
mtws_f3 <- sqrt(comm_f3)
mtws_f4 <- sqrt(comm_f4)
mtws_f5 <- sqrt(comm_f5)

# Retain raw (untransformed) filtered versions for potential downstream usage
mtws_raw_f3 <- comm_f3
mtws_raw_f4 <- comm_f4
mtws_raw_f5 <- comm_f5

# ------------------------------------------------------------------------------
# 6. EXPORT PROCESSED DATA
# ------------------------------------------------------------------------------
write.csv(mtws_f3, "mtwinspan_input_freq3_sqrt.csv", row.names = TRUE)
write.csv(mtws_f4, "mtwinspan_input_freq4_sqrt.csv", row.names = TRUE)
write.csv(mtws_f5, "mtwinspan_input_freq5_sqrt.csv", row.names = TRUE)

write.csv(mtws_raw_f3, "mtwinspan_input_freq3_raw.csv", row.names = TRUE)
write.csv(mtws_raw_f4, "mtwinspan_input_freq4_raw.csv", row.names = TRUE)
write.csv(mtws_raw_f5, "mtwinspan_input_freq5_raw.csv", row.names = TRUE)

# ------------------------------------------------------------------------------
# 7. GENERATE FILTER SUMMARY TABLE
# ------------------------------------------------------------------------------
# Create a comparative summary of the data dimensions across thresholds
summary_table <- data.frame(
  variant = c("clean_raw", "freq3_sqrt", "freq4_sqrt", "freq5_sqrt"),
  plots = c(
    nrow(comm_clean),
    nrow(mtws_f3),
    nrow(mtws_f4),
    nrow(mtws_f5)
  ),
  species = c(
    ncol(comm_clean),
    ncol(mtws_f3),
    ncol(mtws_f4),
    ncol(mtws_f5)
  ),
  sparsity = c(
    sum(comm_clean == 0) / (nrow(comm_clean) * ncol(comm_clean)),
    sum(mtws_f3 == 0) / (nrow(mtws_f3) * ncol(mtws_f3)),
    sum(mtws_f4 == 0) / (nrow(mtws_f4) * ncol(mtws_f4)),
    sum(mtws_f5 == 0) / (nrow(mtws_f5) * ncol(mtws_f5))
  )
)

cat("\n====================================\n")
cat("FILTER SUMMARY TABLE\n")
cat("====================================\n")
print(summary_table)

write.csv(summary_table, "mtwinspan_filter_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. INSPECTION OF THE CHOSEN MATRIX (FREQ_4)
# ------------------------------------------------------------------------------
# Detailed inspection of the optimal threshold (freq >= 4) used in the study
cat("\n====================================\n")
cat("SUMMARY OF mtws_f4 VALUES\n")
cat("====================================\n")
print(summary(as.vector(as.matrix(mtws_f4))))

cat("\nQUANTILES OF NON-ZERO VALUES IN mtws_f4\n")
print(
  quantile(
    as.vector(as.matrix(mtws_f4))[as.vector(as.matrix(mtws_f4)) > 0],
    probs = c(0.05, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
  )
)

# Plot distribution of non-zero cover values for the selected matrix
hist(
  as.vector(as.matrix(mtws_f3))[as.vector(as.matrix(mtws_f3)) > 0],
  breaks = 5,
  main = "Distribution of non-zero values in mtws_f3",
  xlab = "sqrt-transformed cover values"
)



hist(
  as.vector(as.matrix(mtws_f4))[as.vector(as.matrix(mtws_f4)) > 0],
  breaks = 5,
  main = "Distribution of non-zero values in mtws_f4",
  xlab = "sqrt-transformed cover values"
)

hist(
  as.vector(as.matrix(mtws_f5))[as.vector(as.matrix(mtws_f5)) > 0],
  breaks = 5,
  main = "Distribution of non-zero values in mtws_f5",
  xlab = "sqrt-transformed cover values"
)









cat("\nSaved files:\n")
cat("- mtwinspan_input_freq3_sqrt.csv\n")
cat("- mtwinspan_input_freq4_sqrt.csv\n")
cat("- mtwinspan_input_freq5_sqrt.csv\n")
cat("- mtwinspan_input_freq3_raw.csv\n")
cat("- mtwinspan_input_freq4_raw.csv\n")
cat("- mtwinspan_input_freq5_raw.csv\n")
cat("- mtwinspan_filter_summary.csv\n")