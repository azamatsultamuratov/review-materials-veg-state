# ==============================================================================
# Script Name: Evaluation of Optimal Cluster Number for Modified TWINSPAN
# Description: This script determines the optimal number of syntaxonomic clusters 
#              (k) by computing three quantitative ecological metrics across a 
#              range of k values (2 to 15): 
#              1. Number of significant indicator species (multipatt analysis)
#              2. Minimum cluster size (to avoid overly fragmented syntaxa)
#              3. Mean silhouette width (to assess cluster cohesion/separation)
#              The analysis is performed for the total species matrix as well as 
#              separated life-form cohorts (Annuals, Perennials, Woody species).
# ==============================================================================

gc()

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(twinspanR)    
library(dplyr)        
library(tidyr)        
library(ggplot2)      
library(cluster)      
library(vegan)        
library(indicspecies) 

cat("\nProcess started. All packages loaded successfully...\n")

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

# Function 1: Count significant indicator species for a given cluster partition
# Uses permutation tests (nperm = 499) to identify statistically significant 
# diagnostic species (alpha = 0.05).
count_indicators <- function(comm, groups, alpha = 0.05, nperm = 499) {
  groups <- as.factor(groups)
  # Ensure valid groupings before running permutation test
  if(length(unique(groups)) < 2 || ncol(comm) < 2) return(NA)
  
  ind <- tryCatch(
    multipatt(comm, groups, func = "IndVal.g", duleg = TRUE,
              control = how(nperm = nperm)),
    error = function(e) NULL
  )
  
  if(is.null(ind)) return(NA)
  ptab <- ind$sign
  if(is.null(ptab) || nrow(ptab) == 0) return(0)
  
  # Return total count of species with p-value <= alpha
  sum(ptab$p.value <= alpha, na.rm = TRUE)
}

# Function 2: Run Modified TWINSPAN iteratively across different k values
# Compiles evaluation metrics into a formatted data frame for ggplot2.
run_curve_set <- function(comm, dataset_name, k_range, cut_levels, min_group_size) {
  res_list <- list()
  
  for(k in k_range) {
    cat("Evaluating", dataset_name, "- partitions (k) =", k, "\n")
    
    # Safety check: bypass if matrix lacks sufficient variance
    if(nrow(comm) < 2 || ncol(comm) < 2) {
      res_list[[paste0("k", k)]] <- data.frame(
        k = k,
        significant_indicator_species = NA,
        min_cluster_size = NA,
        mean_silhouette_width = NA,
        dataset = dataset_name
      )
      next
    }
    
    # Execute Modified TWINSPAN model
    tw <- twinspan(
      com = comm, modif = TRUE, clusters = k,
      cut.levels = cut_levels, min.group.size = min_group_size,
      quiet = TRUE, show.output.on.console = FALSE
    )
    
    # Extract classification vector
    gr <- as.factor(cut(tw))
    
    # Compute Mean Silhouette Width based on Bray-Curtis dissimilarity
    d <- vegdist(comm, method = "bray")
    sil <- silhouette(as.integer(gr), d)
    mean_sil <- mean(sil[, 3], na.rm = TRUE)
    
    # Determine Minimum Cluster Size
    min_size <- min(table(gr))
    
    # Compute count of Significant Indicator Species
    n_ind <- count_indicators(comm, gr, alpha = 0.05, nperm = 499)
    
    # Append results
    res_list[[paste0("k", k)]] <- data.frame(
      k = k,
      significant_indicator_species = n_ind,
      min_cluster_size = min_size,
      mean_silhouette_width = mean_sil,
      dataset = dataset_name
    )
  }
  
  bind_rows(res_list)
}

# ------------------------------------------------------------------------------
# 2. LOAD PREPROCESSED DATA
# ------------------------------------------------------------------------------
# Load frequency-filtered and sqrt-transformed community matrix from Script 01
comm_all <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)

# Load species metadata (traits/life-forms)
sp_code <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 3. CONSTRUCT LIFE-FORM SPECIFIC MATRICES
# ------------------------------------------------------------------------------
# Match species codes between community matrix and metadata
species_in_matrix <- colnames(comm_all)
sp_match <- sp_code %>% filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

# Isolate columns by functional traits to test sensitivity of the classification
an_cols   <- intersect(sp_match$s_code[sp_match$life_form == "An"], colnames(comm_all)) # Annuals
pr_cols   <- intersect(sp_match$s_code[sp_match$life_form == "Pr"], colnames(comm_all)) # Perennials
wo_cols   <- intersect(sp_match$s_code[sp_match$life_form %in% c("Sr", "Ss")], colnames(comm_all)) # Woody

# Helper function to remove empty rows/columns after subsetting
clean_subset <- function(x) {
  if (is.null(x) || ncol(x) == 0) return(NULL)
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (nrow(x) < 2 || ncol(x) < 2) return(NULL)
  x
}

# Generate cleaned data subsets
comm_all <- clean_subset(comm_all)
comm_an  <- clean_subset(comm_all[, an_cols, drop = FALSE])
comm_pr  <- clean_subset(comm_all[, pr_cols, drop = FALSE])
comm_wo  <- clean_subset(comm_all[, wo_cols, drop = FALSE])

cat("\n====================================\n")
cat("DATASET DIMENSIONS FOR EVALUATION\n")
cat("====================================\n")
cat("All species   :", nrow(comm_all), "plots,", ncol(comm_all), "species\n")
cat("Annuals       :", nrow(comm_an),  "plots,", ncol(comm_an),  "species\n")
cat("Perennials    :", nrow(comm_pr),  "plots,", ncol(comm_pr),  "species\n")
cat("Woody species :", nrow(comm_wo),  "plots,", ncol(comm_wo),  "species\n")

# ------------------------------------------------------------------------------
# 4. COMPUTE CLUSTER METRICS ACROSS K VALUES
# ------------------------------------------------------------------------------
# Standard exploratory range for syntaxonomic classification
k_values <- 2:15 
# Pseudospecies cut levels (standard for TWINSPAN)
cut_use  <- c(0, 0.3, 0.7, 2, 4)
min_group_size_use <- 5

# Execute evaluation iterations
curve_all <- run_curve_set(comm_all, "All species",    k_values, cut_use, min_group_size_use)
curve_an  <- run_curve_set(comm_an,  "Annuals",        k_values, cut_use, min_group_size_use)
curve_pr  <- run_curve_set(comm_pr,  "Perennials",     k_values, cut_use, min_group_size_use)
curve_wo  <- run_curve_set(comm_wo,  "Woody species",  k_values, cut_use, min_group_size_use)

# Combine results into a single data frame for plotting
curve_df <- bind_rows(curve_all, curve_an, curve_pr, curve_wo)
curve_df$dataset <- factor(curve_df$dataset, 
                           levels = c("All species", "Annuals", "Perennials", "Woody species"))

# ------------------------------------------------------------------------------
# 5. GRAPHICAL PARAMETERS (JVS PUBLICATION STANDARDS)
# ------------------------------------------------------------------------------
# Define a clean, high-contrast theme appropriate for academic print
jvs_base_theme <- theme_bw() + 
  theme(
    panel.grid = element_blank(),
    
    # Outer frame border line width
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    
    # X and Y axis typography
    axis.title = element_text(size = 14, colour = "black", family = "sans"), 
    axis.text = element_text(size = 14, colour = "black", family = "sans"),
    
    # Legend aesthetics
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.5, "lines"), 
    legend.position = "inside" 
  )

# Define distinct shapes for accessibility (color-blind friendly approach)
shape_values <- c(15, 1, 2, 5) 

# ------------------------------------------------------------------------------
# 6. GENERATE AND EXPORT HIGH-RESOLUTION PLOTS
# ------------------------------------------------------------------------------
cat("\nRendering and exporting publication-quality figures...\n")

# FIGURE 1: Number of Significant Indicator Species vs. k
p1 <- ggplot(curve_df, aes(x = k, y = significant_indicator_species, 
                           group = dataset, shape = dataset, linetype = dataset)) +
  geom_line(colour = "black", linewidth = 0.5) +
  geom_point(colour = "black", size = 2.5) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = unique(curve_df$k)) +
  labs(x = "Number of clusters (k)", y = "Number of significant indicator species") +
  jvs_base_theme +
  # Place legend at the TOP-RIGHT corner
  theme(legend.position.inside = c(0.85, 0.87)) 

# Export complying with JVS requirements: TIFF format, 600 dpi, LZW compression
ggsave("Figure_1_indicators.tif", plot = p1, width = 180, height = 140, units = "mm", 
       dpi = 600, compression = "lzw")

# FIGURE 2: Minimum Cluster Size vs. k
p2 <- ggplot(curve_df, aes(x = k, y = min_cluster_size, 
                           group = dataset, shape = dataset, linetype = dataset)) +
  geom_line(colour = "black", linewidth = 0.5) +
  geom_point(colour = "black", size = 2.5) +
  # Reference line representing the absolute minimum acceptable plot number per syntaxon
  geom_hline(yintercept = 5, linetype = "dotted", color = "grey40") + 
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = unique(curve_df$k)) +
  labs(x = "Number of clusters (k)", y = "Minimum cluster size") +
  jvs_base_theme +
  theme(legend.position.inside = c(0.85, 0.87))

ggsave("Figure_2_min_size.tif", plot = p2, width = 180, height = 140, units = "mm", 
       dpi = 600, compression = "lzw")

# FIGURE 3: Mean Silhouette Width vs. k
p3 <- ggplot(curve_df, aes(x = k, y = mean_silhouette_width, 
                           group = dataset, shape = dataset, linetype = dataset)) +
  geom_line(colour = "black", linewidth = 0.5) +
  geom_point(colour = "black", size = 2.5) +
  scale_shape_manual(values = shape_values) +
  scale_x_continuous(breaks = unique(curve_df$k)) +
  labs(x = "Number of clusters (k)", y = "Mean silhouette width") +
  jvs_base_theme +
  # Shift legend to TOP-LEFT to avoid occluding data trends
  theme(legend.position.inside = c(0.15, 0.87))

ggsave("Figure_3_silhouette.tif", plot = p3, width = 180, height = 140, units = "mm", 
       dpi = 600, compression = "lzw")

cat("\nProcess successfully completed. Files saved to working directory.\n")