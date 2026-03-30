# ==============================================================================
# Script Name: Species Richness & Vegetation Cover Boxplots (ANOVA)
# Description: Calculates plot-level ecological metrics (Richness and Cover) 
#              stratified by life-form (Total, Annual, Perennial, Woody).
#              Performs One-way ANOVA to test for significant differences across 
#              syntaxonomic clusters and generates multi-panel boxplots.
# ==============================================================================
gc()

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(dplyr)    
library(tidyr)     
library(ggplot2)  
library(patchwork) 

# Generate timestamp to protect exported files from overwriting
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------------------------
# 1. DATA INPUT & PREPARATION
# ------------------------------------------------------------------------------
# Load classification results, raw community data, and species metadata
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
sp_code  <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)

# Match plots strictly between community matrix and classification dataframe
common_plots <- intersect(rownames(comm_raw), class_df$plot_id)
comm_raw <- comm_raw[common_plots, , drop = FALSE]

class_df <- class_df %>% filter(plot_id %in% common_plots)
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]

# Format cluster labels (A to F)
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
class_df$lower_class <- as.character(class_df$lower_class)
class_df$cluster_label <- factor(class_labels[class_df$lower_class], 
                                 levels = c("A", "B", "C", "D1", "D2", "E1", "E2", "F"))

# Match species metadata to columns in the community matrix
species_in_matrix <- colnames(comm_raw)
sp_match <- sp_code %>% filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

# Safety checks
if(any(is.na(sp_match$s_code))) stop("Error: Some species in comm_raw are missing in df_species_code.csv")
if(!"life_form" %in% names(sp_match)) stop("Error: Column 'life_form' not found in df_species_code.csv")

# Define life-form groups (Merging Subshrubs and Shrubs into a single 'Woody' category)
all_cols <- colnames(comm_raw)
an_cols <- intersect(sp_match$s_code[sp_match$life_form == "An"], all_cols)
pr_cols <- intersect(sp_match$s_code[sp_match$life_form == "Pr"], all_cols)
woody_cols <- union(
  intersect(sp_match$s_code[sp_match$life_form == "Sr"], all_cols),
  intersect(sp_match$s_code[sp_match$life_form == "Ss"], all_cols)
)

# ------------------------------------------------------------------------------
# 2. GLOBAL SETTINGS & HELPER FUNCTIONS
# ------------------------------------------------------------------------------
# Define color palette consistent with ordination and dendrograms
cluster_cols <- c("A"="#E76F51", "B"="#DDB7E5", "C"="#A67CA3", "D1"="#F1C40F", 
                  "D2"="#E68A00", "E1"="#A7DCA5", "E2"="#4CAF50", "F"="#0B6623")

# Function to safely calculate row sums (returns 0 if columns list is empty)
safe_row_sums <- function(mat, cols) {
  if(length(cols) == 0) return(rep(0, nrow(mat)))
  rowSums(mat[, cols, drop = FALSE], na.rm = TRUE) 
}

# Function to extract One-way ANOVA p-value
calc_p <- function(df, val_col) {
  f <- as.formula(paste(val_col, "~ cluster_label"))
  fit <- aov(f, data = df)
  p <- summary(fit)[[1]][["Pr(>F)"]][1]
  if(is.na(p)) return("p = NA")
  if(p < 0.001) return("p < 0.001")
  paste0("p = ", format(round(p, 3), nsmall = 3))
}

# Unified function to render standardized boxplots
make_metric_plot <- function(df_sub, y_lab, val_col, show_x_title = TRUE) {
  ggplot(df_sub, aes_string(x = "cluster_label", y = val_col, fill = "cluster_label")) +
    geom_boxplot(width = 0.72, linewidth = 0.85, fatten = 1.6,
                 outlier.shape = 16, outlier.size = 2.5, outlier.alpha = 0.95, outlier.stroke = 0) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3.2, colour = "black") +
    scale_fill_manual(values = cluster_cols, drop = FALSE) +
    labs(x = if (show_x_title) "Clusters" else NULL, y = y_lab) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none",
      axis.line = element_line(colour = "black", linewidth = 0.75),
      axis.ticks = element_line(colour = "black", linewidth = 0.75),
      axis.text = element_text(size = 17, colour = "black"),
      axis.title = element_text(size = 18, colour = "black"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# ==============================================================================
# PART A: SPECIES RICHNESS ANALYSIS (Alpha Diversity)
# ==============================================================================
cat("\nAnalyzing Species Richness...\n")

# Convert cover data to presence/absence matrix (binary)
comm_pa <- (comm_raw > 0) * 1

# Calculate richness per plot stratified by life-form
rich_df <- data.frame(
  plot_id = rownames(comm_pa),
  total_richness     = rowSums(comm_pa > 0, na.rm = TRUE),
  perennial_richness = safe_row_sums(comm_pa, pr_cols),
  annual_richness    = safe_row_sums(comm_pa, an_cols),
  woody_richness     = safe_row_sums(comm_pa, woody_cols),
  stringsAsFactors = FALSE
) %>% left_join(class_df[, c("plot_id", "cluster_label")], by = "plot_id")

# Reshape data to long format for ggplot
plot_df_rich <- rich_df %>%
  pivot_longer(cols = ends_with("richness"), names_to = "metric", values_to = "richness")

# Standardize facet labels
rich_labels <- c("total_richness" = "Total species richness",
                 "perennial_richness" = "Perennial species richness",
                 "annual_richness" = "Annual species richness",
                 "woody_richness" = "Woody species richness")

plot_df_rich$metric <- factor(plot_df_rich$metric, levels = names(rich_labels), labels = rich_labels)

# Calculate Descriptive Statistics & ANOVA p-values
summary_rich <- plot_df_rich %>% group_by(metric, cluster_label) %>%
  summarise(n = n(), mean = mean(richness), sd = sd(richness), min = min(richness), max = max(richness), .groups="drop")

pvals_rich <- plot_df_rich %>% group_by(metric) %>% 
  group_modify(~ data.frame(p_label = calc_p(.x, "richness")))

# Export Results
write.csv(summary_rich, paste0("1_Richness_Summary_", timestamp, ".csv"), row.names = FALSE)
write.csv(pvals_rich, paste0("1_Richness_ANOVA_pvals_", timestamp, ".csv"), row.names = FALSE)

# Generate Individual Panels
met_r <- unique(plot_df_rich$metric)
p_r1 <- make_metric_plot(filter(plot_df_rich, metric == met_r[1]), met_r[1], "richness", FALSE)
p_r2 <- make_metric_plot(filter(plot_df_rich, metric == met_r[2]), met_r[2], "richness", FALSE)
p_r3 <- make_metric_plot(filter(plot_df_rich, metric == met_r[3]), met_r[3], "richness", TRUE)
p_r4 <- make_metric_plot(filter(plot_df_rich, metric == met_r[4]), met_r[4], "richness", TRUE)

# Combine using patchwork and export
p_final_rich <- (p_r1 + p_r2) / (p_r3 + p_r4)
ggsave(paste0("1_Figure_Species_Richness_", timestamp, ".png"), plot = p_final_rich, width = 11, height = 11, dpi = 800)
cat("- Species Richness multi-panel figure saved.\n")


# ==============================================================================
# PART B: VEGETATION COVER ANALYSIS
# ==============================================================================
cat("Analyzing Vegetation Cover...\n")

# Calculate total cover (using raw abundance values) stratified by life-form
cover_df <- data.frame(
  plot_id = rownames(comm_raw),
  total_cover     = rowSums(comm_raw, na.rm = TRUE),
  perennial_cover = safe_row_sums(comm_raw, pr_cols),
  annual_cover    = safe_row_sums(comm_raw, an_cols),
  woody_cover     = safe_row_sums(comm_raw, woody_cols),
  stringsAsFactors = FALSE
) %>% left_join(class_df[, c("plot_id", "cluster_label")], by = "plot_id")

# Reshape data to long format
plot_df_cov <- cover_df %>%
  pivot_longer(cols = ends_with("cover"), names_to = "metric", values_to = "cover_value")

# Standardize facet labels
cov_labels <- c("total_cover" = "Total vegetation cover",
                "perennial_cover" = "Perennial species cover",
                "annual_cover" = "Annual species cover",
                "woody_cover" = "Woody species cover")

plot_df_cov$metric <- factor(plot_df_cov$metric, levels = names(cov_labels), labels = cov_labels)

# Calculate Descriptive Statistics & ANOVA p-values
summary_cov <- plot_df_cov %>% group_by(metric, cluster_label) %>%
  summarise(n = n(), mean = mean(cover_value), sd = sd(cover_value), min = min(cover_value), max = max(cover_value), .groups="drop")

pvals_cov <- plot_df_cov %>% group_by(metric) %>% 
  group_modify(~ data.frame(p_label = calc_p(.x, "cover_value")))

# Export Results
write.csv(summary_cov, paste0("2_Cover_Summary_", timestamp, ".csv"), row.names = FALSE)
write.csv(pvals_cov, paste0("2_Cover_ANOVA_pvals_", timestamp, ".csv"), row.names = FALSE)

# Generate Individual Panels
met_c <- unique(plot_df_cov$metric)
p_c1 <- make_metric_plot(filter(plot_df_cov, metric == met_c[1]), met_c[1], "cover_value", FALSE)
p_c2 <- make_metric_plot(filter(plot_df_cov, metric == met_c[2]), met_c[2], "cover_value", FALSE)
p_c3 <- make_metric_plot(filter(plot_df_cov, metric == met_c[3]), met_c[3], "cover_value", TRUE)
p_c4 <- make_metric_plot(filter(plot_df_cov, metric == met_c[4]), met_c[4], "cover_value", TRUE)

# Combine using patchwork and export
p_final_cov <- (p_c1 + p_c2) / (p_c3 + p_c4)
ggsave(paste0("2_Figure_Vegetation_Cover_", timestamp, ".png"), plot = p_final_cov, width = 11, height = 11, dpi = 800)
cat("- Vegetation Cover multi-panel figure saved.\n")

cat("\nProcess successfully completed. All output files saved.\n")