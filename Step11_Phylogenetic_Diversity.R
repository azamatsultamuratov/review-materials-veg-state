
graphics.off()
gc()

# Load required packages
library(V.PhyloMaker2)
library(ape)
library(picante)
library(ggplot2)
library(patchwork)
library(dplyr)

# -------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -------------------------------------------------------------------
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
sp_code <- read.csv("su_gen_df_species_code.csv", stringsAsFactors = FALSE)
print(class_df)
str(class_df)
print(comm_raw)
str(comm_raw)
print(sp_code)
str(sp_code)

# Convert species codes to full scientific names
col_indices <- match(colnames(comm_raw), sp_code$s_code)
colnames(comm_raw) <- sp_code$sp_name[col_indices]

# Remove empty plots and species (zero-sums)
comm_raw <- comm_raw[rowSums(comm_raw) > 0, colSums(comm_raw) > 0]

# Standardize cluster labels from class_df (A, B, C...)
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
class_df$Cluster_Name <- class_labels[as.character(class_df$lower_class)]

# Retain only common plots present in both datasets
common_plots <- intersect(rownames(comm_raw), class_df$plot_id)
comm_clean <- comm_raw[common_plots, , drop = FALSE]
class_clean <- class_df[match(common_plots, class_df$plot_id), ]

# -------------------------------------------------------------------
# 2. PHYLOGENETIC TREE CONSTRUCTION
# -------------------------------------------------------------------
sp_list <- data.frame(
  species = gsub("_", " ", sp_code$sp_name),
  genus = sapply(strsplit(sp_code$sp_name, "_"), `[`, 1),
  family = NA, 
  stringsAsFactors = FALSE
)

cat("Constructing phylogenetic mega-tree...\n")
tree_res <- phylo.maker(sp.list = sp_list, scenarios = "S3")
full_tree <- tree_res$scenario.3

# -------------------------------------------------------------------
# 3. CALCULATE PLOT-LEVEL PHYLOGENETIC INDICES
# -------------------------------------------------------------------
matched <- match.phylo.comm(full_tree, comm_clean)
phy_matched <- matched$phy
comm_matched <- matched$comm

# A. Faith's Phylogenetic Diversity (PD) per plot
cat("Calculating Faith's PD...\n")
pd_plot <- pd(comm_matched, phy_matched, include.root = TRUE)
pd_plot$plot_id <- rownames(pd_plot)

# B. Standardized Effect Size of Mean Pairwise Distance (SES.MPD) per plot
cat("Calculating SES.MPD (Ecological filtering). This may take a while...\n")
phydist <- cophenetic(phy_matched)
nri_plot <- ses.mpd(comm_matched, phydist, null.model = "taxa.labels", 
                    abundance.weighted = TRUE, runs = 999)
nri_plot$plot_id <- rownames(nri_plot)

# Merge results into a single dataframe
df_metrics <- pd_plot %>%
  left_join(nri_plot %>% select(plot_id, mpd.obs.z, mpd.obs.p), by = "plot_id") %>%
  left_join(class_clean %>% select(plot_id, Cluster_Name), by = "plot_id") %>%
  rename(SES_MPD = mpd.obs.z, P_Value = mpd.obs.p)

# Order clusters sequentially
df_metrics$Cluster_Name <- factor(df_metrics$Cluster_Name, levels = c("A", "B", "C", "D1", "D2", "E1", "E2", "F"))

# -------------------------------------------------------------------
# 4. DRAW BOXPLOT GRAPHICS (MATCHED TO RICHNESS DESIGN)
# -------------------------------------------------------------------
cluster_cols <- c("A"="#E76F51", "B"="#DDB7E5", "C"="#A67CA3", "D1"="#F1C40F", 
                  "D2"="#E68A00", "E1"="#A7DCA5", "E2"="#4CAF50", "F"="#0B6623")

# A. Faith's PD Boxplot (Without jitter points, mean value marked by black diamond)
p_pd <- ggplot(df_metrics, aes(x = Cluster_Name, y = PD, fill = Cluster_Name)) +
  geom_boxplot(width = 0.72, linewidth = 0.85, fatten = 1.6,
               outlier.shape = 16, outlier.size = 2.5, outlier.alpha = 0.95, outlier.stroke = 0) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.2, colour = "black") +
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  labs(title = "A", y = "Faith's PD", x = "") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18),
    axis.line = element_line(colour = "black", linewidth = 0.75),
    axis.ticks = element_line(colour = "black", linewidth = 0.75),
    axis.text = element_text(size = 17, colour = "black"),
    axis.title = element_text(size = 18, colour = "black"),
    plot.margin = margin(8, 8, 8, 8)
  )

# B. SES.MPD Boxplot (Includes null model threshold lines and black diamond mean)
p_nri <- ggplot(df_metrics, aes(x = Cluster_Name, y = SES_MPD, fill = Cluster_Name)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 1) +
  geom_hline(yintercept = -1.96, linetype = "dotted", color = "red", linewidth = 0.8) +
  geom_boxplot(width = 0.72, linewidth = 0.85, fatten = 1.6,
               outlier.shape = 16, outlier.size = 2.5, outlier.alpha = 0.95, outlier.stroke = 0) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.2, colour = "black") +
  scale_fill_manual(values = cluster_cols, drop = FALSE) +
  labs(title = "B", y = "SES.MPD (Z-score)", x = "Clusters") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18),
    axis.line = element_line(colour = "black", linewidth = 0.75),
    axis.ticks = element_line(colour = "black", linewidth = 0.75),
    axis.text = element_text(size = 17, colour = "black"),
    axis.title = element_text(size = 18, colour = "black"),
    plot.margin = margin(8, 8, 8, 8)
  )

# Combine both plots vertically
final_plot <- p_pd / p_nri

# Export final plot
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
ggsave(paste0("Phylogenetic_Metrics_Boxplots_", timestamp, ".png"), 
       plot = final_plot, width = 11, height = 11, dpi = 800, bg = "white")

print(final_plot)
cat("Success! Plot-level phylogenetic boxplots have been saved.\n")

# Save summary mean values for reviewer reference
summary_stats_initial <- df_metrics %>%
  group_by(Cluster_Name) %>%
  summarise(
    Mean_PD = mean(PD, na.rm = TRUE),
    Mean_SES_MPD = mean(SES_MPD, na.rm = TRUE),
    Plots_Clustered_Count = sum(SES_MPD < -1.96, na.rm = TRUE),
    Total_Plots = n()
  )
write.csv(summary_stats_initial, paste0("Phylogenetic_Metrics_Summary_Initial_", timestamp, ".csv"), row.names = FALSE)


# ==============================================================================
# 5. EXPORT RESULTS FOR MANUSCRIPT AND REVIEWERS
# ==============================================================================
cat("\nExporting results to Excel (CSV) format. Please wait...\n")

# -------------------------------------------------------------------
# 1. MAIN PLOT-LEVEL METRICS (Raw Data)
# -------------------------------------------------------------------
export_df_metrics <- df_metrics %>%
  select(plot_id, Cluster_Name, SR, PD, SES_MPD, P_Value) %>%
  arrange(Cluster_Name, plot_id)

write.csv(export_df_metrics, paste0("1_Plot_Level_Phylogenetic_Metrics_", timestamp, ".csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 2. OVERALL STATISTICS BY CLUSTER (Summary Table)
# -------------------------------------------------------------------
summary_stats <- df_metrics %>%
  group_by(Cluster_Name) %>%
  summarise(
    Total_Plots = n(),
    Mean_Species_Richness = round(mean(SR, na.rm = TRUE), 2),
    Mean_Faiths_PD = round(mean(PD, na.rm = TRUE), 2),
    Mean_SES_MPD = round(mean(SES_MPD, na.rm = TRUE), 3),
    Plots_Significantly_Clustered = sum(SES_MPD < -1.96, na.rm = TRUE),
    Percent_Clustered = round((sum(SES_MPD < -1.96, na.rm = TRUE) / n()) * 100, 1)
  )

write.csv(summary_stats, paste0("2_Cluster_Summary_Statistics_", timestamp, ".csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 3. SPECIES MEGA-TREE MATCHING STATUS
# -------------------------------------------------------------------
matching_status <- tree_res$species.list
write.csv(matching_status, paste0("3_Species_Megatree_Matching_Status_", timestamp, ".csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4. EVOLUTIONARY DISTINCTIVENESS (ED)
# -------------------------------------------------------------------
cat("Calculating Evolutionary Distinctiveness (ED)...\n")
ed_results <- evol.distinct(full_tree, type = "fair.proportion")

# Sort species by highest evolutionary distinctiveness
ed_sorted <- ed_results[order(-ed_results$w), ]
names(ed_sorted) <- c("Species_su_name", "Evolutionary_Distinctiveness_Score")

write.csv(ed_sorted, paste0("4_Species_Evolutionary_Distinctiveness_", timestamp, ".csv"), row.names = FALSE)

# ==============================================================================
cat("\nALL REQUIRED DATA HAS BEEN SUCCESSFULLY SAVED!\n")
cat("You can find the exported files in your Working Directory.\n")


# ==============================================================================
# Script: Phylogenetic Tree with Species Frequency Heatmap
# Description: Generates a circular phylogenetic tree paired with an outer 
#              heatmap showing species constancy (frequency) across syntaxa.
# ==============================================================================
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ape) # Required for safe tree pruning

cat("\nAligning phylogenetic tree with heatmap...\n")

# 1. CLUSTER COLORS
cluster_colors <- c(
  "1" = "#E76F51", "2" = "#DDB7E5", "3" = "#A67CA3", "4" = "#F1C40F",
  "5" = "#E68A00", "6" = "#A7DCA5", "7" = "#4CAF50", "8" = "#0B6623"
)

# 2. PREPARE FREQUENCY DATA
synoptic_data <- read.csv(
  "Table2_final_synoptic_table_with_other_species.csv",
  check.names = FALSE
)

synoptic_data$Species_Code <- sp_code$s_code[
  match(synoptic_data$Species, sp_code$sp_name)
]

freq_mat <- synoptic_data %>%
  select(Species_Code, A, B, C, D1, D2, E1, E2, F) %>%
  filter(!is.na(Species_Code)) %>%
  distinct(Species_Code, .keep_all = TRUE) %>%
  tibble::column_to_rownames("Species_Code")

apply_10_interval <- function(x) {
  case_when(
    is.na(x) | x == 0 ~ "0",
    x >= 1  & x <= 10  ~ "1-10 %",
    x >= 11 & x <= 20  ~ "11-20 %",
    x >= 21 & x <= 30  ~ "21-30 %",
    x >= 31 & x <= 40  ~ "31-40 %",
    x >= 41 & x <= 50  ~ "41-50 %",
    x >= 51 & x <= 60  ~ "51-60 %",
    x >= 61 & x <= 70  ~ "61-70 %",
    x >= 71 & x <= 80  ~ "71-80 %",
    x >= 81 & x <= 90  ~ "81-90 %",
    x >= 91 & x <= 100 ~ "91-100 %"
  )
}

freq_data_grad <- as.data.frame(lapply(freq_mat, apply_10_interval))
rownames(freq_data_grad) <- rownames(freq_mat)

levels_order <- c(
  "0", "1-10 %", "11-20 %", "21-30 %", "31-40 %",
  "41-50 %", "51-60 %", "61-70 %", "71-80 %",
  "81-90 %", "91-100 %"
)
freq_data_grad[] <- lapply(freq_data_grad, factor, levels = levels_order)

# 3. COLOR GRADIENT
image_colors <- colorRampPalette(
  c("#F7FCB9", "#ADDD8E", "#41AB5D", "#238443", "#004529", "#002000")
)(10)

# 0 frequency mapped to white (blank)
final_palette <- c("white", image_colors)
names(final_palette) <- levels_order

# 4. PREPARE PHYLOGENETIC TREE (SAFE PRUNING ADDED)
idx_tree <- match(phy_matched$tip.label, sp_code$sp_name)
phy_codes_final <- phy_matched
phy_codes_final$tip.label <- sp_code$s_code[idx_tree]

# =========================================================
# CORE FIX: Safe tree pruning to prevent structural breaks
# =========================================================
# 1. Compare species between tree tips and heatmap matrix
tips_in_tree <- phy_codes_final$tip.label
rows_in_data <- rownames(freq_data_grad)

# 2. Identify tips in the tree that are missing from the data
missing_in_data <- setdiff(tips_in_tree, rows_in_data)

# 3. Officially drop surplus tips using ape::drop.tip
if(length(missing_in_data) > 0) {
  phy_codes_final <- ape::drop.tip(phy_codes_final, missing_in_data)
}
# =========================================================

# Proceed with hierarchical clustering on the safely pruned tree
dist_tree <- cophenetic(phy_codes_final)
groups_codes <- cutree(hclust(as.dist(dist_tree)), k = 8)
split_codes <- split(names(groups_codes), groups_codes)
phy_tree_final <- groupOTU(phy_codes_final, split_codes)

# 5. CRITICAL STEP: ALIGN HEATMAP ROWS TO TREE TIP ORDER
tip_order <- phy_tree_final$tip.label

# Reorder heatmap data rows to match tree tip sequences exactly
freq_data_grad <- freq_data_grad[tip_order, , drop = FALSE]

# Validation check
if (!identical(rownames(freq_data_grad), phy_tree_final$tip.label)) {
  stop("Error: Heatmap row order does not match tree tip order.")
}

# 6. DRAW MAIN TREE
p_main <- ggtree(
  phy_tree_final,
  layout = "circular",
  size = 1.8,
  aes(color = group)
) +
  geom_tiplab(size = 7.5, fontface = "bold", offset = 0.5) +
  scale_color_manual(values = cluster_colors, guide = "none") +
  xlim(-25, NA) +
  theme(legend.position = "none")

# 7. ADD HEATMAP LAYER
q1_final <- gheatmap(
  p_main,
  freq_data_grad,
  offset = 60,
  width = 0.5,
  colnames = TRUE,
  colnames_position = "top",
  colnames_angle = 0,
  colnames_offset_y = 2,
  font.size = 10,
  color = "gray10"
) +
  scale_fill_manual(
    values = final_palette,
    name = "Frequency",
    drop = FALSE,
    guide = guide_legend(
      keywidth = unit(1.5, "cm"),
      keyheight = unit(1.5, "cm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 55),
    legend.text = element_text( size = 55),
    axis.text.x = element_text(face = "bold", size = 55),
    legend.key = element_rect(fill = "white", color = "gray10"),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(1.5, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# THICKEN HEATMAP TILE BORDERS
tile_layers <- which(vapply(q1_final$layers, function(x) inherits(x$geom, "GeomTile"), logical(1)))

if (length(tile_layers) > 0) {
  hm_id <- tile_layers[length(tile_layers)]
  q1_final$layers[[hm_id]]$aes_params$linewidth <- 0.6
  q1_final$layers[[hm_id]]$aes_params$size <- 0.6
  q1_final$layers[[hm_id]]$aes_params$colour <- "gray10"
  q1_final$layers[[hm_id]]$aes_params$color  <- "gray10"
}

# 8. EXPORT GRAPHIC
ggsave(
  "6666Phylogenetic_Heatmap_Final_CustomSpace_fixed.png",
  plot = q1_final,
  width = 40,
  height = 40,
  units = "in",
  dpi = 600,
  bg = "white"
)

cat("\nDone! '6666Phylogenetic_Heatmap_Final_CustomSpace_fixed.png' has been saved.\n")