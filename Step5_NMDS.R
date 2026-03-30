# ==============================================================================
# Script Name: 3D NMDS Ordination with Robust Ellipses
# Description: Performs 3D Non-metric Multidimensional Scaling (NMDS) to visualize 
#              the floristic relationships and spread of syntaxonomic clusters. 
#              Includes dynamically projected robust ellipses (65% confidence) 
#              and spatial centroids for each vegetation state.
# ==============================================================================

# Clear environment graphics and free memory
graphics.off()
gc()

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(vegan)         
library(dplyr)         
library(tibble)        
library(scatterplot3d) 
library(ellipse)       
library(MASS)          

# ------------------------------------------------------------------------------
# 1. DATA INPUT AND PREPARATION
# ------------------------------------------------------------------------------
# Load cluster membership and raw community data
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw_f4 <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)

# Ensure strict alignment of plots between community data and classification
common_plots <- intersect(rownames(comm_raw_f4), class_df$plot_id)
comm_raw <- comm_raw_f4[common_plots, , drop = FALSE]

class_df <- class_df %>% filter(plot_id %in% common_plots)
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]

# Apply square-root transformation to reduce dominance effect prior to ordination
comm_nmds <- sqrt(comm_raw)

# Remove any plots or species that became empty after filtering/transformation
valid_rows <- rowSums(comm_nmds) > 0
valid_cols <- colSums(comm_nmds) > 0
comm_nmds <- comm_nmds[valid_rows, valid_cols, drop = FALSE]

# Re-align classification dataframe to match the cleaned ordination matrix
class_df <- class_df[match(rownames(comm_nmds), class_df$plot_id), ]

# ------------------------------------------------------------------------------
# 2. NMDS ORDINATION (3 AXES)
# ------------------------------------------------------------------------------
set.seed(123) # Set seed for full reproducibility of the ordination configuration

# Execute NMDS with 3 dimensions using Bray-Curtis dissimilarity
nmds_3 <- metaMDS(
  comm_nmds,
  distance = "bray",
  k = 3,
  trymax = 500,           # High number of random starts to find global optimum
  autotransform = FALSE,  # Transformation already applied manually above
  trace = FALSE
)

cat("NMDS Stress value:", round(nmds_3$stress, 3), "\n")

# Extract plot coordinates (scores) and merge with cluster membership
scores_df <- as.data.frame(scores(nmds_3, display = "sites")) %>%
  rownames_to_column("plot_id") %>%
  left_join(class_df, by = "plot_id")

# Convert lower_class to a correctly ordered factor
scores_df$lower_class <- factor(
  scores_df$lower_class,
  levels = as.character(sort(unique(scores_df$lower_class)))
)

# ------------------------------------------------------------------------------
# 3. VISUALIZATION AESTHETICS AND SETTINGS
# ------------------------------------------------------------------------------
# Define standardized syntaxonomic labels
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", 
                  "5"="D2", "6"="E1", "7"="E2", "8"="F")

# Define distinct color palette matching the dendrograms
cluster_cols <- c("1"="#E76F51", "2"="#DDB7E5", "3"="#8F5D99", "4"="#F1C40F", 
                  "5"="#E68A00", "6"="#A7DCA5", "7"="#4CAF50", "8"="#0B6623")

# Function to calculate robust 2D ellipses based on projected coordinates
# Level = 0.65 captures the core 65% density of the cluster, minimizing outlier overlap
get_robust_ellipse <- function(df_sub, level = 0.65) {
  rb <- cov.trob(cbind(df_sub$px, df_sub$py))
  return(as.data.frame(ellipse(rb$cov, centre = rb$center, level = level)))
}

# ------------------------------------------------------------------------------
# 4. EXPORT PARAMETERS (HIGH-RESOLUTION PRINT)
# ------------------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_name <- paste0("Figure_NMDS_3D_Syntaxa_", timestamp, ".png")

# Initialize high-resolution PNG device (1000 DPI for sharp geometric rendering)
png(file_name, width = 220, height = 180, units = "mm", res = 1000, bg = "white")

# Set graphic margins; tcl = -0.5 ensures axis ticks point strictly outward
par(mar = c(6, 5, 2, 4), xpd = TRUE, tcl = -0.5) 

# Calculate coordinate ranges for all axes
xr <- range(scores_df$NMDS1)
yr <- range(scores_df$NMDS2)
zr <- range(scores_df$NMDS3)

# ------------------------------------------------------------------------------
# 5. RENDER 3D SCATTERPLOT BASE
# ------------------------------------------------------------------------------
# Create the 3D projection frame without drawing points yet (pch = NA)
s3d <- scatterplot3d(
  scores_df$NMDS1, scores_df$NMDS2, scores_df$NMDS3,
  xlim = xr, ylim = yr, zlim = zr, 
  pch = NA, angle = 60, box = FALSE, grid = FALSE, 
  axis = TRUE, xlab = "", ylab = "", zlab = "", 
  main = "",              
  cex.axis = 1.5,        # Axis number size
  font.axis = 1,         # Standard (non-bold) font for axis numbers
  col.axis = "black", 
  tick.marks = TRUE
)

# Convert actual 3D coordinates to 2D screen projection coordinates for plotting shapes
xy <- s3d$xyz.convert(scores_df$NMDS1, scores_df$NMDS2, scores_df$NMDS3)
scores_df$px <- xy$x
scores_df$py <- xy$y

# ------------------------------------------------------------------------------
# 6. DRAW ROBUST COVARIANCE ELLIPSES
# ------------------------------------------------------------------------------
# Overlay semi-transparent ellipses to represent the floristic spread of clusters
for (cl in levels(scores_df$lower_class)) {
  sub <- scores_df[scores_df$lower_class == cl, ]
  # Only draw ellipse if the cluster has more than 5 plots to ensure statistical stability
  if (nrow(sub) > 5) {
    ed <- get_robust_ellipse(sub, level = 0.65)
    polygon(ed, col = adjustcolor(cluster_cols[cl], alpha.f = 0.25), 
            border = cluster_cols[cl], lwd = 2.5) 
  }
}

# ------------------------------------------------------------------------------
# 7. PLOT DATA POINTS AND AXIS VECTORS
# ------------------------------------------------------------------------------
# Overlay individual plots (sites) colored by their syntaxonomic cluster
for (cl in levels(scores_df$lower_class)) {
  ii <- scores_df$lower_class == cl
  points(scores_df$px[ii], scores_df$py[ii], pch = 16, 
         col = adjustcolor(cluster_cols[cl], alpha.f = 0.95), cex = 0.8)
}

# Draw directional axis vectors starting from origin (0,0,0)
origin <- s3d$xyz.convert(0, 0, 0)
v1 <- s3d$xyz.convert(max(xr)*1.1, 0, 0)
v2 <- s3d$xyz.convert(0, max(yr)*1.1, 0)
v3 <- s3d$xyz.convert(0, 0, max(zr)*1.1)

segments(origin$x, origin$y, v1$x, v1$y, col = "red", lwd = 2)
segments(origin$x, origin$y, v2$x, v2$y, col = "red", lwd = 2)
segments(origin$x, origin$y, v3$x, v3$y, col = "red", lwd = 2)

# Label the vectors
text(v1$x, v1$y, "NMDS 1", pos = 4, cex = 1.2, font = 2, col = "red")
text(v2$x, v2$y, "NMDS 2", pos = 4, cex = 1.2, font = 2, col = "red")
text(v3$x, v3$y, "NMDS 3", pos = 3, cex = 1.2, font = 2, col = "red")

# ------------------------------------------------------------------------------
# 8. ADD STRESS METRIC AND CENTROID LABELS
# ------------------------------------------------------------------------------
# Display NMDS stress value in the bottom-left corner
stress_val <- round(nmds_3$stress, 3) 
legend("bottomleft", legend = paste("Stress =", stress_val), 
       bty = "n", cex = 1.5, text.font = 1, inset = c(0.05, 0.02))

# Calculate spatial centroids (mean coordinates) for each cluster
centroids_3d <- scores_df %>% 
  group_by(lower_class) %>% 
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), NMDS3 = mean(NMDS3), .groups = "drop")

# Map syntaxonomic labels to centroids
centroids_3d$label <- class_labels[as.character(centroids_3d$lower_class)]
c_xy <- s3d$xyz.convert(centroids_3d$NMDS1, centroids_3d$NMDS2, centroids_3d$NMDS3)

# Draw distinct centroid markers with labels
points(c_xy$x, c_xy$y, pch = 16, cex = 2.5, col = adjustcolor("white", 0.95))
text(c_xy$x, c_xy$y, centroids_3d$label, cex = 1.1, font = 1)

# Add descriptive caption at the bottom
mtext("Projected robust ellipses (65% CI) indicate the floristic spread of lower clusters.", 
      side = 1, line = 4.5, cex = 1.1, font = 2)

# Close device and finalize export
dev.off()
cat("Process completed! 3D NMDS plot saved to:", file_name, "\n")