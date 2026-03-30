# ==============================================================================
# Script Name: Modified TWINSPAN Hierarchy Visualizations
# Description: Generates publication-ready hierarchical diagrams for syntaxa:
#              1. Vertical Cladogram (Ranked branches format)
#              2. Horizontal Dendrogram (Heterogeneity-based format)
#              Both incorporate proportional bands for cluster sizes.
# ==============================================================================

gc()

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(twinspanR)   
library(dplyr)       
library(dendextend)  

# ------------------------------------------------------------------------------
# 1. LOAD DATA & RUN MODIFIED TWINSPAN
# ------------------------------------------------------------------------------
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
class_df$lower_class <- as.character(class_df$lower_class)

cut_use <- c(0, 0.3, 0.7, 2, 4)

tw_lower <- twinspan(
  com = mtws_f4,
  modif = TRUE,
  clusters = 8,
  cut.levels = cut_use,
  min.group.size = 5,
  diss = "bray",
  quiet = TRUE,
  show.output.on.console = FALSE
)

# ------------------------------------------------------------------------------
# 2. EXTRACT CLUSTER SIZES & DEFINE AESTHETICS
# ------------------------------------------------------------------------------
heter_list <- tw_lower$summary$heter
n_final <- tw_lower$summary$clusters

lower_sizes <- class_df %>%
  count(lower_class, name = "n_plots") %>%
  arrange(as.numeric(lower_class))

if (nrow(lower_sizes) != n_final) {
  stop("Error: Number of lower clusters in class_df does not match tw_lower output.")
}

leaf_labels <- c("A", "B", "C", "D1", "D2", "E1", "E2", "F")

cluster_cols <- c(
  "A"  = "#E76F51", "B"  = "#DDB7E5", "C"  = "#A67CA3",
  "D1" = "#F1C40F", "D2" = "#E68A00", "E1" = "#A7DCA5",
  "E2" = "#4CAF50", "F"  = "#0B6623"
)

# ------------------------------------------------------------------------------
# 3. RECONSTRUCT HIERARCHICAL CLUSTERING OBJECT (hclust)
# ------------------------------------------------------------------------------
current <- as.list(seq_len(n_final))
merge_mat <- matrix(NA_integer_, nrow = n_final - 1, ncol = 2)
height_vec <- numeric(n_final - 1)

next_node_id <- 1

to_hclust_id <- function(x) {
  node_id <- attr(x, "node_id")
  if (is.null(node_id)) return(-as.integer(x))
  return(as.integer(node_id))
}

for (m in n_final:2) {
  prev_step <- heter_list[[m - 1]]
  split_idx <- prev_step$which.most.heter
  split_height <- prev_step$cluster.heter[split_idx]
  
  left  <- current[[split_idx]]
  right <- current[[split_idx + 1]]
  
  left_id  <- to_hclust_id(left)
  right_id <- to_hclust_id(right)
  
  row_id <- next_node_id
  merge_mat[row_id, ] <- c(left_id, right_id)
  height_vec[row_id] <- split_height
  
  new_node <- c(unlist(left), unlist(right))
  attr(new_node, "node_id") <- row_id
  
  current <- append(current[-c(split_idx, split_idx + 1)], list(new_node), after = split_idx - 1)
  next_node_id <- next_node_id + 1
}

hc_mtws <- list(
  merge = merge_mat, height = height_vec, order = seq_len(n_final),
  labels = leaf_labels, method = "reconstructed modified TWINSPAN"
)
class(hc_mtws) <- "hclust"

# Proportions
total_n <- sum(lower_sizes$n_plots)
widths <- (lower_sizes$n_plots / total_n) * 100

# ==============================================================================
# 4. FUNCTION: VERTICAL CLADOGRAM (RANKED BRANCHES FORMAT)
# ==============================================================================
plot_mtws_cladogram <- function(file_png, width_cm = 70, height_cm = 50) {
  
  x_right  <- cumsum(widths)
  x_left   <- c(0, head(x_right, -1))
  x_center <- (x_left + x_right) / 2
  
  node_x <- numeric(nrow(hc_mtws$merge))
  for (i in seq_len(nrow(hc_mtws$merge))) {
    left  <- hc_mtws$merge[i, 1]
    right <- hc_mtws$merge[i, 2]
    x1 <- if (left  < 0) x_center[-left]  else node_x[left]
    x2 <- if (right < 0) x_center[-right] else node_x[right]
    node_x[i] <- (x1 + x2) / 2
  }
  
  # --- CORE MODIFICATION FOR CLADOGRAM (RANKED BRANCHES) ---
  # Replaces raw heterogeneity distances with sequential merge ranks (1, 2, 3...)
  node_y_plot <- seq_len(nrow(hc_mtws$merge))
  y_max <- max(node_y_plot)
  
  png(filename = file_png, width = width_cm / 2.54, height = height_cm / 2.54, 
      units = "in", res = 500, bg = "white")
  
  par(mar = c(10, 7, 4, 4), font = 2, xpd = NA)
  
  band_h <- y_max * 0.12
  label_offset1 <- y_max * 0.08
  label_offset2 <- y_max * 0.20
  
  plot(0, type = "n", xlim = c(0, 100), ylim = c(-label_offset2 - band_h * 0.3, y_max * 1.08),
       axes = FALSE, xlab = "", ylab = "Splitting levels (Ranked branches)",
       main = "Modified TWINSPAN Cladogram")
  
  # Y-axis displays ranked levels instead of continuous distances
  axis(2, at = 1:y_max, labels = 1:y_max, las = 1, cex.axis = 1.2, lwd = 2)
  box(lwd = 1.5)
  
  for (i in seq_len(nrow(hc_mtws$merge))) {
    left  <- hc_mtws$merge[i, 1]
    right <- hc_mtws$merge[i, 2]
    
    x1 <- if (left < 0) x_center[-left] else node_x[left]
    y1 <- if (left < 0) 0 else node_y_plot[left]
    x2 <- if (right < 0) x_center[-right] else node_x[right]
    y2 <- if (right < 0) 0 else node_y_plot[right]
    y_mid <- node_y_plot[i]
    
    segments(x0 = x1, y0 = y1, x1 = x1, y1 = y_mid, lwd = 2.2)
    segments(x0 = x2, y0 = y2, x1 = x2, y1 = y_mid, lwd = 2.2)
    segments(x0 = x1, y0 = y_mid, x1 = x2, y1 = y_mid, lwd = 2.2)
  }
  
  rect(xleft = x_left, ybottom = -band_h, xright = x_right, ytop = 0,
       col = cluster_cols[leaf_labels], border = "black", lwd = 1.2)
  text(x = x_center, y = -band_h - label_offset1, labels = leaf_labels, cex = 2.0, font = 2)
  text(x = x_center, y = -band_h - label_offset2, labels = paste0("n=", lower_sizes$n_plots, "\n(", round(widths, 1), "%)"), cex = 1.1, font = 3)
  
  dev.off()
}

# ==============================================================================
# 5. FUNCTION: HORIZONTAL DENDROGRAM (HETEROGENEITY-BASED)
# ==============================================================================
plot_mtws_horizontal <- function(file_png, width_cm = 80, height_cm = 50) {
  
  options(OutDec = ".") 
  
  y_top_vec    <- cumsum(widths)
  y_bottom_vec <- c(0, head(y_top_vec, -1))
  y_center     <- (y_bottom_vec + y_top_vec) / 2
  
  node_y <- numeric(nrow(hc_mtws$merge))
  for (i in seq_len(nrow(hc_mtws$merge))) {
    left  <- hc_mtws$merge[i, 1]
    right <- hc_mtws$merge[i, 2]
    y1 <- if (left  < 0) y_center[-left]  else node_y[left]
    y2 <- if (right < 0) y_center[-right] else node_y[right]
    node_y[i] <- (y1 + y2) / 2
  }
  
  # Uses actual cluster heterogeneity distances
  node_x_plot <- hc_mtws$height
  x_max <- max(node_x_plot)
  
  png(filename = file_png, width = width_cm / 2.54, height = height_cm / 2.54, 
      units = "in", res = 600, bg = "white")
  
  par(mar = c(11, 2, 2, 8), font = 2, xpd = NA)
  
  band_w <- x_max * 0.01         
  label_offset <- x_max * 0.01
  my_ticks <- pretty(c(0, x_max))
  max_tick <- max(my_ticks) 
  
  plot(0, type = "n", xlim = c(max_tick, -band_w - label_offset * 3), 
       ylim = c(100, 0), axes = FALSE, xlab = "", ylab = "", main = "")
  
  axis(1, at = my_ticks, labels = formatC(my_ticks, format = "f", digits = 1), 
       las = 1, cex.axis = 4.0, lwd = 2.5, lwd.ticks = 2.5, tcl = -0.5, padj = 0.8)      
  
  mtext("Cluster heterogeneity", side = 1, line = 8.5, cex = 5.0, font = 1) 
  
  for (i in seq_len(nrow(hc_mtws$merge))) {
    left  <- hc_mtws$merge[i, 1]
    right <- hc_mtws$merge[i, 2]
    y1 <- if (left < 0) y_center[-left] else node_y[left]
    x1 <- if (left < 0) 0 else node_x_plot[left]
    y2 <- if (right < 0) y_center[-right] else node_y[right]
    x2 <- if (right < 0) 0 else node_x_plot[right]
    x_mid <- node_x_plot[i]
    
    segments(x0 = x1, y0 = y1, x1 = x_mid, y1 = y1, lwd = 2.5)
    segments(x0 = x2, y0 = y2, x1 = x_mid, y1 = y2, lwd = 2.5)
    segments(x0 = x_mid, y0 = y1, x1 = x_mid, y1 = y2, lwd = 2.5)
  }
  
  rect(xleft = 0, ybottom = y_bottom_vec, xright = -band_w, ytop = y_top_vec,
       col = cluster_cols[leaf_labels], border = "black", lwd = 1.2)
  text(x = -band_w - label_offset, y = y_center, labels = leaf_labels, cex = 3.5, font = 1, adj = 0)
  
  dev.off()
}

# ==============================================================================
# 6. EXPORT PLOTS
# ==============================================================================
cat("\nGenerating dendrograms...\n")

# 1. Export Vertical Cladogram (Ranked Branches)
plot_mtws_cladogram("1Cladogram_ModifiedTWINSPAN_Ranked.png")
cat("- Vertical cladogram (ranked) saved.\n")

# 2. Export Horizontal Dendrogram (Heterogeneity-based)
plot_mtws_horizontal("2Dendrogram_ModifiedTWINSPAN_Horizontal.png")
cat("- Horizontal dendrogram (heterogeneity) saved.\n")

cat("\nProcess successfully completed!\n")