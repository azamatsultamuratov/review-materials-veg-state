# ==============================================================================
# Script Name: Detrended Correspondence Analysis (DCA)
# Description: Performs DCA on the full community matrix and functional 
#              life-form subsets (Annuals, Perennials, Woody species). 
#              Generates publication-quality ordinations with dynamically 
#              repelled centroid labels and exports axis summaries.
# ==============================================================================

graphics.off() 
gc()           
# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(vegan)
library(dplyr)
library(tibble)

# Generate timestamp to protect exported files from overwriting
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------------------------
# 1. DATA INPUT & PREPARATION
# ------------------------------------------------------------------------------
# Load classification results, raw community data, and species metadata
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw_f4 <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
sp_code <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)

# Match plots between community matrix and classification dataframe
common_plots <- intersect(rownames(comm_raw_f4), class_df$plot_id)
comm_raw <- comm_raw_f4[common_plots, , drop = FALSE]

class_df <- class_df %>% filter(plot_id %in% common_plots)
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]
class_df$lower_class <- as.character(class_df$lower_class)

# Apply square-root transformation to stabilize variances
comm_all <- sqrt(comm_raw)

# Remove any plots or species that became empty
comm_all <- comm_all[rowSums(comm_all) > 0, , drop = FALSE]
comm_all <- comm_all[, colSums(comm_all) > 0, drop = FALSE]

# Re-align classification dataframe
class_df <- class_df[match(rownames(comm_all), class_df$plot_id), ]

# ------------------------------------------------------------------------------
# 2. SUBSETTING BY LIFE-FORM
# ------------------------------------------------------------------------------
species_in_matrix <- colnames(comm_all)
sp_match <- sp_code %>% filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

# Identify columns by life-form category
an_cols <- sp_match$s_code[sp_match$life_form == "An"]
pr_cols <- sp_match$s_code[sp_match$life_form == "Pr"]
woody_cols <- union(sp_match$s_code[sp_match$life_form == "Sr"], 
                    sp_match$s_code[sp_match$life_form == "Ss"])

# Create subset matrices
com_an    <- comm_all[, intersect(colnames(comm_all), an_cols), drop = FALSE]
com_pr    <- comm_all[, intersect(colnames(comm_all), pr_cols), drop = FALSE]
com_woody <- comm_all[, intersect(colnames(comm_all), woody_cols), drop = FALSE]

# Helper function to remove empty rows/columns from subsets
clean_subset <- function(x) {
  if (is.null(x) || ncol(x) == 0) return(NULL)
  x <- x[rowSums(x) > 0, , drop = FALSE]
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (nrow(x) < 2 || ncol(x) < 2) return(NULL)
  return(x)
}

comm_all  <- clean_subset(comm_all)
com_an    <- clean_subset(com_an)
com_pr    <- clean_subset(com_pr)
com_woody <- clean_subset(com_woody)

# ------------------------------------------------------------------------------
# 3. COLOR PALETTE & LABELS
# ------------------------------------------------------------------------------
# Define standardized syntaxonomic labels and colors
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
cluster_cols <- c("1"="#E76F51", "2"="#DDB7E5", "3"="#A67CA3", "4"="#F1C40F", 
                  "5"="#E68A00", "6"="#A7DCA5", "7"="#4CAF50", "8"="#0B6623")

# ------------------------------------------------------------------------------
# 4. HELPER FUNCTIONS FOR DCA PLOTTING
# ------------------------------------------------------------------------------
# Algorithm to repel overlapping centroid labels dynamically
repel_centroid_labels <- function(df, xcol = "DCA1", ycol = "DCA2", pad_x = 0.07, pad_y = 0.095, max_iter = 450) {
  out <- df
  out$lab_x <- out[[xcol]]
  out$lab_y <- out[[ycol]]
  xr <- range(out[[xcol]], na.rm = TRUE)
  yr <- range(out[[ycol]], na.rm = TRUE)
  dx_thr <- diff(xr) * pad_x
  dy_thr <- diff(yr) * pad_y
  
  if (nrow(out) < 2) return(out)
  for (it in 1:max_iter) {
    moved <- FALSE
    for (i in 1:(nrow(out) - 1)) {
      for (j in (i + 1):nrow(out)) {
        dx <- out$lab_x[j] - out$lab_x[i]
        dy <- out$lab_y[j] - out$lab_y[i]
        if (abs(dx) < dx_thr && abs(dy) < dy_thr) {
          if (dx == 0) dx <- runif(1, -1, 1) * dx_thr * 0.25
          if (dy == 0) dy <- runif(1, -1, 1) * dy_thr * 0.25
          shift_x <- sign(dx) * dx_thr * 0.20
          shift_y <- sign(dy) * dy_thr * 0.20
          out$lab_x[i] <- out$lab_x[i] - shift_x
          out$lab_y[i] <- out$lab_y[i] - shift_y
          out$lab_x[j] <- out$lab_x[j] + shift_x
          out$lab_y[j] <- out$lab_y[j] + shift_y
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  return(out)
}

# Main function to run DCA and generate plot panels
run_dca_panel <- function(comm, class_df, panel_title = "", point_cex = 2, centroid_cex = 3.5,
                          label_cex = 1.18, add_legend = FALSE, use_repel = TRUE, axis_cex = 1.35,
                          metric_cex = 1.40, legend_cex = 1.30, legend_pt_cex = 1.60, 
                          legend_position = "bottomleft", inset = c(0.01, 0.01)) {
  
  cdf <- class_df[match(rownames(comm), class_df$plot_id), ]
  dca <- decorana(comm)
  
  scr <- as.data.frame(scores(dca, display = "sites"))
  scr$plot_id <- rownames(scr)
  scr <- scr %>% left_join(cdf[, c("plot_id", "lower_class")], by = "plot_id")
  
  # Extract Eigenvalues and Standard Deviation Units (SDU)
  ev1 <- dca$evals[1]
  ev2 <- dca$evals[2]
  sdu1 <- diff(range(scr$DCA1, na.rm = TRUE))
  sdu2 <- diff(range(scr$DCA2, na.rm = TRUE))
  
  par(bty = "o", lwd = 0.8, mgp = c(3, 1.5, 0))
  plot(scr$DCA1, scr$DCA2, type = "n", xlab = "", ylab = "", 
       main = panel_title, font.main = 1, cex.main = 3.0, font.axis = 1, cex.axis = axis_cex)
  
  # Display axes metrics (SDU and Eigenvalues)
  mtext(paste0("DCA 1: ", round(sdu1, 2), " SDU; EV ", round(ev1, 2)), side = 1, line = 3.2, cex = metric_cex, font = 1)
  mtext(paste0("DCA 2: ", round(sdu2, 1), " SDU; EV ", round(ev2, 2)), side = 2, line = 3.2, cex = metric_cex, font = 1)
  
  # Plot individual site points
  for (cl in sort(unique(scr$lower_class))) {
    ii <- scr$lower_class == cl
    points(scr$DCA1[ii], scr$DCA2[ii], pch = 21, cex = point_cex, bg = "white", col = cluster_cols[cl], lwd = 1.15)
  }
  
  # Calculate and plot cluster centroids
  centroids <- scr %>% group_by(lower_class) %>% 
    summarise(DCA1 = mean(DCA1, na.rm = TRUE), DCA2 = mean(DCA2, na.rm = TRUE), .groups = "drop")
  
  points(centroids$DCA1, centroids$DCA2, pch = 21, cex = centroid_cex, bg = "white", col = "black", lwd = 1.30)
  centroids$label <- class_labels[centroids$lower_class]
  
  # Add labels to centroids (with or without repel algorithm)
  if (use_repel) {
    centroids_lab <- repel_centroid_labels(centroids, pad_x = 0.07, pad_y = 0.095, max_iter = 450)
    segments(centroids_lab$DCA1, centroids_lab$DCA2, centroids_lab$lab_x, centroids_lab$lab_y, col = "black", lwd = 0.95)
    points(centroids_lab$lab_x, centroids_lab$lab_y, pch = 21, bg = "white", col = "white", cex = 1.55)
    text(centroids_lab$lab_x, centroids_lab$lab_y, labels = centroids_lab$label, cex = label_cex, font = 2)
  } else {
    text(centroids$DCA1, centroids$DCA2, labels = centroids$label, cex = label_cex, font = 2)
  }
  
  # Extract plot boundaries to dynamically position text
  usr <- par("usr")
  text(x = usr[2] - 0.02 * (usr[2] - usr[1]), 
       y = usr[3] + 0.04 * (usr[4] - usr[3]),
       labels = paste0("n = ", nrow(comm), ", spp = ", ncol(comm)), 
       adj = c(1, 0), # 1: aligns to the right, 0: aligns to the bottom
       cex = 1.30, font = 1)
  
  # Add legend if requested
  if (add_legend) {
    present_classes <- as.character(sort(unique(scr$lower_class)))
    legend(legend_position, legend = class_labels[present_classes],
           pch = 21, pt.bg = "white", col = cluster_cols[present_classes],
           pt.lwd = 1.15, pt.cex = legend_pt_cex, bty = "n",
           cex = legend_cex, text.font = 1, inset = inset) 
  }
  
  invisible(list(dca = dca, scores = scr, centroids = centroids))
}

# ------------------------------------------------------------------------------
# 5. EXPORT: MAIN ARTICLE FIGURE
# ------------------------------------------------------------------------------
cat("\nGenerating Main DCA Ordination...\n")
main_fig_name <- paste0("Figure_3_Main_DCA_", timestamp, ".png")
png(main_fig_name, width = 180, height = 160, units = "mm", res = 1000, bg = "white")

par(mar = c(5.5, 5.5, 2.5, 1.0), bg = "white", font = 1, font.main = 1, font.lab = 1, font.axis = 1)

res_all <- run_dca_panel(
  comm = comm_all, class_df = class_df, panel_title = "", 
  point_cex = 1.15, centroid_cex = 2.95, label_cex = 1.20,
  add_legend = TRUE, use_repel = FALSE, axis_cex = 1.42,
  metric_cex = 1.2, legend_cex = 1.2, legend_pt_cex = 1.30,
  legend_position = "topright", inset = c(0.01, -0.02)
)

dev.off()

# ------------------------------------------------------------------------------
# 6. EXPORT: APPENDIX FIGURES (Life-Form Subsets)
# ------------------------------------------------------------------------------
cat("Generating Life-Form Subsets DCA (Appendix)...\n")
app_fig_name <- paste0("Appendix_DCA_LifeForms_", timestamp, ".png")
png(app_fig_name, width = 160, height = 360, units = "mm", res = 1000, bg = "white") 

par(mfrow = c(3, 1), mar = c(5.0, 6.0, 4.0, 1.5), oma = c(1.5, 1.0, 1.0, 1.0),
    bg = "white", font = 1, font.main = 1, font.lab = 1, font.axis = 1)

res_an <- run_dca_panel(com_an, class_df, "Annuals", add_legend = TRUE, legend_position = "bottomleft", axis_cex = 2.5, metric_cex = 1.5)
res_pr <- run_dca_panel(com_pr, class_df, "Perennials", add_legend = FALSE, axis_cex = 2.5, metric_cex = 1.5)
res_woody <- run_dca_panel(com_woody, class_df, "Woody species", add_legend = FALSE, axis_cex = 2.5, metric_cex = 1.5)

dev.off()

# ------------------------------------------------------------------------------
# 7. DATA EXPORT: SCORES AND AXIS SUMMARIES
# ------------------------------------------------------------------------------
cat("Exporting scores and summaries to CSV...\n")
write.csv(res_all$scores,   paste0("DCA_scores_all_main_", timestamp, ".csv"), row.names = FALSE)
write.csv(res_an$scores,    paste0("DCA_scores_annual_app_", timestamp, ".csv"), row.names = FALSE)
write.csv(res_pr$scores,    paste0("DCA_scores_perennial_app_", timestamp, ".csv"), row.names = FALSE)
write.csv(res_woody$scores, paste0("DCA_scores_woody_app_", timestamp, ".csv"), row.names = FALSE)

dca_summary <- data.frame(
  dataset = c("All", "Annual", "Perennial", "Woody"),
  axis1_eigen = c(res_all$dca$evals[1], res_an$dca$evals[1], res_pr$dca$evals[1], res_woody$dca$evals[1]),
  axis2_eigen = c(res_all$dca$evals[2], res_an$dca$evals[2], res_pr$dca$evals[2], res_woody$dca$evals[2])
)

write.csv(dca_summary, paste0("DCA_axis_summary_", timestamp, ".csv"), row.names = FALSE)

# Print confirmation
cat("\nSuccess! All DCA plots and CSV summaries exported safely.\n")



