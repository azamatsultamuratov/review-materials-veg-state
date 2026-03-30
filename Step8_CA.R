# ==============================================================================
# Script Name: Interactive CA Biplot Visualization (Plotly)
# Description: Performs Correspondence Analysis (CA) to explore relationships 
#              between syntaxonomic clusters and ecological types (eco_tp). 
#              Generates a highly customizable, interactive HTML biplot and 
#              allows for high-resolution static PNG exports.
# ==============================================================================

gc()

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(dplyr)       
library(tidyr)       
library(vegan)       
library(plotly)      
library(htmlwidgets) 

# Generate timestamp for safe file exports
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# ------------------------------------------------------------------------------
# 1. DATA INPUT & PREPARATION
# ------------------------------------------------------------------------------
# Load cluster membership and ecological traits data
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
eco_df <- read.csv("eco_tp.csv", stringsAsFactors = FALSE)

# Convert IDs and identifiers to character strings to prevent merging errors
class_df$plot_id <- as.character(class_df$plot_id)
class_df$lower_class <- as.character(class_df$lower_class)
eco_df$plot_id <- as.character(eco_df$plot_id)
eco_df$Eco_typ <- as.character(eco_df$Eco_typ)

# Format cluster labels
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
cluster_levels <- c("A", "B", "C", "D1", "D2", "E1", "E2", "F")
class_df$cluster_label <- factor(class_labels[class_df$lower_class], levels = cluster_levels)

# Define logical ordering for Ecological Types
eco_levels <- c("Gp", "PoGa", "PsGp", "GaGp", "Ps", "GpPs", "Ga", "Po", "GaPo", 
                "PsPo", "GaPs", "PoPsGa", "PoPs", "PsGa", "PoGpGa", "GaGpPs", "GaPsPo")
eco_df$Eco_typ <- factor(eco_df$Eco_typ, levels = eco_levels)

# ------------------------------------------------------------------------------
# 2. DATA JOINING AND CONTINGENCY TABLE CREATION
# ------------------------------------------------------------------------------
# Merge datasets
eco_join <- eco_df %>% inner_join(class_df[, c("plot_id", "cluster_label")], by = "plot_id")

# Create a contingency table (cross-tabulation) of Clusters vs. Eco Types
tab_df <- eco_join %>% 
  count(cluster_label, Eco_typ, name = "n") %>% 
  complete(
    cluster_label = factor(cluster_levels, levels = cluster_levels), 
    Eco_typ = factor(eco_levels, levels = eco_levels), 
    fill = list(n = 0)
  )

# Reshape to wide matrix format for vegan processing
tab_wide <- tab_df %>% pivot_wider(names_from = Eco_typ, values_from = n)

tab_mat <- as.data.frame(tab_wide)
rownames(tab_mat) <- tab_mat$cluster_label
tab_mat$cluster_label <- NULL
tab_mat <- as.matrix(tab_mat) # Final matrix for CA

# ------------------------------------------------------------------------------
# 3. CORRESPONDENCE ANALYSIS (CA)
# ------------------------------------------------------------------------------
# Execute CA using vegan's cca() function without constraining variables
ca_mod <- cca(tab_mat)

# Extract eigenvalues to calculate axis variance explanation percentages
eig_vals <- ca_mod$CA$eig
eig_df <- data.frame(
  axis = paste0("CA", seq_along(eig_vals)), 
  eigenvalue = eig_vals, 
  percent = eig_vals / sum(eig_vals) * 100
)

# ------------------------------------------------------------------------------
# 4. EXTRACT ORDINATION SCORES
# ------------------------------------------------------------------------------
# Scaling 2: Emphasizes relationships between species (Eco_types)
row_scores <- as.data.frame(scores(ca_mod, display = "sites", scaling = 2))
row_scores$cluster <- rownames(row_scores)

col_scores <- as.data.frame(scores(ca_mod, display = "species", scaling = 2))
col_scores$eco_type <- rownames(col_scores)

# ------------------------------------------------------------------------------
# 5. VISUALIZATION AESTHETICS (PLOTLY)
# ------------------------------------------------------------------------------
# Define color palette mapping
cluster_cols <- c("A"="#E76F51", "B"="#DDB7E5", "C"="#A67CA3", "D1"="#F1C40F", 
                  "D2"="#E68A00", "E1"="#A7DCA5", "E2"="#4CAF50", "F"="#0B6623")

# Generate dynamic axis labels based on extracted eigenvalues
xlab_txt <- paste0("CA 1 (", round(eig_df$percent[1], 1), "%)")
ylab_txt <- paste0("CA 2 (", round(eig_df$percent[2], 1), "%)")

# Define vector arrows pointing from origin to ecological types
eco_arrows <- lapply(1:nrow(col_scores), function(i) {
  list(
    x = col_scores$CA1[i], y = col_scores$CA2[i],
    xref = "x", yref = "y", axref = "x", ayref = "y",
    ax = 0, ay = 0, text = "", 
    showarrow = TRUE, arrowhead = 2, arrowsize = 1.2, arrowwidth = 1.5, arrowcolor = "grey30"
  )
})

# Define HTML-formatted italicized labels for ecological types
eco_labels <- lapply(1:nrow(col_scores), function(i) {
  list(
    x = col_scores$CA1[i], y = col_scores$CA2[i],
    xref = "x", yref = "y",
    text = paste0("<i>", col_scores$eco_type[i], "</i>"), 
    showarrow = FALSE, 
    # Dynamically anchor text to avoid overlapping the arrows
    xanchor = ifelse(col_scores$CA1[i] > 0, "left", "right"),
    yanchor = ifelse(col_scores$CA2[i] > 0, "bottom", "top"),
    font = list(size = 15, color = "black")
  )
})

# Combine annotations into a single list
all_annotations <- c(eco_arrows, eco_labels)

# ------------------------------------------------------------------------------
# 6. RENDER INTERACTIVE PLOTLY BIPLOT
# ------------------------------------------------------------------------------
export_filename <- paste0("Figure_CA_Plotly_", timestamp)

p <- plot_ly() %>%
  # Add cluster centroids as distinct markers
  add_markers(
    data = row_scores,
    x = ~CA1, 
    y = ~CA2,
    color = ~cluster, 
    colors = cluster_cols,
    marker = list(
      symbol = "star-triangle-down-open", 
      size = 20, 
      line = list(width = 2.5)
    ),
    text = ~paste("Cluster:", cluster), 
    hoverinfo = "text", 
    name = ~cluster
  ) %>%
  # Define layout parameters aligned with academic print standards
  layout(
    xaxis = list(
      title = list(text = xlab_txt, font = list(size = 22, color = "black")), 
      zeroline = TRUE, zerolinecolor = "grey", zerolinewidth = 1,
      showgrid = FALSE, 
      linecolor = "black", 
      linewidth = 1,
      mirror = TRUE, 
      ticks = "outside",       # Direct ticks strictly outwards
      ticklen = 8,             # Tick length
      tickwidth = 1.2,         # Tick thickness
      tickcolor = "black",
      tickfont = list(size = 20, color = "black", family = "Arial")
    ),
    yaxis = list(
      title = list(text = ylab_txt, font = list(size = 22, color = "black")), 
      zeroline = TRUE, zerolinecolor = "grey", zerolinewidth = 1,
      showgrid = FALSE, 
      linecolor = "black", 
      linewidth = 1,
      mirror = TRUE, 
      ticks = "outside",       
      ticklen = 8,             
      tickwidth = 1.2,         
      tickcolor = "black",
      tickfont = list(size = 20, color = "black", family = "Arial")
    ),
    annotations = all_annotations, 
    legend = list(
      title = list(text = "<b>Clusters</b>"), 
      x = 0.98, y = 0.02, 
      xanchor = "right", yanchor = "bottom",
      bgcolor = "rgba(0,0,0,0)", 
      borderwidth = 0,            
      font = list(size = 16, color = "black")
    ),
    plot_bgcolor = "white", 
    paper_bgcolor = "white"
  ) %>%
  # Configure interactive UI options (Camera icon export settings)
  config(
    toImageButtonOptions = list(
      format = "png",
      filename = export_filename,
      height = 700,
      width = 800,
      scale = 10 # High quality export scalar (~1000 DPI equivalent)
    ),
    displayModeBar = TRUE
  )

# Display the interactive plot in the RStudio Viewer pane
print(p)

# ------------------------------------------------------------------------------
# 7. EXPORT FILES
# ------------------------------------------------------------------------------
# Export as a self-contained, interactive HTML widget
html_name <- paste0(export_filename, "_Interactive.html")
saveWidget(p, html_name, selfcontained = TRUE)
cat(sprintf("Interactive HTML saved successfully as: %s\n", html_name))

# NOTE: For direct programmatic PNG export using orca/kaleido engine
# Remove the '#' from the lines below if kaleido is correctly installed:

# tryCatch({
#   png_name <- paste0(export_filename, "_1000DPI.png")
#   save_image(p, png_name, width = 800, height = 700, scale = 14)
#   cat(sprintf("High-resolution static PNG saved successfully as: %s\n", png_name))
# }, error = function(e) {
#   cat("\nNote: Automatic PNG export requires the 'kaleido' package to be configured.\n")
#   cat("You can also export the PNG manually by clicking the camera icon in the RStudio viewer.\n")
# })