# ==============================================================================
# Script Name: Tanglegram - TWINSPAN vs. Phylogenetic Tree
# Description: This script isolates the necessary steps to construct a modified 
#              TWINSPAN hierarchy (hc_mtws) and a phylogenetic tree (phydist, 
#              comm_matched), then face-to-face compares them using a tanglegram.
# ==============================================================================

graphics.off()
gc()

# -------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# -------------------------------------------------------------------
library(V.PhyloMaker2)
library(ape)
library(picante)
library(dendextend)
library(dplyr)
library(twinspanR)

cat("\nPreparing data for Tanglegram. Please wait...\n")

# -------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# -------------------------------------------------------------------
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
sp_code  <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)

# Replace short species codes (s_code) with full names (sp_name) in comm_raw
col_indices <- match(colnames(comm_raw), sp_code$s_code)
colnames(comm_raw) <- sp_code$sp_name[col_indices]

# Group plots by their respective lower_class clusters
comm_raw$plot_id <- rownames(comm_raw)
merged_data <- merge(class_df[, c("plot_id", "lower_class")], comm_raw, by = "plot_id")

# Aggregate species cover (sum) by cluster
comm_clustered <- aggregate(. ~ lower_class, data = merged_data[, -1], FUN = sum)
rownames(comm_clustered) <- paste("Cluster", comm_clustered$lower_class)
comm_clustered$lower_class <- NULL

# -------------------------------------------------------------------
# 2. PHYLOGENETIC TREE AND DISTANCES (phydist & comm_matched)
# -------------------------------------------------------------------
# Prepare species list for V.PhyloMaker2
sp_list <- data.frame(
  species = gsub("_", " ", sp_code$sp_name), 
  genus = sapply(strsplit(sp_code$sp_name, "_"), `[`, 1),
  family = NA, 
  stringsAsFactors = FALSE
)

# Generate the mega-tree
tree_res <- phylo.maker(sp.list = sp_list, scenarios = "S3")
full_tree <- tree_res$scenario.3

# Match the phylogeny with the clustered community matrix
matched <- match.phylo.comm(full_tree, comm_clustered)
phy_matched <- matched$phy
comm_matched <- matched$comm

# Calculate cophenetic phylogenetic distances
phydist <- cophenetic(phy_matched)

# -------------------------------------------------------------------
# 3. RECONSTRUCT TWINSPAN HIERARCHY (hc_mtws)
# -------------------------------------------------------------------
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
cut_use <- c(0, 0.3, 0.7, 2, 4)

# Run Modified TWINSPAN
tw_lower <- twinspan(com = mtws_f4, modif = TRUE, clusters = 8, cut.levels = cut_use, 
                     min.group.size = 5, diss = "bray", quiet = TRUE, show.output.on.console = FALSE)

# Reconstruct hclust object from TWINSPAN splitting history
heter_list <- tw_lower$summary$heter
n_final <- tw_lower$summary$clusters
current <- as.list(seq_len(n_final))
merge_mat <- matrix(NA_integer_, nrow = n_final - 1, ncol = 2)
height_vec <- numeric(n_final - 1)
next_node_id <- 1

to_hclust_id <- function(x) {
  node_id <- attr(x, "node_id")
  if (is.null(node_id)) return(-as.integer(x)) else return(as.integer(node_id))
}

for (m in n_final:2) {
  prev_step <- heter_list[[m - 1]]; split_idx <- prev_step$which.most.heter
  split_height <- prev_step$cluster.heter[split_idx]
  left  <- current[[split_idx]]; right <- current[[split_idx + 1]]
  left_id  <- to_hclust_id(left); right_id <- to_hclust_id(right)
  row_id <- next_node_id; merge_mat[row_id, ] <- c(left_id, right_id)
  height_vec[row_id] <- split_height
  new_node <- c(unlist(left), unlist(right)); attr(new_node, "node_id") <- row_id
  current <- append(current[-c(split_idx, split_idx + 1)], list(new_node), after = split_idx - 1)
  next_node_id <- next_node_id + 1
}

# Finalize hclust object
leaf_labels <- c("A", "B", "C", "D1", "D2", "E1", "E2", "F")
hc_mtws <- list(merge = merge_mat, height = height_vec, order = seq_len(n_final), 
                labels = leaf_labels, method = "reconstructed modified TWINSPAN", 
                dist.method = "cluster heterogeneity")
class(hc_mtws) <- "hclust"

# -------------------------------------------------------------------
# 4. PREPARE DENDROGRAMS FOR TANGLEGRAM
# -------------------------------------------------------------------
cat("\nCombining both trees into a Tanglegram...\n")

# A) TWINSPAN DENDROGRAM
dend_tw <- as.dendrogram(hc_mtws)
dend_tw <- rank_branches(dend_tw) # Force uniform branch spacing (Cladogram style)

# B) PHYLOGENETIC DENDROGRAM
dist_phylo_clusters <- comdist(comm_matched, phydist, abundance.weighted = TRUE)
hc_phylo <- hclust(dist_phylo_clusters, method = "average")
dend_phy <- as.dendrogram(hc_phylo)

# Standardize phylogenetic labels to match TWINSPAN (A, B, C...)
phy_labels <- labels(dend_phy)
phy_labels <- gsub("Cluster 1", "A", phy_labels)
phy_labels <- gsub("Cluster 2", "B", phy_labels)
phy_labels <- gsub("Cluster 3", "C", phy_labels)
phy_labels <- gsub("Cluster 4", "D1", phy_labels)
phy_labels <- gsub("Cluster 5", "D2", phy_labels)
phy_labels <- gsub("Cluster 6", "E1", phy_labels)
phy_labels <- gsub("Cluster 7", "E2", phy_labels)
phy_labels <- gsub("Cluster 8", "F", phy_labels)
labels(dend_phy) <- phy_labels

# Force uniform branch thickness
dend_tw <- set(dend_tw, "branches_lwd", 5.0)
dend_phy <- set(dend_phy, "branches_lwd", 5.0)

# Combine and untangle trees
dl <- dendlist(TWINSPAN = dend_tw, Phylogeny = dend_phy)
dl <- untangle(dl, method = "step1side")

# Color mapping for connecting lines
color_map <- c("A"="#E76F51", "B"="#DDB7E5", "C"="#A67CA3", "D1"="#F1C40F", 
               "D2"="#E68A00", "E1"="#A7DCA5", "E2"="#4CAF50", "F"="#0B6623")
left_labels <- labels(dl[[1]])
line_colors <- color_map[left_labels]

# -------------------------------------------------------------------
# 5. EXPORT HIGH-RESOLUTION TANGLEGRAM
# -------------------------------------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
file_name <- paste0("Figure_Tanglegram_TWINSPAN_vs_Phylogeny_", timestamp, ".png")

png(file_name, width = 22, height = 11, units = "in", res = 1000, bg = "white")

# Plot parameters: thick lines, large axes
par(font = 2, 
    font.main = 2, 
    font.lab = 2, 
    font.axis = 1.5, 
    cex.axis = 3.5,       # Axis number size
    lwd = 7.0,            # Main axis line width
    lwd.ticks = 5.0,      # Axis tick width
    lab = c(10, 10, 7), 
    mgp = c(3, 1.8, 0))   # Offset numbers slightly downward

tanglegram(dl,
           sort = FALSE,
           margin_outer = 6.0,                
           margin_top = 6.0,                  
           margin_bottom = 3.5,               # Shifted distance axes upward
           margin_inner = 6.0,                
           color_lines = line_colors,
           lwd = 5.0,                         # Connecting line width
           edge.lwd = 5.0,                    
           columns_width = c(5, 1.5, 5),      
           highlight_distinct_edges = FALSE,  
           highlight_branches_lwd = FALSE,    
           main_left = "A",                   
           main_right = "B",                  
           lab.cex = 3.5,                     
           cex_main = 4.5,                    
           axes = TRUE,                       
           dLeaf_left = 0,                    
           dLeaf_right = 0)                   

dev.off()

cat(sprintf("\nSuccess! Tanglegram saved as: %s\n", file_name))