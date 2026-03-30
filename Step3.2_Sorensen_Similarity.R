# ==============================================================================
# Script Name: Within-Class Sørensen Similarity Analysis
# Description: Calculates the mean within-class Sørensen similarity (1 - Bray-Curtis 
#              on binary data) to assess the compositional homogeneity of syntaxa.
#              Analysis is performed across hierarchical levels (Upper, Lower, and 
#              Custom Groupings) and stratified by life-form.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(vegan)
library(tibble)

# ------------------------------------------------------------------------------
# 1. DATA INPUT & PREPARATION
# ------------------------------------------------------------------------------
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
sp_code  <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)
raw_use  <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)

class_df$upper_class <- as.character(class_df$upper_class)
class_df$lower_class <- as.character(class_df$lower_class)

# Ensure strict alignment between community matrix and classification
common_plots <- intersect(rownames(raw_use), class_df$plot_id)
raw_use2 <- raw_use[common_plots, , drop = FALSE]
class_df2 <- class_df %>% dplyr::filter(plot_id %in% common_plots)

# Align and filter species metadata
species_in_matrix <- colnames(raw_use2)
sp_match <- sp_code %>% dplyr::filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

# Extract species columns by life-form category
annual_cols    <- sp_match$s_code[sp_match$life_form == "An"]
perennial_cols <- sp_match$s_code[sp_match$life_form == "Pr"]
woody_cols     <- sp_match$s_code[sp_match$life_form %in% c("Sr", "Ss")]

# ------------------------------------------------------------------------------
# 2. RECODE CLASS LABELS TO STANDARDIZED SYNTAXONOMIC NOMENCLATURE
# ------------------------------------------------------------------------------
# Upper hierarchical level mapping (1-6 -> A-F)
upper_map <- c("1" = "A", "2" = "B", "3" = "C", "4" = "D", "5" = "E", "6" = "F")

# Lower hierarchical level mapping (1-8 -> D and E are split)
lower_map <- c("1" = "A", "2" = "B", "3" = "C", "4" = "D1", "5" = "D2", "6" = "E1", "7" = "E2", "8" = "F")

class_df2 <- class_df2 %>%
  dplyr::mutate(
    upper_label = dplyr::recode(upper_class, !!!upper_map),
    lower_label = dplyr::recode(lower_class, !!!lower_map)
  )

# ------------------------------------------------------------------------------
# 3. CORE COMPUTATION FUNCTION: WITHIN-CLASS SØRENSEN SIMILARITY
# ------------------------------------------------------------------------------
calc_within_class_sorensen <- function(comm, plot_ids, species_subset = NULL) {
  
  # Subset matrix to target plots
  sub <- comm[rownames(comm) %in% plot_ids, , drop = FALSE]
  
  # Subset matrix to target life-form if specified
  if (!is.null(species_subset)) {
    species_subset <- intersect(colnames(sub), species_subset)
    sub <- sub[, species_subset, drop = FALSE]
  }
  
  # Remove species with zero occurrences in the current subset
  if (ncol(sub) > 0) sub <- sub[, colSums(sub > 0) > 0, drop = FALSE]
  
  # Remove empty plots after filtering
  if (ncol(sub) > 0) sub <- sub[rowSums(sub > 0) > 0, , drop = FALSE]
  
  # Exception Handling: No species left
  if (ncol(sub) == 0) {
    return(data.frame(n_releves = nrow(sub), n_species = 0, mean_J = NA, median_J = NA, min_J = NA, max_J = NA, status = "n.a."))
  }
  
  # Exception Handling: Fewer than 2 plots (cannot calculate pairwise similarity)
  if (nrow(sub) < 2) {
    return(data.frame(n_releves = nrow(sub), n_species = ncol(sub), mean_J = NA, median_J = NA, min_J = NA, max_J = NA, status = "n.d."))
  }
  
  # Convert to Presence/Absence (Binary matrix)
  pa <- (sub > 0) * 1
  
  # Calculate Sørensen similarity (1 - binary Bray-Curtis)
  d <- vegan::vegdist(pa, method = "bray", binary = TRUE)
  sim <- 1 - as.matrix(d)
  
  # Extract upper triangle of the similarity matrix to avoid duplicate pairwise values
  vals <- sim[upper.tri(sim)]
  
  data.frame(
    n_releves = nrow(pa), n_species = ncol(pa),
    mean_J = mean(vals, na.rm = TRUE),
    median_J = median(vals, na.rm = TRUE),
    min_J = min(vals, na.rm = TRUE),
    max_J = max(vals, na.rm = TRUE),
    status = "ok"
  )
}

# ------------------------------------------------------------------------------
# 4. WRAPPER FUNCTION: AGGREGATE ACROSS ALL CLASSES
# ------------------------------------------------------------------------------
summarise_classes <- function(comm, class_df, class_var, species_subset = NULL, dataset_label = "All") {
  cls <- unique(class_df[[class_var]])
  
  out <- lapply(cls, function(cl) {
    plot_ids <- class_df$plot_id[class_df[[class_var]] == cl]
    res <- calc_within_class_sorensen(comm = comm, plot_ids = plot_ids, species_subset = species_subset)
    res[[class_var]] <- cl
    res$dataset <- dataset_label
    res
  })
  
  dplyr::bind_rows(out) %>%
    dplyr::select(dplyr::all_of(class_var), dataset, n_releves, n_species, mean_J, median_J, min_J, max_J, status)
}

# ------------------------------------------------------------------------------
# 5. DEFINE EXTRA GROUPED SCHEMES (Alternative Hierarchies)
# ------------------------------------------------------------------------------
class_df_extra <- class_df2 %>%
  dplyr::mutate(
    lower_group = dplyr::case_when(
      lower_label %in% c("A", "B") ~ "A,B",
      lower_label %in% c("A", "B", "C") ~ "A,B,C",
      lower_label %in% c("D1", "D2") ~ "D",
      lower_label %in% c("A", "B", "C", "D1", "D2") ~ "A,B,C,D",
      lower_label %in% c("E1", "E2") ~ "E",
      lower_label %in% c("A", "B", "C", "D1", "D2", "E1", "E2") ~ "A,B,C,D,E",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# 6. EXECUTE COMPUTATIONS: LOWER CLASSES
# ------------------------------------------------------------------------------
lower_all       <- summarise_classes(raw_use2, class_df2, "lower_label", NULL,           "All")
lower_annual    <- summarise_classes(raw_use2, class_df2, "lower_label", annual_cols,    "Annual")
lower_perennial <- summarise_classes(raw_use2, class_df2, "lower_label", perennial_cols, "Perennial")
lower_woody     <- summarise_classes(raw_use2, class_df2, "lower_label", woody_cols,     "Woody")

lower_all_out <- dplyr::bind_rows(lower_all, lower_annual, lower_perennial, lower_woody)
write.csv(lower_all_out, "within_lower_sorensen_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. EXECUTE COMPUTATIONS: UPPER CLASSES
# ------------------------------------------------------------------------------
upper_all       <- summarise_classes(raw_use2, class_df2, "upper_label", NULL,           "All")
upper_annual    <- summarise_classes(raw_use2, class_df2, "upper_label", annual_cols,    "Annual")
upper_perennial <- summarise_classes(raw_use2, class_df2, "upper_label", perennial_cols, "Perennial")
upper_woody     <- summarise_classes(raw_use2, class_df2, "upper_label", woody_cols,     "Woody")

upper_all_out <- dplyr::bind_rows(upper_all, upper_annual, upper_perennial, upper_woody)
write.csv(upper_all_out, "within_upper_sorensen_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. EXECUTE COMPUTATIONS: EXTRA GROUPINGS
# ------------------------------------------------------------------------------
valid_extra_df <- class_df_extra %>% dplyr::filter(!is.na(lower_group))

extra_all       <- summarise_classes(raw_use2, valid_extra_df, "lower_group", NULL,           "All")
extra_annual    <- summarise_classes(raw_use2, valid_extra_df, "lower_group", annual_cols,    "Annual")
extra_perennial <- summarise_classes(raw_use2, valid_extra_df, "lower_group", perennial_cols, "Perennial")
extra_woody     <- summarise_classes(raw_use2, valid_extra_df, "lower_group", woody_cols,     "Woody")

extra_all_out <- dplyr::bind_rows(extra_all, extra_annual, extra_perennial, extra_woody)
write.csv(extra_all_out, "within_extra_group_sorensen_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. FORMAT AND EXPORT PUBLICATION-READY WIDE TABLES
# ------------------------------------------------------------------------------
format_wide <- function(df, group_col) {
  df %>%
    dplyr::select(!!sym(group_col), dataset, mean_J, status) %>%
    dplyr::mutate(value = ifelse(status == "ok", sprintf("%.2f", mean_J), status)) %>%
    dplyr::select(-mean_J, -status) %>%
    tidyr::pivot_wider(names_from = dataset, values_from = value)
}

lower_wide <- format_wide(lower_all_out, "lower_label")
upper_wide <- format_wide(upper_all_out, "upper_label")
extra_wide <- format_wide(extra_all_out, "lower_group")

write.csv(lower_wide, "within_lower_sorensen_wide.csv", row.names = FALSE)
write.csv(upper_wide, "within_upper_sorensen_wide.csv", row.names = FALSE)
write.csv(extra_wide, "within_extra_group_sorensen_wide.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 10. GENERATE PLOT-FRIENDLY LONG TABLE
# ------------------------------------------------------------------------------
plot_lower <- lower_all_out %>% dplyr::mutate(level = "Lower") %>% dplyr::rename(class = lower_label)
plot_upper <- upper_all_out %>% dplyr::mutate(level = "Upper") %>% dplyr::rename(class = upper_label)
plot_extra <- extra_all_out %>% dplyr::mutate(level = "Extra_group") %>% dplyr::rename(class = lower_group)

plot_df <- dplyr::bind_rows(plot_lower, plot_upper, plot_extra)
write.csv(plot_df, "within_class_sorensen_plot_long.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 11. FINAL DIAGNOSTICS LOG
# ------------------------------------------------------------------------------
cat("\n--- Data Alignment Diagnostics ---\n")
cat("Original raw plots matched: ", length(common_plots), "\n")
cat("Dimensions of filtered community matrix: ", dim(raw_use2)[1], "plots x", dim(raw_use2)[2], "species\n")
cat("Success: All similarity computations completed and exported.\n")