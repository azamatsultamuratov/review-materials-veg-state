# ==============================================================================
# Script 11a: EXPORT ALL POSSIBLE INTERNAL INFORMATION FROM MODIFIED TWINSPAN
# Description: Dumps all internal structures, summaries, and matrices from the 
#              twinspanR object into a dedicated directory for reproducibility.
# ==============================================================================

library(twinspanR)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)

gc()

# ---------------------------------------------------------
# 1) DATA INPUT & MODEL RUN
# ---------------------------------------------------------
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
comm_raw_f4 <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)

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

common_plots <- intersect(rownames(mtws_f4), rownames(comm_raw_f4))
common_spp   <- intersect(colnames(mtws_f4), colnames(comm_raw_f4))
mtws_use <- mtws_f4[common_plots, common_spp, drop = FALSE]
raw_use  <- comm_raw_f4[common_plots, common_spp, drop = FALSE]

# ---------------------------------------------------------
# 2) EXPORT INTERNAL INFORMATION
# ---------------------------------------------------------
out_dir <- "mtws_internal_exports"
if (!dir.exists(out_dir)) dir.create(out_dir)

# 1) object names
sink(file.path(out_dir, "01_names_tw_lower.txt"))
cat("Names in tw_lower:\n"); print(names(tw_lower))
sink()

# 2) structure
sink(file.path(out_dir, "02_str_tw_lower.txt"))
cat("Structure of tw_lower:\n"); str(tw_lower, max.level = 4)
sink()

# 3) summary
sink(file.path(out_dir, "03_summary_tw_lower.txt"))
cat("tw_lower$summary:\n"); print(tw_lower$summary)
sink()

# 4) heter history
sink(file.path(out_dir, "04_heter_history.txt"))
cat("tw_lower$summary$heter:\n"); print(tw_lower$summary$heter)
sink()

# 5) classif
sink(file.path(out_dir, "05_classif.txt"))
cat("tw_lower$classif:\n"); print(tw_lower$classif)
sink()

# 6) twi
sink(file.path(out_dir, "06_twi.txt"))
cat("tw_lower$twi:\n"); print(tw_lower$twi)
sink()

# 7) spnames
sink(file.path(out_dir, "07_spnames.txt"))
cat("tw_lower$spnames:\n"); print(tw_lower$spnames)
sink()

# 8) dput versions for exact internal parsing
dput(tw_lower$summary, file = file.path(out_dir, "08_summary_dput.txt"))
dput(tw_lower$summary$heter, file = file.path(out_dir, "09_heter_dput.txt"))
dput(tw_lower$classif, file = file.path(out_dir, "10_classif_dput.txt"))
dput(tw_lower$twi, file = file.path(out_dir, "11_twi_dput.txt"))
dput(tw_lower$spnames, file = file.path(out_dir, "12_spnames_dput.txt"))

# 9) final memberships
lower_membership <- data.frame(
  plot_id = rownames(mtws_use),
  lower_class = as.character(cut(tw_lower)),
  stringsAsFactors = FALSE
)
write.csv(lower_membership, file.path(out_dir, "13_lower_membership.csv"), row.names = FALSE)

# 10) save matrix actually used
write.csv(mtws_use, file.path(out_dir, "14_mtws_use_matrix.csv"), row.names = TRUE)

# 11) save raw matrix if available
if (exists("raw_use")) {
  write.csv(raw_use, file.path(out_dir, "15_raw_use_matrix.csv"), row.names = TRUE)
}

# 12) save cut levels and settings
settings_df <- data.frame(object = c("clusters", "min.group.size"), value = c(8, 5))
write.csv(settings_df, file.path(out_dir, "16_basic_settings.csv"), row.names = FALSE)
writeLines("cut.levels = c(0, 0.3, 0.7, 2, 4)", con = file.path(out_dir, "17_cut_levels.txt"))

cat("All internal exports successfully written to folder:", out_dir, "\n")











# ==============================================================================
# Script 11b: MODIFIED TWINSPAN HIERARCHY + DIFFERENTIAL SPECIES EXPORT
# Description: Traverses the clustering tree and calculates differential species 
#              and pseudospecies between left and right branches at each division.
# ==============================================================================

gc()
library(twinspanR)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)

# ---------------------------------------------------------
# 0) INPUT & MODEL RUN
# ---------------------------------------------------------
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
cut_use <- c(0, 0.3, 0.7, 2, 4)

tw_lower <- twinspan(com = mtws_f4, modif = TRUE, clusters = 8, cut.levels = cut_use, 
                     min.group.size = 5, diss = "bray", quiet = TRUE, show.output.on.console = FALSE)

lower_membership <- data.frame(plot_id = rownames(mtws_f4), lower_class = as.character(cut(tw_lower)), stringsAsFactors = FALSE)
lower_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")

# ---------------------------------------------------------
# 1) BUILD TREE FROM HETER HISTORY
# ---------------------------------------------------------
heter_list <- tw_lower$summary$heter
n_final <- tw_lower$summary$clusters

make_leaf <- function(label) { list(type = "leaf", label = as.character(label), labels = as.character(label)) }
make_internal <- function(left, right, height) { list(type = "node", left = left, right = right, height = height, labels = c(left$labels, right$labels)) }

current <- lapply(seq_len(n_final), make_leaf)
for (m in n_final:2) {
  prev_step <- heter_list[[m - 1]]
  split_idx <- prev_step$which.most.heter
  split_height <- prev_step$cluster.heter[split_idx]
  left  <- current[[split_idx]]; right <- current[[split_idx + 1]]
  new_node <- make_internal(left, right, split_height)
  current <- append(current[-c(split_idx, split_idx + 1)], list(new_node), after = split_idx - 1)
}
tree_root <- current[[1]]

# ---------------------------------------------------------
# 2) ASSIGN NODE IDS
# ---------------------------------------------------------
assign_node_ids <- function(node, tree_id = 1) {
  node$tree_id <- tree_id
  if (node$type == "node") {
    node$left  <- assign_node_ids(node$left,  tree_id = tree_id * 2)
    node$right <- assign_node_ids(node$right, tree_id = tree_id * 2 + 1)
  }
  node
}
tree_root <- assign_node_ids(tree_root, 1)

# ---------------------------------------------------------
# 3) PLOT SETS FOR EACH NODE
# ---------------------------------------------------------
plots_by_class <- split(lower_membership$plot_id, lower_membership$lower_class)
get_node_plots <- function(labels_vec) { unlist(plots_by_class[as.character(labels_vec)], use.names = FALSE) }
fmt_lbl <- function(lbls) { paste(ifelse(lbls %in% names(lower_labels), lower_labels[lbls], lbls), collapse = "+") }

# ---------------------------------------------------------
# 4) SPECIES DIFFERENTIAL ANALYSIS
# ---------------------------------------------------------
calc_species_diff <- function(left_plots, right_plots, comm_mat, top_n = 10) {
  left_mat  <- comm_mat[left_plots, , drop = FALSE]; right_mat <- comm_mat[right_plots, , drop = FALSE]
  left_const  <- colMeans(left_mat  > 0, na.rm = TRUE) * 100; right_const <- colMeans(right_mat > 0, na.rm = TRUE) * 100
  diff <- left_const - right_const
  
  out <- data.frame(species = colnames(comm_mat), left_const = as.numeric(left_const), right_const = as.numeric(right_const),
                    diff = as.numeric(diff), abs_diff = abs(as.numeric(diff)), side = ifelse(diff > 0, "LEFT", "RIGHT"), stringsAsFactors = FALSE) %>%
    arrange(desc(abs_diff), desc(abs(left_const + right_const)))
  
  top_left <- out %>% filter(diff > 0) %>% arrange(desc(diff), desc(left_const)) %>% slice_head(n = top_n)
  top_right <- out %>% filter(diff < 0) %>% arrange(diff, desc(right_const)) %>% slice_head(n = top_n)
  
  list(full = out, left = top_left, right = top_right)
}

# ---------------------------------------------------------
# 5) PSEUDOSPECIES DIFFERENTIAL ANALYSIS
# ---------------------------------------------------------
make_pseudo_matrix <- function(comm_mat, cut_levels) {
  pos_cuts <- cut_levels[cut_levels > 0]; pseudo_list <- list()
  for (sp in colnames(comm_mat)) {
    x <- comm_mat[, sp]
    for (ct in pos_cuts) { nm <- paste0(sp, "@", ct); pseudo_list[[nm]] <- as.numeric(x >= ct) }
  }
  pseudo_mat <- as.data.frame(pseudo_list, check.names = FALSE); rownames(pseudo_mat) <- rownames(comm_mat); pseudo_mat
}
pseudo_mat <- make_pseudo_matrix(mtws_f4, cut_use)

calc_pseudo_diff <- function(left_plots, right_plots, pseudo_mat, top_n = 10) {
  left_mat  <- pseudo_mat[left_plots, , drop = FALSE]; right_mat <- pseudo_mat[right_plots, , drop = FALSE]
  left_const  <- colMeans(left_mat,  na.rm = TRUE) * 100; right_const <- colMeans(right_mat, na.rm = TRUE) * 100
  diff <- left_const - right_const
  
  out <- data.frame(pseudospecies = colnames(pseudo_mat), left_const = as.numeric(left_const), right_const = as.numeric(right_const),
                    diff = as.numeric(diff), abs_diff = abs(as.numeric(diff)), side = ifelse(diff > 0, "LEFT", "RIGHT"), stringsAsFactors = FALSE) %>%
    mutate(species = sub("@.*$", "", pseudospecies), cut_level = sub("^.*@", "", pseudospecies)) %>%
    arrange(desc(abs_diff), desc(abs(left_const + right_const)))
  
  top_left <- out %>% filter(diff > 0) %>% group_by(species) %>% arrange(desc(diff), desc(left_const), .by_group = TRUE) %>% slice(1) %>% ungroup() %>% arrange(desc(diff), desc(left_const)) %>% slice_head(n = top_n)
  top_right <- out %>% filter(diff < 0) %>% group_by(species) %>% arrange(diff, desc(right_const), .by_group = TRUE) %>% slice(1) %>% ungroup() %>% arrange(diff, desc(right_const)) %>% slice_head(n = top_n)
  
  list(full = out, left = top_left, right = top_right)
}

# ---------------------------------------------------------
# 6) FORMAT LINES & 7) TRAVERSE TREE AND EXPORT TEXT
# ---------------------------------------------------------
format_species_line <- function(df, plusminus = "+") { if(nrow(df) == 0) return("none"); paste0(plusminus, paste0(df$species, "(", round(abs(df$diff), 1), ")"), collapse = " ") }
format_pseudo_line <- function(df, plusminus = "+") { if(nrow(df) == 0) return("none"); paste0(plusminus, paste0(df$species, "@", df$cut_level, "(", round(abs(df$diff), 1), ")"), collapse = " ") }

node_species_tables <- list()
node_pseudo_tables  <- list()

write_node_text <- function(node, file = "", depth = 0) {
  indent <- paste(rep("  ", depth), collapse = "")
  if(node$type == "node") {
    left_labels <- node$left$labels; right_labels <- node$right$labels
    left_plots <- get_node_plots(left_labels); right_plots <- get_node_plots(right_labels)
    
    sp_res <- calc_species_diff(left_plots, right_plots, mtws_f4, top_n = 8)
    ps_res <- calc_pseudo_diff(left_plots, right_plots, pseudo_mat, top_n = 8)
    
    node_species_tables[[paste0("node_", node$tree_id)]] <<- sp_res$full %>% mutate(node_id = node$tree_id, left_clusters = fmt_lbl(left_labels), right_clusters = fmt_lbl(right_labels), heter = node$height)
    node_pseudo_tables[[paste0("node_", node$tree_id)]] <<- ps_res$full %>% mutate(node_id = node$tree_id, left_clusters = fmt_lbl(left_labels), right_clusters = fmt_lbl(right_labels), heter = node$height)
    
    cat(sprintf("%s%d) heter=%.3f: {%s} | {%s}", indent, node$tree_id, node$height, fmt_lbl(left_labels), fmt_lbl(right_labels)), "\n", file = file, append = TRUE, sep = "")
    cat(sprintf("%s  LEFT differential species:  %s", indent, format_species_line(sp_res$left, "+")), "\n", file = file, append = TRUE, sep = "")
    cat(sprintf("%s  RIGHT differential species: %s", indent, format_species_line(sp_res$right, "-")), "\n", file = file, append = TRUE, sep = "")
    cat(sprintf("%s  LEFT pseudospecies:        %s", indent, format_pseudo_line(ps_res$left, "+")), "\n", file = file, append = TRUE, sep = "")
    cat(sprintf("%s  RIGHT pseudospecies:       %s", indent, format_pseudo_line(ps_res$right, "-")), "\n\n", file = file, append = TRUE, sep = "")
    
    write_node_text(node$left, file = file, depth = depth + 1)
    write_node_text(node$right, file = file, depth = depth + 1)
  } else {
    cl <- node$label; cl_lab <- ifelse(cl %in% names(lower_labels), lower_labels[[cl]], cl); plots <- get_node_plots(cl)
    cat(sprintf("%s%d) lower=%s (%s), N=%d: %s", indent, node$tree_id, cl, cl_lab, length(plots), paste(plots, collapse = " ")), "\n", file = file, append = TRUE, sep = "")
  }
}

# ---------------------------------------------------------
# 8) MAIN TXT EXPORT
# ---------------------------------------------------------
txt_file <- "modified_twinspan_hierarchy_with_differential_species_variant1.txt"
if(file.exists(txt_file)) file.remove(txt_file)

writeLines(c("Modified TWINSPAN hierarchy with differential species export (Variant 1)", paste0("Lower clusters = ", n_final), paste0("Cut levels = ", paste(cut_use, collapse = ", ")), "Differential species = constancy difference between left and right child nodes", "Pseudospecies = strongest single cut-level retained per species per side", ""), con = txt_file)

write_node_text(tree_root, file = txt_file, depth = 0)

# ---------------------------------------------------------
# 9 & 10) EXPORT NODE TABLES
# ---------------------------------------------------------
species_diff_table <- bind_rows(node_species_tables)
pseudo_diff_table  <- bind_rows(node_pseudo_tables)

write.csv(species_diff_table, "node_differential_species_full_variant1.csv", row.names = FALSE)
write.csv(pseudo_diff_table, "node_differential_pseudospecies_full_variant1.csv", row.names = FALSE)

species_top <- species_diff_table %>% group_by(node_id) %>% arrange(desc(abs_diff), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()
pseudo_top <- pseudo_diff_table %>% group_by(node_id, side, species) %>% arrange(desc(abs_diff), .by_group = TRUE) %>% slice(1) %>% ungroup() %>% group_by(node_id) %>% arrange(desc(abs_diff), .by_group = TRUE) %>% slice_head(n = 20) %>% ungroup()

write.csv(species_top, "node_differential_species_top20_variant1.csv", row.names = FALSE)
write.csv(pseudo_top, "node_differential_pseudospecies_top20_variant1.csv", row.names = FALSE)

cat("Success! Files exported for Differential Species Variant.\n")











# ==============================================================================
# Script 11c: ORIGINAL-LIKE TXT HIERARCHY FROM MODIFIED TWINSPAN
# Target Journal: Journal of Vegetation Science (JVS)
# Description: Parses the internal tw_lower$twi object to reconstruct the classic 
#              DOS-era TWINSPAN textual report (Eigenvalues, +/- Preferentials).
# ==============================================================================

gc()
library(twinspanR)
library(dplyr)
library(tidyr)
library(tibble)

# ---------------------------------------------------------
# 0) CHECK INPUTS AND PREPARE MODEL
# ---------------------------------------------------------
mtws_f4 <- read.csv("mtwinspan_input_freq4_sqrt.csv", row.names = 1, check.names = FALSE)
cut_use <- c(0, 0.3, 0.7, 2, 4)
tw_lower <- twinspan(com = mtws_f4, modif = TRUE, clusters = 8, cut.levels = cut_use, min.group.size = 5, diss = "bray", quiet = TRUE, show.output.on.console = FALSE)

lower_membership <- data.frame(plot_id = rownames(mtws_f4), lower_class = as.character(cut(tw_lower)), stringsAsFactors = FALSE)
lower_membership$lower_class <- as.character(lower_membership$lower_class)
lower_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")

if(!exists("tw_lower")) stop("tw_lower object not found in memory.")
if(is.null(tw_lower$twi)) stop("tw_lower$twi is missing.")

# ---------------------------------------------------------
# 1) FLATTEN twi TEXT
# ---------------------------------------------------------
twi_raw <- as.character(unlist(tw_lower$twi, use.names = FALSE))

clean_line <- function(x) { x <- gsub("[[:cntrl:]]", " ", x); x <- gsub("\\s+", " ", x); trimws(x) }
twi_lines <- vapply(twi_raw, clean_line, character(1))
twi_lines <- twi_lines[nchar(twi_lines) > 0]

writeLines(twi_lines, "mtws_full_cleaned_twi_report.txt")

# ---------------------------------------------------------
# 2) SPLIT INTO DIVISION BLOCKS & PARSE
# ---------------------------------------------------------
div_idx <- grep("^DIVISION[[:space:]]+[0-9]+", twi_lines)
if(length(div_idx) == 0) stop("No DIVISION blocks found inside tw_lower$twi.")

blocks <- list()
for(i in seq_along(div_idx)) {
  s <- div_idx[i]
  e <- if(i < length(div_idx)) div_idx[i + 1] - 1 else length(twi_lines)
  blocks[[i]] <- twi_lines[s:e]
}

extract_section <- function(block, start_pattern, stop_patterns) {
  s <- grep(start_pattern, block, ignore.case = TRUE)
  if(length(s) == 0) return(character(0))
  s <- s[1]
  stop_idx <- integer(0)
  for(p in stop_patterns) { z <- grep(p, block, ignore.case = TRUE); z <- z[z > s]; if(length(z) > 0) stop_idx <- c(stop_idx, z[1]) }
  e <- if(length(stop_idx) == 0) length(block) else min(stop_idx) - 1
  out <- block[(s + 1):e]
  out[nchar(out) > 0]
}

collapse_pref <- function(x) { if(length(x) == 0) return(""); paste(trimws(gsub("^[-+ ]+", "", x)), collapse = " ") }

parse_division_block <- function(block) {
  div_line <- block[grep("^DIVISION[[:space:]]+[0-9]+", block)][1]
  div_num <- sub("^DIVISION[[:space:]]+([0-9]+).*", "\\1", div_line)
  
  eig_line <- block[grep("Eigenvalue", block, ignore.case = TRUE)][1]
  eig_val <- if(length(eig_line) > 0 && !is.na(eig_line)) sub(".*Eigenvalue[[:space:]]*=?[[:space:]]*([0-9.]+).*", "\\1", eig_line, ignore.case = TRUE) else NA_character_
  
  thr_line <- block[grep("<[[:space:]]*[-0-9.]+", block)]
  thr_val <- if(length(thr_line) > 0) sub(".*(<[[:space:]]*[-0-9.]+).*", "\\1", thr_line[1]) else ""
  
  neg_pref <- extract_section(block, "^NEGATIVE PREFERENTIALS", c("^POSITIVE PREFERENTIALS", "^NON-PREFERENTIALS", "^DIVISION", "^End of level", "^This is the end"))
  pos_pref <- extract_section(block, "^POSITIVE PREFERENTIALS", c("^NEGATIVE PREFERENTIALS", "^NON-PREFERENTIALS", "^DIVISION", "^End of level", "^This is the end"))
  
  neg_txt <- collapse_pref(neg_pref); pos_txt <- collapse_pref(pos_pref)
  
  summary_txt <- paste0(div_num, ") eig=", eig_val, ": ", if(nchar(pos_txt) > 0) paste0("+", pos_txt, " ") else "", if(nchar(neg_txt) > 0) paste0("-", neg_txt, " ") else "", thr_val)
  summary_txt <- trimws(gsub("\\s+", " ", summary_txt))
  
  list(division = div_num, eig = eig_val, threshold = thr_val, negative_preferentials = neg_txt, positive_preferentials = pos_txt, summary_line = summary_txt, raw_block = block)
}

parsed_blocks <- lapply(blocks, parse_division_block)

# ---------------------------------------------------------
# 4) BUILD TREE STRUCTURE AND MATCH TO DIVISIONS
# ---------------------------------------------------------
heter_list <- tw_lower$summary$heter
n_final <- tw_lower$summary$clusters

make_leaf <- function(label) { list(type = "leaf", label = as.character(label), labels = as.character(label)) }
make_internal <- function(left, right, height) { list(type = "node", left = left, right = right, height = height, labels = c(left$labels, right$labels)) }

current <- lapply(seq_len(n_final), make_leaf)
for (m in n_final:2) {
  prev_step <- heter_list[[m - 1]]; split_idx <- prev_step$which.most.heter
  split_height <- prev_step$cluster.heter[split_idx]
  left <- current[[split_idx]]; right <- current[[split_idx + 1]]
  new_node <- make_internal(left, right, split_height)
  current <- append(current[-c(split_idx, split_idx + 1)], list(new_node), after = split_idx - 1)
}
tree_root <- current[[1]]

assign_node_ids <- function(node, tree_id = 1) {
  node$tree_id <- tree_id
  if(node$type == "node") { node$left <- assign_node_ids(node$left, tree_id * 2); node$right <- assign_node_ids(node$right, tree_id * 2 + 1) }
  node
}
tree_root <- assign_node_ids(tree_root, 1)

plots_by_class <- split(lower_membership$plot_id, lower_membership$lower_class)
fmt_lbl <- function(lbls) { paste(ifelse(lbls %in% names(lower_labels), lower_labels[lbls], lbls), collapse = "+") }

block_map <- parsed_blocks
names(block_map) <- vapply(parsed_blocks, function(x) x$division, character(1))

# ---------------------------------------------------------
# 5) WRITE ORIGINAL-LIKE TXT HIERARCHY
# ---------------------------------------------------------
txt_file <- "modified_twinspan_original_like_hierarchy.txt"
if(file.exists(txt_file)) file.remove(txt_file)

writeLines(c("Modified TWINSPAN original-like hierarchy summary", "Reconstructed from tw_lower$twi + final lower memberships", ""), con = txt_file)

write_tree_with_summary <- function(node, depth = 0, file = "") {
  indent <- paste(rep("  ", depth), collapse = "")
  if(node$type == "node") {
    node_id_chr <- as.character(node$tree_id)
    line <- if(node_id_chr %in% names(block_map)) block_map[[node_id_chr]]$summary_line else paste0(node$tree_id, ") heter=", round(node$height, 3), ": {", fmt_lbl(node$left$labels), "} | {", fmt_lbl(node$right$labels), "}")
    cat(indent, line, "\n", file = file, append = TRUE, sep = "")
    write_tree_with_summary(node$left, depth = depth + 1, file = file)
    write_tree_with_summary(node$right, depth = depth + 1, file = file)
  } else {
    cl <- node$label; cl_lab <- ifelse(cl %in% names(lower_labels), lower_labels[[cl]], cl); plots <- plots_by_class[[cl]]
    cat(indent, node$tree_id, ") N=", length(plots), ": ", paste(plots, collapse = " "), " [lower=", cl, ", ", cl_lab, "]\n", file = file, append = TRUE, sep = "")
  }
}
write_tree_with_summary(tree_root, depth = 0, file = txt_file)

# ---------------------------------------------------------
# 6 & 7) SAVE CSV SUMMARY AND FULL BLOCKS
# ---------------------------------------------------------
division_table <- bind_rows(lapply(parsed_blocks, function(x) { data.frame(division = x$division, eig = x$eig, threshold = x$threshold, positive_preferentials = x$positive_preferentials, negative_preferentials = x$negative_preferentials, summary_line = x$summary_line, stringsAsFactors = FALSE) }))
write.csv(division_table, "modified_twinspan_division_summary_table.csv", row.names = FALSE)

con <- file("modified_twinspan_division_blocks_full.txt", open = "wt")
for(i in seq_along(parsed_blocks)) { writeLines(paste0("===== DIVISION ", parsed_blocks[[i]]$division, " ====="), con); writeLines(parsed_blocks[[i]]$raw_block, con); writeLines("", con) }
close(con)

cat("Success! Original-like TWINSPAN text reports generated.\n")