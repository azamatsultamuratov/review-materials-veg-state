# =========================================================
# 10a: CONSTANCY AND FIDELITY SYNOPTIC TABLE
# =========================================================

gc()
library(dplyr)
library(tidyr)
library(tibble)

# =========================================================
# SETTINGS
# =========================================================
top_n_diag      <- 12
min_const_diag  <- 10
min_all_comp    <- 10

# =========================================================
# 1) INPUT
# =========================================================
class_df <- read.csv("mtws_membership_upper6_lower8.csv", stringsAsFactors = FALSE)
comm_raw_f4 <- read.csv("mtwinspan_input_freq4_raw.csv", row.names = 1, check.names = FALSE)
sp_code <- read.csv("df_species_code.csv", stringsAsFactors = FALSE)

# =========================================================
# 2) PREPARE MATRIX
# =========================================================
common_plots <- intersect(rownames(comm_raw_f4), class_df$plot_id)
comm_raw <- comm_raw_f4[common_plots, , drop = FALSE]
class_df <- class_df %>% filter(plot_id %in% common_plots)
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]
class_df$lower_class <- as.character(class_df$lower_class)

comm_raw <- comm_raw[rowSums(comm_raw) > 0, , drop = FALSE]
comm_raw <- comm_raw[, colSums(comm_raw) > 0, drop = FALSE]
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]

comm_pa <- (comm_raw > 0) * 1

# =========================================================
# 3) LABELS / SPECIES METADATA
# =========================================================
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
species_in_matrix <- colnames(comm_pa)
sp_match <- sp_code %>% filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

if(!"sp_name" %in% names(sp_match)) sp_match$sp_name <- sp_match$s_code
if(!"life_form" %in% names(sp_match)) stop("Column 'life_form' not found in df_species_code.csv")
lf_df <- sp_match %>% select(s_code, sp_name, life_form)

# =========================================================
# 4) CONSTANCY
# =========================================================
constancy_long <- as.data.frame(comm_pa) %>%
  rownames_to_column("plot_id") %>%
  left_join(class_df[, c("plot_id", "lower_class")], by = "plot_id") %>%
  pivot_longer(cols = -c(plot_id, lower_class), names_to = "s_code", values_to = "pa") %>%
  group_by(lower_class, s_code) %>%
  summarise(constancy = mean(pa, na.rm = TRUE) * 100, .groups = "drop")

overall_constancy <- as.data.frame(comm_pa) %>%
  rownames_to_column("plot_id") %>%
  pivot_longer(cols = -plot_id, names_to = "s_code", values_to = "pa") %>%
  group_by(s_code) %>%
  summarise(All = mean(pa, na.rm = TRUE) * 100, n_occ = sum(pa > 0, na.rm = TRUE), .groups = "drop")

constancy_wide <- constancy_long %>%
  mutate(class_label = class_labels[lower_class]) %>%
  select(-lower_class) %>%
  pivot_wider(names_from = class_label, values_from = constancy, values_fill = 0)

# =========================================================
# 5) TPC / COVER SUMMARY
# =========================================================
cover_long <- as.data.frame(comm_raw) %>%
  rownames_to_column("plot_id") %>%
  left_join(class_df[, c("plot_id", "lower_class")], by = "plot_id") %>%
  pivot_longer(cols = -c(plot_id, lower_class), names_to = "s_code", values_to = "cover")

cover_summary <- cover_long %>%
  group_by(s_code) %>%
  summarise(TPC = mean(cover, na.rm = TRUE),
            TPC_occ = ifelse(sum(cover > 0, na.rm = TRUE) > 0, mean(cover[cover > 0], na.rm = TRUE), NA_real_),
            .groups = "drop")

# =========================================================
# 6) SAFE PHI
# =========================================================
phi_one <- function(x, g) {
  x <- as.numeric(x); g <- as.numeric(g)
  a <- as.numeric(sum(x == 1 & g == 1, na.rm = TRUE))
  b <- as.numeric(sum(x == 1 & g == 0, na.rm = TRUE))
  c <- as.numeric(sum(x == 0 & g == 1, na.rm = TRUE))
  d <- as.numeric(sum(x == 0 & g == 0, na.rm = TRUE))
  denom <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  if(is.na(denom) || !is.finite(denom) || denom == 0) return(NA_real_)
  (a * d - b * c) / denom
}

class_levels <- sort(unique(class_df$lower_class))
phi_list <- vector("list", length(class_levels))
names(phi_list) <- class_levels

for(cl in class_levels) {
  g <- ifelse(class_df$lower_class == cl, 1, 0)
  phi_vals <- sapply(colnames(comm_pa), function(sp) { phi_one(comm_pa[, sp], g) })
  phi_list[[cl]] <- data.frame(s_code = names(phi_vals), lower_class = cl, class_label = class_labels[cl], phi = as.numeric(phi_vals), stringsAsFactors = FALSE)
}
phi_table <- bind_rows(phi_list)

# =========================================================
# 7) JOIN
# =========================================================
const_rank <- constancy_long %>%
  left_join(phi_table, by = c("s_code", "lower_class")) %>%
  mutate(class_label = class_labels[lower_class])

syn0 <- constancy_wide %>%
  left_join(overall_constancy, by = "s_code") %>%
  left_join(cover_summary, by = "s_code") %>%
  left_join(lf_df, by = "s_code")

# =========================================================
# 8) ASSIGN EACH SPECIES TO ONE BEST CLASS
# =========================================================
best_class <- const_rank %>%
  filter(!is.na(phi), phi > 0, constancy >= min_const_diag) %>%
  group_by(s_code) %>%
  arrange(desc(phi), desc(constancy), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(s_code, best_block = class_label, best_phi = phi, best_constancy = constancy)

diag_species <- best_class %>%
  group_by(best_block) %>%
  arrange(desc(best_phi), desc(best_constancy), .by_group = TRUE) %>%
  slice_head(n = top_n_diag) %>%
  ungroup()

diag_counts <- diag_species %>%
  count(best_block, name = "n_diag")

# =========================================================
# 9) COMPANION SPECIES
# =========================================================
comp_species <- syn0 %>%
  filter(All >= min_all_comp) %>%
  filter(!(s_code %in% diag_species$s_code)) %>%
  select(s_code)

# =========================================================
# 10) OTHER NON-DIAGNOSTIC SPECIES
# =========================================================
other_species <- syn0 %>%
  filter(!(s_code %in% diag_species$s_code)) %>%
  filter(!(s_code %in% comp_species$s_code)) %>%
  select(s_code)

# =========================================================
# 11) FINAL TABLE CONTENT
# =========================================================
final_species <- bind_rows(
  diag_species %>% select(s_code) %>% distinct(),
  comp_species,
  other_species
) %>% distinct()

syn_final <- syn0 %>%
  filter(s_code %in% final_species$s_code) %>%
  left_join(diag_species %>% select(s_code, best_block, best_phi), by = "s_code") %>%
  mutate(
    block = case_when(
      s_code %in% diag_species$s_code ~ best_block,
      s_code %in% comp_species$s_code ~ "Companion species",
      TRUE ~ "Other non-diagnostic species"
    ),
    best_phi = ifelse(is.na(best_phi), -Inf, best_phi)
  )

block_order <- c("A", "B", "C", "D1", "D2", "E1", "E2", "F", "Companion species", "Other non-diagnostic species")
syn_final$block <- factor(syn_final$block, levels = block_order)

for(nm in c("A", "B", "C", "D1", "D2", "E1", "E2", "F")) {
  if(!nm %in% names(syn_final)) syn_final[[nm]] <- 0
}

syn_final <- syn_final %>% arrange(block, desc(best_phi), desc(All), desc(TPC), sp_name)

species_rows <- syn_final %>%
  transmute(Block = as.character(block), Species = sp_name, LF = life_form, TPC = round(TPC, 2), All = round(All),
            A = round(A), B = round(B), C = round(C), D1 = round(D1), D2 = round(D2), E1 = round(E1), E2 = round(E2), F = round(F))

# =========================================================
# 12) HEADER / META ROWS
# =========================================================
plot_counts <- class_df %>% count(lower_class) %>% mutate(class_label = class_labels[lower_class])

plot_count_row <- data.frame(Block = "Number of plots", Species = "", LF = "", TPC = NA, All = nrow(comm_raw),
                             A = plot_counts$n[match("A", plot_counts$class_label)], B = plot_counts$n[match("B", plot_counts$class_label)], C = plot_counts$n[match("C", plot_counts$class_label)],
                             D1 = plot_counts$n[match("D1", plot_counts$class_label)], D2 = plot_counts$n[match("D2", plot_counts$class_label)],
                             E1 = plot_counts$n[match("E1", plot_counts$class_label)], E2 = plot_counts$n[match("E2", plot_counts$class_label)], F = plot_counts$n[match("F", plot_counts$class_label)])

block_headers <- diag_counts %>%
  transmute(Block = paste0("Diagnostic species ", best_block, " (", n_diag, " taxa)"), Species = "", LF = "", TPC = NA, All = NA,
            A = NA, B = NA, C = NA, D1 = NA, D2 = NA, E1 = NA, E2 = NA, F = NA, order_block = match(best_block, block_order))

comp_header <- data.frame(Block = "Companion species", Species = "", LF = "", TPC = NA, All = NA, A = NA, B = NA, C = NA, D1 = NA, D2 = NA, E1 = NA, E2 = NA, F = NA, order_block = match("Companion species", block_order))
other_header <- data.frame(Block = "Other non-diagnostic species", Species = "", LF = "", TPC = NA, All = NA, A = NA, B = NA, C = NA, D1 = NA, D2 = NA, E1 = NA, E2 = NA, F = NA, order_block = match("Other non-diagnostic species", block_order))

block_headers <- bind_rows(block_headers, comp_header, other_header)
species_rows$order_block <- match(species_rows$Block, block_order)

out_list <- list(plot_count_row)
for(bl in block_order) {
  hdr <- block_headers %>% filter(order_block == match(bl, block_order))
  spp <- species_rows %>% filter(Block == bl)
  if(nrow(hdr) > 0) out_list[[length(out_list) + 1]] <- hdr %>% select(-order_block)
  if(nrow(spp) > 0) out_list[[length(out_list) + 1]] <- spp %>% select(-order_block)
}
table2_final <- bind_rows(out_list)

# =========================================================
# 13) EXTRA OUTPUTS: OTHER SPECIES LIST
# =========================================================
other_species_list <- syn0 %>%
  filter(s_code %in% other_species$s_code) %>%
  transmute(s_code = s_code, Species = sp_name, LF = life_form, TPC = round(TPC, 2), TPC_occ = round(TPC_occ, 2), All = round(All), n_occ = n_occ,
            A = round(A), B = round(B), C = round(C), D1 = round(D1), D2 = round(D2), E1 = round(E1), E2 = round(E2), F = round(F)) %>%
  arrange(desc(All), desc(TPC), Species)

# =========================================================
# 14) EXPORT
# =========================================================
write.csv(table2_final, "Table2_final_synoptic_table_with_other_species.csv", row.names = FALSE)
write.csv(diag_counts, "Table2_diagnostic_counts_by_class.csv", row.names = FALSE)
write.csv(diag_species, "Table2_selected_diagnostic_species.csv", row.names = FALSE)
write.csv(other_species_list, "Other_non_diagnostic_species_list.csv", row.names = FALSE)

cat("\n===== diagnostic counts =====\n"); print(diag_counts)
cat("\n===== companion species n =====\n"); print(nrow(comp_species))
cat("\n===== other non-diagnostic species n =====\n"); print(nrow(other_species))
cat("\n===== preview =====\n"); print(head(table2_final, 100))







# =========================================================
# SETTINGS
# =========================================================
top_n_diag      <- 12
min_const_diag  <- 10
min_all_comp    <- 10


# =========================================================
# 2) PREPARE MATRIX
# =========================================================
common_plots <- intersect(rownames(comm_raw_f4), class_df$plot_id)
comm_raw <- comm_raw_f4[common_plots, , drop = FALSE]
class_df <- class_df %>% filter(plot_id %in% common_plots)
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]
class_df$lower_class <- as.character(class_df$lower_class)

comm_raw <- comm_raw[rowSums(comm_raw) > 0, , drop = FALSE]
comm_raw <- comm_raw[, colSums(comm_raw) > 0, drop = FALSE]
class_df <- class_df[match(rownames(comm_raw), class_df$plot_id), ]

comm_pa <- (comm_raw > 0) * 1

# =========================================================
# 3) LABELS / SPECIES METADATA
# =========================================================
class_labels <- c("1"="A", "2"="B", "3"="C", "4"="D1", "5"="D2", "6"="E1", "7"="E2", "8"="F")
species_in_matrix <- colnames(comm_raw)
sp_match <- sp_code %>% filter(s_code %in% species_in_matrix)
sp_match <- sp_match[match(species_in_matrix, sp_match$s_code), ]

if(!"sp_name" %in% names(sp_match)) sp_match$sp_name <- sp_match$s_code
if(!"life_form" %in% names(sp_match)) stop("Column 'life_form' not found in df_species_code.csv")
lf_df <- sp_match %>% select(s_code, sp_name, life_form)

# =========================================================
# 4) CONSTANCY + PHI FOR GROUPING
# =========================================================
constancy_long <- as.data.frame(comm_pa) %>%
  rownames_to_column("plot_id") %>%
  left_join(class_df[, c("plot_id", "lower_class")], by = "plot_id") %>%
  pivot_longer(cols = -c(plot_id, lower_class), names_to = "s_code", values_to = "pa") %>%
  group_by(lower_class, s_code) %>%
  summarise(constancy = mean(pa, na.rm = TRUE) * 100, .groups = "drop")

overall_constancy <- as.data.frame(comm_pa) %>%
  rownames_to_column("plot_id") %>%
  pivot_longer(cols = -plot_id, names_to = "s_code", values_to = "pa") %>%
  group_by(s_code) %>%
  summarise(All_constancy = mean(pa, na.rm = TRUE) * 100, n_occ = sum(pa > 0, na.rm = TRUE), .groups = "drop")

phi_one <- function(x, g) {
  x <- as.numeric(x); g <- as.numeric(g)
  a <- as.numeric(sum(x == 1 & g == 1, na.rm = TRUE)); b <- as.numeric(sum(x == 1 & g == 0, na.rm = TRUE))
  c <- as.numeric(sum(x == 0 & g == 1, na.rm = TRUE)); d <- as.numeric(sum(x == 0 & g == 0, na.rm = TRUE))
  denom <- sqrt((a + b) * (c + d) * (a + c) * (b + d))
  if(is.na(denom) || !is.finite(denom) || denom == 0) return(NA_real_)
  (a * d - b * c) / denom
}

class_levels <- sort(unique(class_df$lower_class))
phi_list <- vector("list", length(class_levels)); names(phi_list) <- class_levels

for(cl in class_levels) {
  g <- ifelse(class_df$lower_class == cl, 1, 0)
  phi_vals <- sapply(colnames(comm_pa), function(sp) { phi_one(comm_pa[, sp], g) })
  phi_list[[cl]] <- data.frame(s_code = names(phi_vals), lower_class = cl, class_label = class_labels[cl], phi = as.numeric(phi_vals), stringsAsFactors = FALSE)
}
phi_table <- bind_rows(phi_list)

const_rank <- constancy_long %>% left_join(phi_table, by = c("s_code", "lower_class")) %>% mutate(class_label = class_labels[lower_class])

best_class <- const_rank %>% filter(!is.na(phi), phi > 0, constancy >= min_const_diag) %>%
  group_by(s_code) %>% arrange(desc(phi), desc(constancy), .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  select(s_code, best_block = class_label, best_phi = phi, best_constancy = constancy)

diag_species <- best_class %>% group_by(best_block) %>% arrange(desc(best_phi), desc(best_constancy), .by_group = TRUE) %>% slice_head(n = top_n_diag) %>% ungroup()
comp_species <- overall_constancy %>% filter(All_constancy >= min_all_comp, !(s_code %in% diag_species$s_code)) %>% select(s_code)
other_species <- overall_constancy %>% filter(!(s_code %in% diag_species$s_code), !(s_code %in% comp_species$s_code)) %>% select(s_code)

# =========================================================
# 5) COVER SYNOPTIC VALUES
# =========================================================
cover_long <- as.data.frame(comm_raw) %>%
  rownames_to_column("plot_id") %>%
  left_join(class_df[, c("plot_id", "lower_class")], by = "plot_id") %>%
  pivot_longer(cols = -c(plot_id, lower_class), names_to = "s_code", values_to = "cover")

cover_allplots_long <- cover_long %>% group_by(lower_class, s_code) %>% summarise(mean_cover = mean(cover, na.rm = TRUE), .groups = "drop") %>% mutate(class_label = class_labels[lower_class])
cover_allplots_wide <- cover_allplots_long %>% select(s_code, class_label, mean_cover) %>% pivot_wider(names_from = class_label, values_from = mean_cover, values_fill = 0)
overall_cover_allplots <- cover_long %>% group_by(s_code) %>% summarise(All = mean(cover, na.rm = TRUE), .groups = "drop")

cover_occ_long <- cover_long %>% group_by(lower_class, s_code) %>% summarise(mean_cover_occ = ifelse(sum(cover > 0, na.rm = TRUE) > 0, mean(cover[cover > 0], na.rm = TRUE), NA_real_), .groups = "drop") %>% mutate(class_label = class_labels[lower_class])
cover_occ_wide <- cover_occ_long %>% select(s_code, class_label, mean_cover_occ) %>% pivot_wider(names_from = class_label, values_from = mean_cover_occ, values_fill = 0)
overall_cover_occ <- cover_long %>% group_by(s_code) %>% summarise(All = ifelse(sum(cover > 0, na.rm = TRUE) > 0, mean(cover[cover > 0], na.rm = TRUE), NA_real_), .groups = "drop")

for(nm in c("A", "B", "C", "D1", "D2", "E1", "E2", "F")) {
  if(!nm %in% names(cover_allplots_wide)) cover_allplots_wide[[nm]] <- 0
  if(!nm %in% names(cover_occ_wide)) cover_occ_wide[[nm]] <- 0
}

# =========================================================
# 6) FINAL SPECIES SET
# =========================================================
final_species <- bind_rows(diag_species %>% select(s_code) %>% distinct(), comp_species, other_species) %>% distinct()
species_block_df <- bind_rows(diag_species %>% transmute(s_code, block = best_block, best_phi = best_phi), comp_species %>% transmute(s_code, block = "Companion species", best_phi = -Inf), other_species %>% transmute(s_code, block = "Other non-diagnostic species", best_phi = -Inf))
lf_join <- lf_df %>% select(s_code, sp_name, life_form)

# A) COVER TABLE: all plots
cover_synoptic_allplots <- cover_allplots_wide %>% left_join(overall_cover_allplots, by = "s_code") %>% rename(All = All) %>% left_join(lf_join, by = "s_code") %>% left_join(species_block_df, by = "s_code") %>% filter(s_code %in% final_species$s_code)
cover_synoptic_allplots$block <- factor(cover_synoptic_allplots$block, levels = c("A", "B", "C", "D1", "D2", "E1", "E2", "F", "Companion species", "Other non-diagnostic species"))
cover_rows_allplots <- cover_synoptic_allplots %>% arrange(block, desc(best_phi), desc(All), sp_name) %>% transmute(Block = as.character(block), Species = sp_name, LF = life_form, All = round(All, 2), A = round(A, 2), B = round(B, 2), C = round(C, 2), D1 = round(D1, 2), D2 = round(D2, 2), E1 = round(E1, 2), E2 = round(E2, 2), F = round(F, 2))

# B) COVER TABLE: occupied plots
cover_synoptic_occ <- cover_occ_wide %>% left_join(overall_cover_occ, by = "s_code") %>% rename(All = All) %>% left_join(lf_join, by = "s_code") %>% left_join(species_block_df, by = "s_code") %>% filter(s_code %in% final_species$s_code)
cover_synoptic_occ$block <- factor(cover_synoptic_occ$block, levels = c("A", "B", "C", "D1", "D2", "E1", "E2", "F", "Companion species", "Other non-diagnostic species"))
cover_rows_occ <- cover_synoptic_occ %>% arrange(block, desc(best_phi), desc(All), sp_name) %>% transmute(Block = as.character(block), Species = sp_name, LF = life_form, All = round(All, 2), A = round(A, 2), B = round(B, 2), C = round(C, 2), D1 = round(D1, 2), D2 = round(D2, 2), E1 = round(E1, 2), E2 = round(E2, 2), F = round(F, 2))

# =========================================================
# 7) HEADER ROWS & EXPORT
# =========================================================
make_header_row <- function(txt) { data.frame(Block = txt, Species = "", LF = "", All = NA, A = NA, B = NA, C = NA, D1 = NA, D2 = NA, E1 = NA, E2 = NA, F = NA, stringsAsFactors = FALSE) }
blocks <- c("A", "B", "C", "D1", "D2", "E1", "E2", "F", "Companion species", "Other non-diagnostic species")

assemble_tab <- function(rows) { bind_rows(lapply(blocks, function(bl) { hdr <- make_header_row(if(bl %in% c("Companion species", "Other non-diagnostic species")) bl else paste("Diagnostic species", bl)); bind_rows(hdr, rows %>% filter(Block == bl)) })) }

write.csv(assemble_tab(cover_rows_allplots), "Cover_synoptic_table_mean_cover_all_plots.csv", row.names = FALSE)
write.csv(assemble_tab(cover_rows_occ), "Cover_synoptic_table_mean_cover_occupied_plots.csv", row.names = FALSE)

cat("\nSaved:\n- Cover_synoptic_table_mean_cover_all_plots.csv\n- Cover_synoptic_table_mean_cover_occupied_plots.csv\n")
gc()