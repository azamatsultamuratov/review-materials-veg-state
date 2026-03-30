README

Title
Data, metadata, and analysis scripts for the manuscript:
Alternative vegetation states and assembly pathways on the dried Aral Sea bed

Overview
This archive contains the processed vegetation-plot dataset, associated metadata, and R scripts used in the analyses presented in the manuscript. The materials were prepared for vegetation classification, ordination, diversity analyses, ecological-type interpretation, and phylogenetic analyses.

Related manuscript
Alternative vegetation states and assembly pathways on the dried Aral Sea bed

Archive contents

Data and metadata files

1. ready_data_AD_DM.csv
Processed vegetation plot × species matrix used as the main analytical dataset.

Rows represent vegetation plots (relevés).
Columns represent species codes.
The first column contains unique plot identifiers.

Species values in the matrix are given as cover values.
Zero indicates absence of a species in a plot.

2. df_species_code.csv
Primary species metadata table.

This file links the species codes used in the vegetation matrix to full species names and life-form categories.

Main fields:
- s_code: species code used in the vegetation matrix
- sp_name: full species name
- life_form: life-form category

Life-form categories:
- An = annuals
- Pr = perennials
- Sr = shrubs
- Ss = subshrubs

3. gen_df_species_code.csv
Auxiliary taxonomic reference table for phylogenetic analyses.

This file contains older or alternative taxonomic names used to improve species matching in phylogenetic workflows and external phylogenetic resources.

4. eco_tp.csv
Ecological-type metadata table.

This file links plot identifiers to interpreted ecological types.

Main fields:
- plot_id: vegetation plot identifier
- Eco_typ: ecological type assignment

Ecological-type abbreviations:
- Ga = halophyte
- Gp = gypsophyte
- Po = potamophyte
- Ps = psammophyte
- GaGp = halo-gypsophyte
- GaPo = halo-potamophyte
- GaPs = halo-psammophyte
- GpPs = gypso-psammophyte
- PoGa = potamo-halophyte
- PoPs = potamo-psammophyte
- PsGa = psammo-halophyte
- PsGp = psammo-gypsophyte
- PsPo = psammo-potamophyte
- GaGpPs = halo-gypso-psammophyte
- GaPsPo = halo-psammo-potamophyte
- PoGpGa = potamo-gypso-halophyte
- PoPsGa = potamo-psammo-halophyte

Relationships among files
- Plot identifiers in ready_data_AD_DM.csv correspond to plot_id in eco_tp.csv.
- Species columns in ready_data_AD_DM.csv correspond to s_code in df_species_code.csv.
- gen_df_species_code.csv serves as an auxiliary taxonomic reference for phylogenetic analyses.

Purpose of the archived files
These files were used to:
- link species codes in the vegetation matrix to taxon names;
- assign species to life-form groups;
- connect vegetation plots with ecological-type information;
- improve taxonomic matching for phylogenetic analyses;
- preprocess and filter the vegetation matrix for modified TWINSPAN;
- generate classification, ordination, diversity, ecological-type, and phylogenetic outputs.

Scripts included

The archive also contains R scripts used for data preprocessing, classification, ordination, diversity analyses, ecological-type analyses, and phylogenetic analyses. The scripts were used to generate the analytical datasets, derived tables, and figures reported in the manuscript.

Main analytical components

1. Data preprocessing and filtering
The vegetation matrix was cleaned, empty plots and empty species were removed, and frequency-filtered datasets were generated. Matrices filtered at different species-frequency thresholds were created for sensitivity analysis, and square-root transformed versions were exported for modified TWINSPAN classification. :contentReference[oaicite:0]{index=0}

2. Evaluation of classification levels
Alternative numbers of clusters were evaluated using three criteria:
- number of significant indicator species,
- minimum cluster size,
- mean silhouette width.
These evaluations were carried out for the full dataset and for annual, perennial, and woody subsets. :contentReference[oaicite:1]{index=1}

3. Final modified TWINSPAN classification
The scripts produce the final hierarchical classification, including upper and lower cluster levels, cluster membership tables, indicator species tables, dominant species summaries, cluster sizes, and diversity statistics. :contentReference[oaicite:2]{index=2}

4. Within-class compositional homogeneity
Within-class Sørensen similarity was calculated for upper, lower, and grouped classification levels, including separate summaries for annual, perennial, and woody subsets. :contentReference[oaicite:3]{index=3}

5. Ordination analyses
The scripts include:
- three-dimensional NMDS ordination based on species composition,
- DCA ordinations for the full dataset and life-form subsets,
- correspondence analysis of vegetation clusters and ecological types using eco_tp.csv. :contentReference[oaicite:4]{index=4}

6. Richness and vegetation cover analyses
Plot-level species richness and vegetation cover were calculated for total vegetation, annuals, perennials, and woody species, followed by group summaries and one-way ANOVA. :contentReference[oaicite:5]{index=5}

7. Synoptic and diagnostic summaries
The scripts generate synoptic tables based on species constancy, fidelity, mean cover, and diagnostic species assignment for the lower vegetation clusters. :contentReference[oaicite:6]{index=6}

8. Hierarchy reconstruction and visualization
The archive includes scripts used to reconstruct and visualize the modified TWINSPAN hierarchy as cladograms, dendrograms, and node-level differential summaries. :contentReference[oaicite:7]{index=7}

9. Phylogenetic analyses
Phylogenetic workflows include taxonomic matching, tree construction, plot-level phylogenetic diversity metrics, standardized effect sizes of mean pairwise distance, and a phylogenetic heatmap of species occurrence across vegetation clusters. The auxiliary taxonomic file gen_df_species_code.csv was used to improve taxonomic coverage in these analyses. :contentReference[oaicite:8]{index=8}

Notes
- The main vegetation matrix contains plot-level species cover values.
- Zero values represent species absence.
- The primary species metadata are provided in df_species_code.csv.
- The file gen_df_species_code.csv was used specifically for phylogenetic analyses and includes older or alternative names to improve taxonomic coverage.
- The file eco_tp.csv provides ecological-type information for vegetation plots through plot identifiers.
- Legacy plots derived from published sources are retained with source attribution in the associated metadata.

Access conditions
These files are provided for peer review through confidential review access.
A repository-based public archive with a persistent identifier will be provided upon acceptance and cited in the final published version of the manuscript.
Access conditions applying to plot-level vegetation data will be described separately where relevant.
