README

Title
Data, metadata, and analysis scripts for anonymous peer review

Overview
This archive contains the processed vegetation-plot dataset, associated metadata, and R scripts used in the analyses. The archived materials support vegetation classification, ordination, diversity analyses, ecological-type interpretation, synoptic summaries, and phylogenetic analyses.

Related manuscript
Details withheld for anonymous peer review.

Archive contents

A. Data and metadata files

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

B. Scripts

This archive includes 12 main R scripts used to preprocess the data, perform the analyses, and generate derived outputs.

1. Step1_TWINSPAN_Filtering
Reads the main vegetation matrix, removes empty plots and species, applies frequency filtering, and exports raw and square-root transformed matrices for modified TWINSPAN analyses. The script generates filtered datasets for frequency thresholds of 3, 4, and 5 and retains the frequency ≥ 4 matrix for the main analyses.

2. Step2_Optimal_Clusters
Evaluates alternative numbers of clusters for modified TWINSPAN using three criteria:
- number of significant indicator species,
- minimum cluster size,
- mean silhouette width.
The analysis is performed for the full dataset and for annual, perennial, and woody subsets.

3. Step3.2_Sorensen_Similarity
Calculates within-class Sørensen similarity to assess compositional homogeneity of syntaxa. The script summarizes similarity values across upper, lower, and grouped classification levels and across life-form subsets.

4. Step3_TWINSPAN_Classification_and_Diversity
Performs the final modified TWINSPAN classification and defines the hierarchical syntaxonomic structure with 6 upper and 8 lower clusters. The script also calculates indicator species, species summaries, dominant species, cluster sizes, and alpha/beta diversity statistics.

5. Step4_Dendrogram_Visualizations
Reconstructs the modified TWINSPAN hierarchy and exports two diagram types:
- a ranked cladogram,
- a heterogeneity-based dendrogram.

6. Step5_NMDS
Performs three-dimensional NMDS ordination based on species composition, calculates plot scores and centroids, and exports the final ordination graphic with projected ellipses.

7. Step6_DCA
Performs DCA on the full dataset and on annual, perennial, and woody subsets. The script exports main and appendix ordination graphics, together with site scores and axis summaries.

8. Step7_Richness_Cover_ANOVA
Calculates plot-level richness and vegetation cover for total vegetation, annuals, perennials, and woody species. The script performs one-way ANOVA across vegetation clusters and exports summary tables and boxplots.

9. Step8_CA
Performs correspondence analysis using vegetation clusters and ecological-type data from eco_tp.csv. The script generates an interactive CA biplot and exports the corresponding outputs.

10. Step9_synoptic_table
Constructs synoptic tables using species constancy, fidelity, and cover summaries. The script identifies diagnostic species, companion species, and other non-diagnostic species and exports the final synoptic table and related summaries.

11. Step10_TWINSPAN_TXT_Report
Exports internal modified TWINSPAN outputs and textual hierarchy reports. This step includes:
- internal object summaries,
- hierarchy diagnostics,
- differential species and pseudospecies exports,
- reconstructed original-style text hierarchy reports.

12. Step11_Phylogenetic_Diversity
Performs phylogenetic analyses, including taxonomic matching, tree construction, Faith’s phylogenetic diversity, SES.MPD, and phylogenetic heatmap generation. The script also exports plot-level phylogenetic metrics and cluster-level summaries.

13. Step12_Tanglegram
Constructs a tanglegram comparing the modified TWINSPAN hierarchy with the phylogenetic tree structure and exports the final comparison graphic.

Purpose of the archived files
These files were used to:
- link species codes in the vegetation matrix to taxon names;
- assign species to life-form groups;
- connect vegetation plots with ecological-type information;
- improve taxonomic matching for phylogenetic analyses;
- preprocess and filter the vegetation matrix for modified TWINSPAN;
- evaluate alternative classification levels;
- generate classification, ordination, diversity, synoptic, ecological-type, and phylogenetic outputs.

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
