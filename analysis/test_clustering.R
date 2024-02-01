# Random seeds
read_downsample_seed <- "0" # random seqtk seed for downsampling reads - 0 for no downsampling, integer >1 for downsampled seed
R_random_seed <- 100 # including for cell downsampling

# Downsample fractions
seu_read_fraction_after_downsampling <- "1.0" # fraction of reads after downsampling - any number from (0,1.0]
scan_read_fraction_after_downsampling <- "1.0" # fraction of reads after downsampling - any number from (0,1.0]

seu_cell_fraction_after_downsampling <- "1.0" # fraction of cells after downsampling - any number from (0,1.0]
scan_cell_fraction_after_downsampling <- "1.0" # fraction of cells after downsampling - any number from (0,1.0]

# Matrix generation methods
seu_matrix_generation_method <- "kb" # str["kb", "cellranger"]
seu_matrix_generation_method_version <- "0.28.0"
seu_matrix_include_unspliced <- FALSE # kb only; if TRUE --> .total.mtx (kb SUA); if FALSE --> .cell.mtx (KB SA) or .mtx (older versions of KB)
seu_matrix_source <- "generated" # "generated" (kb/cellranger count) or "downloaded" (eg 10x website)
seu_matrix_qc <- "raw" # "raw" or "filtered"
seu_used_batch <- FALSE

scan_matrix_generation_method <- "kb" # str["kb", "cellranger"]
scan_matrix_generation_method_version <- "0.28.0"
scan_matrix_include_unspliced <- FALSE # kb only; if TRUE --> .total.mtx (kb SUA); if FALSE --> .cell.mtx (KB SA) or .mtx (older versions of KB)
scan_matrix_source <- "generated" # "generated" (kb/cellranger count) or "downloaded" (eg 10x website)
scan_matrix_qc <- "raw" # "raw" or "filtered"
scan_used_batch <- FALSE

# Function argument and input settings
analysis_methods <- "scanpy_like" # str["default", "seurat_like", "scanpy_like"]
data_input <- "scanpy" # str["default", "seurat", "scanpy"]

# Package versions
seurat_version <- "4.3.0"
scanpy_version <- "1.9.5"

# Custom parameters
inflection_UMI_manual <- NULL # number >=0; or NULL to have automatic selection, especially necessary for lower fracs (e.g., 30 for frac=0.02, 20 for frac=0.01)
min_cells <- 3
min_features <- 200
max_pct_mct <- 20 # default 5
max_n_genes_by_counts_scanpy <- 12000 # default 2500

seu_num_pcs <- 50 # number 1-50; or NULL to select after elbow plot visualization
scan_num_pcs <- 50
umap_knn_k <- 50
umap_leiden_clustering_resolution <- 0.8

# Output path specifications
project_base_path <- "/workspace/analysis"
data_name <- "5k_pbmc_v3"   #!!!
include_seeds_in_file_paths <- FALSE
save_data <- FALSE
dpi <- 350



set.seed(R_random_seed)
options(max.print = 130)

seu_read_fraction_after_downsampling <- gsub("\\.", "_", as.character(seu_read_fraction_after_downsampling)) # fraction of reads after downsampling, as string representation and using underscores in place of decimal
scan_read_fraction_after_downsampling <- gsub("\\.", "_", as.character(scan_read_fraction_after_downsampling))
seu_cell_fraction_after_downsampling <- gsub("\\.", "_", as.character(seu_cell_fraction_after_downsampling))
scan_cell_fraction_after_downsampling <- gsub("\\.", "_", as.character(scan_cell_fraction_after_downsampling))

seu_matrix_generation_method_version <- gsub("\\.", "_", seu_matrix_generation_method_version)
scan_matrix_generation_method_version <- gsub("\\.", "_", scan_matrix_generation_method_version)

seu_included_transcripts <- ifelse(seu_matrix_include_unspliced == TRUE, "sua", "sa")
scan_included_transcripts <- ifelse(scan_matrix_include_unspliced == TRUE, "sua", "sa")

seu_matrix_generation_method_full <- glue::glue("{seu_matrix_generation_method}{seu_matrix_generation_method_version}_{seu_matrix_qc}_{seu_matrix_source}")
scan_matrix_generation_method_full <- glue::glue("{scan_matrix_generation_method}{scan_matrix_generation_method_version}_{scan_matrix_qc}_{scan_matrix_source}")

seu_kb_major_version <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", seu_matrix_generation_method_version))
scan_kb_major_version <- as.integer(sub("^[0-9]+_([0-9]+)_.*$", "\\1", scan_matrix_generation_method_version))

# Seurat vs Scanpy
group1_color <- "#D55E00"
group2_color <- "#56B4E9"



seu_data_path <- glue::glue("{project_base_path}/count_matrix_collection/{data_name}/{seu_matrix_generation_method_full}/frac{seu_read_fraction_after_downsampling}_seed{read_downsample_seed}")
scan_data_path <- glue::glue("{project_base_path}/count_matrix_collection/{data_name}/{scan_matrix_generation_method_full}/frac{scan_read_fraction_after_downsampling}_seed{read_downsample_seed}")



if (seu_matrix_generation_method == scan_matrix_generation_method &&
    seu_matrix_generation_method_version == scan_matrix_generation_method_version &&
    seu_included_transcripts == scan_included_transcripts &&
    seu_matrix_source == scan_matrix_source &&
    seu_matrix_qc == scan_matrix_qc) {
    matrix_generation_method_full <- glue::glue("{seu_matrix_generation_method_full}_{seu_included_transcripts}")
} else {
    matrix_generation_method_full <- glue::glue("seu_{seu_matrix_generation_method_full}_{seu_included_transcripts}__scan_{scan_matrix_generation_method_full}_{scan_included_transcripts}")
}

read_fraction_after_downsampling <- ifelse(seu_read_fraction_after_downsampling == scan_read_fraction_after_downsampling, seu_read_fraction_after_downsampling,
                                           paste("seu", seu_read_fraction_after_downsampling, "vs", "scan", scan_read_fraction_after_downsampling, sep = "_")
)

cell_fraction_after_downsampling <- ifelse(seu_cell_fraction_after_downsampling == scan_cell_fraction_after_downsampling, seu_cell_fraction_after_downsampling,
                                           paste("seu", seu_cell_fraction_after_downsampling, "vs", "scan", scan_cell_fraction_after_downsampling, sep = "_")
)

if (include_seeds_in_file_paths) {
    read_fraction_after_downsampling <- glue::glue("{read_fraction_after_downsampling}_seed{read_downsample_seed}")
    cell_fraction_after_downsampling <- glue::glue("{cell_fraction_after_downsampling}_seed{R_random_seed}")
}

output_base_path <- glue::glue("{project_base_path}/output/{data_name}/seuratv{seurat_version}_vs_scanpyv{scanpy_version}/methods_{analysis_methods}/input_{data_input}/{matrix_generation_method_full}/cell_fraction_{cell_fraction_after_downsampling}/read_fraction_{read_fraction_after_downsampling}")

output_data_file_paths <- list(
    markers = glue::glue("{output_base_path}/data_files/markers.rds"),
    results_scan = glue::glue("{output_base_path}/data_files/results_scan.rds"),
    markers2 = glue::glue("{output_base_path}/data_files/markers2.rds"),
    seu_object = glue::glue("{output_base_path}/data_files/seu.rds"),
    scan_adata = glue::glue("{output_base_path}/data_files/adata.h5ad")
)

# FALSE to have no save
file_paths <- list(
    euler_stats_before_QC_file = FALSE, # glue::glue("{output_base_path}/stats/euler_stats_beforeQC.txt"),
    euler_stats_after_QC_file = glue::glue("{output_base_path}/stats/euler_stats_afterQC.txt"),
    pca_knn_clustering_umap_file = glue::glue("{output_base_path}/stats/pca_knn_clustering_umap_stats.txt"),
    de_stats_file = glue::glue("{output_base_path}/stats/de_stats.txt"),
    pre_filtering_upset_cell = FALSE, # glue::glue("{output_base_path}/plots/pre_filtering_upset_cell.tiff"),
    pre_filtering_upset_gene = FALSE, # glue::glue("{output_base_path}/plots/pre_filtering_upset_gene.tiff"),
    
    knee_plot = FALSE, # glue::glue("{output_base_path}/plots/knee_plot.tiff"),
    umi_scatterplot = FALSE, # glue::glue("{output_base_path}/plots/umi_scatterplot.tiff"),
    
    violin_counts_comparison <- FALSE, # glue::glue("{output_base_path}/plots/violin_counts_comparison.tiff"),
    seu_violin_file_path = FALSE, # glue::glue("{output_base_path}/plots/seu_violin_plot.tiff"),
    scan_violin_file_path_genes = FALSE, # glue::glue("{output_base_path}/plots/scan_violin_plot_genes.tiff"),
    scan_violin_file_path_counts = FALSE, # glue::glue("{output_base_path}/plots/scan_violin_plot_counts.tiff"),
    scan_violin_file_path_mt = FALSE, # glue::glue("{output_base_path}/plots/scan_violin_plot_mt.tiff"),
    
    upset_cells = glue::glue("{output_base_path}/plots/upset_cells.tiff"),
    upset_genes = glue::glue("{output_base_path}/plots/upset_genes.tiff"),
    upset_hvgs = glue::glue("{output_base_path}/plots/upset_hvgs.tiff"),
    upset_markers_genes_only = glue::glue("{output_base_path}/plots/upset_marker_genes_only.tiff"),
    upset_markers = glue::glue("{output_base_path}/plots/upset_markers.tiff"),
    euler_before_qc_cell_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_cells_beforeQC.tiff"),
    euler_before_qc_gene_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_genes_beforeQC.tiff"),
    
    euler_after_qc_cell_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_cells_afterQC.tiff"),
    euler_after_qc_gene_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_genes_afterQC.tiff"),
    euler_after_qc_hvg_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_hvgs_afterQC.tiff"),
    euler_after_qc_marker_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_markers.tiff"),
    euler_after_qc_marker_manual_bonferroni_file_path = FALSE, # glue::glue("{output_base_path}/plots/euler_markers_manual_bonferroni.tiff"),
    euler_after_qc_marker_genes_only = FALSE, # glue::glue("{output_base_path}/plots/euler_markers_genes.tiff"),
    
    pca_elbow_filepath_combined = FALSE, # glue::glue("{output_base_path}/plots/pca_elbow_combined.tiff"),
    pca_12_overlay_filepath = glue::glue("{output_base_path}/plots/pca_scatterplot_12.tiff"),
    pca_34_overlay_filepath = FALSE, # glue::glue("{output_base_path}/plots/pca_scatterplot_34.tiff"),
    pca_loading_diffs = FALSE, # glue::glue("{output_base_path}/plots/pc_loading_diffs.tiff"),
    pca_eigs_diff = FALSE, # glue::glue("{output_base_path}/plots/pc_eig_diff.tiff"),
    pca_cluster_filepath_seu = FALSE, # glue::glue("{output_base_path}/plots/pca_scatterplot_clusters_seu.tiff"),
    pca_cluster_filepath_scan = FALSE, # glue::glue("{output_base_path}/plots/pca_scatterplot_clusters_scan.tiff"),
    combined_pc_variance_loadings_plot = glue::glue("{output_base_path}/plots/combined_pc_variance_loadings_plot.tiff"),
    jaccards = FALSE, # glue::glue("{output_base_path}/plots/jaccards.tiff"),
    knn_scatterplot = FALSE, # glue::glue("{output_base_path}/plots/knn_scatterplot.tiff"),
    jaccard_degree_scatterplot = glue::glue("{output_base_path}/plots/jaccard_degree_scatterplot.tiff"),
    pheatmap = FALSE, # glue::glue("{output_base_path}/plots/cluster_pheatmap.tiff"),
    alluvial = glue::glue("{output_base_path}/plots/cluster_alluvial.tiff"),
    alluvial_legend = glue::glue("{output_base_path}/plots/cluster_alluvial_legend.tiff"),
    alluvial_legend_high_alpha = glue::glue("{output_base_path}/plots/cluster_alluvial_legend_high_alpha.tiff"),
    umap_seu = glue::glue("{output_base_path}/plots/umap_seu.tiff"),
    umap_scan = glue::glue("{output_base_path}/plots/umap_scan.tiff"),
    umap_seu_clusters_scan = glue::glue("{output_base_path}/plots/umap_seu_clusters_scan.tiff"),
    umap_scan_clusters_seu = glue::glue("{output_base_path}/plots/umap_scan_clusters_seu.tiff"),
    umap_jaccard_degree_scatterplot = glue::glue("{output_base_path}/plots/umap_jaccard_degree_scatterplot.tiff"),
    umap_jaccard_knn_density = glue::glue("{output_base_path}/plots/umap_jaccard_knn_density.tiff"),
    umap_jaccard_knn_density_seu_facet = glue::glue("{output_base_path}/plots/umap_jaccard_knn_density_seu_facet.tiff"),
    umap_jaccard_knn_density_scan_facet = glue::glue("{output_base_path}/plots/umap_jaccard_knn_density_scan_facet.tiff"),
    umap_alluvial = glue::glue("{output_base_path}/plots/umap_alluvial.tiff"),
    umap_alluvial_legend = glue::glue("{output_base_path}/plots/umap_alluvial_legend.tiff"),
    umap_umap_leiden_seu = glue::glue("{output_base_path}/plots/umap_umap_leiden_seu.tiff"),
    umap_umap_leiden_scan = glue::glue("{output_base_path}/plots/umap_umap_leiden_scan.tiff"),
    logFC_histogram_magnitude_file_path = FALSE, # glue::glue("{output_base_path}/plots/logFC_histogram_magnitude.tiff"),
    logFC_histogram_signed_file_path = FALSE, # glue::glue("{output_base_path}/plots/logFC_histogram_signed.tiff"),
    wilcoxon_histogram_magnitude_file_path = FALSE, # glue::glue("{output_base_path}/plots/wilcoxon_histogram_magnitude.tiff"),
    wilcoxon_histogram_signed_file_path = FALSE, # glue::glue("{output_base_path}/plots/wilcoxon_histogram_signed.tiff"),
    
    logFC_scatterplot_file_path = glue::glue("{output_base_path}/plots/logFC_scatterplot.tiff"),
    wilcoxon_scatterplot_file_path = glue::glue("{output_base_path}/plots/wilcoxon_scatterplot.tiff"),
    logFC_scatterplot_file_path_with_legend = glue::glue("{output_base_path}/plots/logFC_scatterplot_with_legend.tiff"),
    wilcoxon_scatterplot_file_path_with_legend = glue::glue("{output_base_path}/plots/wilcoxon_scatterplot_with_legend.tiff"),
    logFC_scatterplot_outliers_removed_file_path = FALSE, # glue::glue("{output_base_path}/plots/logFC_scatterplot_no_outliers.tiff"),
    wilcoxon_scatterplot_outliers_removed_file_path = FALSE, # glue::glue("{output_base_path}/plots/wilcoxon_scatterplot_no_outliers.tiff"),
    
    logFC_boxplot_magnitude_file_path = FALSE, # glue::glue("{output_base_path}/plots/logFC_boxplot_magnitude.tiff"),
    logFC_boxplot_signed_file_path = FALSE, # glue::glue("{output_base_path}/plots/logFC_boxplot_signed.tiff"),
    wilcoxon_boxplot_magnitude_file_path = FALSE, # glue::glue("{output_base_path}/plots/wilcoxon_boxplot_magnitude.tiff"),
    wilcoxon_boxplot_signed_file_path = FALSE, # glue::glue("{output_base_path}/plots/wilcoxon_boxplot_signed.tiff"),
    
    FC_histogram_magnitude_file_path = FALSE, # glue::glue("{output_base_path}/plots/FC_histogram_magnitude.tiff"),
    FC_histogram_signed_file_path = FALSE # glue::glue("{output_base_path}/plots/FC_histogram_signed.tiff")
)

for (path in output_data_file_paths) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

if (save_data) {
    for (path in file_paths) {
        if (is.character(path)) {
            # Extract the directory part of the path
            specific_output_path <- dirname(path)
            
            # Create the directory if it does not exist
            if (!dir.exists(specific_output_path)) {
                dir.create(specific_output_path, recursive = TRUE, showWarnings = FALSE)
            }
        }
    }
} else {
    for (i in seq_along(file_paths)) {
        file_paths[[i]] <- FALSE
    }
}

for (file in c(file_paths$euler_stats_after_QC_file, file_paths$pca_knn_clustering_umap_file, file_paths$de_stats_file)) {
    if (is.character(file)) {
        sink(file = file, append = FALSE)
        sink()
    }
}



library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(eulerr)
library(scattermore)
library(DropletUtils)
library(glue)
library(bluster)
library(ggforce)
library(ggplotify)
library(grid)
library(gtable)
library(ggalluvial)
theme_set(theme_bw())

source(glue("{project_base_path}/scripts/data_analysis_helper.R"))
source(glue("{project_base_path}/scripts/plotting_and_stats.R"))



if (analysis_methods == "default" || analysis_methods == "scanpy_like") {
    scanpy_hvg_flavor <- "seurat"
    n_top_genes <- NULL
    scanpy_scale_max <- NULL
    scanpy_pca_zero_center <- TRUE
    scan_n_neighbors <- 15
    scanpy_clustering_algorithm <- "leiden"
    scanpy_resolution <- 1
    scanpy_cluster_iters <- -1
    scanpy_umap_min_dist <- 0.5
    scanpy_correction_method <- "benjamini-hochberg"
} else if (analysis_methods == "seurat_like") {
    scanpy_hvg_flavor <- "seurat_v3"
    n_top_genes <- 2000
    scanpy_scale_max <- 10
    scanpy_pca_zero_center <- FALSE
    scan_n_neighbors <- 20
    scanpy_clustering_algorithm <- "louvain"
    scanpy_resolution <- 0.8
    scanpy_cluster_iters <- 10
    scanpy_umap_min_dist <- 0.3
    scanpy_correction_method <- "bonferroni"
} else {
    paste(analysis_methods, "is not a valid input for analysis_methods. Please choose from 'default', 'seurat_like', or 'scanpy_like'.")
}

if (analysis_methods == "default" || analysis_methods == "seurat_like") {
    seurat_hvg_flavor <- "vst"
    seu_mean_cutoff <- c(0.1, 8)
    seu_dispersion_cutoff <- c(1, Inf)
    seu_vars_to_regress <- NULL
    seurat_scale_max <- 10
    seu_n_neighbors <- 20
    seurat_clustering_algorithm <- "louvain"
    seu_resolution <- 0.8
    seu_umap_method <- "uwot"
    seu_umap_min_dist <- 0.3
    seu_umap_metric <- "cosine"
    # correction method = bonferroni
} else if (analysis_methods == "scanpy_like") {
    seurat_hvg_flavor <- "mean.var.plot"
    seu_mean_cutoff <- c(0.0125, 3)
    seu_dispersion_cutoff <- c(0.5, Inf)
    seu_vars_to_regress <- c("nCount_RNA", "pct_mt")
    seurat_scale_max <- Inf
    seu_n_neighbors <- 15
    seurat_clustering_algorithm <- "leiden"
    seu_resolution <- 1
    seu_umap_method <- "umap-learn"
    seu_umap_min_dist <- 0.5
    seu_umap_metric <- "correlation"
    # correction method = benjamini-hochberg
} else {
    paste(analysis_methods, "is not a valid input for analysis_methods. Please choose from 'default', 'seurat_like', or 'scanpy_like'.")
}



seu <- readRDS(output_data_file_paths$seu_object)


seu <- FindClusters(seu, verbose = FALSE, algorithm = 4, resolution = seu_resolution)


saveRDS(seu, file = output_data_file_paths$seu_object)