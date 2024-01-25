group1_name <- "Seurat"
group2_name <- "Scanpy"
rds_path <- "/workspace/analysis/output/5k_pbmc_v3/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/markers2.rds"

library(tidyverse)
library(Seurat)
library(ggforce)
library(ggplotify)
library(ggalluvial)
library(glue)
theme_set(theme_bw())

source("/workspace/analysis/scripts/plotting_and_stats.R")

file_paths <- list()
file_paths$logFC_scatterplot_file_path <- glue::glue("/workspace/analysis/logFC_scatterplot.tiff")
file_paths$wilcoxon_scatterplot_file_path <- glue::glue("/workspace/analysis/wilcoxon_scatterplot.tiff")
file_paths$logFC_scatterplot_file_path_with_legend <- glue::glue("/workspace/analysis/logFC_scatterplot_with_legend.tiff")

markers2 <- readRDS(rds_path)

if ((group1_name == "Seurat" && group2_name == "Scanpy") || (group2_name == "Seurat" && group1_name == "Scanpy")) {
    markers2$p_val_adj_r[markers2$p_val_adj_r == 0] <- .Machine$double.xmin
    markers2$p_val_adj_py[markers2$p_val_adj_py == 0] <- .Machine$double.xmin
} else {
    markers2[[glue("p_val_adj.{group1_name}")]][markers2[[glue("p_val_adj.{group1_name}")]] == 0] <- .Machine$double.xmin
    markers2[[glue("p_val_adj.{group2_name}")]][markers2[[glue("p_val_adj.{group2_name}")]] == 0] <- .Machine$double.xmin
}

logFC_scatterplot <- plot_scatterplot_de_logfc(markers2, group1_name = group1_name, group2_name = group2_name, save = file_paths$logFC_scatterplot_file_path, outliers_excluded = FALSE)
pvaladj_scatterplot <- plot_scatterplot_de_wilcoxon(markers2, group1_name = group1_name, group2_name = group2_name, save = file_paths$wilcoxon_scatterplot_file_path, outliers_excluded = FALSE)

logFC_scatterplot_with_legend <- plot_scatterplot_de_logfc(markers2, group1_name = group1_name, group2_name = group2_name, save = file_paths$logFC_scatterplot_file_path_with_legend, outliers_excluded = FALSE, show_legend = TRUE)
