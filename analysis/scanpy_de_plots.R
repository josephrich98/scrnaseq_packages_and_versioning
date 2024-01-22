library(tidyverse)
library(Seurat)
library(ggforce)
library(ggplotify)
library(ggalluvial)
library(glue)
theme_set(theme_bw())

source("/workspace/analysis/scripts/plotting_and_stats.R")

markers2 <- readRDS("/workspace/analysis/output/5k_pbmc_v3/seuratv4.3.0_vs_scanpyv1.9.5/methods_scanpy_like/input_scanpy/kb0_28_0_raw_generated_sa/cell_fraction_1_0/read_fraction_1_0/data_files/markers2.rds")

markers2$p_val_adj_r[markers2$p_val_adj_r == 0] <- .Machine$double.xmin
markers2$p_val_adj_py[markers2$p_val_adj_py == 0] <- .Machine$double.xmin

file_paths <- list()
file_paths$logFC_scatterplot_file_path = glue::glue("/workspace/analysis/logFC_scatterplot.tiff")
file_paths$wilcoxon_scatterplot_file_path = glue::glue("/workspace/analysis/wilcoxon_scatterplot.tiff")
file_paths$logFC_scatterplot_file_path_with_legend = glue::glue("/workspace/analysis/logFC_scatterplot_with_legend.tiff")
file_paths$wilcoxon_scatterplot_file_path_with_legend = glue::glue("/workspace/analysis/wilcoxon_scatterplot_with_legend.tiff")

seurat_vs_scanpy_logFC_scatterplot <- plot_scatterplot_de(markers2, metric="logFC", save = file_paths$logFC_scatterplot_file_path, outliers_excluded = FALSE)
seurat_vs_scanpy_pvaladj_scatterplot <- plot_scatterplot_de(markers2, metric="pvaladj", save = file_paths$wilcoxon_scatterplot_file_path, outliers_excluded = FALSE)

seurat_vs_scanpy_logFC_scatterplot_with_legend <- plot_scatterplot_de(markers2, metric="logFC", save = file_paths$logFC_scatterplot_file_path_with_legend, outliers_excluded = FALSE, show_legend = TRUE)
seurat_vs_scanpy_pvaladj_scatterplot_with_legend <- plot_scatterplot_de(markers2, metric="pvaladj", save = file_paths$wilcoxon_scatterplot_file_path_with_legend, outliers_excluded = FALSE, show_legend = TRUE)

