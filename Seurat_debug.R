library(Seurat)
library(Matrix)
library(tidyverse)
library(scattermore)
library(DropletUtils)
library(glue)

read_count_output_modified <- function(dir, name, unspliced = FALSE, batch = FALSE, tcc = FALSE) {
    dir <- normalizePath(dir, mustWork = TRUE)
    if (unspliced) {
        m <- readMM(paste0(dir, "/", name, ".total.mtx"))
    } else {
        if (file.exists(paste0(dir, "/", name, ".cell.mtx"))) {
            # If the first file exists, read from it
            m <- readMM(paste0(dir, "/", name, ".cell.mtx"))
        } else {
            m <- readMM(paste0(dir, "/", name, ".mtx"))
        }
    }
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # The matrix read has cells in rows
    ge <- if (tcc) ".ec.txt" else ".genes.txt"
    con_genes <- file(paste0(dir, "/", name, ge))
    if (batch) {
        con_bcs <- file(paste0(dir, "/", name, ".barcodes.combined.txt"))
    } else {
        con_bcs <- file(paste0(dir, "/", name, ".barcodes.txt"))
    }
    genes <- readLines(con_genes)
    barcodes <- readLines(con_bcs)
    colnames(m) <- barcodes
    rownames(m) <- genes
    close(con_genes)
    close(con_bcs)
    return(m)
}


seu1_data_path <- "/workspace/analysis/count_matrix_collection/SC3_v3_NextGem_SI_PBMC_10K/kb0_28_0/frac1_0_seed0"

res_mat1 <- read_count_output_modified(seu1_data_path, name = "cells_x_genes", tcc = FALSE)
tot_counts1 <- Matrix::colSums(res_mat1)
bc_rank1 <- barcodeRanks(res_mat1)

UMI_cutoff1 <- metadata(bc_rank1)$inflection

res_mat_filtered1 <- res_mat1[, tot_counts1 > UMI_cutoff1]
res_mat_filtered1 <- res_mat_filtered1[Matrix::rowSums(res_mat_filtered1) > 0, ]

seu1 <- CreateSeuratObject(counts = res_mat_filtered1, min.cells = 3, min.features = 200)

mt_genes <- data.frame(ensembl_gene_id = c("ENSG00000210049", "ENSG00000211459", "ENSG00000210077", "ENSG00000210082", "ENSG00000209082", "ENSG00000198888", "ENSG00000210100", "ENSG00000210107", "ENSG00000210112", "ENSG00000198763", "ENSG00000210117", "ENSG00000210127", "ENSG00000210135", "ENSG00000210140", "ENSG00000210144", "ENSG00000198804", "ENSG00000210151", "ENSG00000210154", "ENSG00000198712", "ENSG00000210156", "ENSG00000228253", "ENSG00000198899", "ENSG00000198938", "ENSG00000210164", "ENSG00000198840", "ENSG00000210174", "ENSG00000212907", "ENSG00000198886", "ENSG00000210176", "ENSG00000210184", "ENSG00000210191", "ENSG00000198786", "ENSG00000198695", "ENSG00000210194", "ENSG00000198727", "ENSG00000210195", "ENSG00000210196"))

assay_gene_names1 <- rownames(seu1[["RNA"]])
assay_gene_names_trimmed1 <- gsub("\\..*", "", assay_gene_names1)
common_genes1 <- intersect(mt_genes$ensembl_gene_id, assay_gene_names_trimmed1)
common_genes_with_version1 <- assay_gene_names1[match(common_genes1, assay_gene_names_trimmed1)]
seu1[["pct_mt"]] <- PercentageFeatureSet(seu1, features = common_genes_with_version1)

seu1 <- subset(seu1, pct_mt < 20)
seu1 <- NormalizeData(seu1, verbose = FALSE)

seu1 <- FindVariableFeatures(seu1, verbose = FALSE)

seu1 <- ScaleData(seu1, verbose = FALSE)

seu1 <- RunPCA(seu1, npcs = 50, verbose = FALSE, seed.use = 42)

# debug(Seurat::FindNeighbors)

seu1 <- FindNeighbors(seu1, reduction = "pca")

# undebug(Seurat::FindNeighbors)

# debug(Seurat::FindClusters)
# debug(FindClusters.Seurat)
# debug(FindClusters.default)
# undebug(Seurat::FindClusters)

seu1 <- FindClusters(seu1, verbose = FALSE, algorithm = 1)

# debug(Seurat::RunUMAP)
# undebug(Seurat::RunUMAP)

# seu1_umap_info <- RunUMAP(seu1, dims = 1:50)

debug(Seurat::FindAllMarkers)

markers_seu1 <- FindAllMarkers(seu1)

