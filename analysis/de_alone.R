seu <- readRDS("/workspace/analysis/seu_before_de_scanlike.rds")
markers <- Seurat::FindAllMarkers(seu, logfc.threshold=0, min.pct=0, return.thresh=1.0001)
markers$p_val_adj = p.adjust(markers$p_val, method='BH')
saveRDS(markers, file = "/workspace/analysis/seu_markers_after_de_scanlike.rds")