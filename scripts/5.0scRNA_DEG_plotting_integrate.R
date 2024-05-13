library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringi)
library(stringr)

# Plot DEGs
args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # sample_path <- "sample/BA-1_209_NormR/BA-1_209_NormR_final.rds" "BA-1_209_NormR_Object_Annotated.rds" # BA-1_209_NormR_final.rds
df_cluster_path <- args[2]
df_celltype_path <- args[3]
df_singler1_path <- args[4] # "sample/BA-1_209_NormR/BA-1_209_NormR_SingleR_encode_main.csv"
df_singler2_path <- args[5]

sample <- head(tail(strsplit(sample_path, "/", fixed = TRUE)[[1]], n=2), n=1)
analysis_out <- sapply(strsplit(sample_path, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

# sample_path <- "Merge_all/Merge_all_final.rds"
rnaAggr <- readRDS(sample_path)

# # write NA values as unknown
# rnaAggr@meta.data$SingleR_hpca_main

rnaAggr@meta.data <- rnaAggr@meta.data %>%
  mutate(SingleR_hpca_main = ifelse(is.na(SingleR_hpca_main), 'unknown', SingleR_hpca_main))

rnaAggr@meta.data <- rnaAggr@meta.data %>%
  mutate(SingleR_encode_main = ifelse(is.na(SingleR_encode_main), 'unknown', SingleR_encode_main))

DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- JoinLayers(rnaAggr)
rnaAggr <- ScaleData(object = rnaAggr, features = rownames(rnaAggr))
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.4, reduction = "pca")

plot_DEG_markers <- function(rnaAggr, var_col, analysis_out, num_top=10){
  # var_col = "SingleR_encode_main"
  DefaultAssay(rnaAggr) <-"RNA"
  rnaAggr <- SetIdent(rnaAggr, value = var_col)
  
  markers <- FindAllMarkers(rnaAggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  names(markers)[names(markers) == "avg_log2FC"] <- "avglogFC"
  
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(avglogFC > 2) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    slice_head(n = num_top) %>%
    ungroup() -> top10
  
  pdf(file.path(analysis_out, paste0(sample, "_DEG_heatmap_", var_col, ".pdf")),width=16, height=16)
  p <- DoHeatmap(rnaAggr, features = top10$gene) + NoLegend()
  print(p)
  dev.off()
  
  png(file.path(analysis_out, paste0(sample, "_DEG_heatmap_", var_col, ".png")),width=16, height=16, units='in', res=300)
  p <- DoHeatmap(rnaAggr, features = top10$gene) + NoLegend()
  print(p)
  dev.off()
  
  markers.to.plot <- unique(top10$gene)
  pdf(file.path(analysis_out, paste0(sample,"_DEG_dotplot_", var_col, ".pdf")),width=14, height=16)
  p <- DotPlot(object = rnaAggr, features = markers.to.plot) & coord_flip() & theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
  
  png(file.path(analysis_out, paste0(sample,"_DEG_dotplot_", var_col, ".png")),width=16, height=16, units='in', res=300)
  p <- DotPlot(object = rnaAggr, features = markers.to.plot) & coord_flip() & theme(axis.text.x = element_text(angle = 90))
  print(p)
  dev.off()
  
  return(markers)
}

markers_cluster <- plot_DEG_markers(rnaAggr, "seurat_clusters", analysis_out)
write.csv(markers_cluster, df_cluster_path)

markers_cluster <- plot_DEG_markers(rnaAggr, "assign.ident.breast", analysis_out)
write.csv(markers_cluster, df_celltype_path)

markers_cluster <- plot_DEG_markers(rnaAggr, "SingleR_hpca_main", analysis_out)
write.csv(markers_cluster, df_singler1_path)

markers_cluster <- plot_DEG_markers(rnaAggr, "SingleR_encode_main", analysis_out)
write.csv(markers_cluster, df_singler2_path)



