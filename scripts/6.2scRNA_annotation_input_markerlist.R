library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)
library(reshape2)
library(tidyverse)

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # sample_path <- "sample/final.rds"
signatures_path <- args[2] # anno1_dir <- "sample/BA-1_209_NormR/BA-1_209_NormR_SingleR_hpca_main.csv"
output_dir <- args[3]

analysis_out <- sapply(strsplit(output_dir, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))
sample <- head(tail(strsplit(output_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

signatures <- read.csv(signatures_path)
signatures <- stack(signatures)
# head(signatures)
colnames(signatures) <- c("genes", "celltype")
signatures <-signatures %>% select("celltype",  "genes"  )
# save(signatures, file = "signatures.Rdata")
# write.csv(df.sig, "breast_signatures.csv", row.names = F)

head(signatures)
celltypes <- as.character(unique(signatures[,1]))
signature_list <- sapply(1:length(celltypes),function(x){
  return((list(signatures[which(signatures[,1]==celltypes[x]),2])))})
names(signature_list) <- celltypes
# saveRDS(signature_list, "breast_signature_list.rds")

print(signature_list)

rnaAggr <- readRDS( file.path(sample_path))
rnaAggr <- JoinLayers(rnaAggr)
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.4, reduction = "pca")
Idents(rnaAggr) <- rnaAggr$seurat_clusters
markers <- FindAllMarkers(rnaAggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
names(markers)[names(markers) == "avg_log2FC"] <- "avglogFC"


RNA <- rnaAggr
genes = markers
min.score = 0.5

cluster_celltype_score = sapply(as.integer(unique(RNA@meta.data$seurat_clusters))-1, function(x){
  idx = genes$cluster==x
  avglogFC = genes$avglogFC[idx]
  names(avglogFC) = toupper(genes$gene[idx])
  score_cluster = sapply(signature_list, function(y){
    score = sum(avglogFC[y], na.rm = TRUE) / log2(length(y))
    return(score)
  })
})

colnames(cluster_celltype_score) = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
cellscore_max = apply(cluster_celltype_score, 2, max, na.rm = TRUE)
cellscore_max_celltype = apply(cluster_celltype_score, 2, function(x){
  if (max(x) < min.score){
    return("Others")
  }else{
    return(rownames(cluster_celltype_score)[which.max(x)])
  }
})

RNA@meta.data$custom.marker.assign.ident = as.character(as.integer(RNA@meta.data$seurat_clusters)-1)
current.cluster.ids = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
new.cluster.ids = cellscore_max_celltype

RNA@meta.data$custom.marker.assign.score = cellscore_max[RNA@meta.data$custom.marker.assign.ident]
RNA@meta.data$custom.marker.assign.ident = plyr::mapvalues(x = RNA@meta.data$custom.marker.assign.ident,
                                             from = current.cluster.ids, to = new.cluster.ids)

pdf(file.path( analysis_out, paste0(sample, "_UMAP_by_celltype_custom_marker_list_based.pdf")))
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "custom.marker.assign.ident", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

png(file.path( analysis_out, paste0(sample, "_UMAP_by_celltype_custom_marker_list_based.png")),  width=10, height=12, units='in', res=300)
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "custom.marker.assign.ident", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

write.csv(RNA@meta.data, file.path(analysis_out, paste0(sample, "_custom_marker_annotation.csv")))
saveRDS(RNA , file.path(output_dir))

