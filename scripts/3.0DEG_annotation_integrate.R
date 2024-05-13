library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)

args <- commandArgs( trailingOnly = TRUE )
sample_path_tmp <- args[1]
sample_path = strsplit(sample_path_tmp,' ')[[1]]

# output_dir <- args[2]
output_dir <- strsplit(args[2],' ')[[1]]
analysis_out <- sapply(strsplit(output_dir, "/", fixed = TRUE),
                   function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

rnaAggr <- readRDS( file.path(sample_path))

rnaAggr <- JoinLayers(rnaAggr)
DefaultAssay(rnaAggr) <- "RNA"
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 0.4, reduction = "pca")
Idents(rnaAggr) <- rnaAggr$seurat_clusters
markers <- FindAllMarkers(rnaAggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
names(markers)[names(markers) == "avg_log2FC"] <- "avglogFC"

saveRDS(markers, file.path(analysis_out,paste0(integrate_sample_name, "_DEG_markers_by_cluster.rds")))
write.csv(markers, file.path(analysis_out,paste0(integrate_sample_name,"_DEG_markers_by_cluster.csv")))

# signature_list <- readRDS( "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/breast_signature_list.rds")
signature_list <- readRDS( "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/new_breast_marker_list.rds")


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

RNA@meta.data$assign.ident = as.character(as.integer(RNA@meta.data$seurat_clusters)-1)
current.cluster.ids = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
new.cluster.ids = cellscore_max_celltype

RNA@meta.data$assign.score = cellscore_max[RNA@meta.data$assign.ident]
RNA@meta.data$assign.ident = plyr::mapvalues(x = RNA@meta.data$assign.ident,
                                             from = current.cluster.ids, to = new.cluster.ids)

# head(RNA@meta.data)
pdf(file.path( analysis_out, paste0(integrate_sample_name,"_UMAP_by_celltype_breast_markerbased.pdf")))
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "assign.ident", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

png(file.path( analysis_out, paste0(integrate_sample_name,"_UMAP_by_celltype_breast_markerbased.png")), width=10, height=12, units='in', res=300)
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "assign.ident", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

# saveRDS(RNA , file.path(output_dir))

markers2 <- c( "FOXA1", "ESR1", "PGR", "PRLR", "SPDEF", # HS
               "ELF5", "MFGE8", "KIT", "ALDH1A3", "PROM1", #AV
               "OXTR", "ACTA2",  "MYLK","MYH11", "MYL9", # Basal
               "COL1A1","COL1A2","COL3A1", "COL6A2","LUM",# Fibroblasts
               "VWF","CDH5","RASIP1","RGS5","MSN", # Endothelial
               "ACTA2","MCAM","CSPG4","PDGFRB" # pericytes
               
)
png(file.path(analysis_out, paste0(integrate_sample_name, "_UMAP_feature_plot_breast_markers.png")), width=16, height=4*ceiling(length(markers2)/5.0), units='in', res=300)
p <- FeaturePlot(obj = RNA, features = markers2, ncol  = 5,min.cutoff = 0)
print(p)
dev.off()

markers2 <- c( "CSF1R", "CD68", "CD163", # MACROPHAGE
               "CD3D", "CD8A", "CD4",
               "CD79A",  "MS4A1","JCHAIN" # B cells
)
png(file.path(analysis_out, paste0(integrate_sample_name, "_UMAP_feature_plot_immune_markers.png")), width=16, height=4*ceiling(length(markers2)/3.0), units='in', res=300)
p <- FeaturePlot(obj = RNA, features = markers2, ncol  = 3,min.cutoff = 0)
print(p)
dev.off()

rnaAggr$assign.ident.breast <- RNA$assign.ident

### todo: add CIBERSORT marker based annotation ###
# and save annotation as csv: barcode + annotation
print("Perform Annotation using CIBERSORT markers")
# load("/mnt/cfce-rcsm/projects/yijiajiang_analysis/data/human.immune.CIBERSORT.RData")
# signatures <- human.immune.CIBERSORT
# celltypes <- as.character(unique(signatures[,1]))
# signature_list <- sapply(1:length(celltypes),function(x){
#   return(toupper(as.character(signatures[which(signatures[,1]==celltypes[x]),2])))})
# names(signature_list) <- celltypes
# # saveRDS(signature_list, "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/cibersort_marker_list.rds")

# data(human.immune.CIBERSORT)

RNA <- rnaAggr
genes = markers
min.score = 0.5
signature_list <- readRDS("/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/cibersort_marker_list.rds")

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

RNA@meta.data$assign.ident.cibersort = as.character(as.integer(RNA@meta.data$seurat_clusters)-1)
current.cluster.ids = as.character(as.integer(unique(RNA@meta.data$seurat_clusters))-1)
new.cluster.ids = cellscore_max_celltype

RNA@meta.data$assign.score.cibersort = cellscore_max[RNA@meta.data$assign.ident.cibersort]
RNA@meta.data$assign.ident.cibersort = plyr::mapvalues(x = RNA@meta.data$assign.ident.cibersort,
                                                       from = current.cluster.ids, to = new.cluster.ids)

head(RNA@meta.data)
pdf(file.path( analysis_out, paste0(integrate_sample_name, "_UMAP_by_celltype_cibersort_markerbased.pdf")))
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "assign.ident.cibersort", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

png(file.path( analysis_out, paste0(integrate_sample_name, "_UMAP_by_celltype_cibersort_markerbased.png")),  width=10, height=12, units='in', res=300)
p1 <- DimPlot(RNA, reduction = "umap",  group.by = "assign.ident.cibersort", label=TRUE) + ggtitle("scRNA Annotation")
print(p1)
dev.off()

saveRDS(RNA , file.path(output_dir))








