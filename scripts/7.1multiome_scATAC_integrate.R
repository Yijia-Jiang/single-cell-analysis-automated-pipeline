library(Signac)
library(Seurat)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
library(patchwork)
set.seed(101)
library(future)
library(tidyverse)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)

args <- commandArgs( trailingOnly = TRUE )
input_file <- args[1] #
annotFile <- args[2] # 
output_data_dir <- args[3]
umap1_dir <- args[4]

# input_file <- "/home/yj976/scpipeline/yj_version4/analysis/scMultiome/scATAC/merge_1/merge_1_combine.rds"
# annotFile <- "metasheet.csv"
# output_data_dir <- "analysis/scMultiome/scATAC/merge_1/merge_1_combine_integrated.rds"

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)


tmp_ann <- read.csv(annotFile, sep=",", header=T, row.names=1, 
                    stringsAsFactors=FALSE, check.names=F)

rowNames <- rownames(tmp_ann)
colNames <- colnames(tmp_ann)

atacAggr <- readRDS(file.path(input_file)) # 
samples <- intersect(unique(atacAggr$dataset), rownames(tmp_ann))
# print(samples)
cat("samples finally used: ", samples, "\n")

#select samples common for metasheet and tag-count file
tmp_ann <- as.data.frame(tmp_ann[samples,])
colnames(tmp_ann) <- colNames
rownames(tmp_ann) <- samples

if(length(unique(tmp_ann[,"batch"])) == 1){
  print("All samples belong to the same batch; No Need to performing harmony")
  atacAggr <- FindNeighbors(atacAggr, reduction = 'lsi', dims = 2:30)
  atacAggr <- FindClusters(atacAggr, resolution = 0.4)
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.pdf")))
  p1 <- DimPlot(atacAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(atacAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.pdf")))
  p1 <- DimPlot(atacAggr, reduction = "umap",  group.by = "orig.ident")+ggtitle("scATAC Sample Integration")
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(atacAggr, reduction = "umap",  group.by = "orig.ident") +ggtitle("scATAC Sample Integration")
  print(p1)
  dev.off()
  
  saveRDS(atacAggr, output_data_dir)
} else{
  new_names <- as.list(tmp_ann[rownames(tmp_ann), "batch"])
  print(new_names)
  names(new_names) <- rownames(tmp_ann)
  Idents(object = atacAggr) <- atacAggr$dataset
  atacAggr <- RenameIdents(object = atacAggr, new_names)
  atacAggr$batch <- Idents(atacAggr)
  
  # ifnb <- SplitObject(atacAggr, split.by = "dataset")
  ifnb <- SplitObject(atacAggr, split.by = "batch")
  
  # find integration anchors
  integration.anchors <- FindIntegrationAnchors(
    object.list = ifnb,
    reduction = "rlsi",
    dims = 2:30
  )
  
  # integrate LSI embeddings
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = atacAggr[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:30
  )
  
  # create a new UMAP using the integrated embeddings
  integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
  
  # # newly added
  gene.activities <- GeneActivity(integrated)
  integrated[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
  DefaultAssay(integrated) <- "ACTIVITY"
  integrated <- NormalizeData(integrated)
  integrated <- ScaleData(integrated , features = rownames(integrated ))
  
  # integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
  # integrated <- NormalizeData(
  #   object = integrated,
  #   assay = 'RNA',
  #   normalization.method = 'LogNormalize',
  #   scale.factor = median(integratedc$nCount_RNA)
  # )
  
  pdf(umap1_dir)
  p2 <- DimPlot(integrated, group.by = "dataset")
  print(p2)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")))
  p2 <- DimPlot(integrated, group.by = "dataset")
  print(p2)
  dev.off()
  
  saveRDS(integrated, file = output_data_dir) 
  
}










