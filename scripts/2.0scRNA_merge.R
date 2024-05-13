library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)

# data_folder <- "."

args <- commandArgs( trailingOnly = TRUE )

input <- args[1]
inputFiles = strsplit(input,' ')[[1]]
samples <- strsplit(args[2],',')[[1]]
output_data_dir <- args[3]
umap_cluster_dir <- args[4]
umap_sample_dir <- args[5]

inputfile_samples <- unlist(lapply(inputFiles, function(s) head(tail(strsplit(s, "/", fixed = TRUE)[[1]], n=2), n=1)))
names(inputFiles) <- inputfile_samples

inputFiles <- inputFiles[names(inputFiles) %in% samples]
print(inputFiles)

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

ifnb.list = c()
for (i in c(1:length(inputFiles))){
  
  rds_folder <- inputFiles[i]
  sample <- names(rds_folder)
  print(rds_folder)
  print(sample)
  
  subset <- readRDS(file.path(rds_folder))
  
  ifnb.list <- c(ifnb.list, subset)
}

# names(ifnb.list) <- samples
names(ifnb.list) <- names(inputFiles)

# Output merged r object
rnaAggr <- merge(ifnb.list[[1]], y = unlist(ifnb.list)[2:length(unlist(ifnb.list))], add.cell.ids =names(ifnb.list) )
unique((rnaAggr@meta.data$orig.ident))
head(rnaAggr@meta.data)

# rnaAggr <- SCTransform(rnaAggr, verbose = TRUE)
rnaAggr <- NormalizeData(rnaAggr)
rnaAggr <- FindVariableFeatures(rnaAggr)
rnaAggr <- ScaleData(rnaAggr)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:50, verbose = TRUE, reduction = "pca")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 1, reduction = "pca")
rnaAggr <- RunUMAP(rnaAggr, dims = 1:50, verbose = TRUE, reduction = "pca")

# saveRDS(rnaAggr, file.path(output_dir, "Object_Filtered_Merged.rds"))
saveRDS(rnaAggr, file.path(output_data_dir))

# rnaAggr <- readRDS(file.path(output_dir, "Object_Filtered_Merged.rds"))

# pdf(file.path(output_dir, "UMAP_by_cluster_before_integration.pdf"))
pdf(file.path(umap_cluster_dir))
p1 <- DimPlot(rnaAggr, reduction = "umap",  label = TRUE) + NoLegend()+ ggtitle("scRNA Seurat Clustering")
print(p1)
dev.off()

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_before_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(rnaAggr, reduction = "umap",  label = TRUE) + NoLegend()+ ggtitle("scRNA Seurat Clustering")
print(p1)
dev.off()

# pdf(file.path(output_dir, "UMAP_by_sample_before_integration.pdf"))
pdf(file.path(umap_sample_dir))
p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident") + ggtitle("scRNA Sample Integration")
print(p1)
dev.off()

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_before_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident") + ggtitle("scRNA Sample Integration")
print(p1)
dev.off()
