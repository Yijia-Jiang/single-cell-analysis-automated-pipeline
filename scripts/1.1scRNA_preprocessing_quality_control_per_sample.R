suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
# library(tidyverse)
suppressMessages(library(stringi))
suppressMessages(library(stringr))
library(gridExtra)

# path <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/cellranger/Ellisen_BCA/per_sample_outs/BI54/count/sample_filtered_feature_bc_matrix.h5"
# sample <- "BI54"

# print(sample)

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # 
output_dir <- args[2] # 
nCount_RNA_min = as.numeric(args[3]) 
nCount_RNA_max = as.numeric(args[4])
mito_rate_max = as.numeric(args[5])


print(paste0("nCount_RNA_min: ", nCount_RNA_min))
print(paste0("mitochondria qc: " ,mito_rate_max))
print(paste0("nCount_RNA_max qc: " ,nCount_RNA_max))

# sample_path <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/cellranger/Ellisen_BCA/per_sample_outs/BI54/count/sample_filtered_feature_bc_matrix.h5"
# output_dir <- "analysis/single_sample/BI54/BI54_QC.rds"
# nCount_RNA_min = 200
# nCount_RNA_max = 8000
# mito_rate_max = 5

plot_dir <- sapply(strsplit(output_dir, "/", fixed = TRUE),function(i) paste(head(i, -1), collapse = "/"))
sample <- head(tail(strsplit(output_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

data_folder <- file.path(sample_path)
data <- Read10X_h5(data_folder)


# If multiome structure, read data$gene expression
if(length(data)==2){
  data <- data$`Gene Expression`
} else{
  # else, read data directly
  data <- data
}

seurat = CreateSeuratObject(counts = data, project = sample)
ncell <- ncol(seurat)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
print(seurat)

pdf(file.path(plot_dir, paste0(sample, "_QC_plot_before.pdf")))
p <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()
print("before QC")
print(seurat)

png(file.path(plot_dir, paste0(sample, "_QC_plot_before_original.png")), width = 6, height = 6, units='in', res=300)
p <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()
print("before QC")
print(seurat)

# seurat <- readRDS("analysis/single_sample/BI54/BI54_QC.rds")
# seurat$assay <- "RNA"
seurat <- JoinLayers(seurat)
Idents(seurat) <- seurat$orig.ident
# pdf(file.path( plot_dir, paste0(sample, "_QC_plot_after_test.pdf")))
# p <- VlnPlot(seurat, layer="counts", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# print(p)
# dev.off()

r1<-VlnPlot(seurat, features = c("nFeature_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r2<-VlnPlot(seurat, features = c("nCount_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r3<-VlnPlot(seurat, features = c("percent.mt"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
png(file.path(plot_dir, paste0(sample, "_QC_plot_before_adjust.png")),width = 3, height = 3, units='in', res=300)
p <- grid.arrange(r1,r2,r3,ncol = 4)
print(p)
dev.off()

rnaAggr <- subset(seurat, subset = nCount_RNA > nCount_RNA_min  & nCount_RNA < nCount_RNA_max & percent.mt < mito_rate_max)
# rnaAggr <- subset(seurat, subset = nCount_RNA > 200  & nCount_RNA < 8000 & percent.mt < 5)
# rnaAggr <- SCTransform(rnaAggr, verbose = TRUE)
rnaAggr <- NormalizeData(rnaAggr)
rnaAggr <- FindVariableFeatures(rnaAggr)
rnaAggr <- ScaleData(rnaAggr)

print("after QC")
print(rnaAggr )
pdf(file.path( plot_dir, paste0(sample, "_QC_plot_after.pdf")))
p <- VlnPlot(rnaAggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()

png(file.path(plot_dir, paste0(sample, "_QC_plot_after_original.png")), width = 6, height = 6, units='in', res=300)
p <- VlnPlot(rnaAggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(p)
dev.off()
print("before QC")


png(file.path(plot_dir, paste0(sample, "_QC_plot_after_adjusted.png")),width = 5, height = 5, units='in', res=300)
r1<-VlnPlot(rnaAggr, features = c("nFeature_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r2<-VlnPlot(rnaAggr, features = c("nCount_RNA"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r3<-VlnPlot(rnaAggr, features = c("percent.mt"),ncol=3,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))

p <- grid.arrange(r1,r2,r3,ncol = 4)
print(p)
dev.off()


rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)

# pdf(file.path(output_dir, paste0("PCA_ElbowPlot.pdf")))
# ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
# dev.off()

rnaAggr <- FindNeighbors(rnaAggr, dims = 1:50, verbose = TRUE, reduction = "pca")
rnaAggr <- FindClusters(rnaAggr, verbose = TRUE, resolution = 1, reduction = "pca")
rnaAggr <- RunUMAP(rnaAggr, dims = 1:50, verbose = TRUE, reduction = "pca")


pdf(file.path(plot_dir, paste0(sample, "_UMAP_by_cluster_number.pdf")))
p1 <- DimPlot(rnaAggr, reduction = "umap",  label = TRUE) + NoLegend()+ ggtitle("scRNA Seurat Clustering")
print(p1)
dev.off()

png(file.path(plot_dir, paste0(sample, "_UMAP_by_cluster_number_original.png")), width = 6, height = 6, units='in', res=300)
p1 <- DimPlot(rnaAggr, reduction = "umap",  label = TRUE) + NoLegend()+ ggtitle("scRNA Seurat Clustering")
print(p1)
dev.off()

png(file.path(plot_dir, paste0(sample, "_UMAP_by_cluster_number_adjust.png")), width = 3, height = 3, units='in', res=300)
p1 <- DimPlot(rnaAggr, reduction = "umap",label = TRUE,pt.size =0.001)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6),axis.title.y = element_text(color = "grey20", size = 4))
print(p1)
dev.off()


saveRDS(rnaAggr, file.path(output_dir))

