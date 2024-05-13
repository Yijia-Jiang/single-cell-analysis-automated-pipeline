library(Signac)
library(Seurat)
library(SeuratWrappers)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86) #---GRCh38 (hg38)
# library(ggploT1)
library(patchwork)
set.seed(101)
library(future)
library(tidyverse)
# plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) #50G
library(GenomicRanges)
# devtools::install_github("kassambara/ggpubr")
# library(ggpubr)
# source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")

# parent_dir <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/"
# output_dir <- "analysis/integration/merge_all/unified.peaks.rds"

args <- commandArgs( trailingOnly = TRUE )
input_file <- args[1] #input_file <- readRDS("analysis/scATAC/integration/merge_1/merge_1_combine.rds")
annotFile <- args[2]
output_data_dir <- args[3]
umap1_dir <- args[4]

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)


tmp_ann <- read.csv(annotFile, sep=",", header=T, row.names=1, 
                    stringsAsFactors=FALSE, check.names=F)

# print(tmp_ann)
## Remove all cols that have compare information
#if(any(grepl("comp_", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                 ("comp_"))) }

#if(any(grepl("comp_", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                 ("comp_"))) }

# ## Remove all cols that have fastq info
# if(any(grepl("Treat_fastq", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                        ("Treat_fastq"))) }
# 
# if(any(grepl("Ctrl_fastq", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                       ("Ctrl_fastq"))) }

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
  
  # gene.activities <- GeneActivity(integrated)
  # DefaultAssay(integrated) <- "ACTIVITY"
  # integrated <- NormalizeData(integrated)
  # integrated <- ScaleData(integrated , features = rownames(integrated ))
  
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.pdf")))
  p1 <- DimPlot(atacAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(atacAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.pdf")))
  p1 <- DimPlot(atacAggr, reduction = "umap",  group.by = "dataset")+ggtitle("scATAC Sample Integration")
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(atacAggr, reduction = "umap",  group.by = "dataset") +ggtitle("scATAC Sample Integration")
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
  
  # gene.activities <- GeneActivity(integrated)
  # DefaultAssay(integrated) <- "ACTIVITY"
  # integrated <- NormalizeData(integrated)
  # integrated <- ScaleData(integrated , features = rownames(integrated ))
  
  # create a new UMAP using the integrated embeddings
  integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
  
  pdf(umap1_dir)
  p2 <- DimPlot(integrated, group.by = "dataset")
  print(p2)
  dev.off()
  
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.pdf")))
  p1 <- DimPlot(integrated, reduction = "umap",  group.by = "dataset")+ggtitle("scATAC Sample Integration")
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")))
  p2 <- DimPlot(integrated, group.by = "dataset")+ggtitle("scATAC Sample Integration")
  print(p2)
  dev.off()
  
  saveRDS(integrated, file = output_data_dir) 
  
}










