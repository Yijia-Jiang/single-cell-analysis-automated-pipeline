
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stringi)
library(stringr)

args <- commandArgs( trailingOnly = TRUE )
input_file <- args[1]
annotFile <- args[2]
output_data_dir <- args[3]

# analysis_out <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

# input_file <- "Merge_exp2/Merge_exp2.rds"
# annotFile <- "/home/yj976/scpipeline/cobra_simple/metasheet.csv"
# output_data_dir <- "Merge_exp2/Merge_exp2_harmony_integrated.rds"
# analysis_out <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

tmp_ann <- read.csv(annotFile, sep=",", header=T, row.names=1, 
                    stringsAsFactors=FALSE, check.names=F)

## Remove all cols that have compare information
if(any(grepl("comp_", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
                                                                                 ("comp_"))) }

if(any(grepl("comp_", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
                                                                                 ("comp_"))) }

# ## Remove all cols that have fastq info
# if(any(grepl("Treat_fastq", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                        ("Treat_fastq"))) }
# 
# if(any(grepl("Ctrl_fastq", colnames(tmp_ann)))) { tmp_ann <- dplyr::select(tmp_ann, -(starts_with
#                                                                                       ("Ctrl_fastq"))) }

rowNames <- rownames(tmp_ann)
colNames <- colnames(tmp_ann)

rnaAggr <- readRDS(file.path(input_file))
samples <- intersect(unique(rnaAggr$orig.ident), rownames(tmp_ann))

cat("samples finally used: ", samples, "\n")

#select samples common for metasheet and tag-count file
tmp_ann <- as.data.frame(tmp_ann[samples,])
colnames(tmp_ann) <- colNames
rownames(tmp_ann) <- samples

if(length(unique(tmp_ann[,"batch"])) == 1){
  print("All samples belong to the same batch; No Need to performing harmony")
  rnaAggr <- FindClusters(rnaAggr, resolution = 0.4)
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.pdf")))
  p1 <- DimPlot(rnaAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(rnaAggr, group.by = "seurat_clusters" ) 
  print(p1)
  dev.off()
  
  pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.pdf")))
  p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident")+ggtitle("scRNA Sample Integration")
  print(p1)
  dev.off()
  
  png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")), width=6, height=6, units='in', res=300)
  p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident") +ggtitle("scRNA Sample Integration")
  print(p1)
  dev.off()
  saveRDS(rnaAggr, output_data_dir)
  
} else{

new_names <- as.list(tmp_ann[rownames(tmp_ann), "batch"])
names(new_names) <- rownames(tmp_ann)
Idents(object = rnaAggr) <- rnaAggr$orig.ident
rnaAggr <- RenameIdents(object = rnaAggr, new_names)
rnaAggr$batch <- Idents(rnaAggr)

# Remove batch effects using harmony
rnaAggr <- IntegrateLayers(object = rnaAggr, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony",
                           verbose = FALSE)
rnaAggr <- FindNeighbors(rnaAggr, reduction = "harmony", dims = 1:30)
rnaAggr <- RunUMAP(rnaAggr, reduction = "harmony", dims = 1:30, reduction.name = "umap")
rnaAggr <- FindClusters(rnaAggr, resolution = 0.4)

saveRDS(rnaAggr, output_data_dir)

rnaAggr <- FindClusters(rnaAggr, resolution = 0.4)

pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.pdf")))
p1 <- DimPlot(rnaAggr, group.by = "seurat_clusters" ) 
print(p1)
dev.off()

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_after_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(rnaAggr, group.by = "seurat_clusters" ) 
print(p1)
dev.off()

pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.pdf")))
p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident")+ggtitle("scRNA Sample Integration")
print(p1)
dev.off()

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_after_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(rnaAggr, reduction = "umap",  group.by = "orig.ident") +ggtitle("scRNA Sample Integration")
print(p1)
dev.off()
}




