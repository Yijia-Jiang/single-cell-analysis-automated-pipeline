# Integration

library(Seurat) # 3.02
library(tidyverse)
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
# library(DoubletFinder)
# library(openxlsx)
# library(here)
library(dplyr)
library(scDblFinder)
library(stringr)
library(tibble)
# data_folders <- list.files("../cellranger/two_sequence_library_per_sample_outs")
# sample_path <- "./sample/MGH23039R/MGH23039R_Object_Filtered.rds"

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1]
output_dir <- args[2]
output_df <- args[3]

plot_dir <- sapply(strsplit(output_dir, "/", fixed = TRUE),
                   function(i) paste(head(i, -1), collapse = "/"))
sample <- head(tail(strsplit(output_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

rnaAggr <- readRDS(sample_path)

object <- rnaAggr
library(SingleCellExperiment)
Idents(object) <- object$seurat_clusters
object[["RNA"]] <- as(object[["RNA"]], Class="Assay")
sce <- as.SingleCellExperiment(object)

set.seed(10010101)
# sce.dbl <- scDblFinder(sce, clusters=colLabels(sce))
sce.dbl <- scDblFinder(sce, clusters=sce$seurat_clusters)
df <- colData(sce.dbl)
df_new <- df[c("scDblFinder.class", "scDblFinder.cluster", "scDblFinder.score")]
df_new <- as.data.frame(df_new)

rnaAggr <- AddMetaData(rnaAggr, df_new)

pdf(file.path(plot_dir, paste0(sample, "_UMAP_doublet.pdf")))
p2<-DimPlot(rnaAggr, group.by="scDblFinder.class")
print(p2)
dev.off()

png(file.path(plot_dir, paste0(sample, "_UMAP_doublet.png")), width = 6, height = 6, units='in', res=300)
p2<-DimPlot(rnaAggr, group.by="scDblFinder.class")
print(p2)
dev.off()

write.csv(df_new, output_df)

# Filter out singlet cells 
rnaAggr_orig <- readRDS(sample_path)
rnaAggr_orig <- AddMetaData(rnaAggr_orig, df_new)
rnaAggr_filt <- subset(rnaAggr_orig, subset=scDblFinder.class=="singlet")
saveRDS(rnaAggr_filt, output_dir)

