library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)


args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # sample_path <- "sample/BA-1_209_NormR/BA-1_209_NormR_Object_Annotated.rds"
anno1_dir <- args[2] # anno1_dir <- "sample/BA-1_209_NormR/BA-1_209_NormR_SingleR_hpca_main.csv"
anno2_dir <- args[3] # /mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/analysis/sample/BI54/BI54_SingleR_encode_main.csv
metadata_dir <- args[4]
final_dir <- args[5]

seurat <- readRDS(sample_path)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# # view cell cycle scores and phase assignments
# pdf(ridgeplot_dir)
# RidgePlot(seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
# dev.off()
# 
# pdf(pcaplot_dir)
# seurat <- RunPCA(seurat, features = c(s.genes, g2m.genes))
# DimPlot(seurat)
# dev.off()

anno1 <- read.csv(anno1_dir, row.names = 1)
anno2 <- read.csv(anno2_dir, row.names = 1)

seurat <- AddMetaData(seurat, anno1)
seurat <- AddMetaData(seurat, anno2)

write.csv(seurat@meta.data, metadata_dir)
saveRDS(seurat, final_dir)

