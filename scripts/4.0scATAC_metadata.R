library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)
library(Signac)

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # sample_path <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/M2_atac_QC_processed.rds"
anno1_dir <- args[2] # anno1_dir <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/scATAnno_PBMC/M2_atac_query_annotated.csv"
metadata_dir <- args[3]
final_dir <- args[4]

sample <- head(tail(strsplit(sample_path, "/", fixed = TRUE)[[1]], n=2), n=1)
print(sample)

seurat <- readRDS(sample_path)

meta <- read.csv(anno1_dir, check.names=FALSE, row.names=1)
meta$barcodes <- rownames(meta)

for(i in unique(meta$dataset)){
  rownames(meta) <- str_replace_all(meta$barcodes, paste0("_", i), "")
}
seurat <- AddMetaData(seurat, meta)

write.csv(seurat@meta.data, metadata_dir)
saveRDS(seurat, final_dir)

