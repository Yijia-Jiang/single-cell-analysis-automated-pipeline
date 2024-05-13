library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)
library(EnhancedVolcano)
library(stringr)
library(Matrix)

args <- commandArgs( trailingOnly = TRUE )
inputdata <- args[1] # ./merge_all_combine_integrated.rds # "home/path/merge_unHealthy_final.rds"
h5ad_mtx <- args[2]
h5ad_barcode <- args[3]
h5ad_feature <- args[4]
h5ad_metadata <- args[5]

analysis_out <- sapply(strsplit(h5ad_mtx, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

# dir.create(file.path(analysis_out), showWarnings = F)

rnaAggr <- readRDS(file.path(inputdata))
DefaultAssay(rnaAggr) <-"RNA"
rnaAggr <- JoinLayers(rnaAggr)
counts <- rnaAggr[["RNA"]]$counts

convert_mtx <- function(mat, h5ad_mtx, h5ad_feature, h5ad_barcode) {
  # dir.create(file.path("./", "merged_sample_final"), showWarnings = F)
  # model_dir <- file.path("./", "merged_sample_final")
  # writeMM(mat, file.path(analysis_out, paste0('matrix.mtx')))
  writeMM(mat, file.path(h5ad_mtx))
  write(rownames(mat), file.path(h5ad_feature))
  write(colnames(mat), file.path(h5ad_barcode))
  # write(rownames(mat), file.path(analysis_out, paste0( 'features.tsv')))
  # write(colnames(mat), file.path(analysis_out, paste0( 'barcodes.tsv')))
}

convert_mtx(counts, h5ad_mtx, h5ad_feature, h5ad_barcode)
write.csv(rnaAggr@meta.data, file.path(h5ad_metadata))
