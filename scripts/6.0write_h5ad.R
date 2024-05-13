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
outputdata <- args[2] # ./merge_all_combine_integrated/merge_all_combine_integrated.mtx look for mtx # home/path/merge_unHealthy_final/mtx

inputdata  <- "./merge_all_combine_integrated.rds"
outputdata <- "./merge_all_combine_integrated/merge_all_combine_integrated.mtx"

analysis_out <- sapply(strsplit(outputdata, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

dir.create(file.path(analysis_out), showWarnings = F)

rnaAggr <- readRDS(file.path(inputdata))
DefaultAssay(rnaAggr) <-"RNA"
rnaAggr <- JoinLayers(rnaAggr)
counts <- rnaAggr[["RNA"]]$counts

convert_mtx <- function(mat, analysis_out) {
  # dir.create(file.path("./", "merged_sample_final"), showWarnings = F)
  # model_dir <- file.path("./", "merged_sample_final")
  # writeMM(mat, file.path(analysis_out, paste0('matrix.mtx')))
  writeMM(mat, file.path(outputdata))
  write(rownames(mat), file.path(analysis_out, paste0( 'features.tsv')))
  write(colnames(mat), file.path(analysis_out, paste0( 'barcodes.tsv')))
}

convert_mtx(counts, analysis_out)

write.csv(rnaAggr@meta.data, file.path(analysis_out, paste0( 'metadata_final.csv')))
