library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)
library(Signac)
library(tidyverse)
args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # sample_path <- "analysis/scATAC/integration/merge_1/merge_1_combine_integrated.rds"
anno1_dir <- args[2] # anno1_dir <- "analysis/scATAC/integration/merge_1/scATAnno_PBMC/query_annotated.csv"
metadata_dir <- args[3]
final_dir <- args[4]

sample <- head(tail(strsplit(sample_path, "/", fixed = TRUE)[[1]], n=2), n=1)
print(sample)

seurat <- readRDS(sample_path)

meta <- read.csv(anno1_dir, check.names=FALSE, row.names=1)
meta$barcodes <- rownames(meta)


for(i in unique(meta$dataset)){
meta$barcodes <- str_replace_all(meta$barcodes, paste0("_", i), "")
}

rownames(meta) <- paste0(meta$dataset, "_",meta$barcodes)

# meta <- meta %>% select(-c("celltypes", "tissue", "dataset"))

seurat <- AddMetaData(seurat, meta)

write.csv(seurat@meta.data, metadata_dir)
saveRDS(seurat, final_dir)

