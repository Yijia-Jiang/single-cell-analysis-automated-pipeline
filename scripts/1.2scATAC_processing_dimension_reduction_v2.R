library(Signac)
library(Seurat)
library(ggplot2)
library(stringr)

library(Signac) #version 0.2.1
library(Seurat) #version 3.0.2
library(GenomeInfoDb)
# BiocManager::install("Bioconductor/GenomeInfoDb")
library(harmony) #version 1.0
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(here)
library(tibble)
library(dplyr)
library(gridExtra)
set.seed(1234)

# sessionInfo()
args <- commandArgs( trailingOnly = TRUE )

# Input CSV file; or path
filepath <- args[1] # QC.rds

output_dir <- args[2] # sample_QC_processed.rds
print(output_dir)
# get file path
# inputFiles <- list.files(file.path(filepath , list.files(filepath)), pattern = ".tsv.gz$", full.names = TRUE)
# names(inputFiles) <-  lapply(inputFiles, function(x) str_split(x,"/")[[1]][3] )    ## extract names of the samples
# print(inputFiles)

sample <- head(tail(strsplit(output_dir, "/", fixed = TRUE)[[1]], n=2), n=1)
plot_dir <- sapply(strsplit(output_dir, "/", fixed = TRUE),function(i) paste(head(i, -1), collapse = "/"))

# filepath <- '/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/M1/outs/filtered_peak_bc_matrix.h5'
# fragpath <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/M1/outs/fragments.tsv.gz"
# metapath <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/M1/outs/singlecell.csv"
# output_dir <- "analysis/single_sample/M1/M1_QC.rds"

pbmc <- readRDS(file.path(filepath))

# dir.create(file.path("analysis"), showWarnings = F)
# dir.create(file.path("analysis/single_sample"), showWarnings = F)
# dir.create(file.path(plot_dir), showWarnings = F)

# Dim reduction
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
gene.activities <- GeneActivity(pbmc)

pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
# pbmc <- NormalizeData(
#   object = pbmc,
#   assay = 'ACTIVITY',
#   normalization.method = 'LogNormalize',
#   scale.factor = median(pbmc$nCount_RNA)
# )

# newly added
DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc , features = rownames(pbmc ))


DefaultAssay(pbmc) <- "peaks"
pdf(file.path(plot_dir, paste0(sample, "_UMAP_by_cluster.pdf")))
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()

png(file.path(plot_dir, paste0(sample, "_UMAP_by_cluster.png")), width = 12, height = 4, units='in', res=300)
##Fragment length periodicity:
DimPlot(object = pbmc, label = TRUE) + NoLegend()
dev.off()


saveRDS(pbmc, file = file.path(output_dir))






