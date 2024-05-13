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


library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
# devtools::install_github("kassambara/ggpubr")
# library(ggpubr)
# source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")

# parent_dir <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/"
# output_dir <- "analysis/integration/merge_all/unified.peaks.rds"
# args <- commandArgs( trailingOnly = TRUE )

#inputfile = args[1] # final.rds
#output_data_dir = args[2]

inputfile <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/M2_atac_QC_processed.rds"
# plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
#                       function(i) paste(head(i, -1), collapse = "/"))
# integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

# TODO: 
## output should within scATAnno_atlas: like scATAnno_PBMC/diff_motif/plot.pdf
atac <- readRDS(inputfile)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

DefaultAssay(atac) <- "peaks"

metadata <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/scATAnno_PBMC/M2_atac_query_annotated.csv"
sample <- "M2_atac"

meta <- read.csv(metadata, row.names = 1)
meta$barcodes <- rownames(meta)
rownames(meta)

library(stringr)
rownames(meta) <- str_replace_all(rownames(meta), paste0("_", sample), "")
atac <- AddMetaData(atac, meta)

# Idents(atac) <- atac$seurat_clusters
Idents(atac) <- atac$cluster_annotation
da_peaks <- FindAllMarkers(
  object = atac,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

saveRDS(top.da.peak, "top.da.peak.rds")

# test enrichment
enriched.motifs <- FindMotifs(
  object = atac,
  features = top.da.peak
)

pdf("test_motif.pdf")
MotifPlot(
  object = atac,
  motifs = head(rownames(enriched.motifs))
)
dev.off()


# mouse_brain <- RunChromVAR(
#   object = mouse_brain,
#   genome = BSgenome.Mmusculus.UCSC.mm10
# )
# 
# DefaultAssay(mouse_brain) <- 'chromvar'
# 
# # look at the activity of Mef2c
# p2 <- FeaturePlot(
#   object = mouse_brain,
#   features = "MA0497.1",
#   min.cutoff = 'q10',
#   max.cutoff = 'q90',
#   pt.size = 0.1
# )
# p1 + p2
# 
# 
# differential.activity <- FindMarkers(
#   object = mouse_brain,
#   ident.1 = 'Pvalb',
#   ident.2 = 'Sst',
#   only.pos = TRUE,
#   mean.fxn = rowMeans,
#   fc.name = "avg_diff"
# )
# 
# MotifPlot(
#   object = mouse_brain,
#   motifs = head(rownames(differential.activity)),
#   assay = 'peaks'
# )

# saveRDS(count_list,  "./analysis/count_list.rds")
saveRDS(atac,  "test_atac.rds")

