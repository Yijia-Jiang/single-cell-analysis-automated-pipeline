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
args <- commandArgs( trailingOnly = TRUE )

inputfile = args[1] # final.rds
# metadata = args[2]
# sample=args[3]
da_peaks_dir = args[2]
motif_dir = args[3]
output_data_dir = args[4]

# inputfile <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/M2_atac_final.rds"
# metadata <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/scATAnno_PBMC/M2_atac_query_annotated.csv"
# sample <- "M2_atac"

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
sample <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

# TODO: 
## output should within scATAnno_atlas: like scATAnno_PBMC/diff_motif/plot.pdf
atac <- readRDS(inputfile)

# # Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# # saveRDS(pfm, "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/scatac_motif_pfm.rds")

# pfm <- readRDS( "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/scatac_motif_pfm.rds")

# seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "NCBI"

# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

DefaultAssay(atac) <- "peaks"

meta <- read.csv(metadata, row.names = 1)
meta$barcodes <- rownames(meta)
rownames(meta)

library(stringr)
rownames(meta) <- str_replace_all(rownames(meta), paste0("_", sample), "")
atac <- AddMetaData(atac, meta)

# Idents(atac) <- atac$seurat_clusters
Idents(atac) <- atac$cluster_annotation

# da_peaks <- FindMarkers(
#   object = atac,
#   ident.1="Memory CD4 T",
#   ident.2="Central memory CD8 T",
#   only.pos = TRUE,
#   test.use = 'LR',
#   min.pct = 0.05,
#   latent.vars = 'nCount_peaks'
# )

atac <- subset(atac, cluster_annotation!= "unknown")


# atac <- subset(atac, cluster_annotation%in%c("Memory CD4 T" ,"Central memory CD8 T", "pDC" ))

da_peaks_all <- FindAllMarkers(
  object = atac,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

write.csv(da_peaks_all, file.path(plot_folder, "da_peaks_all.csv"))
saveRDS(da_peaks_all, da_peaks_dir)

filtered_da_peaks <- da_peaks_all[grepl("^chr", row.names(da_peaks_all)), ]

plot_folder <- "test_motif"
dir.create(plot_folder)
for(celltype in unique(atac$cluster_annotation)){
  print(celltype)
# get top differentially accessible peaks
  tmp <- filtered_da_peaks[filtered_da_peaks$cluster == celltype, ]
  print(head(tmp))
  top.da.peak <- rownames(tmp[tmp$p_val < 0.005 & tmp$pct.1 > 0.2, ])
# top.da.peak <- rownames(filtered_da_peaks[filtered_da_peaks$p_val < 0.005 & filtered_da_peaks$pct.1 > 0.2, ])
# saveRDS(top.da.peak, "top.da.peak.rds")

# test enrichment
enriched.motifs <- FindMotifs(
  object = atac,
  features = top.da.peak
)

write.csv(enriched.motifs, file.path(plot_folder, paste0("enriched_motifs_", celltype, ".csv")))

pdf(file.path(plot_folder, paste0("enriched_motifs_", celltype, ".pdf")))
# pdf("test_motif.pdf")
MotifPlot(
  object = atac,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

}

saveRDS(enriched.motifs, motif_dir)


pdf(output_data_dir)
# pdf("test_motif.pdf")
MotifPlot(
  object = atac,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

# library(chromVAR)
# 
# atac <- RunChromVAR(
#   object = atac,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
# )
# 
# DefaultAssay(atac) <- 'chromvar'
# 
# # # look at the activity of Mef2c
# # p2 <- FeaturePlot(
# #   object = mouse_brain,
# #   features = "MA0497.1",
# #   min.cutoff = 'q10',
# #   max.cutoff = 'q90',
# #   pt.size = 0.1
# # )
# # p1 + p2
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
# saveRDS(atac,  "test_atac.rds")

