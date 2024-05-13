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
sample=args[2]
da_peaks_dir = args[3]
# motif_dir = args[3]
# output_data_dir = args[4]

# inputfile <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/M2_atac_final.rds"
# metadata <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample/M2_atac/scATAnno_PBMC/M2_atac_query_annotated.csv"
# sample <- "M2_atac"

plot_folder <- sapply(strsplit(da_peaks_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
# sample <- head(tail(strsplit(da_peaks_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

# TODO: 
## output should within scATAnno_atlas: like scATAnno_PBMC/motif/plot.pdf
atac <- readRDS(inputfile)

# # Get a list of motif position frequency matrices from the JASPAR database
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

Idents(atac) <- atac$cluster_annotation

if(length(unique(atac$cluster_annotation)) == 1) {
  print("There is only one celltype, unable to perform TF analysis")
  message <- "There is only one celltype, unable to perform TF analysis"
  saveRDS(message, da_peaks_dir)
  
} else {

atac <- subset(atac, cluster_annotation!= "unknown")

# if all celltypes are unknown, then stop executing

# try:
# textPlot <- function(plotname, string){
#   par(mar=c(0,0,0,0))
#   pdf(paste0(plotname, ".pdf"))
#   plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
#   text(x = 0.5, y = 0.5, paste(string), cex = 4, col = "black", family="serif", font=2, adj=0.5)
#   dev.off()
# }
# 
# textPlot("test", "no motif")


# tryCatch(da_peaks_all <- FindAllMarkers(
#   object = atac,
#   only.pos = TRUE,
#   test.use = 'LR',
#   min.pct = 0.05,
#   latent.vars = 'nCount_peaks'
# ),
#          error = function(e) {
#            print(e)
#            print("There are no differential motifs")
#            stop()
#          })


da_peaks_all <- FindAllMarkers(
  object = atac,
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

print(head(da_peaks_all))

write.csv(da_peaks_all, file.path(plot_folder, "da_peaks_all.csv"))


filtered_da_peaks <- da_peaks_all[grepl("^chr", row.names(da_peaks_all)), ]

if(nrow(filtered_da_peaks) <1 ){
  print("No differential motifs found")
  msg <- "No differential motifs found"
  saveRDS(msg, da_peaks_dir)
  
} else {
  saveRDS(da_peaks_all, da_peaks_dir)
# plot_folder <- "test_motif"
# dir.create(plot_folder)
for(celltype in unique(atac$cluster_annotation)){
  print(celltype)
# get top differentially accessible peaks
  tmp <- filtered_da_peaks[filtered_da_peaks$cluster == celltype, ]
  print(head(tmp))
  top.da.peak <- rownames(tmp[tmp$p_val < 0.005 & tmp$pct.1 > 0.2, ])

  # test enrichment
  enriched.motifs <- FindMotifs(
    object = atac,
    features = top.da.peak
  )
  
  write.csv(enriched.motifs, file.path(plot_folder, paste0("enriched_motifs_", celltype, ".csv")))
  
  pdf(file.path(plot_folder, paste0("enriched_motifs_", celltype, ".pdf")))
  # pdf("test_motif.pdf")
  p <- MotifPlot(
    object = atac,
    motifs = head(rownames(enriched.motifs))
  )
  print(p)
  dev.off()
  
  
  png(file.path(plot_folder, paste0("enriched_motifs_", celltype, ".png")), width = 6, height = 6, units='in', res=300)
  # pdf("test_motif.pdf")
  p <- MotifPlot(
    object = atac,
    motifs = head(rownames(enriched.motifs))
  )
  print(p)
  dev.off()

}

}
}

