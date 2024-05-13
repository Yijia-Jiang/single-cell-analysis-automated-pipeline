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
filepath <- args[1] # 
fragpath <- args[2]
metapath <- args[3]
# ref <- args[4] # hg38
output_dir <- args[4] # QC.rds
df_output_dir <- args[5]
df_output_dir_barcodes <- args[6]
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

counts <- Read10X_h5(filename = file.path(filepath))
metadata <- read.csv(
  file = metapath,
  header = TRUE,
  row.names = 1)


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = fragpath)

seurat <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata,
  project=sample
)

# if (ref=="hg38"){
# extract gene annotations from EnsDb
annotations <- readRDS("/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/annotations.rds")
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(seurat) <- annotations
# } 

# Quality Control
# compute nucleosome signal score per cell
combined <- seurat
combined <- NucleosomeSignal(object = combined)
# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE, assay = 'peaks', verbose = T)
# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

# dir.create(file.path("analysis"), showWarnings = F)
# dir.create(file.path("analysis/single_sample"), showWarnings = F)
# dir.create(file.path(plot_dir), showWarnings = F)


##Fragment length periodicity:
# FragmentHistogram(object = combined, group.by = 'nucleosome_group')
print(file.path(plot_dir, paste0(sample, "_QC_plot_before.pdf")))
pdf(file.path(plot_dir, paste0(sample, "_QC_plot_before.pdf")))
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 4
)
dev.off()

png(file.path(plot_dir, paste0(sample, "_QC_plot_before_original.png")), width = 12, height = 4, units='in', res=300)
##Fragment length periodicity:
VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 4
)
dev.off()

png(file.path(plot_dir, paste0(sample, "_QC_plot_before_adjust.png")), width = 10, height = 5, units='in', res=300)
r1<-VlnPlot(combined, features = c('pct_reads_in_peaks'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r2<-VlnPlot(combined, features = c('peak_region_fragments'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r3<-VlnPlot(combined, features = c('TSS.enrichment'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r4<-VlnPlot(combined, features = c('nucleosome_signal'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
p <- grid.arrange(r1,r2,r3,r4,ncol = 4)
print(p)
dev.off()


## remove outlier cells based on QC metrics:
combined.pro <- subset(
  x = combined,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

# todo: output cell barcode list

# todo: plot scatterplot

pdf(file.path(plot_dir, paste0(sample, "_QC_plot_after.pdf")))
VlnPlot(
  object = combined.pro,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 5
)
dev.off()

png(file.path(plot_dir, paste0(sample, "_QC_plot_after_original.png")), width = 6, height = 6, units='in', res=300)
##Fragment length periodicity:
VlnPlot(
  object = combined.pro,
  features = c('pct_reads_in_peaks', 'peak_region_fragments','TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 4
)
dev.off()

png(file.path(plot_dir, paste0(sample, "_QC_plot_after_adjust.png")), width = 10, height = 5, units='in', res=300)
r1<-VlnPlot(combined.pro, features = c('pct_reads_in_peaks'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r2<-VlnPlot(combined.pro, features = c('peak_region_fragments'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r3<-VlnPlot(combined.pro, features = c('TSS.enrichment'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
r4<-VlnPlot(combined.pro, features = c('nucleosome_signal'),ncol=4,pt.size = 0)+ theme(legend.position = "none", axis.text.y = element_text(color = "grey20", size = 4),axis.text.x = element_text(color = "grey20", size = 4),axis.title.x = element_text(color = "grey20", size = 4),title =element_text(size=6))
p <- grid.arrange(r1,r2,r3,r4,ncol = 4)
print(p)
dev.off()
# todo: output cell barcode list
# write.csv(combined.pro@meta.data, file.path(plot_dir, paste0(sample, "_cell_filtered.csv")))
write.csv(combined.pro@meta.data, file.path(df_output_dir))

# save barcodes only
# data <- readRDS("/mnt/cfce-rcsm/projects/development2/scpipe/yj_version4/analysis/scATAC/single_sample/M1_atac/M1_atac_QC.rds")
barcode_df <- rownames(combined.pro@meta.data)
write.table(barcode_df, file.path(df_output_dir_barcodes), row.names=FALSE, col.names=FALSE,  quote = F)

# process orig.ident information
# combined.pro@meta.data$orig.ident <- gsub("_.*", "", rownames(combined.pro@meta.data))
saveRDS(combined.pro, file = file.path(output_dir))




