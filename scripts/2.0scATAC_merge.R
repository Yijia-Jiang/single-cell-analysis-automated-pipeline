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
# devtools::install_github("kassambara/ggpubr")
# library(ggpubr)
# source(file = "/home/longzhilin/Analysis_Code/Plot_colorPaletters.R")

# parent_dir <- "/mnt/cfce-rcsm/projects/yijiajiang_analysis/240226_X202SC24012007-Z01-F002_Sonsoles_scATAC_human/cellranger/"
# output_dir <- "analysis/integration/merge_all/unified.peaks.rds"
args <- commandArgs( trailingOnly = TRUE )
print(args[1])
bed_inputFiles = strsplit(args[1],' ')[[1]]
print(bed_inputFiles)
meta_inputFiles = strsplit(args[2],' ')[[1]]
frag_inputFiles= strsplit(args[3],' ')[[1]]
# samplenames <- strsplit(args[4],' ')[[1]]
# merge_samples <- strsplit(args[5],',')[[1]]
# output_data_dir <- args[6] # "merge_sample_combine.rds"
# combined_peaks_dir <- args[7] # combined_peaks.rds
# concat_robject_dir <- args[8] # count_list.rds
# umap1_dir <- args[9]
# samples <- strsplit(args[4],',')[[1]]
samplenames <- strsplit(args[4],',')[[1]]
print(samplenames)
output_data_dir <- args[5] # "merge_sample_combine.rds"
combined_peaks_dir <- args[6] # combined_peaks.rds
concat_robject_dir <- args[7] # count_list.rds
umap1_dir <- args[8]
# print(umap1_dir)
# umap2_dir <- args[9]
# bed_inputFiles <- c(file.path(parent_dir,"M1", "/outs/peaks.bed"),
#                     file.path(parent_dir,"B1", "/outs/peaks.bed")
#                 # file.path(parent_dir,"M5", "/outs/peaks.bed")
#                 )
# 
# meta_inputFiles <- c(file.path("analysis/single_sample/M1/", "M1_cell_filtered.csv"),
#                      file.path("analysis/single_sample/B1/", "B1_cell_filtered.csv")
#                      # file.path("analysis/single_sample/M5/", "M5_cell_filtered.csv")
# )
# 
# frag_inputFiles <- c(file.path(parent_dir,"M1", "/outs/fragments.tsv.gz"),
#                      # file.path(parent_dir,"M5", "/outs/fragments.tsv.gz")
#                      file.path(parent_dir, "B1", "/outs/fragments.tsv.gz")
# )
# 
# # samplenames <- c("M1", "M5")
# samplenames <- c("M1", "B1")

# inputfile_samples <- unlist(lapply(inputFiles, function(s) head(tail(strsplit(s, "/", fixed = TRUE)[[1]], n=2), n=1)))
# names(inputFiles) <- inputfile_samples
# 
# inputFiles <- inputFiles[names(inputFiles) %in% samples]
# print(inputFiles)

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

### Creating a common peak set by merging:
gr_peak_list <- c()
for(bed_file in bed_inputFiles){
  print(bed_file)
peaks <- read.table(
  file = file.path(bed_file),
  col.names = c("chr", "start", "end")
)
print(head(peaks))
# # # convert to genomic ranges
gr.peaks <- makeGRangesFromDataFrame(peaks)
print(gr.peaks)
# 
gr_peak_list <- c(gr_peak_list, gr.peaks)
}

combined.peaks <- UnifyPeaks(object.list = unlist(gr_peak_list), mode = "reduce")
combined.peaks


# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20] 
print(combined_peaks_dir)

saveRDS(combined.peaks, combined_peaks_dir)
print("Done with obtaining combined peaks")

count_list <- c()
# # load metadata
for(i in seq(length(bed_inputFiles))){
  print(i)
  print(samplenames[i])
  md.T1 <- read.table(
    file = file.path(meta_inputFiles[i]),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ]
  print(meta_inputFiles[i])
  frags.T1 <- CreateFragmentObject(
      path = file.path(frag_inputFiles[i]),
      cells = rownames(md.T1)
    )
  print(frag_inputFiles[i])
  T1.counts <- FeatureMatrix(
    fragments = frags.T1,
    features = combined.peaks,
    cells = rownames(md.T1)
  )
  
  # saveRDS(T1.counts, paste0("./analysis/", samplenames[i], "_new_count.rds"))
  # count_list <- c(count_list, T1.counts)
  T1_assay <- CreateChromatinAssay(T1.counts, fragments = frags.T1)
  T1 <- CreateSeuratObject(T1_assay, assay = "peaks", meta.data = md.T1, project=samplenames[i])

  T1$dataset <- samplenames[i]

  count_list <- c(count_list, T1)

}

# saveRDS(combined.peaks,"./analysis/combined.peaks.rds")
# saveRDS(count_list,  "./analysis/count_list.rds")
saveRDS(count_list,  concat_robject_dir)

# combined.peaks <- readRDS("./v3_analysis_240409/combined.peaks.rds")
# count_list <- readRDS("./v3_analysis_240409/count_list.rds")


# # Rename cells if provided
list_seurat <- count_list
combined <- merge(list_seurat[[1]], y = unlist(list_seurat)[2:length(unlist(list_seurat))], add.cell.ids = samplenames )

library(SeuratObject)
# extract gene annotations from EnsDb
annotations <- readRDS("/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/annotations.rds")
# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(combined) <- annotations

atacAggr<-combined
atacAggr <- FindTopFeatures(atacAggr, min.cutoff = 'q1')
atacAggr <- RunTFIDF(atacAggr)
atacAggr <- RunSVD(
  object = atacAggr,
  assay = 'peaks')
atacAggr <- RunUMAP(atacAggr, reduction = "lsi", dims = 2:30)
atacAggr <- FindNeighbors(atacAggr, reduction = "lsi", dims = 2:30)
atacAggr <- FindClusters(atacAggr, reduction = "lsi", dims = 2:30)

saveRDS(atacAggr, file = output_data_dir) 

### save PDF

pdf(umap1_dir) # _UMAP_by_sample_before_integration.pdf
p1 <- DimPlot(atacAggr, group.by = "dataset")
print(p1)
dev.off()


pdf(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_before_integration.pdf")))
p1 <- DimPlot(atacAggr, reduction = "umap",  label = TRUE, group.by = "seurat_clusters") +  ggtitle("scATAC Seurat Clustering")
print(p1)
dev.off()

### save png

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_cluster_before_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(atacAggr, reduction = "umap",  label = TRUE, group.by = "seurat_clusters") +  ggtitle("scATAC Seurat Clustering")
print(p1)
dev.off()

png(file.path(plot_folder, paste0(integrate_sample_name, "_UMAP_by_sample_before_integration.png")), width=6, height=6, units='in', res=300)
p1 <- DimPlot(atacAggr, reduction = "umap",  group.by = "dataset") + ggtitle("scATAC Sample Integration")
print(p1)
dev.off()

# ifnb <- SplitObject(atacAggr, split.by = "dataset")
# 
# # find integration anchors
# integration.anchors <- FindIntegrationAnchors(
#   object.list = ifnb,
#   # anchor.features = rownames(pbmc.multi),
#   reduction = "rlsi",
#   dims = 2:30
# )
# 
# # integrate LSI embeddings
# integrated <- IntegrateEmbeddings(
#   anchorset = integration.anchors,
#   reductions = atacAggr[["lsi"]],
#   new.reduction.name = "integrated_lsi",
#   dims.to.integrate = 1:30
# )
# 
# # create a new UMAP using the integrated embeddings
# integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
# # pdf("v3_analysis_240409/UMAP_by_dataset_after_integration.pdf")
# pdf(umap2_dir)
# p2 <- DimPlot(integrated, group.by = "dataset")
# print(p2)
# dev.off()
# 
# saveRDS(integrated, file = output_data_dir) 

# pdf("analysis/integration/merge_all/UMAP_by_dataset_after_integration_by_sample.pdf")
# pdf(umap1_dir)
# p2 <- DimPlot(integrated, group.by = "orig.ident")
# print(p2)
# dev.off()

# saveRDS(integrated, file = "scATAC_QC_integrated.rds") 

