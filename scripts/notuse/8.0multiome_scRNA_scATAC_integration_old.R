# this script will eliminate doublets from an aggregated snRNA object prior to preprocessing an aggregated snRNA dataset
library(Seurat) # 3.02
library(ggplot2)
library(harmony) # 1.0
library(Rcpp)
library(DoubletFinder)
library(openxlsx)
library(here)
library(dplyr)
library(stringr)
library(tibble)
set.seed(1234)

library(Signac)

args <- commandArgs( trailingOnly = TRUE )

rna_inputfile = args[1] # combine_integrated_annotated.rds
atac_inputfile = args[2] # combine_integrated.rds
output_data_dir = args[3] 

plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

# rna_inputfile <- ""
# atac_inputfile <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/integration/merge_1/merge_1_combine_integrated.rds"

#### atac
atac <- readRDS(atac_inputfile)
DefaultAssay(atac) <- "peaks"
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pdf(file.path(plot_folder, "ATAC_UMAP_only.pdf"))
DimPlot(atac, group.by = "orig.ident", label = FALSE) +  ggtitle("peaks")
dev.off()

#### rna
# Adding ScaleData will introduce errors in FindTransferAnchors since gene activities are not scaled 
rna <- readRDS(rna_inputfile)
DefaultAssay(rna) <- "RNA"
rna <- SCTransform(rna,  verbose = FALSE)
# rna <- RunPCA(rna, verbose = FALSE)
# rna <- RunUMAP(rna, dims = 1:30, verbose = FALSE)

### integrate
ifnb.list <- list(rna, atac)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", anchor.features = features,
                                        verbose = FALSE)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
                            verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)

pdf(file.path(plot_folder,"RNA_ATAC_integration.pdf"))
DimPlot(integrated,group.by = "orig.ident")
dev.off()

predicted.labels <- TransferData(
  anchorset = anchors,
  refdata = rna$assign.ident.breast
  #weight.reduction = atac[["lsi"]], # this is actually an lsi reduction
  #dims = 1:30
)

write.csv(predicted.labels, file.path(plot_folder, "atac_predicted_labels.csv"))

atac <- AddMetaData(atac, metadata = predicted.labels)
pdf(file.path(plot_folder,"scATAC_predicted_celltypes_gene_activity_anchor.pdf"))
DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
  ggtitle("scATAC-seq Predicted Celltypes") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
dev.off()

write.csv(atac@meta.data, file.path(plot_folder, "atac_metadata_annotated.csv"))

saveRDS(integrated, output_data_dir)


# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = anchors, refdata = refdata, weight.reduction = atac[["lsi"]], dims = 2:30)

# this line adds the imputed data matrix to the pbmc.atac object
atac[["RNA_impute"]] <- imputation
coembed <- merge(x = rna, y = atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)


head(coembed@meta.data)
tail(coembed@meta.data)
dim(coembed)
coembed$celltype <- ifelse(!is.na(coembed$assign.ident.breast), coembed$assign.ident.breast, coembed$predicted.id)

pdf(file.path(plot_folder,"RNA_ATAC_coembed_celltype.pdf"))
Idents(coembed) <- coembed$celltype
p1 <- DimPlot(coembed)
print(p1)
dev.off()

pdf(file.path(plot_folder,"RNA_ATAC_coembed_assay.pdf"))
Idents(coembed) <- coembed$orig.ident
p1 <- DimPlot(coembed)
print(p1)
dev.off()

unique(coembed$orig.ident)
pdf(file.path(plot_folder,"RNA_coembed.pdf"))
DimPlot(subset(coembed, subset = orig.ident == "scRNA"))
dev.off()

pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
DimPlot(subset(coembed, subset = orig.ident != "scRNA"))
dev.off()

saveRDS(coembed, file.path(plot_folder, "rna_atac_coembed.rds"))
# # # # identify anchors to transfer cell labels from snRNAseq to snATACseq "RNA" gene activity scores
# # transfer.anchors <- FindTransferAnchors(
# #   reference = rna,
# #   query = atac,
# #   reference.assay = 'RNA',
# #   query.assay = 'RNA',
# #   reduction = 'cca',
# # )
# # 
# # # save transfer anchors for plotting coembedded snRNA and snATAC
# # # saveRDS(transfer.anchors, file = file.path("cellranger_atac_prep/transfer_anchors.rds"))
# # transfer.anchors <- readRDS( file = file.path("cellranger_atac_prep/transfer_anchors.rds"))
# 
# # this is high-resolution celltype prediction which is great for predicting celltypes
# # but may not be useful for thresholding snATAC data
# # predicted.labels <- TransferData(
# #   anchorset = transfer.anchors,
# #   refdata = rna$major_annotation,
# #   weight.reduction = atac[["lsi"]], # this is actually an lsi reduction
# #   dims = 2:30
# # )
# # 
# # # add predicted cell types to the snATACseq object
# # atac <- AddMetaData(atac, metadata = predicted.labels)
# # head(atac@meta.data)
# # pdf("scATAC_predicted_celltypes_new.pdf")
# # DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
# #   ggtitle("scATAC-seq Predicted Celltypes") + 
# #   NoLegend() + scale_colour_hue(drop = FALSE)
# # dev.off()
# # 
# # p2 <- DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
# #   ggtitle("scATAC-seq Predicted Celltypes") + 
# #   NoLegend() + scale_colour_hue(drop = FALSE)
# # p3 <- DimPlot(rna, reduction = "umap", assay = "RNA", label = TRUE, repel = TRUE) +
# #   ggtitle("scRNA-seq Annotated Celltypes") + 
# #   NoLegend()
# # CombinePlots(plots = list(p2, p3))
# 
# 
# # atac <- readRDS(file = file.path("/Users/jiang/Dropbox (Partners HealthCare)/CFCE_Collaboration_Labs/20230907_NovaSeq_scATAC_breast/scATAC-cellranger/aggr_noNorm/analysis/cellranger_atac_prep/scATAC_transfer_anchors.rds"))
# 
# # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# # full transcriptome if we wanted to
# genes.use <- VariableFeatures(rna)
# refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]
# 
# # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["pca"]], dims = 2:20)
# 
# # this line adds the imputed data matrix to the pbmc.atac object
# atac[["RNA_impute"]] <- imputation
# coembed <- merge(x = rna, y = atac)
# 
# # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# # datasets
# coembed
# coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
# coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
# coembed <- RunUMAP(coembed, dims = 1:30)
# 
# 
# head(coembed@meta.data)
# tail(coembed@meta.data)
# dim(coembed)
# coembed$celltype <- ifelse(!is.na(coembed$major_annotation), coembed$major_annotation, coembed$predicted.id)
# 
# pdf(file.path(plot_folder,"RNA_ATAC_coembed.pdf"))
# Idents(coembed) <- coembed$celltype
# DimPlot(coembed)
# dev.off()
# 
# unique(coembed$orig.ident)
# pdf(file.path(plot_folder,"RNA_coembed.pdf"))
# DimPlot(subset(coembed, subset = orig.ident == "10x_RNA"))
# dev.off()
# 
# pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
# DimPlot(subset(coembed, subset = orig.ident != "10x_RNA"))
# dev.off()