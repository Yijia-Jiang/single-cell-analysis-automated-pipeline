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

# rna_inputfile <- "analysis/scMultiome/scRNA/merge_1/merge_1_combine_integrated_annotated.rds"
# atac_inputfile <- "analysis/scMultiome/scATAC/merge_2/merge_2_combine_integrated.rds"
# output_data_dir <- "analysis/scMultiome/integration/merge_1/merge_1_scrna_scatac_integrated.rds"


plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)

print(unique(atac$dataset))

if (unique(atac$dataset) == 1){
  atac <- readRDS(atac_inputfile)
  DefaultAssay(atac) <- "peaks"

  pdf(file.path(plot_folder, "ATAC_UMAP_only.pdf"))
  DimPlot(atac, group.by = "dataset") +  ggtitle("peaks")
  dev.off()
  
  #### rna
  # Adding ScaleData will introduce errors in FindTransferAnchors since gene activities are not scaled 
  rna <- readRDS(rna_inputfile)
  rna <- JoinLayers(rna)
  rna$assay <- "RNA"
  DefaultAssay(rna) <- "RNA"
  # rna <- SCTransform(rna,  verbose = FALSE)
  # rna <- RunPCA(rna, verbose = FALSE)
  # rna <- RunUMAP(rna, dims = 1:30, verbose = FALSE)
  
  # Identify anchors
  transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                          reference.assay = "RNA", query.assay = "RNA", reduction = "cca")
  
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                       weight.reduction = atac[["integrated_lsi"]], dims = 2:30)
  
  celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                       weight.reduction = atac[["lsi"]], dims = 2:30)
  
  atac <- AddMetaData(atac, metadata = celltype.predictions)
  
  write.csv(celltype.predictions , file.path(plot_folder, "ATAC_predicted_labels.csv"))
  
  pdf(file.path(plot_folder,"ATAC_predicted_celltypes_gene_activity_anchor.pdf"))
  DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
    ggtitle("scATAC-seq Predicted Celltypes") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
  dev.off()
  
  write.csv(atac@meta.data, file.path(plot_folder, "ATAC_metadata_annotated.csv"))
  
  # saveRDS(integrated, output_data_dir)
  
  # full transcriptome if we wanted to
  genes.use <- VariableFeatures(rna)
  refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]
  
  # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
  # imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
  #                            dims = 2:30)
  
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["integrated_lsi"]],
                             dims = 2:30)
  atac[["RNA"]] <- imputation
  
  coembed <- merge(x = rna, y = atac)
  # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
  # datasets
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  
  head(coembed@meta.data)
  dim(coembed)
  coembed$celltype <- ifelse(!is.na(coembed$assign.ident.breast), coembed$assign.ident.breast, coembed$predicted.id)
  
  pdf(file.path(plot_folder,"RNA_ATAC_coembed_celltype.pdf"))
  Idents(coembed) <- coembed$celltype
  p1 <- DimPlot(coembed, label = TRUE, repel = TRUE)
  print(p1)
  dev.off()
  
  pdf(file.path(plot_folder,"RNA_ATAC_coembed_assay.pdf"))
  Idents(coembed) <- coembed$assay
  p1 <- DimPlot(coembed)
  print(p1)
  dev.off()
  
  
  pdf(file.path(plot_folder,"RNA_coembed.pdf"))
  DimPlot(subset(coembed, subset = assay == "RNA"))
  dev.off()
  
  pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
  DimPlot(subset(coembed, subset = assay != "RNA"))
  dev.off()
  
  saveRDS(coembed, file.path(plot_folder, "RNA_ATAC_coembed.rds"))
  
} else {
#### atac
atac <- readRDS(atac_inputfile)
DefaultAssay(atac) <- "peaks"

gene.activities <- GeneActivity(atac)
atac[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
    object = atac,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(atac$nCount_RNA)
  )
DefaultAssay(atac) <- 'ACTIVITY'
atac <- ScaleData(atac, features = rownames(atac))
# atac$assay <- "ATAC"
DefaultAssay(atac) <- "peaks"
# atac <- RunSVD(atac)
# atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# DefaultAssay(atac) <- "ATAC"
pdf(file.path(plot_folder, "ATAC_UMAP_only.pdf"))
DimPlot(atac, group.by = "dataset") +  ggtitle("peaks")
dev.off()

#### rna
# Adding ScaleData will introduce errors in FindTransferAnchors since gene activities are not scaled 
rna <- readRDS(rna_inputfile)
rna <- JoinLayers(rna)
rna$assay <- "RNA"
DefaultAssay(rna) <- "RNA"
# rna <- SCTransform(rna,  verbose = FALSE)
# rna <- RunPCA(rna, verbose = FALSE)
# rna <- RunUMAP(rna, dims = 1:30, verbose = FALSE)

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                     weight.reduction = atac[["integrated_lsi"]], dims = 2:30)

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                     weight.reduction = atac[["lsi"]], dims = 2:30)

atac <- AddMetaData(atac, metadata = celltype.predictions)

write.csv(celltype.predictions , file.path(plot_folder, "ATAC_predicted_labels.csv"))

pdf(file.path(plot_folder,"ATAC_predicted_celltypes_gene_activity_anchor.pdf"))
DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
  ggtitle("scATAC-seq Predicted Celltypes") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
dev.off()

write.csv(atac@meta.data, file.path(plot_folder, "ATAC_metadata_annotated.csv"))

# saveRDS(integrated, output_data_dir)

# full transcriptome if we wanted to
genes.use <- VariableFeatures(rna)
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
#                            dims = 2:30)

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["integrated_lsi"]],
                           dims = 2:30)
atac[["RNA"]] <- imputation

coembed <- merge(x = rna, y = atac)
# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

head(coembed@meta.data)
dim(coembed)
coembed$celltype <- ifelse(!is.na(coembed$assign.ident.breast), coembed$assign.ident.breast, coembed$predicted.id)

pdf(file.path(plot_folder,"RNA_ATAC_coembed_celltype.pdf"))
Idents(coembed) <- coembed$celltype
p1 <- DimPlot(coembed, label = TRUE, repel = TRUE)
print(p1)
dev.off()

pdf(file.path(plot_folder,"RNA_ATAC_coembed_assay.pdf"))
Idents(coembed) <- coembed$assay
p1 <- DimPlot(coembed)
print(p1)
dev.off()


pdf(file.path(plot_folder,"RNA_coembed.pdf"))
DimPlot(subset(coembed, subset = assay == "RNA"))
dev.off()

pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
DimPlot(subset(coembed, subset = assay != "RNA"))
dev.off()

saveRDS(coembed, file.path(plot_folder, "RNA_ATAC_coembed.rds"))
}

# ### integrate
# ifnb.list <- list(rna, atac)
# features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
# anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", anchor.features = features,
#                                         verbose = FALSE)
# integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", 
#                             verbose = FALSE)
# integrated <- RunPCA(integrated, verbose = FALSE)
# integrated <- RunUMAP(integrated, dims = 1:30, verbose = FALSE)
# 
# pdf(file.path(plot_folder,"RNA_ATAC_integration.pdf"))
# DimPlot(integrated,group.by = "orig.ident")
# dev.off()
# 
# predicted.labels <- TransferData(
#   anchorset = anchors,
#   refdata = rna$assign.ident.breast
#   #weight.reduction = atac[["lsi"]], # this is actually an lsi reduction
#   #dims = 1:30
# )
# 
# write.csv(predicted.labels, file.path(plot_folder, "atac_predicted_labels.csv"))
# 
# atac <- AddMetaData(atac, metadata = predicted.labels)
# pdf(file.path(plot_folder,"scATAC_predicted_celltypes_gene_activity_anchor.pdf"))
# DimPlot(atac, group.by = "predicted.id", label = TRUE, repel = TRUE) +
#   ggtitle("scATAC-seq Predicted Celltypes") + 
#   NoLegend() + scale_colour_hue(drop = FALSE)
# dev.off()
# 
# write.csv(atac@meta.data, file.path(plot_folder, "atac_metadata_annotated.csv"))
# 
# saveRDS(integrated, output_data_dir)
# 
# 
# # note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# # full transcriptome if we wanted to
# 
# # refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# # (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
# imputation <- TransferData(anchorset = anchors, refdata = refdata, weight.reduction = atac[["lsi"]], dims = 2:30)
# 
# # this line adds the imputed data matrix to the pbmc.atac object
# atac[["RNA_impute"]] <- imputation
# coembed <- merge(x = rna, y = atac)
# 
# # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# # datasets
# coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
# coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
# coembed <- RunUMAP(coembed, dims = 1:30)
# 
# 
# head(coembed@meta.data)
# tail(coembed@meta.data)
# dim(coembed)
# coembed$celltype <- ifelse(!is.na(coembed$assign.ident.breast), coembed$assign.ident.breast, coembed$predicted.id)
# 
# pdf(file.path(plot_folder,"RNA_ATAC_coembed_celltype.pdf"))
# Idents(coembed) <- coembed$celltype
# p1 <- DimPlot(coembed)
# print(p1)
# dev.off()
# 
# pdf(file.path(plot_folder,"RNA_ATAC_coembed_assay.pdf"))
# Idents(coembed) <- coembed$orig.ident
# p1 <- DimPlot(coembed)
# print(p1)
# dev.off()
# 
# unique(coembed$orig.ident)
# pdf(file.path(plot_folder,"RNA_coembed.pdf"))
# DimPlot(subset(coembed, subset = orig.ident == "scRNA"))
# dev.off()
# 
# pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
# DimPlot(subset(coembed, subset = orig.ident != "scRNA"))
# dev.off()
# 
# saveRDS(coembed, file.path(plot_folder, "rna_atac_coembed.rds"))
