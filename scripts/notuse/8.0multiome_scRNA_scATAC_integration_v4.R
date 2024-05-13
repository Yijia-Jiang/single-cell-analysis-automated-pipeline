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

rna_inputfile <- "analysis/scMultiome/scRNA/merge_2/merge_2_combine_integrated_annotated.rds"
atac_inputfile <- "analysis/scMultiome/scATAC/merge_2/merge_2_combine_integrated.rds"
output_data_dir <- "analysis/scMultiome/integration/merge_2/merge_2_scrna_scatac_integrated.rds"


plot_folder <- sapply(strsplit(output_data_dir, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))
integrate_sample_name <- head(tail(strsplit(output_data_dir, "/", fixed = TRUE)[[1]], n=2), n=1)



### atac
atac <- readRDS(atac_inputfile)
atac$assay <- "ATAC"
print(atac)
print(unique(atac$dataset))

DefaultAssay(atac) <- "peaks"
pdf(file.path(plot_folder, "ATAC_UMAP_only.pdf"))
DimPlot(atac, group.by = "dataset") +  ggtitle("peaks")
dev.off()

#### rna
# Adding ScaleData will introduce errors in FindTransferAnchors since gene activities are not scaled 
rna <- readRDS(rna_inputfile)
rna <- JoinLayers(rna)
rna$assay <- "RNA"

print(rna)
DefaultAssay(rna) <- "RNA"
# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = VariableFeatures(object = rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")




tryCatch(celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                              weight.reduction = atac[["integrated_lsi"]], dims = 2:30), 
         error = function(e) { 
           print(e)
           # skip_to_next <-  TRUE
         })

# if(skip_to_next==TRUE) { next } 

tryCatch(celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
                                              weight.reduction = atac[["lsi"]], dims = 2:30), 
         error = function(e) { 
           print(e)
           # skip_to_next <-  TRUE
         })

# if(skip_to_next==TRUE) { next } 

# if(skip_to_next==TRUE) { next } 

### test if "lsi" or integrated_lsi is present
# celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
#                                      weight.reduction = atac[["integrated_lsi"]], dims = 2:30)
# 
# celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$assign.ident.breast,
#                                      weight.reduction = atac[["lsi"]], dims = 2:30)

print(head(celltype.predictions))
atac <- AddMetaData(atac, metadata = celltype.predictions)

print(head(atac))
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

# tryCatch(imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["integrated_lsi"]],
#                                     dims = 2:30),
#          error = function(e) { 
#            # print(e)
#            imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
#                                       dims = 2:30)
#            
#            # skip_to_next <-  TRUE
#          })


tryCatch(imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["integrated_lsi"]],
                                    dims = 2:30), 
         error = function(e) { 
           print(e)
           # skip_to_next <-  TRUE
         })

# if(skip_to_next==TRUE) { next } 

tryCatch(imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
                                    dims = 2:30), 
         error = function(e) { 
           print(e)
           # skip_to_next <-  TRUE
         })


# imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["integrated_lsi"]],
#                            dims = 2:30)
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

pdf(file.path(plot_folder,"RNA_coembed_celltype.pdf"))
DimPlot(subset(coembed, subset = assay == "RNA"), group.by = "celltype")
dev.off()

pdf(file.path(plot_folder,"ATAC_coembed.pdf"))
DimPlot(subset(coembed, subset = assay != "RNA"))
dev.off()

pdf(file.path(plot_folder,"ATAC_coembed_celltype.pdf"))
DimPlot(subset(coembed, subset = assay != "RNA"), group.by = "celltype")
dev.off()

saveRDS(coembed, file.path(output_data_dir))

# saveRDS(coembed, file.path(plot_folder, "RNA_ATAC_coembed.rds"))
