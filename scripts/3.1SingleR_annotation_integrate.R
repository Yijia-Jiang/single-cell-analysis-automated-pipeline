library(Seurat)
library(ggplot2)
# library(tidyverse)
library(stringi)
library(stringr)
library(SingleR)
library(Matrix)
library(tidyverse)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("celldex")

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # object filtered.rds: /mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/analysis/sample/BI54/BI54_Object_Filtered.rds
anno1_dir <- args[2] # /mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/analysis/sample/BI54/BI54_SingleR_hpca_main.csv
anno2_dir <- args[3] # /mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/analysis/sample/BI54/BI54_SingleR_encode_main.csv

# sample_path <- "/home/yj976/scpipeline/cobra_simple/Merge_all/Merge_all_Integrated.rds"
# anno1_dir <- "/home/yj976/scpipeline/cobra_simple/Merge_all/Merge_all_SingleR_hpca_main.csv"
# anno2_dir <- "/home/yj976/scpipeline/cobra_simple/Merge_all/Merge_all_SingleR_encode_main.csv"
integrate_sample_name <- head(tail(strsplit(anno1_dir, "/", fixed = TRUE)[[1]], n=2), n=1)
analysis_out <- sapply(strsplit(anno1_dir, "/", fixed = TRUE),
                   function(i) paste(head(i, -1), collapse = "/"))

# rnaAggr <- readRDS( file.path("/mnt/cfce-rcsm/projects/yijiajiang_analysis/20231207_Ellisen_fixed_scRNA/analysis/sample/BI54/BI54_Object_Filtered.rds"))

rnaAggr <- readRDS( file.path(sample_path))
### todo: add singleR annotation ###
# and save annotation as csv: barcode + annotation
rnaAggr <- JoinLayers(rnaAggr)

blueprintEncode.ref <- BlueprintEncodeData() 
#get reference atlas 
hpca.ref <- HumanPrimaryCellAtlasData()

#### subset reference ####
# library(remotes)
# remotes::install_version("matrixStats", version="1.1.0") # sessionInfo()

# https://github.com/LTLA/SingleR/issues/256
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE.
data <- as(rnaAggr[["RNA"]]$counts, Class = "dgCMatrix")
blueEncode.pred <- SingleR(test=data, assay.type.test=1, 
                           ref=blueprintEncode.ref, labels=blueprintEncode.ref$label.main)

hpca.main <- SingleR(test=data, assay.type.test=1, 
                     ref=hpca.ref, labels=hpca.ref$label.main)

hpca.main_anno <- as.data.frame(hpca.main) %>% select(labels)
colnames(hpca.main_anno) <- "SingleR_hpca_main"
df1 <- hpca.main_anno
df1$Barcode <- rownames(df1)
df1 <- df1 %>% select("Barcode", "SingleR_hpca_main")
write.csv(df1, file.path(anno1_dir), row.names = F)

rnaAggr_orig <- readRDS( file.path(sample_path))
rnaAggr_orig <- AddMetaData(rnaAggr_orig, df1)
pdf(file.path( analysis_out, paste0(integrate_sample_name,"_UMAP_by_celltype_singleR_HPCA.pdf")), height=10, width=12)
p1 <- DimPlot(rnaAggr_orig, reduction = "umap",  group.by = "SingleR_hpca_main", label=TRUE) + ggtitle("scRNA Annotation SingleR HPCA") 
print(p1)
dev.off()

png(file.path( analysis_out, paste0(integrate_sample_name,"_UMAP_by_celltype_singleR_HPCA.png")),width=10, height=12, units='in', res=300)
p1 <- DimPlot(rnaAggr_orig, reduction = "umap",  group.by = "SingleR_hpca_main", label=TRUE) + ggtitle("scRNA Annotation SingleR HPCA") 
print(p1)
dev.off()

#### Encode annotation ####
encode_anno <- as.data.frame(blueEncode.pred) %>% select(labels)
colnames(encode_anno) <- "SingleR_encode_main"
df1 <- encode_anno %>% select("SingleR_encode_main")
df1$Barcode <- rownames(df1)
df1 <- df1 %>% select("Barcode", "SingleR_encode_main")
write.csv(df1, file.path(anno2_dir), row.names = F)

rnaAggr_orig <- readRDS( file.path(sample_path))
rnaAggr_orig <- AddMetaData(rnaAggr_orig, df1)
pdf(file.path( analysis_out, paste0(integrate_sample_name, "_UMAP_by_celltype_singleR_ENCODE.pdf")), height=12, width=16)
p1 <- DimPlot(rnaAggr_orig, reduction = "umap",  group.by = "SingleR_encode_main", label=TRUE) + ggtitle("scRNA Annotation SingleR ENCODE")
print(p1)
dev.off()

png(file.path( analysis_out, paste0(integrate_sample_name, "_UMAP_by_celltype_singleR_ENCODE.png")), width=10, height=12, units='in', res=300)
p1 <- DimPlot(rnaAggr_orig, reduction = "umap",  group.by = "SingleR_encode_main", label=TRUE) + ggtitle("scRNA Annotation SingleR ENCODE")
print(p1)
dev.off()



