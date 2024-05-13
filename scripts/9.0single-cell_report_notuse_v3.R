# Function to create multiple tabs
# https://stackoverflow.com/questions/46877451/generating-tabs-in-r-markdown-with-a-loop

args <- commandArgs( trailingOnly = TRUE )

single_path = args[1] # "/home/yj976/scpipeline/yj_version4/analysis/scRNA/single_sample"
integrate_path = args[2] # combine_integrated.rds
atlas = args[3]
single_path_atac <- args[4]
integrate_path_atac <- args[5]
integrate_path_scrna_scatac <- args[6]
output <- args[7]

outputname <- strsplit(output, "/")[[1]][length(strsplit(output, "/")[[1]])]
plot_folder <- sapply(strsplit(output, "/", fixed = TRUE),
                      function(i) paste(head(i, -1), collapse = "/"))

print(outputname)

project_dir <- here::here()

# ### scRNA
# single_path <- "/home/yj976/scpipeline/yj_version4/analysis/scRNA/single_sample"
# integrate_path <- "/home/yj976/scpipeline/yj_version4/analysis/scRNA/integration"

single_sample_names <- list.files(single_path)
merge_sample_names <- list.files(integrate_path )

# ### scATAC
# atlas = "PBMC"
# single_path_atac <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/single_sample"
# integrate_path_atac <- "/home/yj976/scpipeline/yj_version4/analysis/scATAC/integration"

single_sample_names_atac <- list.files(single_path_atac)
merge_sample_names_atac <- list.files(integrate_path_atac )

# ### integration
# integrate_path_scrna_scatac <- "/home/yj976/scpipeline/yj_version4/analysis/scMultiome/integration"
rna_atac_merge_sample_names <- list.files(integrate_path_scrna_scatac)

make.tabs <- function(samplename, images, titles){
  res <- NULL
  res <- c(res, '#### ', samplename, '\n')
  for(i in seq_along(images)){
    
    res <- c(res, 
             "```{r, echo=FALSE, fig.cap='", titles[i], "', out.width = '60%', fig.topcaption=TRUE}", '\n',
             "knitr::include_graphics('", images[i], "')", '\n',
             '```', '\n\n')
    
  }
  return(res)
}

make.df <- function(samplename, dfs){
  library(tidyverse)
  res <- NULL
  for(i in seq_along(dfs)){
    summary<-read.csv(df, row.names = 1)
    summary<-as.data.frame((summary))
    res <- c(res, 
             "```{r, echo=FALSE, fig.cap='", titles[i], "', out.width = '50%'}", '\n',
             "print(summary %>% head(5))", '\n',
             '```', '\n\n')
    
  }
  return(res)
  
}

# make.tabs.all <- function()
res2 <- NULL
# Create the Rmd to knit
res2 <- c(res2, 
  '---
title: "Single-cell Pipeline"
author: Yijia Jiang
institution: CFCE -Dana Darber
output: html_document
---
## scRNA-seq data
### Result {.tabset}

')

for (i in c(single_sample_names)){
  res2 <- c(
    res2,
  make.tabs(samplename = i, 
            images = c(file.path(project_dir, single_path , i, paste0(i, '_QC_plot_before_original.png')),
                       file.path(project_dir, single_path , i, paste0(i, "_QC_plot_after_original.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_UMAP_doublet.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_UMAP_by_cluster_number_original.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_UMAP_by_celltype_breast_markerbased.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_DEG_heatmap_assign.ident.breast.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_UMAP_feature_plot_breast_markers.png")),
                       file.path(project_dir,single_path , i, paste0(i, "_UMAP_feature_plot_immune_markers.png")),
                       file.path(project_dir,single_path , i, "GSEA", paste0(i, "_per_cluster"), paste0(i, "_per_cluster 1_HALLMARK_terms.png"))
                       ),
            titles = c("QC metrics of cells before QC filtering", 
                       "QC metrics of cells after QC filtering",
                       "Doublet detection",
                       "UMAP Clustering by cluster numbers",
                       "Celltype annotation based on breast markers",
                       "Heatmap of differentially expressed genes across annotated celltypes",
                       "Gene expression levels of breast atlas markers",
                       "Gene expression levels of immune cell markers",
                       "GSEA hallmark by cluster"
                       ))
  )
  
  # make.df(samplename = i,
  #         df = c(file.path(single_path , i, paste0(i, '_DEG_diffgene_per_cluster.csv')))
  # )
  
}

for (i in c(merge_sample_names)){
  res2 <- c(
    res2,
    make.tabs(samplename = i,
              images = c(file.path(project_dir,integrate_path , i, paste0(i, '_UMAP_by_sample_before_integration.png')),
                         file.path(project_dir,integrate_path , i, paste0(i, "_UMAP_by_sample_after_integration.png")),
                         file.path(project_dir,integrate_path , i, paste0(i, "_UMAP_by_celltype_breast_markerbased.png")),
                         file.path(project_dir,integrate_path , i, paste0(i, "_DEG_heatmap_assign.ident.breast.png")),
                         file.path(project_dir,integrate_path , i, paste0(i, "_UMAP_feature_plot_breast_markers.png")),
                         file.path(project_dir,integrate_path , i, paste0(i, "_UMAP_feature_plot_immune_markers.png")),
                         file.path(project_dir,integrate_path , i, "GSEA", paste0(i, "_per_cluster"), paste0(i, "_per_cluster 1_HALLMARK_terms.png"))
              ),
              titles = c("UMAP by samples before integration",
                         "UMAP by samples after batch correction",
                         "Celltype annotation based on breast markers",
                         "Heatmap of differentially expressed genes across annotated celltypes",
                         "Heatmap of differentially expressed genes across annotated celltypes",
                         "Gene expression levels of breast atlas markers",
                         "Gene expression levels of immune cell markers",
                         "GSEA hallmark by cluster"
              ))

  )

}


res2 <- c(res2, 
          '
## scATAC-seq data

### Result {.tabset}
')

for (i in c(single_sample_names_atac)){
  res2 <- c(
    res2,
    make.tabs(samplename = i, 
              images = c(file.path(project_dir,single_path_atac , i, paste0(i, '_QC_plot_before_original.png')),
                         file.path(project_dir,single_path_atac , i, paste0(i, "_QC_plot_after_original.png")),
                         file.path(project_dir,single_path_atac , i, paste0(i, "_UMAP_by_cluster.png")),
                         file.path(project_dir,single_path_atac , i, paste0("scATAnno_", atlas), "celltype_assignment/query_assignment/Legend_side", paste0( "umap_3.query_cluster_annotation.png"))

                         ),
              titles = c("QC metrics of cells before QC filtering", 
                         "QC metrics of cells after QC filtering",
                         "UMAP Clustering by cluster numbers",
                         "UMAP Clustering by scATAnno Annotation"
              ))
  )
  
  # make.df(samplename = i,
  #         df = c(file.path(single_path , i, paste0(i, '_DEG_diffgene_per_cluster.csv')))
  # )
  
}

for (i in c(merge_sample_names_atac)){
  res2 <- c(
    res2,
    make.tabs(samplename = i,
              images = c(file.path(project_dir,integrate_path_atac , i, paste0(i, '_UMAP_by_sample_before_integration.png')),
                         file.path(project_dir,integrate_path_atac , i, paste0(i, "_UMAP_by_sample_after_integration.png")),
                         file.path(project_dir,integrate_path_atac , i, paste0("scATAnno_", atlas), "celltype_assignment/query_assignment/Legend_side", paste0( "umap_3.query_cluster_annotation.png"))
              ),
              titles = c("UMAP by scATAC samples before integration",
                         "UMAP by scATAC samples after batch correction",
                         "UMAP Clustering by scATAnno Annotation"
              ))

  )

}

res2 <- c(res2, 
          '
## scRNA-scATAC integration 

### Result {.tabset}
')

for (i in c(rna_atac_merge_sample_names)){
  res2 <- c(
    res2,
    make.tabs(samplename = i,
              images = c(file.path(project_dir,integrate_path_scrna_scatac , i, paste0(i, '_RNA_ATAC_coembed_assay.png')),
                         file.path(project_dir,integrate_path_scrna_scatac , i, paste0(i, '_RNA_ATAC_coembed_sample.png')),
                         file.path(project_dir,integrate_path_scrna_scatac , i, paste0(i, "_RNA_ATAC_coembed_celltype.png"))

              ),
              titles = c("UMAP integration by assays",
                         "UMAP integration by samples",
                         "UMAP integration by celltypes"
              ))
  )


}


res2 <- cat(res2,
            sep = "",
            file = output)

rmarkdown::render(output)
