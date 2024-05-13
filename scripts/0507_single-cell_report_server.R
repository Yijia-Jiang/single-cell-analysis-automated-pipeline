# Function to create multiple tabs
# https://stackoverflow.com/questions/46877451/generating-tabs-in-r-markdown-with-a-loop

single_path <- "/home/yj976/scpipeline/yj_version4/analysis/scRNA/single_sample"
integrate_path <- "/home/yj976/scpipeline/yj_version4/analysis/scRNA/integration"

single_sample_names <- list.files(single_path)
merge_sample_names <- list.files(integrate_path )

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
  # return(res)
  
}

# make.tabs.all <- function()
res2 <- NULL
# Create the Rmd to knit
res2 <- c(res2, 
  '---
title: "single-cell pipeline"
output: html_document
---
## scRNA-seq data
### Result {.tabset}

```{r}
# library(dplyr)
```
')

for (i in c(single_sample_names)){
  res2 <- c(
    res2,
  make.tabs(samplename = i, 
            images = c(file.path(single_path , i, paste0(i, '_QC_plot_before_original.png')),
                       file.path(single_path , i, paste0(i, "_QC_plot_after_original.png")),
                       file.path(single_path , i, paste0(i, "_UMAP_doublet.png")),
                       file.path(single_path , i, paste0(i, "_UMAP_by_cluster_number_original.png")),
                       file.path(single_path , i, paste0(i, "_UMAP_by_celltype_breast_markerbased.png")),
                       file.path(single_path , i, paste0(i, "_DEG_heatmap_assign.ident.breast.png")),
                       file.path(single_path , i, paste0(i, "_UMAP_feature_plot_breast_markers.png")),
                       file.path(single_path , i, paste0(i, "_UMAP_feature_plot_immune_markers.png")),
                       file.path(single_path , i, "GSEA", paste0(i, "_per_cluster"), paste0(i, "_per_cluster 1_HALLMARK_terms.png")),
                       file.path(single_path , i, "GSEA", paste0(i, "_per_celltype_marker_annot"), paste0(i, "_per_celltype_marker_annot Macrophages_HALLMARK_terms.pdf"))
                       ),
            titles = c("QC metrics of cells before QC filtering", 
                       "QC metrics of cells after QC filtering",
                       "Doublet detection",
                       "UMAP Clustering by cluster numbers",
                       "Celltype annotation based on breast markers",
                       "Heatmap of differentially expressed genes across annotated celltypes",
                       "Gene expression levels of breast atlas markers",
                       "Gene expression levels of immune cell markers",
                       "GSEA hallmark by cluster",
                       "GSEA hallmark by celltype"
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
              images = c(file.path(integrate_path , i, paste0(i, '_UMAP_by_sample_before_integration.png')),
                         file.path(integrate_path , i, paste0(i, "_UMAP_by_sample_after_integration.png")),
                         file.path(integrate_path , i, paste0(i, "_UMAP_by_celltype_breast_markerbased.png")),
                         file.path(integrate_path , i, paste0(i, "_DEG_heatmap_assign.ident.breast.png")),
                         file.path(integrate_path , i, paste0(i, "_UMAP_feature_plot_breast_markers.png")),
                         file.path(integrate_path , i, paste0(i, "_UMAP_feature_plot_immune_markers.png")),
                         file.path(integrate_path , i, "GSEA", paste0(i, "_per_cluster"), paste0(i, "_per_cluster 1_HALLMARK_terms.png"))
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

res2 <- cat(res2,
  sep = "",
  file = "file_summary.Rmd")

# Render the Rmd created into html here
rmarkdown::render("file_summary.Rmd")
