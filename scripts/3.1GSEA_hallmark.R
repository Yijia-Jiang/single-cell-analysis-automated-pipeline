suppressMessages(library(magrittr))
suppressMessages(library(org.Hs.eg.db))

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")

suppressMessages(library(clusterProfiler)) #install 
suppressMessages(library(tidyverse))
suppressMessages(library(ggh4x))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(stringi))
suppressMessages(library(stringr))
suppressMessages(library(EnhancedVolcano))


args <- commandArgs( trailingOnly = TRUE )
diff1 <- args[1] # cluster
diff2 <- args[2] # marker based
diff3 <- args[3] # singleR HPCA
diff4 <- args[4] # singler ENCODE
output_dir <- args[5] # singler ENCODE
dir.create(output_dir, showWarnings = F)
########################## hallmark plot #############################
pcut = 10
minGSSize = 3
nPerm = 1000

########################## hallmark plot #############################
hallmarkplot <- function(Outdir, deseq2_mat, h_file, Condition, Treatment, 
                         Control, pcut = 10, minGSSize = 3, nPerm = 1000, 
                         plotname)
{
  ####read in data and convert to entrez ID
  data <- deseq2_mat
  # data <- read.csv(deseq2_mat, header = TRUE) %>% as.data.frame()
  data <- na.omit(data)
  colnames(data)[1] <- 'gene_id'
  # geneList <- sign(data$avg_log2FC) * (-log10(data$p_val_adj))
  geneList <- data$avglogFC
  geneList <- ifelse(geneList == Inf, 200, geneList)
  geneList <- ifelse(geneList == -Inf, -200, geneList)
  GeneIDSymbol <- toTable(org.Hs.egSYMBOL)
  names(geneList) <- GeneIDSymbol[match(data[,1],GeneIDSymbol$symbol),'gene_id']
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[!duplicated(names(geneList))]
  
  
  ####function of KEGG pathway enrichment
  GSEAKEGG <- function(geneList,minGSSize,nPerm,pcut){
    set.seed(1234)
    gsea.kegg <- gseKEGG(geneList = geneList,
                         organism = "hsa",
                         keyType = "kegg",
                         nPerm = nPerm, # number permutations
                         minGSSize = minGSSize,
                         pAdjustMethod = "BH",
                         pvalueCutoff = pcut, # padj cutoff value
                         verbose = FALSE,
                         seed = TRUE)
    return (gsea.kegg)
  }
  
  ####function of mSigDB hallmark geneset enrichment
  GSEAHALLMARK <- function(geneList,minGSSize,pcut,hallmark){
    set.seed(1234)
    gsea.hallmark <- GSEA(geneList = geneList,
                          TERM2GENE = hallmark, 
                          minGSSize = minGSSize,
                          pAdjustMethod = "BH",
                          pvalueCutoff = pcut, # padj cutoff value
                          verbose = FALSE,
                          seed = TRUE)
    return (gsea.hallmark)
  }
  
  
  format <- function(gseaRes){  
    enrich.res <- gseaRes@result %>%
      mutate(group = ifelse(NES > 0,paste("Enrich in",Treatment) ,paste("Enrich in",Control))) %>%
      mutate(significance = cut(p.adjust,breaks=c(-Inf, 0.01, 0.05, 0.1, Inf), label=c("***", "**", "*", ""))) %>%
      dplyr::filter(p.adjust <= pcut)
    
    res <- NULL
    if(dim(subset(enrich.res, NES > 0))[1] >= 10){
      res <- rbind(res,enrich.res[1:10,])
    }
    if(dim(subset(enrich.res, NES > 0))[1] < 10){
      res <- rbind(res,subset(enrich.res, NES > 0))
    }
    if(dim(subset(enrich.res, NES < 0))[1] >= 10){
      total <- dim(enrich.res)[1]
      res <- rbind(res,enrich.res[(total-9):total,])
    }
    if(dim(subset(enrich.res, NES < 0))[1] < 10){
      res <- rbind(res,subset(enrich.res, NES < 0))
    }
    res$Description <- factor(res$Description,levels = levels(reorder(res$Description,res$NES)))
    res.new <- res[match(levels(res$Description),res$Description),]
    return(res.new)
  }
  
  #new gesa plot
  bubble_plot <- function(data, treatment, control) {
    
    final_up <- data[data$NES > 0,]
    final_down <- data[data$NES < 0,]
    
    final_up <- final_up[order(final_up$pvalue, decreasing = FALSE),]
    final_down <- final_down[order(final_down$pvalue, decreasing = FALSE),]
    
    final <- rbind(head(final_up, 15), head(final_down, 15))
    final <- final[order(final$NES, decreasing = TRUE),]
    # final$qvalues <- ifelse(final$qvalue <= 0.15, final$qvalue, limit)
    final$qvalues <- final$qvalue
    final <- final[final$qvalues < 0.25, ]
    final$Description <- factor(final$Description, levels = rev(unique(final$Description)))
    
    final$condition <- ifelse(final$NES > 0, paste0("Enriched in ", treatment), paste0("Enriched in ", control))
    final$condition <- factor(final$condition, levels = c(paste0("Enriched in ", control), paste0("Enriched in ", treatment)))
    p <- ggplot(final, aes(x= NES, y= factor(Description), label = qvalues, color=qvalues,size=setSize)) + geom_point() +
      facet_grid(~condition, scales = "free_x") + theme_bw()
    
    p <- p + scale_color_gradient(low = "red", high = "blue", limits = c(0,0.25)) +
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = 18),
            strip.text = element_text(size = 16),
            axis.title.x = element_text(size = 18))
    
    write.csv(final, file = file.path(Outdir, paste(plotname,"_terms_table.csv", sep = "")))
    
    return(p)
  }
  
  
  writeoutput <- function(result,onto){
    p <- bubble_plot(result@result, Treatment, Control)
    pdf(file.path(Outdir, paste(plotname,"_",onto,"_terms.pdf", sep = "")),width = 30, height = 15)
    print(p)
    dev.off()
    
    png(file.path(Outdir, paste(plotname,"_",onto,"_terms.png", sep = "")),width = 30, height = 15, units='in', res=300)
    print(p)
    dev.off()
    
  }
  
  
  
  ###For HALLMARK
  hallmark_gene <- read.gmt(h_file)
  hallmark <- GSEAHALLMARK(geneList,minGSSize,pcut,hallmark_gene)
  hallmark@result$ID %<>% str_replace_all("GOBP_", "")
  hallmark@result$Description %<>% str_replace_all("GOBP_", "")
  
  # writeoutput(hallmark,"HALLMARK")
  writeoutput(hallmark,"GOBP")
  
}


draw_each_diff <- function(f){
  tmp <- head(tail(strsplit(f, "/", fixed = TRUE)[[1]], n=1), n=1)
  tmp <- str_replace(tmp, "_DEG_diffgene", "")
  sample <- str_replace(tmp, ".csv", "")
  print(sample)
  Outdir = file.path(output_dir, sample)
  print(Outdir)
  # dir.create("GSEA")
  dir.create(Outdir, showWarnings = F)
  deseq2_mat = file.path(f)
  
  h_file = "/mnt/cfce-rcsm/projects/yijiajiang_analysis/scpipeline/common_data/c5.go.bp.v2023.2.Hs.entrez.gmt.txt"
  Condition = "Comparison"

  data <- read.csv(deseq2_mat, header = TRUE) %>% as.data.frame()
  
  for(i in unique(data$cluster)){
    print(i)
    tmp <- data %>% filter(cluster == i)
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(hallmarkplot(Outdir = Outdir,
                          deseq2_mat = tmp,
                          h_file = h_file,
                          Condition = "Comparison",
                          Treatment = paste0("cluster ", i),
                          Control = "Others",
                          plotname = paste(sample, i)), 
             error = function(e) { 
               print(e)
               skip_to_next <-  TRUE
             })
    
    if(skip_to_next==TRUE) { next } 
    
    # hallmarkplot(Outdir = Outdir,
    #              deseq2_mat = tmp,
    #              h_file = h_file,
    #              Condition = "Comparison",
    #              Treatment = paste0("cluster ", i),
    #              Control = "Others",
    #              plotname = paste(sample, i))
  }
  
}


draw_each_diff(diff1) # sample/sample_diffgene.csv
draw_each_diff(diff2)
draw_each_diff(diff3)
draw_each_diff(diff4)


# diff1 <- "BA-1_209_NormR_DEG_diffgene_per_celltype_marker_annot.csv"
# f <- diff1
# # sample <- str_replace(f, "_DEG.csv", "")
# sample <- str_replace(f, "_DEG_diffgene", "")
# sample <- str_replace(sample, ".csv", "")
# Outdir = file.path("./GSEA", sample)
# dir.create("./GSEA")
# dir.create(Outdir)
# deseq2_mat = file.path("./", f)
# h_file = "/Users/jiang/Dropbox (Partners HealthCare)/CFCE_Collaboration_Labs/20230823_fixed_scRNA_11278/analysis/Merge_RNA_PBMC/c5.go.bp.v2023.2.Hs.entrez.gmt.txt"
# Condition = "Comparison"
# # Treatment = "Responder"
# # Control = "NonResponder"
# # plotname = paste("new", sample)
# data <- read.csv(deseq2_mat, header = TRUE) %>% as.data.frame()
# 
# for(i in unique(data$cluster)){
#   print(i)
#   tmp <- data %>% filter(cluster == i)
#   hallmarkplot(Outdir = Outdir,
#                deseq2_mat = tmp,
#                h_file = "/Users/jiang/Dropbox (Partners HealthCare)/CFCE_Collaboration_Labs/20230823_fixed_scRNA_11278/analysis/Merge_RNA_PBMC/c5.go.bp.v2023.2.Hs.entrez.gmt.txt",
#                Condition = "Comparison",
#                Treatment = paste0("cluster ", i),
#                Control = "Others",
#                plotname = paste(sample, i))
# }

# data %>% filter(cluster == 0)
#
# hallmarkplot(Outdir = file.path("DEG",celltype),
#              deseq2_mat = file.path("DEG", celltype, f),
#              h_file = "/Users/jiang/Dropbox (Partners HealthCare)/CFCE_Collaboration_Labs/20230823_fixed_scRNA_11278/analysis/Merge_RNA_PBMC/c5.go.bp.v2023.2.Hs.entrez.gmt.txt",
#              Condition = "Comparison",
#              Treatment = "Responder",
#              Control = "NonResponder",
#              plotname = paste("new", sample))

# Treatment = gsub("_over.*", "", "Bcells_stim_Res_over_NonRes_DEG")
# Control = gsub(".*over_", "", "Bcells_stim_Res_over_NonRes_DEG")
# plotname = paste("new", "Bcells_stim_Res_over_NonRes_DEG")




