suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

args <- commandArgs( trailingOnly = TRUE )
input <- args[1]

# input <- "analysis/scRNA/single_sample/BI54_rna/BI54_rna_stats.csv analysis/scRNA/single_sample/MGH23039R_rna/MGH23039R_rna_stats.csv analysis/scRNA/single_sample/BA-1_209_NormR_rna/BA-1_209_NormR_rna_stats.csv analysis/scRNA/single_sample/BA-8_374_NormL_rna/BA-8_374_NormL_rna_stats.csv"
inputFiles = strsplit(input,' ')[[1]]
df_stat <- args[2]


# inputFiles <- c("/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/BA-8_374_NormL_rna_stats.csv",
#               "/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/BA-8_374_NormL_rna_stats.csv")

library(plyr)
df_final = data.frame()
for(file in inputFiles){
  print(file)
  df <- read.csv(file, check.names = F)
  # print(head(df))
  df_final <- rbind.fill(df_final, df)
  
}

print(df_final)

write.csv(df_final, df_stat, row.names = F)

