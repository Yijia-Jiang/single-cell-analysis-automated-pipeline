suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # metadata.csv
df_stat <- args[2]

analysis_out <- sapply(strsplit(sample_path, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

# sample_path <- "/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/merge_1_final_metadata.csv"
data <- read.csv(sample_path)
col2 <- colorRampPalette(brewer.pal(9,"Paired"))


stat1 <- data %>% group_by(orig.ident, assign.ident.breast) %>%
  summarise(count = n()) %>% mutate(freq = count / sum(count)) %>% as.data.frame()
library(reshape)
stat1 <- cast(stat1, orig.ident~assign.ident.breast, value = "freq") %>% as.data.frame()

stat2 <- data %>% group_by(orig.ident) %>%
  summarise(number_of_cells = n()) %>% as.data.frame()

stat3 <- data %>% group_by(orig.ident) %>%
  summarise(median_genes_per_cell = median(nFeature_RNA)) %>% as.data.frame()

stat4 <- data %>% group_by(orig.ident) %>%
  summarise(median_mito_percent = median(percent.mt)) %>% as.data.frame()

df1 <- merge(x = stat2, y = stat3, by = "orig.ident", all.x = TRUE)
df2 <- merge(x = df1, y = stat4, by = "orig.ident", all.x = TRUE) 
df3 <- merge(x = df2, y = stat1, by = "orig.ident", all.x = TRUE)
head(df3)

# write.csv(df3, "/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/merge_1_stats.csv", row.names = F)
write.csv(df3, df_stat, row.names = F)

