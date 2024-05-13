suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # metadata.csv
df_stat <- args[2]

analysis_out <- sapply(strsplit(sample_path, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

# sample_path <- "/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/M2_atac_final_metadata.csv"
data <- read.csv(sample_path)
col2 <- colorRampPalette(brewer.pal(9,"Paired"))


stat1 <- data %>% group_by(dataset, cluster_annotation) %>%
  summarise(count = n()) %>% mutate(freq = round( (count / sum(count)),digits = 3)) %>% as.data.frame()
library(reshape)
stat1 <- cast(stat1, dataset~cluster_annotation, value = "freq") %>% as.data.frame()


stat0 <- data %>% group_by(dataset) %>%
  summarise(total_molecules = sum(nCount_peaks)) %>% as.data.frame()

stat2 <- data %>% group_by(dataset) %>%
  summarise(number_of_cells = n()) %>% as.data.frame()

stat3 <- data %>% group_by(dataset) %>%
  summarise(median_peaks_per_cell = median(nFeature_peaks)) %>% as.data.frame()

# stat4 <- data %>% group_by(orig.ident) %>%
#   summarise(median_mito_percent = median(percent.mt)) %>% as.data.frame()

df0 <- merge(x = stat0, y = stat2, by = "dataset", all.x = TRUE)
df1 <- merge(x = df0, y = stat3, by = "dataset", all.x = TRUE)
df3 <- merge(x = df1, y = stat1, by = "dataset", all.x = TRUE)
head(df3)

# write.csv(df3, "/Users/jiang/Dropbox (Partners HealthCare)/Single-cell_pipeline/scpipeline/scpipe_yijia_version4/BA-8_374_NormL_rna_stats.csv", row.names = F)
write.csv(df3, df_stat, row.names = F)

