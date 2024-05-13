suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # metadata.csv
df_barplot1 <- args[2]

# df_barplot2 <- args[3]
# df_barplot3 <- args[4]

analysis_out <- sapply(strsplit(sample_path, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

data <- read.csv(sample_path)
col2 <- colorRampPalette(brewer.pal(9,"Paired"))

# Plot Barplot for marker based annotation
pdf(df_barplot1)

stat1 <- data %>% group_by(dataset, cluster_annotation) %>%
  summarise(count = n()) %>% as.data.frame()
p1 <- ggplot(stat1, aes(fill=cluster_annotation, y=count, x=dataset)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$cluster_annotation)))) + 
  guides(fill=guide_legend(title="proportion of Celltype by samples")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

print(p1)

dev.off()

print("Done Plotting Barplot celltype by samples")

