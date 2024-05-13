suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))

args <- commandArgs( trailingOnly = TRUE )
sample_path <- args[1] # metadata.csv
df_barplot1 <- args[2]
# df_barplot2 <- args[3] # hpca
# df_barplot3 <- args[4] # "sample/BA-1_209_NormR/BA-1_209_NormR_SingleR_encode_main.csv"
# df_barplot4 <- args[5]
df_barplot2 <- args[3]
df_barplot3 <- args[4]

analysis_out <- sapply(strsplit(sample_path, "/", fixed = TRUE),
                       function(i) paste(head(i, -1), collapse = "/"))

data <- read.csv(sample_path)
col2 <- colorRampPalette(brewer.pal(9,"Paired"))

# Plot Barplot for marker based annotation
pdf(df_barplot1)

stat1 <- data %>% group_by(orig.ident, assign.ident.breast) %>%
  summarise(count = n()) %>% as.data.frame()
p0 <- ggplot(stat1, aes(fill=assign.ident.breast, y=count, x=orig.ident)) +
  geom_bar(position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values = col2(length(unique(stat1$assign.ident.breast)))) +
  guides(fill=guide_legend(title="proportion of Breast MarkerBased Celltype by samples")) +
  ylab("Proportion") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

# data <- read.csv(sample_path)

stat1 <- data %>% group_by(orig.ident, assign.ident.cibersort) %>%
  summarise(count = n()) %>% as.data.frame()
p1 <- ggplot(stat1, aes(fill=assign.ident.cibersort, y=count, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$assign.ident.cibersort)))) + 
  guides(fill=guide_legend(title="proportion of CIBERSORT MarkerBased Celltype by samples")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

# data <- read.csv(sample_path)
stat1 <- data %>% group_by(orig.ident, SingleR_hpca_main) %>% 
  summarise(count = n()) %>% as.data.frame()
p2 <- ggplot(stat1, aes(fill=SingleR_hpca_main, y=count, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$SingleR_hpca_main)))) + 
  guides(fill=guide_legend(title="proportion of SingleR HPCA celltype by samples")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

# data <- read.csv(sample_path)
stat1 <- data %>% group_by(orig.ident, SingleR_encode_main) %>% 
  summarise(count = n()) %>% as.data.frame()
p3<-ggplot(stat1, aes(fill=SingleR_encode_main, y=count, x=orig.ident)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$SingleR_encode_main)))) + 
  guides(fill=guide_legend(title="proportion of SingleR ENCODE celltype by samples")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())
print(p0)
print(p1)
print(p2)
print(p3)
dev.off()

print("Done Plotting Barplot celltype by samples")

# Plot barplot for cell cycle
pdf(df_barplot2 )
stat1 <- data %>% group_by(seurat_clusters, Phase) %>% 
  summarise(count = n()) %>% as.data.frame()
p4 <- ggplot(stat1, aes(fill=Phase, y=count, x=factor(seurat_clusters))) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$Phase)))) + 
  guides(fill=guide_legend(title="proportion of cellcycle per cluster")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())
print(p4)
dev.off()
print("Done Plotting Barplot cellcycle by cluster")

pdf(df_barplot3)
stat1 <- data %>% group_by(assign.ident.breast, Phase) %>%
  summarise(count = n()) %>% as.data.frame()
p50<- ggplot(stat1, aes(fill=Phase, y=count, x=assign.ident.breast)) +
  geom_bar(position="fill", stat="identity", width = 0.5) +
  scale_fill_manual(values = col2(length(unique(stat1$Phase)))) +
  guides(fill=guide_legend(title="proportion of cellcycle per celltype marker based")) +
  ylab("Proportion") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

stat1 <- data %>% group_by(assign.ident.cibersort, Phase) %>% 
  summarise(count = n()) %>% as.data.frame()
p5<- ggplot(stat1, aes(fill=Phase, y=count, x=assign.ident.cibersort)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$Phase)))) + 
  guides(fill=guide_legend(title="proportion of cellcycle per celltype marker based")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

stat1 <- data %>% group_by(SingleR_hpca_main, Phase) %>% 
  summarise(count = n()) %>% as.data.frame()
p6<- ggplot(stat1, aes(fill=Phase, y=count, x=SingleR_hpca_main)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$Phase)))) + 
  guides(fill=guide_legend(title="proportion of cellcycle per celltype SingleR HPCA")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

stat1 <- data %>% group_by(SingleR_encode_main, Phase) %>% 
  summarise(count = n()) %>% as.data.frame()
p7<- ggplot(stat1, aes(fill=Phase, y=count, x=SingleR_encode_main)) + 
  geom_bar(position="fill", stat="identity", width = 0.5) + 
  scale_fill_manual(values = col2(length(unique(stat1$Phase)))) + 
  guides(fill=guide_legend(title="proportion of cellcycle per celltype SingleR HPCA")) + 
  ylab("Proportion") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90),
        axis.title.x=element_blank())

print(p50)
print(p5)
print(p6)
print(p7)
dev.off()

print("Done Plotting Barplot cellcycle by celltype")

# # Plot Barplot for singleR based annotation
# pdf( df_barplot2 )
# data <- read.csv(sample_path)
# stat1 <- data %>% group_by(orig.ident, SingleR_hpca_main) %>% 
#   summarise(count = n()) %>% as.data.frame()
# ggplot(stat1, aes(fill=SingleR_hpca_main, y=count, x=orig.ident)) + 
#   geom_bar(position="fill", stat="identity", width = 0.5) + 
#   scale_fill_manual(values = col2(length(unique(stat1$SingleR_hpca_main)))) + 
#   guides(fill=guide_legend(title="Celltype Proportion SingleR HPCA Based")) + 
#   ylab("Proportion") + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
#         axis.text.x = element_text(angle = 90),
#         axis.title.x=element_blank())
# dev.off()
# print("Done Plotting Barplot2")
# 
# # Plot Barplot for singleR based annotation2
# pdf(df_barplot3 )
# data <- read.csv(sample_path)
# stat1 <- data %>% group_by(orig.ident, SingleR_encode_main) %>% 
#   summarise(count = n()) %>% as.data.frame()
# ggplot(stat1, aes(fill=SingleR_encode_main, y=count, x=orig.ident)) + 
#   geom_bar(position="fill", stat="identity", width = 0.5) + 
#   scale_fill_manual(values = col2(length(unique(stat1$SingleR_encode_main)))) + 
#   guides(fill=guide_legend(title="Celltype Proportion SingleR ENCODE Based")) + 
#   ylab("Proportion") + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"), 
#         axis.text.x = element_text(angle = 90),
#         axis.title.x=element_blank())
# dev.off()
# print("Done Plotting Barplot3")


