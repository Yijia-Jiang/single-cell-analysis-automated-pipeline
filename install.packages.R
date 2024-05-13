if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("org.Hs.eg.db")

library(remotes)
remotes::install_version("matrixStats", version="1.1.0") # sessionInfo()


install.packages('irlba')

library(devtools)
devtools::install_github("teunbrand/ggh4x")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")


# scATAC

# install signac 1.18 new versions
install.packages("Signac")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("plyranges")

BiocManager::install("chromVAR")
install.packages('ggseqlogo')

install.packages("openxlsx")
install.packages("reshape")

# python
# pip install pyparsing