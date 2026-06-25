if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("Biobase")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("GSEABase")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
install.packages(dplyr)
BiocManager::install('seandavi/GEOquery')

if (!require('remotes')) install.packages('remotes') # make sure you have Rtools installed first! if not, then run:
#install.packages('installr')
#install.Rtools()
remotes::install_github('talgalili/installr')

install(gplots)
install(ggrepel)

BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("fgsea")
BiocManager::install("biomaRt")
install("DT")
