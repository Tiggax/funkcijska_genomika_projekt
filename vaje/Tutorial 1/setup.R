# Installing R
# Windows:
#   https://cloud.r-project.org/bin/windows/base/
# Mac:
#   Go to https://cloud.r-project.org/bin/macosx/
# Installing RStudio
# Go to https://www.rstudio.com/products/rstudio/download/#download.
# Installing tidyverse
# Open RStudio. Go to Tools > Install Packages. Enter tidyverse, then select Install.
# You are ready for Data Science Essentials for R!

#recall 
#
# ctrl+enter execute line
# ctrl+shfit+enter execute whole file
# ctrl+shift+c comment selected

#installing some packages

installandload<-function(
  normalpackages=c("rlang","tidyverse","readxl","writexl","lubridate","stringr","esquisse","miniUI","pheatmap","RColorBrewer","hexbin","BiocManager","reshape2","ggnewscale","ggbeeswarm"),
  biocpackages=c("DESeq2","SummarizedExperiment","GEOquery","biomaRt","org.Hs.eg.db","org.Mm.eg.db","pathview","clusterProfiler","EnhancedVolcano","apeglm","genefilter","tximport","tximeta","GenomicFeatures","edgeR")){
  new.packages <- normalpackages[!(normalpackages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) {install.packages(new.packages,dependencies =T,quiet = T)}
  new.bpackages <- biocpackages[!(biocpackages %in% installed.packages()[,"Package"])]
  if(length(new.bpackages)) BiocManager::install(new.bpackages,update=T,ask=F)
  all<-c(normalpackages,biocpackages)
  invisible(lapply(all,function(x){
    if(require(x,character.only=T,quietly=T)) library(x, quietly=T,character.only=T)}))
  message("the packages")
  print(all[all %in% installed.packages()[,"Package"]])
  message("are loaded")
}
installandload()

#
#BiocManager::install("DESeq2")

# Tutorial overview
# 1. Installation of R and packages.
# 2. Recall R, recall R markdown (Submit projects in R markdown)
# 3. Download data and metadata from GEO omnibus. (GEOquerry)
# 4. Import into R/make it ready for DESeq2 (note how data was obtained(Salmon, featureCounts,...))
# 5. Exploratory data analysis
# 6. RNA-seq: Differential expression analysis (DE) of genes with DESeq2
# 7. Import annotations (various packages: annotation.dbi,  biomart, org.Hs.eg.db, org.Mm.eg.db )
# 8. GSEA - gene set enrichment analysis / ORA - over representation analysis: for identification KEGG, GO, REACTOME categories (clusterProfiler,...)
# 9. optional: transcripts to genes (packages tximport,tximeta, GenomicFeatures)
# 10. optional: alternative analysis with EBSeq (edgeR) and compare results with DESeq2


#GEOquery for importing from GEO OMNIBUS website. Imports a ?? object. Suitable for microarrays. Not suitable for DESeq2
#SummarizedExperiment is a class in R. consist of metadata and a dataframe? You can transform GEOquery object into SummarizedExperiment
#tximport or tximeta are for importing dataset (like Salmon) when you are given trainscipt IDS instead of Gene IDs