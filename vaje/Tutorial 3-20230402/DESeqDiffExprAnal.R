# Differential expression analysis
# chapter 5 in the url
# browseURL("http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html")

#load our DESeq2 data sets from before
ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`characteristics_ch1.1+genotype.ch1`

#5.1 Running the differential expression pipeline
dds <- DESeq(dds)
#5.2 Building the results table
res <- results(dds)
res

#you can only compare 2 factors at a time
#its up to you pick 2 and focus on these
res <- results(dds, contrast=c("genotype.ch1","wild.type","Mutant.DNAJC6"))
mcols(res, use.names = TRUE)
summary(res)

#recall p-values and adjusted p-values: https://www.youtube.com/watch?v=K8LQSvtjcEo

res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
#5.3 Other comparisons
dds$characteristics_ch1.1
results(dds, contrast = c("characteristics_ch1.1", "timepoint..DIV.0", "timepoint..DIV.4"))

#5.4 Multiple testing
# how many genes are differentially expressed based on p-value?
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))

# 5 % of 20514 = 1025,7
# 1025.7/1193 = 85 %
# If we just considered the list of genes with a p value below 0.05 as differentially expressed, this list should therefore be expected to contain up to 1580 / 5170 = 31% false positives.

# DESeq2 uses the Benjamini-Hochberg (BH) adjustment (Benjamini and Hochberg 1995) as implemented in the base R p.adjust function

# Hence, if we consider a fraction of 10% false positives acceptable, we can consider all genes with an adjusted p value below 10% = 0.1 as significant. How many such genes are there?

sum(res$padj < 0.1, na.rm=TRUE)
# We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation:
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

# …and with the strongest up-regulation:
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])


#6 Plotting results
#6.1Counts plot

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("genotype.ch1"))
?plotCounts

geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("genotype.ch1","characteristics_ch1.1"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = genotype.ch1, y = count, color = characteristics_ch1.1)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = genotype.ch1, y = count, color = characteristics_ch1.1, group = characteristics_ch1.1)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


#MA plot

# BiocManager::install("apeglm")
# library("apeglm")

resultsNames(dds)
res <- lfcShrink(dds, coef="genotype.ch1_wild.type_vs_Mutant.DNAJC6", type="apeglm")
DESeq2::plotMA(res, ylim = c(-5, 5))



# gene clustering

# library("genefilter")
vsd <- vst(dds, blind = FALSE)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[,c("characteristics_ch1.1","genotype.ch1")])
pheatmap(mat = mat, fontsize_row = 5)
colData(vsd)

