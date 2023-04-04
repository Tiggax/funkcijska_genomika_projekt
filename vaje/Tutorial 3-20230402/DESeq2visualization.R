#run setup
#recall we have prepared DESeq data sets already
# Restore the list of DESeq data sets
ddslist<-readRDS(file = "dds.rds")
#example of accesing it:
ddslist$genotype.ch1@metadata
ddslist$genotype.ch1@colData
counts(ddslist$genotype.ch1)



# Exploratory analysis and visualization
# browseVignettes("DESeq2")
# browseURL("http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html")
# DESeq2 has implemented two methods for transforming the data into a form which is more suitable for visualization of data.
# variance-stabilizing transformation (VST)
# regularized logarithm (rlog))
# 
# Methods for visualizing the data
# MDS
# PCA
# Dendrogram
# heatmap of the sample-to-sample distances


dds<-ddslist$`characteristics_ch1.1+genotype.ch1`
#$characteristics_ch1.1+genotype.ch1` # ddslist$genotype.ch1, ddslist$`characteristics_ch1.1+genotype.ch1`, ddslist$source_name_ch1
nrow(dds)

# filtering (other criteria can be used as well)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

# transformation / Why? What are other options?
# transformation option 1-vsd
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

# transformation option 2-rlog
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

## show the effect of transformation
dds2 <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds2, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


## sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists

# visualize the distances in a heatmap
vsd$title
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <-  paste( vsd$genotype.ch1,vsd$characteristics_ch1.1, sep = " - " )
#rownames(sampleDistMatrix) <- vsd$title
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)




# PCA plot

#plotPCA(vsd, intgroup = "title")
plotPCA(vsd, intgroup = c("genotype.ch1","characteristics_ch1.1"))

## MDS
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))


ggplot(mds, aes(x = `1`, y = `2`, color = genotype.ch1, shape = characteristics_ch1.1)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

ggplot(mds, aes(x = `1`, y = `2`, color = title)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# check how you can show two variables are better for the analysis

  
  
