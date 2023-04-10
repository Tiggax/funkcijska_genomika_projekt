#browseVignettes("edgeR")
#edgeRUsersGuide() which is same as: https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

#recall we have prepared DESeq data sets already
# Restore the list of DESeq data sets
ddslist<-readRDS(file = "dds.rds")
#example of accesing it:
ddslist$genotype.ch1@metadata
ddslist$genotype.ch1@colData
counts(ddslist$genotype.ch1)

#creating edgeR data set
count_data<-as.data.frame(counts(ddslist$genotype.ch1))
groups<-ddslist$genotype.ch1@colData$genotype.ch1

d <- DGEList(counts=count_data, group=groups)


#-----------


#Filtration
cpm <- cpm(d)
keep <- rowSums(cpm>1) >= ncol(cpm)/2 # adjust for your dataset
# Here, a CPM value of 1 means that a gene is expressed if it has at least 20 counts (with a total of 20 million reads)
d <- d[keep, ]
d$samples$lib.size <- colSums(d$counts)


#Normalizing the data

# normalisation by the method of trimmed mean of M-values (TMM) is performed using the calcNormFactors
d <- calcNormFactors(d, method= "TMM")
d$samples$norm.factors




#Expression distribution of samples for unnormalised and normalised data (link)
# https://bioconductor.github.io/BiocWorkshops/rna-seq-analysis-is-easy-as-1-2-3-with-limma-glimma-and-edger.html


#Data Exploration
plotMDS(d, method="bcv")
col=as.numeric(d$samples$group)
legend("bottomleft",
       as.character(unique(d$samples$group)), col=1:3, pch=20)


#Estimating the Dispersion
d1 <- estimateCommonDisp(d, verbose=T)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)


#Differential Expression
#The exact test is only applicable to experiments with a single factor.
unique(groups)
et12 <- exactTest(d1, pair=unique(groups))
topTags(et12, n=10)
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)



design <- model.matrix(~groups) 

y <- estimateDisp(d,design) 

#To perform quasi-likelihood F-tests: 
  fit <- glmQLFit(y,design) 
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf) 


#To perform likelihood ratio tests: 
  fit <- glmFit(y,design) 
lrt <- glmLRT(fit,coef=2) 
topTags(lrt)

