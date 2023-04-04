# 5. Gene Ontology Analysis  

# Gene Ontology defines concepts/classes used to describe gene function, and relationships between these concepts. It classifies functions along three aspects:
#   
# MF: Molecular Function
#     molecular activities of gene products
# CC: Cellular Component
#     where gene products are active
# BP: Biological Process
#     pathways and larger processes made up of the activities ofmultiple gene products
# http://amigo.geneontology.org/amigo/dd_browse
# http://geneontology.org/
# GO terms are organized in a directed acyclic graph, where edge between the terms represent parent-child relationship.

#browseURL("http://yulab-smu.top/clusterProfiler-book/")
#5.2


database<-org.Hs.eg.db #org.Mm.eg.db, org.Hs.eg.db
ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`characteristics_ch1.1+genotype.ch1`

#from DiffExprAnal we found these genes to be differentially expressed. 
dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype.ch1","wild.type","Mutant.DNAJC6"))

# dfannotation2<-mapIds(x=database,
#                       keys = rownames(counts(dds)),
#                       column = "ENTREZID",
#                       keytype = "ALIAS",
#                       multiVals = "first")




res<-res[!is.na(res$padj),]
SigExpGenes<-res[(res$padj < 0.05) & (abs(res$log2FoldChange)>1),] 
rownames(SigExpGenes)%>%length #significally expressed genes
rownames(counts(dds))%>%length #all tested genes

ggo1 <- groupGO(gene     = rownames(SigExpGenes),
               OrgDb    = database,
               keyType = "SYMBOL",
               ont      = "BP",
               level    = 1,
               readable = F)
#SMYBOL requires  readable = F

ggo1[,1:4]

ggo2 <- groupGO(gene     = rownames(SigExpGenes),
                OrgDb    = database,
                keyType = "SYMBOL",
                ont      = "BP",
                level    = 2,
                readable = F)

ggo2[,1:4]

#5.3

#?enrichGO

ego <- enrichGO(gene          = rownames(SigExpGenes),
                universe      = rownames(counts(dds)),
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = F)
ego[,-8]%>%head
dfego<-as.data.frame(ego[,-8])
dfego%>%head

#which BPs were found expressed? i) pvalueCutoff > p-value, ii) pvalueCutoff > adjusted p-value and iii) qvalueCutoff > q-value.
dfego[dfego$p.adjust<ego@pvalueCutoff & dfego$pvalue<ego@pvalueCutoff & dfego$qvalue<ego@qvalueCutoff,2]

barplot(ego, showCategory = 20)

goplot(ego)

#library(ggnewscale)

#In order to consider the potentially biological complexities in which a gene may belong to multiple annotation categories and provide information of numeric changes if available, we developed cnetplot function to extract the complex association.

cnetplot(ego, categorySize = "pvalue", foldChange = rownames(SigExpGenes))


