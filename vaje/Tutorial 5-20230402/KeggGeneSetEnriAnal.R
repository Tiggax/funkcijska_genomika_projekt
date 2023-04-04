



#browseURL("http://yulab-smu.top/clusterProfiler-book/")

database<-org.Hs.eg.db #org.Mm.eg.db, org.Hs.eg.db
ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`characteristics_ch1.1+genotype.ch1`

#from DiffExprAnal we found these genes to be differentially expressed. 
dds <- DESeq(dds)
dds$genotype.ch1
res <- results(dds, contrast=c("genotype.ch1","wild.type","Mutant.DNAJC6"))
res

res<-res[!is.na(res$padj),]
SigExpGenes<-res[(res$padj < 0.05) & (abs(res$log2FoldChange)>1),] 
rownames(SigExpGenes)%>%length #significally expressed genes
rownames(counts(dds))%>%length #all tested genes including NAs
res%>%nrow #all tested genes without NAs

#6. KEGG analysis https://www.genome.jp/kegg/kegg3a.html

search_kegg_organism('mmu', by='kegg_code')
search_kegg_organism('sapiens', by='scientific_name')

#?enrichKEGG

SigExpGenes$ENTREZID<-AnnotationDbi::mapIds(x=database,
                       keys = rownames(SigExpGenes),
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
 
SigExpGenes%>%head

kk <- enrichKEGG(gene         = SigExpGenes$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05)

kkres<-as.data.frame(kk@result)
kkres%>%head(10)
significantkkres<-kkres[kkres$pvalue<kk@pvalueCutoff & kkres$p.adjust<kk@pvalueCutoff & kkres$qvalue < kk@qvalueCutoff,2,drop=F]
significantkkres



#same as above
head(kk)[,2,drop=F]


# visualize pathways

#opens up web page with photo of pathway for our data
#colors represent expression
#browseKEGG(kk, "hsa04512")
#links<-sapply(rownames(head(kk)),function(x){browseKEGG(kk, x)})
#links

# picture generated in R
#library("pathview")
geneList<-SigExpGenes$log2FoldChange
names(geneList)<-SigExpGenes$ENTREZID

download.pathway<-function(keggid){
folder<-"pics"
if (!("pics" %in% list.files())) {dir.create(folder)}
setwd(folder)
pic <- pathview(gene.data  = geneList,
                     pathway.id = keggid,
                     species    = stringr::str_extract(keggid,"[A-z]+"),
                     limit      = list(gene=max(abs(geneList)), cpd=1))
setwd("..")
}
rownames(significantkkres)
lapply(rownames(significantkkres),download.pathway)

#download.pathway("hsa04110")

#
#hsa04110 <- pathview(gene.data  = geneList,
# pathway.id = "hsa04512",
# species    = "hsa",
# limit      = list(gene=max(abs(geneList)), cpd=1))


#---------------

#6.3 KEGG Module over-representation test

mkk <- enrichMKEGG(gene         = SigExpGenes$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 pvalueCutoff = 0.05)
head(mkk)
mkkres<-as.data.frame(mkk@result)

