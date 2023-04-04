
ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`characteristics_ch1.1+genotype.ch1`

#library(AnnotationDbi) #part of biomaRt

database<-org.Hs.eg.db
#use appropriate database
#org.Mm.eg.db ->Mouse 
#org.Hs.eg.db ->Human
# GO.db, KEGG.db
#more databases at https://bioconductor.org/packages/3.12/data/annotation/

#check what GENE IDs can look like in your database
columns(database)
sampleIDS<-sapply(columns(database),function(x){head(keys(database,keytype=x),3)})
sampleIDS

#compare with your GENE ID and find the match
rownames(counts(dds))%>%head

#match(es) found

head(keys(database,keytype="ALIAS"),20)
head(keys(database,keytype="SYMBOL"),20)
rownames(counts(dds)) %>%length()
sum(rownames(counts(dds)) %in% keys(database,keytype="ALIAS"))
sum(rownames(counts(dds)) %in% keys(database,keytype="SYMBOL"))

#functions: mapIds, select

#select returns a dataframe
dfannotation1<-AnnotationDbi::select(x=database,
       keys = rownames(counts(dds)),
       column = c("ENSEMBL","ENTREZID"),
       keytype = "SYMBOL",
       multiVals = "first")
dfannotation1%>%head(10)

#mapIds returns a vector
dfannotation2<-mapIds(x=database,
                      keys = rownames(counts(dds)),
                      column = "ENSEMBL",
                      keytype = "ALIAS",
                      multiVals = "first")

dfannotation2%>%head(10)
#if you get many NAs, change to appropriate database!


#carefull! one-to-many relationships may and will arise 
dfannotation2[["ZNF446"]]
duplicates<-dfannotation2[duplicated(dfannotation2)]
duplicates[!is.na(duplicates)]

dfannotation2[grep("ENSG00000083838",dfannotation2)]





#  Volcano plot

#browseVignettes("EnhancedVolcano")

#annotate and change rownames to annotation

# recall
# dds2fresh <- DESeqDataSetFromMatrix(countData = df2,
#                                     colData = se2@colData,
#                                     design = ~ source_name_ch1,
#                                     metadata=metadata(se2))

dds<-ddslist$`source_name_ch1`
df2<-as.data.frame(counts(dds))
colData<-colData(dds)
metadata<-metadata(dds)
database<-org.Hs.eg.db

dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('source_name_ch1','CAFs.from.genetically.stable.oral.cancers','normal.human.oral.fibroblasts'))
res <- lfcShrink(dds,
                 contrast = c('source_name_ch1','CAFs.from.genetically.stable.oral.cancers','normal.human.oral.fibroblasts'), res=res, type = 'normal')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

#What if we want different GENEIDS on the volcano plot?

rownames(df2)%>%head
symbols <- AnnotationDbi::select(database, keys = rownames(df2),
                  column = c("SYMBOL"), keytype = "ENTREZID")
#df2joined<-cbind(df2,symbols) #may not work properly if order is switched at some point
df2%>%mutate(RowNames=rownames(.))%>%right_join(symbols,by=c("RowNames"="ENTREZID"))%>%filter(!is.na(SYMBOL))->dfjoined
dfjoined%>%head

dfjoined%>%select_if(is.numeric)->dfannotated
rownames(dfannotated)<-dfjoined[,"SYMBOL"]

dds2annotated <- DESeqDataSetFromMatrix(countData = dfannotated,
                                    colData = colData,
                                    design = ~ source_name_ch1,
                                    metadata=metadata)


dds <- dds2annotated
dds <- DESeq(dds, betaPrior=FALSE)

levels(dds$source_name_ch1)
res <- results(dds,
               contrast = c('source_name_ch1','CAFs.from.genetically.stable.oral.cancers','normal.human.oral.fibroblasts'))
res <- lfcShrink(dds,
                 contrast = c('source_name_ch1','CAFs.from.genetically.stable.oral.cancers','normal.human.oral.fibroblasts'), res=res, type = 'normal')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')




