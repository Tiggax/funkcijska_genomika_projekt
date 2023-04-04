ddslist<-readRDS(file = "dds.rds")
dds<-ddslist$`characteristics_ch1.1+genotype.ch1`
database<-org.Hs.eg.db

browseVignettes("biomaRt")
# https://www.ensembl.org/biomart/martview/
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
datasets[grep("[mM]ouse",datasets$description),]
datasets[grep("[Hh]uman",datasets$description),]


ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
filters[1:10,]
attributes = listAttributes(ensembl)
attributes[1:50,1]
attributes[grep("symbol",attributes$name),]
attributes[grep("entrez",attributes$name),]

annotation1<-AnnotationDbi::select(x=database,
                      keys = rownames(counts(dds)),
                      column = c("ENTREZID","SYMBOL"),
                      keytype = "SYMBOL",
                      multiVals = "first")
annotation1%>%dim()
length(rownames(counts(dds)))

annotation2<-getBM(attributes=c("ensembl_gene_id","description","chromosome_name", 'entrezgene_id',"uniprot_gn_symbol","hgnc_symbol"), 
                   filters = "entrezgene_id", 
                  values =  annotation1$ENTREZID, 
                  mart = ensembl)

?getBM

annotation2%>%head(20)
class(annotation1$ENTREZID)
class(annotation2$entrezgene_id)
annotation2$entrezgene_id<-as.character(annotation2$entrezgene_id)
annotation3<-annotation1%>%inner_join(annotation2,by=c("ENTREZID" = "entrezgene_id"))
annotation3%>%head(10)
#----------------------------
#Examine boxplot of the counts for each sample and of the normalized counts
counts(dds)%>%head
df<-as.data.frame(counts(dds))
boxplot(x = as.list(df))

#library(reshape2)

#changing names to fit plot
newnames<-paste(str_extract(dds$genotype.ch1,"^[^\\.]+"),str_sub(dds$timepoint.ch1,1,1),str_extract(dds$timepoint.ch1,"[^ ]+$"),sep="")

colnames(df)<-newnames

#non-normalized boxplots arent that great
ggplot(data = melt(df), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

#vsd normalization and boxplot
df2<-as.data.frame(assay(vsd))
colnames(df2)<-newnames

ggplot(data = melt(df2), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

#vsd transformation and boxplot
vsd <- vst(dds)
assay(vsd)%>%head
df2<-as.data.frame(assay(vsd))
colnames(df2)<-newnames

ggplot(data = melt(df2), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))

#rlog transformation and boxplot
rld <- rlog(dds)
assay(rld)%>%head
df3<-as.data.frame(assay(vsd))
colnames(df3)<-newnames

ggplot(data = melt(df3), aes(x=variable, y=value)) + geom_boxplot(aes(fill=variable))
