#browseVignettes("SummarizedExperiment")
#browseVignettes("GEOquery")

#GSE165111 - this is a microarray, not the focus of this course
GEOid="GSE165111"
gse = getGEO(GEOid)[[1]]

utils::browseURL(getwd()) #opens the download directory which is based on your project directory #opens the download directory which is based on your project directory

se = as(gse, "SummarizedExperiment")
se
assays(se)$exprs%>%head
colData(se)
metadata(se)

getwd() #setwd() to change it
if(GEOid %in% list.files()){"you already have it"}else{getGEOSuppFiles(GEOid)} #get a tar file full of subfiles, look for it in your working directory: 

#GSE151190 - this is RNA sequencing, a.k.a. high throughput sequencing, the focus of the course
GEOid2="GSE151190"
gse<-getGEO(GEOid2,GSEMatrix =TRUE)[[1]]
se = as(gse, "SummarizedExperiment")
se
assays(se)$exprs     #does not work, we get sample list instead of raw data
metadata(se)
assays(se)[[1]]

if(GEOid2 %in% list.files()){"you already have it"}else{getGEOSuppFiles(GEOid2)} #now import it the old fashioned way into R (read.table), no efficent way to do it :((

utils::browseURL(paste(getwd(),GEOid2,sep="\\")) #opens the download directory 


#this should work if there is only file and the separator is \t, another common option is a comma
files<-list.files(GEOid2,full.names=T)
length(files)
df<-read.table(files[1], skip=0, sep = "\t", header=TRUE, row.names = 1)
df%>%rownames
df%>%head
#this file is same the ftp download from geo omnibus