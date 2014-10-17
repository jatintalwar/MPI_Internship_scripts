######script for DeSeq analysis for C.elegans data

library(stringr)
source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
library (DESeq2)
library(DESeq)
biocLite("genefilter")
library(genefilter)
library("RColorBrewer")
library(ReportingTools)  ## no such package as "ReportingTools"
biocLite("biomaRt")
library(biomaRt)
install.packages('gplots')
library(gplots)
biocLite("GOstats")
library(GOstats)
install.packages('hwriter')
library(hwriter)

###############    reading the count table : ###############

count_table <- read.table("count_table_complete_HTSeq_tophat.txt",header=T,row.names=1)
View(count_table)
count_table_metadata <- data.frame(row.names=colnames(count_table),condition=factor(c("control","control","control","Gld4_KO","Gld4_KO","Gld4_KO"),levels=c("control","Gld4_KO")),libType=c("Single-end","Single-end","Single-end","Single-end","Single-end","Single-end"))
condition
count_table_metadata

dds <- DESeqDataSetFromMatrix(countData = count_table,colData = count_table_metadata, design = ~condition)
dds
dds <- DESeq(dds)
res <- results(dds)
View(res)
write.table(res,file= "C.elegans_control_vs_GLD4_Deseq_analysis.txt")
resOrdered <- res[order(res$log2FoldChange)]
head(resOrdered)

plotMA(res, main="DESeq2_Analysis_C.elegans_Data", ylim=c(-2,2))
colnames(res)


###################      heatmap     #########################

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
hmcol
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale="row",dendrogram="none", trace="none", margin=c(10,6),labCol=NULL)

?heatmap.2


############### to extract transfrormed values ##############
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)


#### pca analysis #########
print(plotPCA(rld, intgroup=c("condition")))
