library(stringr)
library (DESeq2)
library(genefilter)
library("RColorBrewer")
library(ReportingTools)
library(biomaRt)
library(gplots)
library(GOstats)
library(hwriter)

################################### Creating the raw Count table ########################################
# read in the count table
# be aware that DESeq2 analyzes two conditions against each other! So split your count table to include the two conditions you need
# this is what is called here FirstComp ??? a table of multiple columns which looks like that:
featureName  salm1	salm2	salm3	ctrl1	ctrl2	ctrl3
128up-RA_1	27	12	12	14	9	2
128up-RA_2	63	26	47	46	12	28
128up-RA_3	239	99	244	261	45	299
128up-RA_4	95	21	53	61	12	66
14-3-3epsilon-RA_1	721	272	216	511	77	184
14-3-3epsilon-RA_2	1917	1230	1286	2548	363	1649
???


################################### DESeq2 Analysis ########################################
#trying the DESEq2 package with no replica
#################### from here, each comparison will be ran separately. ########################

colData <- colData1 <- data.frame(row.names=names(FirstComp), condition = c("salm", "leg"))

# the metadData should  look like that:
sample	condition
name	salm1	salm
name	salm2	salm
name	salm3	salm
name	ctrl1	ctrl
name	ctrl2	ctrl
name	ctrl3	ctrl
###!####to change when changing comparison
Comp 
<- tenthComp ###!####to change when changing comparison

cds <- DESeqDataSetFromMatrix (
  countData = Comp,
  colData   = colData,  
  design    = ~
    condition
)

fit = DESeq(cds)

saveRDS( fit, paste("DEgenes_fitObject_", as.character(colData(cds)[[1]][1]), "vs", as.character(colData(cds)[[1]][3]) , ".RDS", sep="") )

#extracting results
res = results(fit)

#independent filtering
# create bins using the quantile function
qs<-c(0,quantile( res$baseMean[res$baseMean>0],0:7/7) )
# "cut" the genes into the bins
bins<-cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins)<-paste0("~",round(.5*qs[-1]+.5*qs[-length(qs)]))
# calculate the ratio of $p$ values less than .01 for each bin
ratios<-tapply( res$pvalue, bins,function(p)mean( p<.01,na.rm=TRUE) )
# plot these ratios
#barplot(ratios,xlab="mean normalized count",ylab="ratio of small $p$ values")

cat("\n\n creating diagnostic plots \n\n")
pdf(paste("diagnosticPlots_", as.character(colData(cds)[[1]][1]), "_", as.character(colData(cds)[[1]][3]), "_featureCounts.pdf", sep=""), onefile=T)
DESeq2::plotMA( fit , ylim=c(-1,1), main="MA plot")
plotDispEsts( fit , main="estimating dispersion factors")
hist( res$pvalue,breaks=20,col="grey", main="histogram of results p-value distribution", xlab="p-values")
barplot(ratios,xlab="mean normalized count",ylab="ratio of small $p$ values")
plot(attr(res,"filterNumRej"),type="b",xlab="quantiles of'baseMean'",ylab="number of rejections")
dev.off()

cat("\n\n doing biomaRt analysis\n\n ")
res$FBgen<-sapply(strsplit(rownames(res),split="\\+"),"[",1)
ensembl = useMart( "ensembl", dataset = "dmelanogaster_gene_ensembl" )
genemap <- getBM( attributes = c("flybase_gene_id", "ensembl_gene_id", "entrezgene", "flybasename_gene", "external_gene_id"), 
                  filters = "flybase_gene_id", 
                  values = as.character(res$FBgen), 
                  mart = ensembl )
idx <- match( res$FBgen, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$gene_symbol <- genemap$flybasename_gene[ idx ]

resSig<-res[c(which(res$log2FoldChange>=1) ,which(res$log2FoldChange <= -1)), ] #significant genes according to the FC value (as no significance can be calculated!)

cat("\n\n calculating regularized-logarithm transformation \n\n")
rld<-rlog( fit )
rownames(rld) <- res$gene_symbol
##euclDist
sampleDists <- dist( t( assay(rld) ) )
sampleDists
##DistHeatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$condition, sep="-" )
colnames(sampleDistMatrix) <- NULL   
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste("clusteringPlots_rlogtransform_", as.character(colData(cds)[[1]][1]), "_", as.character(colData(cds)[[1]][3]), ".pdf", sep=""), onefile=TRUE)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
##PCA
ramp <- 1:2/2
cols <- c(rgb(ramp, 0, 0))
print( plotPCA( rld, intgroup = c( "condition"), col=cols ) )
##topVarGenes - geneHeatmap
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

# exporting the tables with raw and normalized data
RawCounts <- counts(fit)
resMerged <- merge(res, RawCounts, by=0, all.x=T, all.y=F)
rownames(resMerged) <- resMerged[,1]
#resMerged <- resMerged[,c(8,9,10,3,11,12)]
resMerged <- resMerged[,-c(1,2,4,5,6,7,8)] 
resMerged <- resMerged[, c(2,3,1,4:dim(resMerged)[2])]

write.table(as.data.frame(resMerged), file=paste('DESeq2_analysisResults_', as.character(colData(cds)[[1]][1]), "_", as.character(colData(cds)[[1]][3]), "_allGenes.txt", sep=""), sep="\t")