source("http://bioconductor.org/biocLite.R")
# biocLite("reactome.db")
library(reactome.db)
library(RcppArmadillo)
library(colorspace)
library(lattice)
library(RODBC)
library(Matrix)
library(survival)
library(Rcpp)
library(genefilter)
library(BiocStyle)
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(DESeq2)
library(RColorBrewer)
library(stringr)
library(biomaRt)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(gage)
library(gageData)
library(Biobase)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(SummarizedExperiment)
#######





################################
## 2 Preparing count matrices ##
################################
# obtain the count table of the experiment directly from a pre-saved file. The RNA-seq was aligned to Hg38 by STAR
# read RNA-seq read data. 
#setwd("~/Documents/Bioinformatics analysis/RNA-seq/Testrun/STAR/results")
testseq <- read.csv("~/Documents/Bioinformatics analysis/RNA-seq/Testrun/STAR/results/testseq.csv")
# Use the column one (Ensemble names) as columnn names. 
testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
# Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
testseq <- data.frame(testseq[c(-1,-2,-3,-4),])
par(mar=c(3,12,2,1))
boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)

###############################
###############################
## Compare ASO-Neg and ASO-4 ##
###############################
###############################
## generate dataset for siRNA treatment
testASO <- data.frame(testseq[,c(13,14,15,16,21,22,23,24)])

par(mar=c(3,12,2,1))
boxplot(testASO, outline=FALSE, horizontal=TRUE, las=1)

## Prefiltering: by removing rows in which there are no reads or nearly no reads
guideDataASO <- testASO[rowSums(testASO)>1,]
head(guideDataASO)
dim(guideDataASO)
## Lets see how the data looks in a box plot
par(mar=c(3,12,2,1))
boxplot(guideDataASO, outline=FALSE, horizontal=TRUE, las=1)

# Time to create a design for our "modelling" 
guideDesignASO <- data.frame(row.names = colnames(guideDataASO),
                             condition = c(rep("ASO-Neg",4),rep("ASO-4",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddsASO <- DESeqDataSetFromMatrix(countData = guideDataASO,colData = guideDesignASO,design = ~ condition)
ddsASO

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEASO <- DESeq(ddsASO)
resASO<-results(ddsDEASO)
resASO<-resASO[order(resASO$log2FoldChange),]
## An MA-plot21 provides a useful overview for an experiment with a two-group comparison
## The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis 
## (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene is represented with a dot. 
## Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
plotMA(ddsDEASO, main="DESeq2", ylim=c(-2,2))
## The red points indicate genes for which the log2 fold change was significantly higher than 0.5 or less than -0.5 
# (treatment resulting in more than doubling or less than halving of the normalized counts) with adjusted p value less than 0.1.
resLFC1 <- results(ddsDEASO, lfcThreshold=0.5)
table(resLFC1$padj < 0.1)
plotMA(resLFC1,main="DESeq2", ylim=c(-2,2))

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
resASO <- data.frame(resASO)
columns(org.Hs.eg.db)
resASO$symbol = mapIds(org.Hs.eg.db,
                       keys=row.names(resASO), 
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")
resASO$entrez = mapIds(org.Hs.eg.db,
                       keys=row.names(resASO), 
                       column="ENTREZID",
                       keytype="ENSEMBL",
                       multiVals="first")
resASO$name =   mapIds(org.Hs.eg.db,
                       keys=row.names(resASO), 
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
summary(resASO)
head(resASO, 10)


########################################################
## Gene-set enrichment analysis for ASO-4 and ASO-neg ##
########################################################

## This database works with Entrez IDs, 
## so we will need the entrezid column that we added earlier to the res object.
## First, we subset the results table, res, to only those genes for which the Reactome database has data 
## (i.e, whose Entrez ID we find in the respective key column of reactome.db and 
## for which the DESeq2 test gave an adjusted p value that was not NA.
## resASO$entrez <- convertIDs( row.names(resASO), "ENSEMBL", "ENTREZID", org.Hs.eg.db )
res2 <- resASO[ resASO$entrez %in% keys( reactome.db, "ENTREZID" ) &
               !is.na( resASO$pvalue) , ]
head(res2)
suppressWarnings( reactomeTable <- AnnotationDbi::select( reactome.db, 
                                                          keys=res2$entrez,
                                                          keytype="ENTREZID", 
                                                          columns=c("ENTREZID","REACTOMEID") ) )
head(reactomeTable)

incm <- do.call( rbind, with(reactomeTable, 
                             tapply(ENTREZID, 
                                    factor(REACTOMEID), 
                                    function(x) res2$entrez %in% x ) ))
colnames(incm) <- res2$entrez
str(incm)

incm <- incm[ rowSums(incm) >= 5, ]

testCategory <- function( reactomeID ) {
  isMember <- incm[ reactomeID, ]
  data.frame(
    reactomeID = reactomeID,
    numGenes = sum( isMember ),
    avgLFC = mean( res2$log2FoldChange[isMember] ),
    strength = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
    pvalue = t.test( res2$log2FoldChange[ isMember ] )$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]] ) }

testCategory("109581")

reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )












## Using select, a function from AnnotationDbi for querying database objects, 
## we get a table with the mapping from Entrez IDs to Reactome Path IDs :
reactomeTable <- AnnotationDbi::select( reactome.db, 
                                        keys=as.character(res2$entrez), keytype="ENTREZID", 
                                        columns=c("ENTREZID","REACTOMEID") )
head(reactomeTable)
## The next code chunk transforms this table into an incidence matrix. 
## This is a Boolean matrix with one row for each Reactome Path and one column for each unique gene in res2, 
## which tells us which genes are members of which Reactome Paths.
incm <- do.call( rbind, with(reactomeTable, tapply( 
  ENTREZID, factor(REACTOMEID), function(x) res2$entrez %in% x ) ))
colnames(incm) <- res2$entrez
str(incm)
## We remove all rows corresponding to Reactome Paths with less than 20 or more than 80 assigned genes.
within <- function(x, lower, upper) (x>=lower & x<=upper)
incm <- incm[ within(rowSums(incm), lower=20, upper=80), ]

## To test whether the genes in a Reactome Path behave in a special way in our experiment, 
## we calculate a number of statistics, including a t-statistic to see whether the average of the genes’ log2 fold change values in the gene set is different from zero. 
## To facilitate the computations, we define a little helper function:
testCategory <- function( reactomeID ) {
  cat("reactomeId: ", reactomeID, "\n")
  cat("length of incm: ", length(incm), "\n")
  cat("nrow of incm: ", nrow(incm), "\n")
  cat("ncol of incm: ", ncol(incm), "\n")
  cat("incm[0,0]: ", incm[0,0], "\n")
  cat("incm[1,1]: ", incm[1,1], "\n")
  #cat("incm[x,x+1]: ", incm[reactomeID, reactomeID+1], "\n")
  cat("incm[x,]: ", incm[reactomeID,], "\n")
  #isMember <- incm[ reactomeID, ]
  isMember <- FALSE
  data.frame( 
    reactomeID  = reactomeID,
    numGenes    = sum( isMember ),
    avgLFC      = mean( res2$log2FoldChange[isMember] ),
    sdLFC       = sd( res2$log2FoldChange[isMember] ),
    zValue      = mean( res2$log2FoldChange[isMember] ) /sd( res2$log2FoldChange[isMember] ),
    strength    = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
    pvalue      = t.test( res2$log2FoldChange[ isMember ] )$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]],
    stringsAsFactors = FALSE ) }
## The function can be called with a Reactome Path ID:
testCategory("109606")

reactomeResult <- do.call( rbind, lapply( rownames(incm), testCategory ) )
reactomeResult$padjust <- p.adjust( reactomeResult$pvalue, "BH" )
reactomeResultSignif <- reactomeResult[ reactomeResult$padjust < 0.05, ]
head( reactomeResultSignif[ order(-reactomeResultSignif$strength), ] )
