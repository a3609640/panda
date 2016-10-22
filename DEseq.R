####################
## 1 Preparations ##
####################
# set global chunk options and load the neccessary packages
source("https://bioconductor.org/biocLite.R")
library(genefilter)
biocLite("BiocStyle")
library(BiocStyle)
biocLite("rmarkdown")
library(rmarkdown)
library(geneplotter)
library(ggplot2)
library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Hs.eg.db)
library(AnnotationDbi)
biocLite("pathview")
library(pathview)
biocLite("gage")
library(gage)
biocLite("gageData")
library(gageData)

################################
## 2 Preparing count matrices ##
################################
# obtain the count table of the experiment directly from a pre-saved file. The RNA-seq was aligned to Hg38 by STAR
# read RNA-seq read data. 
setwd("~/Documents/Su Wu/Documents/Bioinformatics/RNA-seq")
testseq <- read.csv("testseq.csv")
# Use the column one (Ensemble names) as columnn names. 
testseq <- data.frame(testseq[,-1], row.names=testseq[,1])
# Remove the first four rows (N_unmapped,N_multimapping,N_noFeature and N_ambiguous)
testseq <- data.frame(testseq[c(-1,-2,-3,-4),])
par(mar=c(3,12,2,1))
boxplot(testseq, outline=FALSE, horizontal=TRUE, las=1)

#######################################
## Quality control of the count data ##
#######################################
## Remove rows in which there are no reads or nearly no reads
guideData <- testseq[rowSums(testseq)>1,]
head(guideData)
dim(guideData)
## Lets see how the data looks in a box plot
par(mar=c(3,12,2,1))
boxplot(guideData, outline=FALSE, horizontal=TRUE, las=1)
## create a design for our "modelling" 
guideDesign <- data.frame(row.names = colnames(guideData),
                          condition = c(rep("Mock",4),rep("siNeg",4),rep("siBF1",4),rep("ASO-neg",4),rep("ASO-1",4),rep("ASO-4",4)))
## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
dds <- DESeqDataSetFromMatrix(countData = guideData,colData = guideDesign,design = ~ condition)
dds
head(assay(dds))

##############################################
## Normalization: estimation of size factor ##
##############################################
## defines a virtual reference sample by taking the median of each gene’s values across samples 
## then computes size factors as the median of ratios of each sample to the reference sample
DESeq2Table <- estimateSizeFactors(dds)
sizeFactors(DESeq2Table)

################################################################
## Differential expression analysis: estimation of dispersion ##
################################################################
# use a Negative–Binomial-distribution (NB) to obtain an estimate of the dispersion parameter for each gene
DESeq2Table <- estimateDispersions(DESeq2Table)
plotDispEsts(DESeq2Table)

####################################################
## Statistical testing of Differential expression ##
#####################################################
# perform the statistical testing for differential expression and extract its results
DESeq2Table <-  nbinomWaldTest(DESeq2Table)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
table(DESeq2Res$padj < 0.1)

#######################
## standard analysis ##
#######################
# DESeq function performs a default analysis through the steps:
# (1) estimation of size factor: estimateSizeFactors
# (2) estimation of dispersion: estimateDispersions
# (3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
ddsDE <- DESeq(dds)
ddsres <- results(ddsDE)
summary(ddsres)

###############################################
## Extracting differentially expressed genes ## 
###############################################
## to make MA-plot from base means
table(ddsres$padj < 0.1)
plotMA(ddsres, main="DESeq2", ylim=c(-2,2))

###########################################################
## Regularized log transformation for clustering and PCA ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# Regularized log transform is to stabilize the variance of the data and to make its distribution roughly symmetric
# Note that an important argument for this function is blind (TRUE by default). 
# The default “blinds” the normalization to the design. 
# This is very important so as to not bias the analyses (e.g. class discovery)
rld=rlog(dds,blind=TRUE)

# Hierarchical clustering using rlog transformation
dists=dist(t(assay(rld)))
plot(hclust(dists), labels=guideDesign$condition)

## Finally lets do some more data separation
## Principal Component analysis
plotPCA(rld, intgroup = c("condition"))
data <- plotPCA(rld, intgroup = c( "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
# We can then use this data to build up the plot, specifying that the color of the points should reflect dexamethasone treatment 
# and the shape should reflect the cell line.
qplot(PC1, PC2, color=condition, shape=condition, data=data) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

############################################################
##Plot of normalized counts for a single gene on log scale##
############################################################
SREBF1 <- plotCounts(dds, "ENSG00000072310", intgroup = "condition", normalized = TRUE, transform = TRUE, "SREBF1", xlab = "group")
SREBF2 <- plotCounts(dds, "ENSG00000198911", intgroup = "condition", normalized = TRUE, transform = TRUE, "SREBF2", xlab = "group")
SCD <- plotCounts(dds, "ENSG00000099194", intgroup = "condition", normalized = TRUE, transform = TRUE, "SCD", xlab = "group")
FASN <- plotCounts(dds, "ENSG00000169710", intgroup = "condition", normalized = TRUE, transform = TRUE, "FASN", xlab = "group")
YY1 <- plotCounts(dds, "ENSG00000100811", intgroup = "condition", normalized = TRUE, transform = TRUE, "YY1", xlab = "group")
JUN <- plotCounts(dds, "ENSG00000177606", intgroup = "condition", normalized = TRUE, transform = TRUE, "JUN", xlab = "group")

#######################################
## Gene ontology enrichment analysis ##
#######################################
## get average gene expressions for each of the genes and then find non–DE genes that show a similar expression as the DE–genes. 
## These genes are then our background.
overallBaseMean <- as.matrix(ddsres[, "baseMean", drop = F])
sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))
backG <- c()
for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
}
backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]
## remove DE genes from background and the get the total number of genes in the background
backG <- setdiff(backG,  anSig$ENSEMBL)
length(backG)
## Plotting the density of the average expressions, 
## shows that the background matching has worked reasonably well.
multidensity( list( 
  all= log2(DESeq2Res[,"baseMean"]) ,
  foreground =log2(DESeq2Res[anSig$ENSEMBL, "baseMean"]), 
  background =log2(DESeq2Res[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

###################
## Running topGO ##
###################

## extract our significant genes and their annotation
sigGenes <- rownames(subset(ddsres, padj < 0.1))
anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys=rownames(ddsres), 
                              columns=c("SYMBOL","SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

head(anSig, 5)

onts = c( "MF", "BP", "CC" )
geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(anSig$ENSEMBL,  backG) 
inSelection =  geneIDs %in% anSig$ENSEMBL 
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]


tab = as.list(onts)
names(tab) = onts
for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results, “GenTable” produces a table of significant GO categories
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = 200)
  
}

topGOResults <- rbind.fill(tab)
write.csv(topGOResults, file = "topGOResults.csv")

#############################
#############################
## Compare siNeg and siBF1 ##
#############################
#############################
## generate dataset for siRNA treatment
testsiRNA <- data.frame(testseq[,c(5,6,7,8,9,10,11,12)])

par(mar=c(3,12,2,1))
boxplot(testsiRNA, outline=FALSE, horizontal=TRUE, las=1)

## Prefiltering: by removing rows in which there are no reads or nearly no reads
guideDatasiRNA <- testsiRNA[rowSums(testsiRNA)>1,]
head(guideDatasiRNA)
dim(guideDatasiRNA)
## Lets see how the data looks in a box plot
par(mar=c(3,12,2,1))
boxplot(guideDatasiRNA, outline=FALSE, horizontal=TRUE, las=1)

# Time to create a design for our "modelling" 
guideDesignsiRNA <- data.frame(row.names = colnames(guideDatasiRNA),
                          condition = c(rep("siNeg",4),rep("siBF1",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddssiRNA <- DESeqDataSetFromMatrix(countData = guideDatasiRNA,colData = guideDesignsiRNA,design = ~ condition)
ddssiRNA

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEsiRNA <- DESeq(ddssiRNA)
plotMA(ddsDEsiRNA, main="DESeq2", ylim=c(-2,2))
ddsressiRNA <- results(ddsDEsiRNA)
plotMA(ddsressiRNA,main="DESeq2", ylim=c(-2,2))

res = ddsressiRNA[order(ddsressiRNA$pvalue),]
res <- data.frame(res)
summary(res)


###########################################################
## Regularized log transformation for clustering and PCA ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# Note that an important argument for this function is blind (TRUE by default). 
# The default “blinds” the normalization to the design. 
# This is very important so as to not bias the analyses (e.g. class discovery)
rldsiRNA=rlog(ddssiRNA,blind=TRUE)
topVarGenes <- head( order( rowVars( assay(rldsiRNA) ), decreasing=TRUE ), 50 )
heatmap.2( assay(rldsiRNA)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", margins=c(2,10),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rldsiRNA)$condition ] )

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
columns(org.Hs.eg.db)
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")
summary(res)
head(res, 10)

############################
## KEGG pathways analysis ##
############################
## kegg.sets.hs is a named list of 229 elements
## Each element is a character vector of member gene Entrez IDs for a single KEGG pathway
data(kegg.sets.hs)
## sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
data(sigmet.idx.hs)
## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of sinaling and metabolic pathways only.
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)

## Generate a named vector of fold changes, where the names of the values are the Entrez gene IDs, for the gage() function
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=5) %>%
  .$id %>%
  as.character()

keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
# plot multiple pathways (plots saved to disk and returns a throwaway list object)
detach("package:dplyr",unload=TRUE)
pv.out.list <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
str(pv.out.list)

## map metabolites on a large scale using pathway.id= '01100'
## Cancer  - Homo sapiens (human) pathway.id= '05200'
## Melanoma - Homo sapiens (human) pathway.id= '05218'
pv.out <- pathview(gene.data = foldchanges, pathway.id= 'hsa05218', species = "hsa", kegg.native = T, same.layer = T)


##############################

deseq2.res <- results(ddsDEsiRNA)
deseq2.fc=deseq2.res$log2FoldChange
names(deseq2.fc)=rownames(deseq2.res)
exp.fc=deseq2.fc
out.suffix="deseq2"

require(gage)
data(kegg.gs)
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
#view first 3 pathways as demo
pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(gene.data = exp.fc, pathway.id = pid,species = "hsa", out.suffix=out.suffix))



dev.print(pv.out.list)
