####################
## 1 Preparations ##
####################
# set global chunk options and load the neccessary packages
chooseCRANmirror()

source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("BiocStyle")
biocLite("rmarkdown")
biocLite("DESeq2")
biocLite("pathview")
biocLite("gage")
biocLite("gageData")
biocLite("RcppArmadillo")
install.packages("https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.6.100.0.0.tar.gz", repos=NULL, type="source")

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
setwd("~/Documents/Bioinformatics analysis/RNA-seq/Testrun/STAR/results")
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
                          condition = c(rep("Mock",4),rep("siNegative",4),rep("siSREBF1",4),rep("ASO-Neg",4),rep("ASO-1",4),rep("ASO-4",4)))
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
# DESeq2Table <- estimateSizeFactors(dds)
# sizeFactors(DESeq2Table)

################################################################
## Differential expression analysis: estimation of dispersion ##
################################################################
# use a Negative–Binomial-distribution (NB) to obtain an estimate of the dispersion parameter for each gene
# DESeq2Table <- estimateDispersions(DESeq2Table)
# plotDispEsts(DESeq2Table)

####################################################
## Statistical testing of Differential expression ##
#####################################################
# perform the statistical testing for differential expression and extract its results
# DESeq2Table <-  nbinomWaldTest(DESeq2Table)
# DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH")
# table(DESeq2Res$padj < 0.1)

#######################
## standard analysis ##
#######################
# DESeq function performs a default analysis through the steps:
# (1) estimation of size factor: estimateSizeFactors
# (2) estimation of dispersion: estimateDispersions
# (3) Negative Binomial GLM fitting and Wald statistics: nbinomWaldTest
ddsDE <- DESeq(dds) ##一步到位，不需要太多步骤
ddsres <- results(ddsDE) ##得到结果
summary(ddsres)
res <- data.frame(ddsres)

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

## Principal Component analysis
plotPCA(rld, intgroup = c("condition"),ntop = 50000)
data <- plotPCA(rld, intgroup = c( "condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(data=data, aes_string(x="PC1", y="PC2", color="condition")) + 
      geom_point(size=3) + 
      theme_bw() + 
      xlim(-6, 6) + 
      ylim(-6, 6) +
      theme(text = black.bold.18.text, 
            axis.text = black.bold.18.text,
            axis.line.x = element_line(color="black", size=1),
            axis.line.y = element_line(color="black", size=1),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size=1),
            panel.background = element_blank(),
            legend.position=c(1,0),
            legend.justification=c(1.02,-0.02),) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) 

############################################################
## Print 3D PCA plot
############################################################
############################################################
plotPCA3D <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE){
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  PC3 = pca$x[, 3],
                  group = group,
                  intgroup.df,
                  name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  message("Generating plotly plot")
  p <- plotly::plot_ly(data = d,
                       x = ~PC1,
                       y = ~PC2,
                       z = ~PC3,
                       color = group,
                       mode = "markers",
                       type = "scatter3d")
  return(p)
}


plotPCA3D(rld, intgroup = "condition", ntop = 500, returnData = FALSE)
############################################################



############################################################
########## Find the genes enriched in each component########
############################################################

biocLite("FactoMineR")
library("FactoMineR")
biocLite("factoextra")
library("factoextra")
# res.pca <- PCA(t(assay(rld)[select, ]), scale.unit = TRUE,graph = TRUE)
# print(res.pca)

### remove scientific notation in printing with this code!
options(scipen=999)

eigenvalues <- res.pca$eig
head(eigenvalues[, 1:2])
eigen <- eigenvalues[1:10,]
# Make a scree plot using base graphics : A scree plot is a graph of the eigenvalues/variances associated with components.
barplot(eigen[, 2], names.arg=1:nrow(eigen), 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
lines(x = 1:nrow(eigen), eigen[, 2], 
      type="b", pch=19, col = "red")
# Make the scree plot using the package factoextra
fviz_screeplot(res.pca, ncp=10)
##  identify the most correlated variables with a given principal component
res.desc <- dimdesc(res.pca, axes = c(1,2),proba = 0.05)
# Description of dimension 1
head(res.desc)
res.desc$Dim.1
dim1 <-data.frame(res.desc$Dim.1)
columns(org.Hs.eg.db)
dim1$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
dim1$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
dim1$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(dim1), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
summary(dim1)

res.desc$Dim.2
dim2 <-data.frame(res.desc$Dim.2)
columns(org.Hs.eg.db)
dim2$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
dim2$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
dim2$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(dim2), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
summary(dim2)

############################################################
##Plot of normalized counts for a single gene on log scale##
############################################################
# plotcount: "normalized" whether the counts should be normalized by size factor (default is TRUE)
# plotcount: "transform" whether to present log2 counts (TRUE) or to present the counts on the log scale (FALSE, default)
# re-arrange x-ase according to the following order: "Mock","siNeg","siBF1","ASO-neg","ASO-1","ASO-4"
# theme_bw() removes background color in the graph, guides(fill=FALSE) removes legends
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)

SREBF1 <- plotCounts(dds, gene="ENSG00000072310", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
SREBF1$condition <- factor(SREBF1$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SREBF1, aes(x=condition, y=log2(count), fill=condition)) + 
      geom_boxplot()+ 
      ylim(8, 10.5)+ 
      theme_bw()+ 
      guides(fill=FALSE)+
      theme(text = black.bold.18.text, 
            axis.text = black.bold.18.text,
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.line.x = element_line(color="black", size=1),
            axis.line.y = element_line(color="black", size=1),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(.25, "cm"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black",size=1),
            panel.background = element_blank())+
      labs(title = "SREBF1",x=" ", y= "log2(read counts)")

## *************************************************************
SCD <- plotCounts(dds, gene="ENSG00000099194", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
SCD$condition <- factor(SCD$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SCD, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(13, 14.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "SCD",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
FASN <- plotCounts(dds, gene="ENSG00000169710", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
FASN$condition <- factor(FASN$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(FASN, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(12, 13.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "FASN",x=" ", y= "log2(read counts)")
## *************************************************************
## *************************************************************
ACACA <- plotCounts(dds, gene="ENSG00000278540", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
ACACA$condition <- factor(ACACA$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(ACACA, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(9.7, 10.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "ACACA",x=" ", y= "log2(read counts)")
## *************************************************************
## *************************************************************
ACSL1 <- plotCounts(dds, gene="ENSG00000169710", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
ACSL1$condition <- factor(ACSL1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(ACSL1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(11.5, 13.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "ACSL1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
BRAF <- plotCounts(dds, gene="ENSG00000157764", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
BRAF$condition <- factor(SCD$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(BRAF, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(6, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "BRAF",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
SPIRE1 <- plotCounts(dds, gene="ENSG00000134278", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
SPIRE1$condition <- factor(SPIRE1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(SPIRE1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(4, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "SPIRE1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
INSIG1 <- plotCounts(dds, gene="ENSG00000186480", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
INSIG1$condition <- factor(INSIG1$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(INSIG1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "INSIG1",x=" ", y= "log2(read counts)")
## *************************************************************
## *************************************************************
YY1 <- plotCounts(dds, gene="ENSG00000100811", 
                     intgroup="condition", 
                     normalized = TRUE, 
                     transform = TRUE,
                     returnData=TRUE)
YY1$condition <- factor(YY1$condition, 
                           levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(YY1, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(8, 10.5)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "YY1",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
FABP7 <- plotCounts(dds, gene="ENSG00000164434", 
                  intgroup="condition", 
                  normalized = TRUE, 
                  transform = TRUE,
                  returnData=TRUE)
FABP7$condition <- factor(FABP7$condition, 
                        levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(FABP7, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(6, 9)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "FABP7",x=" ", y= "log2(read counts)")
## *************************************************************

## *************************************************************
MITF <- plotCounts(dds, gene="ENSG00000187098", 
                   intgroup="condition", 
                   normalized = TRUE, 
                   transform = TRUE,
                   returnData=TRUE)
MITF$condition <- factor(MITF$condition, 
                         levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))

ggplot(MITF, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  ylim(10, 12)+ 
  theme_bw()+ 
  guides(fill=FALSE)+
  theme(text = black.bold.18.text, 
        axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(color="black", size=1),
        axis.line.y = element_line(color="black", size=1),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(.25, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",size=1),
        panel.background = element_blank())+
  labs(title = "MITF",x=" ", y= "log2(read counts)")
## *************************************************************

SREBF2 <- plotCounts(dds, gene="ENSG00000198911", intgroup="condition", normalized = TRUE, transform = TRUE,returnData=TRUE)
SREBF2$condition <- factor(SREBF2$condition, levels = c("Mock","siNegative","siSREBF1","ASO-neg","ASO-1","ASO-4"))
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(SREBF2, aes(x=condition, y=count, fill=condition)) + 
  geom_boxplot()+ 
  theme_bw()+ guides(fill=FALSE)+
  theme(text = black.bold.18.text, axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(colour="black"),
        axis.line.y = element_line(colour="black"))+
  labs(title = "SREBF2",x=" ", y= "log2(readcount)")

HMGCR <- plotCounts(dds, gene="ENSG00000113161", intgroup="condition", normalized = TRUE, transform = TRUE,returnData=TRUE)
HMGCR$condition <- factor(HMGCR$condition, levels = c("Mock","siNegative","siSREBF1","ASO-neg","ASO-1","ASO-4"))
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(HMGCR, aes(x=condition, y=count, fill=condition)) + 
  geom_boxplot()+ 
  theme_bw()+ guides(fill=FALSE)+
  theme(text = black.bold.18.text, axis.text = black.bold.18.text,
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(colour="black"),
        axis.line.y = element_line(colour="black"))+
  labs(title = "HMGCR",x=" ", y= "log2(readcount)")


FABP7 <- plotCounts(dds, gene="ENSG00000164434", intgroup="condition", normalized = TRUE, transform = TRUE,returnData=TRUE)
FABP7$condition <- factor(FABP7$condition, levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(FABP7, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
   theme_bw()+ guides(fill=FALSE)+
  theme(text = black.bold.18.text, axis.text = black.bold.18.text,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(colour="black"),
        axis.line.y = element_line(colour="black"))+
  labs(title = "FABP7",x=" ", y= "log2(readcount)")

NRAS <- plotCounts(dds, gene="ENSG00000213281", intgroup="condition", normalized = TRUE, transform = TRUE,returnData=TRUE)
NRAS$condition <- factor(NRAS$condition, levels = c("Mock","siNegative","siSREBF1","ASO-Neg","ASO-1","ASO-4"))
black.bold.18.text <- element_text(face = "bold", color = "black", size = 18)
ggplot(NRAS, aes(x=condition, y=log2(count), fill=condition)) + 
  geom_boxplot()+ 
  theme_bw()+ guides(fill=FALSE)+
  theme(text = black.bold.18.text, axis.text = black.bold.18.text,
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line.x = element_line(colour="black"),
        axis.line.y = element_line(colour="black"))+
  labs(title = "NRAS",x=" ", y= "log2(readcount)")


########################################################
###### Enrichment Analysis using clusterPofiler ########
########################################################
# install exprAnalysis package from github
install.packages("devtools")
library(devtools)
#  the latest stable release that passed TRAVIS CI check
devtools::install_github("ShirinG/exprAnalysis", build_vignettes=TRUE, ref = "stable.version0.1.0")
library(exprAnalysis)
# data <- testseq[rowSums(DESeq2::counts(testseq)) > 1, ]
# data_DESeq <- DESeq2::DESeq(data)
expmatrix_DESeq <- DESeq2::rlog(data_DESeq, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)
library("pcaExplorer")
pcaExplorer(data_DESeq, expmatrix_DESeq)
groups <- as.factor(c(rep("Ctrl",4), rep("TolLPS",4), rep("TolS100A8",4), rep("ActLPS",4)))
pca_plot(expmatrix, groups)
pca_plot_enrich(expmatrix, groups)
DESeq2::resultsNames(data_DESeq)
res <- DESeq2::results(data_DESeq, contrast=list("treatmentActLPS", "treatmentCtrl"), cooksCutoff = 0.99, independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH")
summary(res)
mcols(res)$description
# order results table by the smallest adjusted p value:
res <- res[order(res$padj),]
results = as.data.frame(dplyr::mutate(as.data.frame(res), sig=ifelse(res$padj<0.05, "FDR<0.05", "Not Sig")), row.names=rownames(res))
head(results)
results_anno <- geneAnnotations(input=results, keys=row.names(results), column=c("ENTREZID", "ENSEMBL"), keytype="SYMBOL", organism = "human")
head(results_anno)



## Loading required package: DOSE
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
biocLite("fgsea")
library(clusterProfiler)
## Loading required package: AnnotationDbi
library(org.Hs.eg.db)
OrgDb <- org.Hs.eg.db # can also be other organisms

geneList <- as.vector(results_anno$log2FoldChange)
names(geneList) <- results_anno$ENTREZID
gene <- na.omit(results_anno$ENTREZID)


# Group GO
ggo <- clusterProfiler::groupGO(gene     = gene,
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
head(summary(ggo)[,-5])
barplot(ggo, drop=TRUE, showCategory=12)
# GO over-representation test
ego <- clusterProfiler::enrichGO(gene          = gene,
                                 OrgDb         = OrgDb,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05, 
                                 readable      = TRUE)
head(summary(ego)[,-8])
barplot(ego, showCategory=25)
clusterProfiler::dotplot(ego, showCategory=25)
#clusterProfiler::plotGOgraph(ego)
## KEGG over-representation test
kk <- clusterProfiler::enrichKEGG(gene         = gene,
                                  organism     = 'hsa',
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff  = 0.05)
head(summary(kk)[,-8])
cnetplot(kk, categorySize="geneNum", foldChange=geneList)




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
ressiRNA<-results(ddsDEsiRNA)
ressiRNA<-ressiRNA[order(ressiRNA$log2FoldChange),]
## An MA-plot21 provides a useful overview for an experiment with a two-group comparison
## The log2 fold change for a particular comparison is plotted on the y-axis and the average of the counts normalized by size factor is shown on the x-axis 
## (“M” for minus, because a log ratio is equal to log minus log, and “A” for average). Each gene is represented with a dot. 
## Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
plotMA(ddsDEsiRNA, main="DESeq2", ylim=c(-2,2))
## The red points indicate genes for which the log2 fold change was significantly higher than 0.5 or less than -0.5 
# (treatment resulting in more than doubling or less than halving of the normalized counts) with adjusted p value less than 0.1.
resLFC1 <- results(ddsDEsiRNA, lfcThreshold=0.5)
table(resLFC1$padj < 0.1)
plotMA(resLFC1,main="DESeq2", ylim=c(-2,2))

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
ressiRNA <- data.frame(ressiRNA)
columns(org.Hs.eg.db)
ressiRNA$symbol = mapIds(org.Hs.eg.db,
                         keys=row.names(ressiRNA), 
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")
ressiRNA$entrez = mapIds(org.Hs.eg.db,
                         keys=row.names(ressiRNA), 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
ressiRNA$name =   mapIds(org.Hs.eg.db,
                         keys=row.names(ressiRNA), 
                         column="GENENAME",
                         keytype="ENSEMBL",
                         multiVals="first")
summary(ressiRNA)
head(ressiRNA, 10)

###########################################################
## Heatmap for top variable genes across siNeg and siBF1 ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# The default “blinds” the normalization to the design. 
# topVarGenes looks at the row variance of the transformed values regardless of which samples come from which group. 
# Differential expression testing asks whether the difference across group is large relative to the within-group variance. 
# So these are different ways of ranking genes.
# calculate the variance for each gene, # select the top 50 genes by variance
rldsiRNA <- rlog(ddsDEsiRNA,blind=TRUE)
vdsiRNA <- varianceStabilizingTransformation(ddsDEsiRNA, blind=TRUE)
topVarGenes <- head(order(rowVars(assay(vdsiRNA)), decreasing=TRUE), 500)
mat <- assay(rldsiRNA)[topVarGenes, ]
mat <- mat - rowMeans(mat)
heatmap.2(mat,labRow = NA,labCol=NA, scale="row",lhei = c(2, 8),
          trace="none", dendrogram="none", margins=c(2,10),
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          ColSideColors = c(Control="gray", DPN="darkgreen", OHT="orange")[colData(rldsiRNA)$condition ])


############################
## KEGG pathways analysis ##
############################

data(kegg.gs)
data(go.gs)
data(carta.gs)
lapply(kegg.gs[1:3],head)
lapply(go.gs[1:3],head)
## kegg.sets.hs is a named list of 229 elements
## Each element is a character vector of member gene Entrez IDs for a single KEGG pathway
data(kegg.sets.hs)
## sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
data(sigmet.idx.hs)
## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of sinaling and metabolic pathways only.
kegg.sets.hs.sigmet <-  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs.sigmet, 3)

## Generate a named vector of fold changes, where the names of the values are the Entrez gene IDs, for the gage() function
foldchanges = -ressiRNA$log2FoldChange
names(foldchanges) = ressiRNA$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
keggres.sigmet = gage(foldchanges, gsets=kegg.sets.hs.sigmet, same.dir=TRUE)
gores = gage(foldchanges, gsets=go.gs, same.dir=TRUE)
cartares = gage(foldchanges, gsets=carta.gs, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres.sigmet, head,10)
lapply(gores, head,10)
lapply(cartares, head,10)

# Get the pathways
library(dplyr)
keggrespathways <- data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
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

write.table(keggres.sigmet$greater, file = "keggres.sigmet.greater.txt",sep = "\t")
write.table(keggres.sigmet$less, file = "keggres.sigmet.less.txt",sep = "\t")

## map metabolites on a large scale using pathway.id= '01100'
## Cancer  - Homo sapiens (human) pathway.id= '05200'
## Melanoma - Homo sapiens (human) pathway.id= '05218'
## hsa04010 MAPK signaling pathway   
## hsa04910 Insulin signaling pathway 
## hsa00061 fatty acid biosynthesis pathway 
pv.out <- pathview(gene.data = foldchanges, pathway.id= 'hsa04210', species = "hsa", kegg.native = T, same.layer = F)


######################################################
## derive the non-redundant signficant gene set lists
kegg.p <- gage(foldchanges, gsets=kegg.sets.hs.sigmet,same.dir=TRUE)
kegg.esg.up <- esset.grp(kegg.p$greater,
                         foldchanges, gsets=kegg.sets.hs.sigmet,
                         test4up = TRUE, output = TRUE, outname = "kegg.up", make.plot = FALSE)

names(kegg.esg.up)
head(kegg.p,4)
head(kegg.esg.up$essentialSets, 4)
head(kegg.esg.up$setGroups, 4)
head(kegg.esg.up$coreGeneSets, 40)
######################################################
## kegg test for 2-directional changes
kegg.2d.p <- gage(foldchanges, gsets = kegg.sets.hs.sigmet,same.dir = FALSE)
head(kegg.2d.p$greater)
head(kegg.2d.p$stats)
rownames(kegg.2d.p$greater)[1:3]

#################################
## Off-target analysis of ASO  ##
#################################
## generate dataset for siRNA treatment
testASO <- data.frame(testseq[,c(9,10,11,12,17,18,19,20)])

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
                             condition = c(rep("siBF1",4),rep("ASO-BF1",4)))

## object construction
## Construct DESeqDataSet with the count matrix, countData, and the sample information, colData
ddsASO <- DESeqDataSetFromMatrix(countData = guideDataASO,colData = guideDesignASO,design = ~ condition)
ddsASO

###########################################################################
## standard analysis to make MA-plot from base means and log fold changes##
###########################################################################
ddsDEASO <- DESeq(ddsASO)
plotMA(ddsDEASO, main="DESeq2", ylim=c(-2,2))
ddsresASO <- results(ddsDEASO)
plotMA(ddsresASO,main="DESeq2", ylim=c(-2,2))

# results extracts a result table from a DESeq analysis giving base means across samples, 
# log2 fold changes, standard errors, test statistics, p-values and adjusted p-values.
# lfcThreshold:a non-negative value which specifies a log2 fold change threshold. 
# The default value is 0, corresponding to a test that the log2 fold changes are equal to zero. 
# The user can specify the alternative hypothesis using the altHypothesis argument, which defaults to testing for log2 fold changes greater in absolute value than a given threshold. 

resLFC1ASO <- results(ddsDEASO, lfcThreshold=1)
table(resLFC1ASO$padj < 0.1)
plotMA(resLFC1ASO,main="DESeq2", ylim=c(-2,2))

# a scatter plot of log2 fold changes (on the y-axis) versus the mean of normalized counts (on the x-axis)
plotMA(resLFC1ASO, main="siSREBF1 vs ASO-1",ylim=c(-3,3), xlab="Mean of Normalized Counts", ylab="Log2 Fold Change")
topGene <- rownames(resLFC1ASO)[which.max(resLFC1ASO$log2FoldChange)]
with(resLFC1ASO[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})



resASO = ddsresASO[order(ddsresASO$pvalue),]
resASO <- data.frame(resASO)
summary(resASO)

###########################################################
## Regularized log transformation for clustering and PCA ##
###########################################################
# The regularized log transform can be obtained using the rlog() function. 
# Note that an important argument for this function is blind (TRUE by default). 
# The default “blinds” the normalization to the design. 
# This is very important so as to not bias the analyses (e.g. class discovery)
rldASO=rlog(ddsASO,blind=TRUE)
topVarGenes <- head( order( rowVars( assay(rldASO) ), decreasing=TRUE ), 50 )
heatmap.2( assay(rldASO)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", margins=c(2,10),
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rldASO)$condition ] )

#######################################################
## Add Entrez IDs, gene symbols, and full gene names ##
#######################################################
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






#############################
#############################
## Compare ASO-Neg and ASO-4 ##
#############################
#############################
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
source("http://bioconductor.org/biocLite.R")
biocLite("reactome.db")
library( "reactome.db" )
## This database works with Entrez IDs, 
## so we will need the entrezid column that we added earlier to the res object.
## First, we subset the results table, res, to only those genes for which the Reactome database has data 
## (i.e, whose Entrez ID we find in the respective key column of reactome.db and 
## for which the DESeq2 test gave an adjusted p value that was not NA.
res2 <- resASO[ resASO$entrez %in% keys( reactome.db, "ENTREZID" ) & !is.na( resASO$padj ) , ]
head(res2)
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

##################################################
## KEGG pathways analysis for ASO-4 and ASO-neg ##
##################################################

data(kegg.gs)
data(go.gs)
data(carta.gs)
lapply(kegg.gs[1:3],head)
lapply(go.gs[1:3],head)
## kegg.sets.hs is a named list of 229 elements
## Each element is a character vector of member gene Entrez IDs for a single KEGG pathway
data(kegg.sets.hs)
## sigmet.idx.hs is an index of numbers of sinaling and metabolic pathways in kegg.set.gs.
data(sigmet.idx.hs)
## kegg.sets.hs[sigmet.idx.hs] gives you the “cleaner” gene sets of signaling and metabolic pathways only.
kegg.sets.hs.sigmet <-  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs.sigmet, 3)


## Generate a named vector of fold changes, where the names of the values are the Entrez gene IDs, 
## for the gage() function
foldchanges = -resASO$log2FoldChange
names(foldchanges) = resASO$entrez
head(foldchanges)

# Get the results
# Run enrichment analysis on all log fc
keggres = gage(foldchanges, 
               gsets=kegg.sets.hs, 
               same.dir=F)
keggres.sigmet = gage(foldchanges, 
                      gsets=kegg.sets.hs.sigmet, 
                      same.dir=F)
gores = gage(foldchanges, 
             gsets=go.gs, 
             same.dir=TRUE)
cartares = gage(foldchanges, 
                gsets=carta.gs, 
                same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head,10)
lapply(keggres.sigmet, head,10)
lapply(keggres.gs, head,10)

lapply(gores, head,10)
lapply(cartares, head,10)

# Get the pathways
library(dplyr)
keggrespathways <- data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
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

write.table(keggres.sigmet$greater, 
            file = "keggres.sigmet.greater.txt",
            sep = "\t")
write.table(keggres.sigmet$less, 
            file = "keggres.sigmet.less.txt",
            sep = "\t")

## map metabolites on a large scale using pathway.id= '01100'
## Cancer  - Homo sapiens (human) pathway.id= '05200'
## Melanoma - Homo sapiens (human) pathway.id= '05218'
## hsa04010 MAPK signaling pathway   
## hsa04910 Insulin signaling pathway 
## hsa00061 fatty acid biosynthesis pathway 
pv.out <- pathview(gene.data = foldchanges, 
                   pathway.id= 'hsa04210', 
                   species = "hsa", 
                   kegg.native = T, 
                   same.layer = F)


############################################################################################################
############################################################################################################
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(exprAnalysis)
library(clusterProfiler)
library(DOSE)

# order results table by the smallest adjusted p value:
resASO <- resASO[order(resASO$padj),]
results = as.data.frame(dplyr::mutate(as.data.frame(resASO), 
                                      sig=ifelse(resASO$padj<0.05, 
                                                 "FDR<0.05", "Not Sig")), 
                        row.names=rownames(resASO))
head(results)


DEgenes_DESeq <- results[which(abs(results$log2FoldChange) > log2(1.5) & results$padj < 0.05),]

p = ggplot2::ggplot(results, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig)) +
  ggplot2::scale_color_manual(values = c("red", "black")) +
  ggplot2::ggtitle("Volcano Plot of DESeq2 analysis")

p + ggrepel::geom_text_repel(data=results[1:10, ], ggplot2::aes(label=rownames(results[1:10, ])))



OrgDb <- org.Hs.eg.db # can also be other organisms

geneList <- as.vector(resASO$log2FoldChange)
names(geneList) <- rownames(resASO)
gene <- na.omit(resASO$entrez)


# Group GO
ggo <- clusterProfiler::groupGO(gene     = gene,
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
head(summary(ggo)[,-5])
barplot(ggo, drop=TRUE, showCategory=12)

# GO over-representation test
ego <- clusterProfiler::enrichGO(gene          = gene,
                                 OrgDb         = OrgDb,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.05, 
                                 readable      = TRUE)
head(summary(ego))
barplot(ego, showCategory=25)
clusterProfiler::dotplot(ego, showCategory=25)
#enrichMap(ego)
cnetplot(ego, categorySize="pvalue", foldChange=geneList)
clusterProfiler::plotGOgraph(ego)

# enrichGO test the whole GO corpus and enriched result may contains very general terms. With dropGO function, user can remove specific GO terms or GO level from results obtained from both enrichGO and compareCluster.


## KEGG over-representation test
kk <- clusterProfiler::enrichKEGG(gene         = gene,
                                  organism     = 'hsa',
                                  pvalueCutoff = 0.05)
head(summary(kk))
barplot(kk, showCategory=8)
clusterProfiler::dotplot(kk)
cnetplot(kk, categorySize="geneNum", foldChange=geneList)
