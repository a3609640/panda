library(plyr)
library(data.table)

setwd("~/Single cell RNA-seq")
singleRNAseq <- fread("GSE72056_melanoma_single_cell_revised_v2.txt")

# -----------------------------------------------------------------------------------------------
#select the genes of interest
SREBF1 <- singleRNAseq[singleRNAseq$Cell=="SREBF1",] 
FASN <- singleRNAseq[singleRNAseq$Cell=="FASN",] 
SCD <- singleRNAseq[singleRNAseq$Cell=="SCD",] 
ACACA <- singleRNAseq[singleRNAseq$Cell=="ACACA",]
CCND3 <- singleRNAseq[singleRNAseq$Cell=="CCND3",]
KDM5B <- singleRNAseq[singleRNAseq$Cell=="KDM5B",]
MITF <- singleRNAseq[singleRNAseq$Cell=="MITF",]
AXL <- singleRNAseq[singleRNAseq$Cell=="AXL",]
PPARGC1A <- singleRNAseq[singleRNAseq$Cell=="PPARGC1A",]
PPARGC1B <- singleRNAseq[singleRNAseq$Cell=="PPARGC1B",]
tumor <- singleRNAseq[singleRNAseq$Cell=="tumor",]
malignant <- singleRNAseq[singleRNAseq$Cell=="malignant(1=no,2=yes,0=unresolved)",]

# -----------------------------------------------------------------------------------------------
#combine them into a new table "geneset"
geneset <- rbind(tumor,malignant,SREBF1,FASN,SCD,ACACA,CCND3,KDM5B,MITF,AXL,PPARGC1A,PPARGC1B)
## transpose the table
# first remember the names
n <- geneset$Cell
# transpose all but the first column (name)
geneset1 <- as.data.frame(t(geneset[,-1]))
colnames(geneset1) <- n

# ?  geneset1$myfactor <- factor(row.names(geneset1))
str(geneset1) # Check the column types

#rename column name `malignant(1=no,2=yes,0=unresolved)' in geneset1 to 'malignant' in geneset2
geneset2 <- rename(geneset1,c('malignant(1=no,2=yes,0=unresolved)'='malignant'))

# select RNA-seq data from malignant or non-malignant cells
geneset3 <- subset(geneset2, malignant> 0) # total geneset
geneset4 <- subset(geneset2, malignant> 1) # malignant geneset
totalgeneset <- subset(geneset2, malignant> 0)
malignantgeneset <- subset(geneset2, malignant== 2)
nonmalignantgeneset <- subset(geneset2, malignant== 1)

# bad transpose geneset2 <- t(geneset)
# sort table based first on malignant (2 is malignant, 1 is normal), then second on tumor groups
# sortedgeneset1 <- geneset1[order(-geneset1$`malignant(1=no,2=yes,0=unresolved)`,geneset1$tumor), ]

# -----------------------------------------------------------------------------------------------
# load package ggplot2
library(ggplot2) 

# histogram to check the distribution of SREBF1 expression data, histogram is not a very good graphic tool here.
hist(geneset4$SREBF1)
hist(log2(geneset4$SREBF1))

# use geo_histogram to check the distribution of SREBF1 and FASN expression data in malignant cells, non-malignant cell only
ggplot(malignantgeneset, aes(x=log2(SREBF1))) + geom_histogram(fill="white", colour="black")
ggplot(nonmalignantgeneset, aes(x=log2(SREBF1))) + geom_histogram(fill="white", colour="black")
ggplot(malignantgeneset, aes(x=log2(FASN))) + geom_histogram(fill="white", colour="black")
ggplot(nonmalignantgeneset, aes(x=log2(FASN))) + geom_histogram(fill="white", colour="black")
ggplot(malignantgeneset, aes(x=log2(SCD))) + geom_histogram(fill="white", colour="black")
ggplot(nonmalignantgeneset, aes(x=log2(SCD))) + geom_histogram(fill="white", colour="black")
ggplot(malignantgeneset, aes(x=log2(ACACA))) + geom_histogram(fill="white", colour="black")
ggplot(nonmalignantgeneset, aes(x=log2(ACACA))) + geom_histogram(fill="white", colour="black")
ggplot(malignantgeneset, aes(x=log2(MITF))) + geom_histogram(fill="white", colour="black")
ggplot(nonmalignantgeneset, aes(x=log2(MITF))) + geom_histogram(fill="white", colour="black")

# use geo_histogram to check the distribution of SREBF1 and FASN expression data in both non-malignant and malignant cells
geneset5 <- geneset3
geneset5$malignant <- factor(geneset5$malignant)
levels(geneset5$malignant)
library(plyr)
geneset5$malignant <- revalue(geneset5$malignant, c("1"="Non-malignant", "2"="Malignant"))
# histogram over malignant or nonmalignant cells
ggplot(geneset5, aes(x=log2(SREBF1))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
ggplot(geneset5, aes(x=log2(FASN))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
ggplot(geneset5, aes(x=log2(SCD))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
ggplot(geneset5, aes(x=log2(ACACA))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
ggplot(geneset5, aes(x=log2(MITF))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
ggplot(geneset5, aes(x=log2(AXL))) + geom_histogram(fill="white", colour="black") + facet_grid(malignant~ .)
# overlay of the histogram on the same graph
ggplot(geneset5, aes(x=log2(SREBF1), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)
ggplot(geneset5, aes(x=log2(FASN), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)
ggplot(geneset5, aes(x=log2(SCD), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)
ggplot(geneset5, aes(x=log2(ACACA), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)
ggplot(geneset5, aes(x=log2(MITF), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)
ggplot(geneset5, aes(x=log2(AXL), fill=malignant)) + geom_histogram(position="identity", alpha=0.4)

ggplot(geneset5, aes(x=log2(SREBF1), fill=malignant)) + geom_density(alpha=.3)
ggplot(geneset5, aes(x=log2(FASN), fill=malignant)) + geom_density(alpha=.3)
ggplot(geneset5, aes(x=log2(SCD), fill=malignant)) + geom_density(alpha=.3)
ggplot(geneset5, aes(x=log2(ACACA), fill=malignant)) + geom_density(alpha=.3)
ggplot(geneset5, aes(x=log2(MITF), fill=malignant)) + geom_density(alpha=.3)
ggplot(geneset5, aes(x=log2(AXL), fill=malignant)) + geom_density(alpha=.3)

# use ggvis to graph data distribution
library(ggvis)
geneset5 %>% 
  ggvis(x=~log2(FASN), fill=~malignant) %>% 
  layer_densities(fill := "green") 

#compare averge SREBF1 expression level in malignant and non-malignant cells using boxplot graph
ggplot(geneset5, aes(x=factor(malignant), y=SREBF1)) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=FASN)) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=SCD)) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=ACACA)) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=MITF)) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=AXL)) + geom_boxplot()

ggplot(geneset5, aes(x=factor(malignant), y=log2(SREBF1))) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=log2(FASN))) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=log2(SCD))) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=log2(ACACA))) + geom_boxplot()
ggplot(geneset5, aes(x=factor(malignant), y=log2(MITF))) + geom_boxplot()

#compare averge SREBF1 expression level in various tumors using boxplot graph
ggplot(malignantgeneset, aes(x=factor(tumor), y=SREBF1)) + geom_boxplot()
ggplot(malignantgeneset, aes(x=factor(tumor), y=FASN)) + geom_boxplot()
ggplot(malignantgeneset, aes(x=factor(tumor), y=SCD)) + geom_boxplot()
ggplot(malignantgeneset, aes(x=factor(tumor), y=ACACA)) + geom_boxplot()

# scatter plot
ggplot(malignantgeneset, aes(x=log2(SREBF1), y=log2(SCD)), colour=cyl) + geom_point()
ggplot(malignantgeneset, aes(log2(SREBF1), log2(FASN))) + geom_point()

# put different colors on malignant vs nonmalignant cells
# ggplot(geneset5, aes(x=AXL, y=MITF, colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(AXL), y=log2(MITF), colour=factor(malignant))) + geom_point()
# ggplot(geneset5, aes(x=KDM5B, y=CCND3, colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(KDM5B), y=log2(CCND3), colour=factor(malignant))) + geom_point()
# ggplot(geneset5, aes(x=FASN, y=SREBF1, colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(SCD), y=log2(FASN), colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(SCD), y=log2(MITF), colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(CCND3), y=log2(MITF), colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(SCD), y=log2(SREBF1), colour=factor(malignant))) + geom_point()
ggplot(geneset5, aes(x=log2(SCD), y=log2(CCND3), colour=factor(malignant))) + geom_point()

# inter-tumor distritution of different gene expression
ggplot(malignantgeneset, aes(x=log2(AXL), y=log2(MITF), colour=factor(tumor))) + geom_point()
ggplot(malignantgeneset, aes(x=log2(KDM5B), y=log2(CCND3), colour=factor(tumor))) + geom_point()
ggplot(malignantgeneset, aes(x=log2(SCD), y=log2(FASN), colour=factor(tumor))) + geom_point()
ggplot(malignantgeneset, aes(x=log2(SCD), y=log2(SREBF1), colour=factor(tumor))) + geom_point()

# does SCD define any different gene expression? did not find any
ggplot(malignantgeneset, aes(x=log2(AXL), y=log2(MITF), colour=SCD)) + geom_point()
ggplot(malignantgeneset, aes(x=log2(KDM5B), y=log2(CCND3), colour=SCD)) + geom_point()
ggplot(malignantgeneset, aes(x=log2(SCD), y=log2(FASN), colour=SCD)) + geom_point()

ggplot(totalgeneset, aes(x=log2(AXL), y=log2(MITF), colour=SCD)) + geom_point()
ggplot(totalgeneset, aes(x=log2(KDM5B), y=log2(CCND3), colour=SCD)) + geom_point()
ggplot(totalgeneset, aes(x=log2(SCD), y=log2(FASN), colour=SCD)) + geom_point()

ggplot(geneset4, aes(x=log2(AXL), y=log2(MITF), colour=SREBF1)) + geom_point()
ggplot(geneset4, aes(x=log2(KDM5B), y=log2(CCND3), colour=SREBF1)) + geom_point()
ggplot(geneset4, aes(x=log2(SCD), y=log2(FASN), colour=SREBF1)) + geom_point()


# -----------------------------------------------------------------------------------------------
#t-test
t.test(x=nonmalignantgeneset$MITF, y=malignantgeneset$MITF)
t.test(x=nonmalignantgeneset$AXL, y=malignantgeneset$AXL)
t.test(x=nonmalignantgeneset$FASN, y=malignantgeneset$FASN)
t.test(x=nonmalignantgeneset$SCD, y=malignantgeneset$SCD)
t.test(x=nonmalignantgeneset$SREBF1, y=malignantgeneset$SREBF1)
t.test(x=nonmalignantgeneset$ACACA, y=malignantgeneset$ACACA)
t.test(x=nonmalignantgeneset$CCND3, y=malignantgeneset$CCND3)

# Creating a Three-Dimensional Scatter Plot 
library(rgl)
open3d()
plot3d(x=malignantgeneset$MITF, y=malignantgeneset$FASN, z=malignantgeneset$SCD, col=rainbow(1000))

# -----------------------------------------------------------------------------------------------
# heatmap
nba_matrix <- data.matrix(totalgeneset)
nba_matrix1 <- subset(nba_matrix, select = -tumor )
nba_matrix2 <- subset(nba_matrix1, select = -malignant )
heatmap(nba_matrix2, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))

library(gplots)
heatmap.2(nba_matrix2,dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none')

distancem <- dist(nba_matrix2)
hclust_completem <- hclust(distancem, method = "complete")
dendcompletem <- as.dendrogram(hclust_completem)
heatmap(nba_matrix2, Rowv=dendcompletem, Colv=NA, scale="column")
