library(data.table)
library(ggplot2)
setwd("~/GTEx")

gene_median <- fread("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.csv",showProgress = T)
gene <- fread("GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.csv",showProgress = T)
Annotations <- fread("GTEx_Data_V6_Annotations_SampleAttributesDS.txt",showProgress = T)
Annotations <-Annotations[,c(1,7)]
## showProgress = T is necessary, otherwise "Error: isLOGICAL(showProgress) is not TRUE"

######################################################################################################
plotGene <- function(goi) {
  go <- gene[gene$Description==goi,]
  go<-t(go)
  ## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
  go <- data.frame(go)
  setDT(go, keep.rownames = TRUE)[]
  colnames(go) <- c("SAMPID", "goi")
  #A one line option is: df$names<-rownames(df)
  goexpression <- merge(go,Annotations,by='SAMPID')
  ## somehow the numbers of SREBF1 columns are all changed into character 
  goexpression$SMTSD <- as.factor(goexpression$SMTSD)
  #  In particular, as.numeric applied to a factor is meaningless, and may happen by implicit coercion. 
  #  To transform a factor f to approximately its original numeric values, as.numeric(levels(f))[f] is recommended. 
  goexpression$goi <- as.numeric(levels(goexpression$goi))[goexpression$goi]
  goexpression <- as.data.frame(goexpression)
  ## draw boxplot for FASN expression across different tissues
  mean <- within(goexpression, SMTSD <-reorder(SMTSD, log2(goi), median))
  ggplot(mean, aes(x=SMTSD, y=log2(goi))) + geom_boxplot()+ theme_bw()+ 
    labs(x = "Tissue types (GTEx)",
         y = paste("log2(", goi, "RNA counts)")) +
    theme(axis.title=element_text(face="bold",size=12,color="black"),
          axis.text=element_text(size=12,angle = 90, hjust = 1, face="bold",color="black"),
          axis.line.x = element_line(color="black"),
          axis.line.y = element_line(color="black"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12, colour = "black"),
          legend.text = element_text(size=12,face="bold",colour = 'black'),
          legend.title = element_text(face = "bold", size = 12, colour = "black"),
          legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))
}

plotGene(goi="FASN")
plotGene(goi="SCD")
plotGene(goi="SREBF1")
plotGene(goi="HMGCR")
plotGene(goi="HMGCS1")
plotGene(goi="SREBF2")
plotGene(goi="MITF")

######################################################################################################
SREBF1 <- gene[gene$Description=="SREBF1",]
SREBF1<-t(SREBF1)
## SREBP1 is generated as a matrix somehow, and it needs to be converted into data frame
SREBF1 <- as.data.frame(SREBF1)
setDT(SREBF1, keep.rownames = TRUE)[]
colnames(SREBF1) <- c("SAMPID","SREBF1")
#A one line option is: df$names<-rownames(df)

SREBF1expression <- merge(SREBF1,Annotations,by='SAMPID')
class(SREBF1expression$SMTSD)
class(SREBF1expression$SREBF1)
## somehow the numbers of SREBF1 columns are all changed into character 

SREBF1expression$SMTSD <- as.factor(SREBF1expression$SMTSD)
levels(SREBF1expression$SMTSD )
#  In particular, as.numeric applied to a factor is meaningless, and may happen by implicit coercion. 
#  To transform a factor f to approximately its original numeric values, as.numeric(levels(f))[f] is recommended. 
SREBF1expression$SREBF1 <- as.numeric(levels(SREBF1expression$SREBF1))[SREBF1expression$SREBF1]
class(SREBF1expression)
SREBF1expression <- as.data.frame(SREBF1expression)
write.csv(SREBF1expression,"SREBF1expression")

## draw boxplot for SREBF1 expression across different tissues
mean <- within(SREBF1expression, SMTSD<-reorder(SMTSD, log2(SREBF1), median))
ggplot(mean, aes(x=SMTSD, y=log2(SREBF1))) + geom_boxplot()+ theme_bw()+ 
  labs(x = "Tissue types",
       y = paste("log2(SREBF1 RNA counts)")) +
        theme(axis.title=element_text(face="bold",size=12,color="black"),
        axis.text=element_text(size=12,angle = 90, hjust = 1, face="bold",color="black"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 12, colour = "black"),
        legend.text = element_text(size=12,face="bold",colour = 'black'),
        legend.title = element_text(face = "bold", size = 12, colour = "black"),
        legend.justification=c(1,1))+ guides(fill=guide_legend(title=NULL))

SREBF1group <-SREBF1expression[,c(2,3)]
class(SREBF1group$SREBF1)
SREBF1.list <- with(SREBF1group, split(SREBF1,SMTSD))

n.obs <- sapply(SREBF1.list, length)
seq.max <- seq_len(max(n.obs))
SREBF1.matrix <- (sapply(SREBF1.list, "[", i = seq.max))
write.csv(SREBF1.matrix,"SREBF1.matrix.csv")


