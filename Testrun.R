# To convert the Ensembl IDs in the rownames of res to gene symbols and add them as a new column
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library("org.Hs.eg.db")

setwd("/Volumes/G-DRIVE\ mobile\ USB-C/Analysis/Testrun/STAR/results")

# convert the values in first column of data frame into row names
cols<-c("unstranded","1st strand","2nd strand")

inputFiles<-c(
  "test1-1ReadsPerGene.out.tab",
  "test1-2ReadsPerGene.out.tab",
  "test1-3ReadsPerGene.out.tab",
  "test1-4ReadsPerGene.out.tab",
  "test2-1ReadsPerGene.out.tab",
  "test2-2ReadsPerGene.out.tab",
  "test2-3ReadsPerGene.out.tab",
  "test2-4ReadsPerGene.out.tab",
  "test3-1ReadsPerGene.out.tab",
  "test3-2ReadsPerGene.out.tab",
  "test3-3ReadsPerGene.out.tab",
  "test3-4ReadsPerGene.out.tab",
  "test4-1ReadsPerGene.out.tab",
  "test4-2ReadsPerGene.out.tab",
  "test4-3ReadsPerGene.out.tab",
  "test4-4ReadsPerGene.out.tab",
  "test5-1ReadsPerGene.out.tab",
  "test5-2ReadsPerGene.out.tab",
  "test5-3ReadsPerGene.out.tab",
  "test5-4ReadsPerGene.out.tab",
  "test6-1ReadsPerGene.out.tab",
  "test6-2ReadsPerGene.out.tab",
  "test6-3ReadsPerGene.out.tab",
  "test6-4ReadsPerGene.out.tab"
)

# Sample format of each file:
# N_unmapped      105984  105984 105984
# N_multimapping  171224  171224 171224
# N_noFeature     177533 3433150 277677
# N_ambiguous     319796    9239 136511
# ENSG00000223972      0       0      0
# ENSG00000227232     10       0     10
# ENSG00000278267      0       0      0
# ...
# 
processFile <- function(fileName) {
  table <- read.table(fileName, stringsAsFactors=T)
  
  # STAR generates data files that row names (for metadata like the number of
  # unmapped segments, and for Ensembl gene designations) in the first column.
  # We'd like to use those as the row names.
  table.with.rownames <- data.frame(table[,-1], row.names=table[,1])
  
  # Furthermore, we'd like to rename the columns, using the test designator.
  # For example, for data from the file 'test6-4ReadsPerGene.out.tab', we'd like
  # the column headers to be 'test6.4.<N>', where N is the column number.
  colnames(table.with.rownames) <-
    paste0(
      # This inner 'paste' will generate 'test6.4.' from
      # 'test6-4ReadsPerGene.out.tab'.
      paste(
        substring(fileName, 1, 5),
        substring(fileName, 7, 7),
        '',  # gives us a trailing '.'
        sep='.'),
      1:3)
  return(table.with.rownames)
}


# Merge the gene count tables.  The desired format of the merged data is:
#
#                 test1.1.1 test1.1.2 test1.1.3 test1.1.4 test1.1.5 ...
# N_unmapped      105984       105984    105984
# N_multimapping  171224       171224    171224
# N_noFeature     177533      3433150    277677
# N_ambiguous     319796         9239    136511
# ENSG00000223972      0            0        0
# ENSG00000227232     10            0       10
# ENSG00000278267      0            0        0
# ...
#
processedData <- lapply(inputFiles, processFile)
test <- do.call(cbind.data.frame, processedData)

# check the basic properties of SRP070081
dim(test)
nrow(test)
length(test)
colnames(test)

test$symbol <- mapIds(org.Hs.eg.db,
                      keys=row.names(test),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")
