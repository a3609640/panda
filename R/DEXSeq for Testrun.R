source("http://bioconductor.org/biocLite.R")
biocLite("DEXSeq")
library(DEXSeq)
suppressPackageStartupMessages( library( "DEXSeq" ) )
library(BiocParallel)



## check the two Python scripts, dexseq_prepare_annotation.py and dexseq_count.py, that come with the DEXSeq package.
pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)

## read files into system, notice wildcard here: only files starting with test1 or test5
#countFiles1 = list.files(path = "/PHShome/sw542/HTSeq/ASO1", pattern="*.text", full.names=TRUE)
#flattenedFile = list.files(path = "/PHShome/sw542/HTSeq", pattern="*.chr.gff", full.names=TRUE)

countFiles1 = list.files(  path = "/home/data/su/data/Analysis/Testrun/HTSeq/ASO1", pattern="*.text", full.names=TRUE)
flattenedFile = list.files(path = "/home/data/su/data/Analysis/Testrun/HTSeq", pattern="*.chr.gff", full.names=TRUE)

basename(countFiles1)
basename(flattenedFile)


sampleTable1 = data.frame(
  row.names = c("test1-1","test1-2","test1-3","test1-4","test5-1", "test5-2","test5-3","test5-4"),
  condition = c("Mock", "Mock", "Mock","Mock","ASO-1", "ASO-1","ASO-1", "ASO-1"),
  libType = c( "paired-end", "paired-end", "paired-end", "paired-end"))
sampleTable1

dxd1 = DEXSeqDataSetFromHTSeq(
  countFiles1,
  sampleData=sampleTable1,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

colData(dxd1)
head( counts(dxd1), 5 )
split( seq_len(ncol(dxd1)), colData(dxd)$exon )
head( featureCounts(dxd1), 5 )
head( rowRanges(dxd1), 3 )
sampleAnnotation(dxd1)


BPPARAM = MulticoreParam(workers=4)
dxd1 = estimateSizeFactors(dxd1)
dxd1 = estimateDispersions(dxd1, BPPARAM=BPPARAM)
dxd1 = testForDEU(dxd1, BPPARAM=BPPARAM)
dxd1 = estimateExonFoldChanges(dxd1, fitExpToVar="condition", BPPARAM=BPPARAM)

## plot Exon usage graph
dxrASO1 = DEXSeqResults( dxd1 )
table( dxr2$padj < 0.1 )
# SREBF1
plotDEXSeq( dxrASO1, "ENSG00000072310", legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
plotDEXSeq( dxrASO1, "ENSG00000072310",expression=FALSE, splicing=TRUE,legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
# FASN
plotDEXSeq( dxrASO1, "ENSG00000169710", legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
# SCD
plotDEXSeq( dxrASO1, "ENSG00000099194", legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )

# SPIRE1
plotDEXSeq( dxrASO1, "ENSG00000134278", legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
plotDEXSeq( dxrASO1, "ENSG00000134278", expression=FALSE, splicing=TRUE,legend=TRUE, cex.axis=1.2, cex=1.3,
            lwd=2 )
