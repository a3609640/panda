source("https://bioconductor.org/biocLite.R")

# install the dependent files for DiffBind
biocLite("GenomicFeatures")
biocLite("GenomicRanges")
biocLite("SummarizedExperiment")
biocLite("VariantAnnotation")
biocLite("RMySQL")
biocLite("systemPipeR")
biocLite("codetools")
biocLite("DiffBind")

library(GenomicFeatures)
library(GenomicRanges)
library(SummarizedExperiment)
library(VariantAnnotation)
library(RMySQL)
library(systemPipeR)
library(codetools)
library(lattice)
library(Matrix)
library(spatial)
library(DiffBind)
# setwd(system.file("extra", package="DiffBind"))
setwd("~/Documents/Bioinformatics analysis/ChIP analysis/SREBP1 four cell lines")

# Reading in the peaksets
samples <- read.csv("Diffbinddatarep.csv")
names(samples)
tamoxifen <- dba(sampleSheet="Diffbinddatarep.csv")
tamoxifen
# a correlation heatmap can be generated which gives an initial
# clustering of the samples using the cross-correlations of each row of the binding matrix
plot(tamoxifen) 

# Counting reads
tamoxifen <- dba.count(tamoxifen)
tamoxifen
# create a correlation plot showing the overlap of the peaks for all the samples.
plot(tamoxifen) 
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)

# Counting reads: Once a consensus peakset has been derived, DiffBind can use the
# supplied sequence read files to count how many reads overlap each interval for each
# unique sample. The peaks in the consensus peakset may be re-centered and trimmed
# based on calculating their summits (point of greatest read overlap) in order to provide
# more standardized peak intervals. The final result of counting is a binding affinity matrix
# containing a (normalized) read count for each sample at every potential binding site.
# With this matrix, the samples can be re-clustered using affinity, rather than occupancy,
# data. The binding affinity matrix is used for QC plotting as well as for subsequent
# differential analysis
# This keeps the peaks at a consistent width (in this case,
# with summits=250, the peaks will be 500bp, extending 250bp up- and down- stream of the
# summit):
tamoxifen <- dba.count(tamoxifen, summits=250)
tamoxifen
plot(tamoxifen) 
dba.plotHeatmap(tamoxifen)
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_ID)

# Establishing a contrast
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_TISSUE)
tamoxifen <- dba.analyze(tamoxifen)
plot(tamoxifen, contrast=1)

dba.plotPCA(tamoxifen, contrast=1,label=DBA_TISSUE)

tamoxifen <- dba.count(tamoxifen)
tamoxifen <- dba.contrast(tamoxifen)
tamoxifen <- dba.analyze(tamoxifen)
tamoxifen.DB <- dba.report(tamoxifen)
