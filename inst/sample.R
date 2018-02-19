## Base library
library(ChipPComp)
## To generate fasta subsets for motif analysis
library(BSgenome.Hsapiens.UCSC.hg19)

## Replace this root to wherever the ChipPCompData libary is stored
root <- "/Users/agriffith/masters/r-workspace/ChipPCompData/inst/extdata/"
peakList <- sapply(c("cd34", "cem_1","erythroid","jurkat"),
                   function(cellType){paste0(root,"tal1_hg_", cellType, "_no_mock_peaks.xls")})
rawData <- sapply(c("cd34", "cem","erythroid","jurkat"),
                  function(cellType){paste0(root,"tal1_hg_", cellType, "_aligned.bed")})
categories <- c("cd34","cem_1","erythroid","jurkat")
ccca<-makeCCCA(rawData,peakList,categories,10)

## visualize the principal components
plotPCA(ccca)

## Based on the visualized principal components select subsections
nm<-normalizePRC(ccca$prc)
## Name leukemia, where the normalized 1st pc is less than 1 sd of the mean
ccca<-addRegion(ccca,"Leukemia",nm[,1]<(mean(nm[,1])-1*sd(nm[,1])))
## Name Erythroid, where the normalized 1st pc is greater than 1 sd of the mean
ccca<-addRegion(ccca,"Erythroid",nm[,1]>(mean(nm[,1])+1*sd(nm[,1])))
## Name CD34, where the normalized 2nd pc is greater than 1 sd of the mean
ccca<-addRegion(ccca,"CD34",nm[,2]>(mean(nm[,2])+1*sd(nm[,2])))

genome<-BSgenome.Hsapiens.UCSC.hg19
ccca<-addFasta(ccca,genome)

## The bed files for CD34
ccca$afs[ccca$reg[,"CD34"],c("chr","start","end")]
## Output the basepairs for CD34
ccca$fasta[ccca$reg[,"CD34"],]



