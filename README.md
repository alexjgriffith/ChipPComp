# ChipPComp

## Example
These examples use the data sets that will be in ChipPCompData, change the root directory to wherever you have this package saved.

Load the libraries
``` R
## Base library
library(ChipPComp)
## To generate fasta subsets for motif analysis
library(BSgenome.Hsapiens.UCSC.hg19)
```

Set the root to wherever you have saved the ChipPCompData package. Alternativly you could provide your own peaks and aligned sequences. Note that the peaks must be in the XLS format output by MACS 1.3 or 2.0+ and the aligned sequences must be in a Bed3 format (i.e. <chr><start><end>).

``` R
root <- "/Users/agriffith/masters/r-workspace/ChipPCompData/inst/extdata/"
```

Generate file names
``` R
peakList <- sapply(c("cd34", "cem_1","erythroid","jurkat"),
                   function(cellType){paste0(root,"tal1_hg_", cellType, "_no_mock_peaks.xls")})
rawData <- sapply(c("cd34", "cem","erythroid","jurkat"),
                  function(cellType){paste0(root,"tal1_hg_", cellType, "_aligned.bed")})
categories <- c("cd34","cem_1","erythroid","jurkat")
```

`makeCCCA` will generate the AFS, UDM and principal components for the input sequences and peaks. The final variable is a cut off for the peaks. (pvalue = 10^-value)
``` R
ccca<-makeCCCA(rawData,peakList,categories,10)
```

If you have a genome installed you can now add it to the `CCCA` object to resolve the base pairs under individual peaks.
``` R
genome<-BSgenome.Hsapiens.UCSC.hg19
ccca<-addFasta(ccca,genome)
```

Simply visualize the principal components (PCs). The second arguments is a tuple of the PCs that you want to visualize.
```
plotPCA(ccca,c(1,2))
```

Based on the visualized PCs select subsections that seperate categories. For example, PC1 seperates Leukemic from Erythroid cell types and PC2 isolates CD34.
```R
## normalize the principal components
nm<-normalizePRC(ccca$prc)
## Name leukemia, where the normalized 1st pc is less than 1 sd of the mean
ccca<-addRegion(ccca,"Leukemia",nm[,1]<(mean(nm[,1])-1*sd(nm[,1])))
## Name Erythroid, where the normalized 1st pc is greater than 1 sd of the mean
ccca<-addRegion(ccca,"Erythroid",nm[,1]>(mean(nm[,1])+1*sd(nm[,1])))
## Name CD34, where the normalized 2nd pc is greater than 1 sd of the mean
ccca<-addRegion(ccca,"CD34",nm[,2]>(mean(nm[,2])+1*sd(nm[,2])))
```

You can now ouput the peaks important to a specific category of cells
```R
## The bed files for CD34
ccca$afs[ccca$reg[,"CD34"],c("chr","start","end")]
```

You can also output the basepairs if you wish to do motif analysis
``` R
## Output the basepairs for CD34
ccca$fasta[ccca$reg[,"CD34"],]
```
