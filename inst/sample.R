

tal1_hg_cd34_no_mock_peaks.xls

root <- "/Users/agriffith/masters/r-workspace/ChipPCompData/inst/extdata/"

bedFiles <- sapply(c("cd34","cem","erythroid","jurkat"),function(cellType){paste0(root,"tal1_hg_", cellType, "_aligned.bed")})

peakFiles <- sapply(c("cd34", "cem_1","erythroid","jurkat"),function(cellType){paste0(root,"tal1_hg_", cellType, "_no_mock_peaks.xls")})

categories <- c("cd34","cem_1","erythroid","jurkat")
makeAFS(peakFiles,)


biocLite("Rsamtools")

library("Rsamtools")
library("ShortRead")

library (ShortRead)

bamFile="/Users/agriffith/masters/r-workspace/ChipPCompData/inst/extdata/tal1_hg_cd34_aligned.bam"

which <- IRangesList(seq1=IRanges(1000,2000))
what <- c("rname", "strand", "pos", "qwidth","seq")
index <- indexBam(file)


param <- ScanBamParam(which=which, what=what)

bam <- scanBam (bamfile,bf)

bf <- BamFile(bamfile)




bf <- BamFile(bamfile)

seqnames(seqinfo(bf))

bf<-scanBam(bamFile,ScanBamParam(what=c('rname','pos', 'qwidth','strand','qname')))

data.frame(char=bf$rname,start=bf$pos,end=bf$pos+bf$qwidth)

loadChar<-function(char,bamFile){
    param <- ScanBamParam(what=c('pos', 'qwidth'),                      
                          which=GRanges(char, IRanges(1, 1e7)),
                          flag=scanBamFlag(isUnmappedQuery=FALSE))
    x <- scanBam(bamFile,  param=param)[[1]]
    print(char)
    str()
    data.frame(char=char,start=x$pos,
               end= x$pos+x$qwidth)
}


do.call(rbind,Filter(function(x){dim(x)==2},lapply(seqnames(seqinfo(bf)),loadChar,bamfile)))

data.frame(char=bf[[1]]$rname,start=bf[[1]]$pos,end=bf[[1]]$pos+bf[[1]]$qwidth)


lapply(seqnames(seqinfo(bf)),loadChar,bamfile)

scanBamWhat()
