#' Compare multiple ChIP data sets across experimental conditions
#'
#' @docType package
#' @name ChipPComp
#' @author "Alexander Griffith <griffitaj@@gmail.com>"
NULL

#' @useDynLib ChipPComp, unityOutput, pileup
#' @import parallel
#' @import Biostrings
#' @import graphics
#' @importFrom stats var
#' @importFrom utils read.table str write.table
NULL


#' Read Peaks XLS format
#'
#' Read in peaks that MACs has output in its XLS format. Works with MACS
#' 2.1 and 1.3.1
#' @param file The MACS xls output to be imported
#' @param name A string as a moniker for this data set
#' @param pvalue log10 pvalue used as a cut off, default 0
#' @return a data.frame <chr><summit><name> | NULL if file has no entries
#' @examples
#' # need to create sample.xls
#' filename<-system.file("extdata","sample1.xls", package = "CCCA")
#' readPeaksXLS(filename)
#' @export
readPeaksXLS<-function(file,name=file,pvalue=0){
    bedData<-read.table(file,header=TRUE,skip="#")
    bedData<-bedData[bedData$X.log10.pvalue>pvalue,]
    if(dim(bedData)[1]>0){
        ret<-data.frame(bedData$chr,bedData$abs_summit,name)
        colnames(ret)<-c("chr","summit","name")
    }
    else
        ret<-NULL    
    ret
}


#' unity output
#'
#' Wrapper for the unityOutput C function.
#' 
#' @param peaks A 2xN matrix representing overlaping peak regions
#' @param intChr N length list of chomosomes represented as integers
#' @param intSummit The summit location 
#' @param intname The data set or name that the peak come froms in int form
._unityOutput<-function(peaks,intChr,intSummit,intname){
    nchr<-length(unique(intname))
    lpeaks<-dim(peaks)[1]
    retChr<-integer(lpeaks)
    retMatrix<-integer(lpeaks*nchr)
    retSummit<-integer(lpeaks)
    data<-.C(unityOutput,as.integer(intChr),
             as.integer(intSummit),
             as.integer(intname),
             as.integer(peaks[,1]),
             as.integer(peaks[,2]),
             as.integer(lpeaks),
             as.integer(nchr),chr=retChr,summit=retSummit,matrix=retMatrix)
    list(chro=data$chr,summit=data$summit,
         matrix=t(matrix(data$matrix,nrow=nchr)))
}

#' Get Pile Up
#'
#' A minimal wrapper for the pileup function from hg19Height.c
#' returns a integer vector of computed read pileups
#' @param file The raw data file of interest
#' @param bed The preloaded bed infromation including bed$start and bed$end
#' @param chroms a list of chromosomes whos string values have
#' been repalced with ranks
#' @param peakLength The length of the bed data provided
._getPileUp<-function(file,bed,chroms,peakLength){
    file<-normalizePath(file)
    if(!file.exists(file)){
        stop(paste0("Can't find file: ",file))
    }
    if(length(bed$start)!=length(chroms))
        stop("Chromosomes and start length differ")
    start<-as.integer(as.character(bed$start))
    end<-as.integer(as.character(bed$end))
    chroms<-as.character(chroms)
    peaknum<-as.integer(peakLength)
    score<-as.integer(rep(0,peakLength))
    results<-.C(pileup,file,chrom=chroms,start=start,end=end,
                peaknum=peaknum,score=score)
    results$score
}

#' Make AFS
#'
#' Take a list of peaks and create a unified set. The unified set is
#' created by grouping summits that fall within \code{width} of one
#' another. The mean of the grouped summits is used as a shared summit
#' between the grouped peaks. The start and end of the new peak is defined
#' as the new sumit +/ 1/2 \code{width}.
#' @param peakList A list of strings refering to file locations containing
#' Peak data
#' @param categories The monikers these peak files will be refered to ass
#' @param format The format to be imported, either bed or xls
#' @param pvalue The pvalue cut off to be applied, only valid with xls format
#' @param width The min width summits must be from one another in order to
#' cluster
#' @return The AFS generated from the input peak lists <chr><summit><catagory>
#' @examples
#' peakList<-sapply(c("sample1.xls","sample2.xls","sample3.xls"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' categories<-c("s1","s2","s3")
#' sampleafs<-makeAFS(peakList,categories,pvalue=20)
#' @export
makeAFS<-function(peakList,categories,format="xls",pvalue=NULL,width=700){
    if(format!="xls" & !is.null(pvalue)){
        stop("to have a pvalue cut off the file format must be in
MACS XLS output")
    }
    if(! format %in% c("xls") ){
        stop("format options are (xls)")
    }
    ## loadBedFun<-switch(format,
    ##                    xls=function(a,b=a){readPeaksXLS(a,b,pvalue)},
    ##                    bed=readPeaksBed
    ##                    )
    ## Fix to include pvlaue
    loadBedFun<-function(a,b=a) {
        readPeaksXLS(a,b,ifelse(is.null(pvalue),0,pvalue))
    }
    ## Load each of th files
    files<-lapply(._mapziplist(peakList,categories),
                  function(x) do.call(loadBedFun, as.list(x)))
    ## Unify the peaks based on width (default 700)
    ## Fix issue with readPeaksXLS
    testBed<-._unifyBedFile(
        ._sortDataFrame(
            do.call(rbind,Filter(function(x) ! is.null(x), files)),
            "chr","summit"),width)
    ## Make sure none of the peaks have values less than 0

    testBedSE<-._shiftFromZero(testBed$summit)
    ## Return a data frame of form <chr><start><end><h1>...<hn>
    width<-dim(testBed)[2]
    retData<-data.frame(chr=testBed[,1],
                        start=testBedSE[,1],
                        end=testBedSE[,2])
    ret<-cbind(retData,testBed[,3:width])
    ret
}

#' Write AFS
#' Write the AFS to a table
#' @param data AFS object to be saved
#' @param fname Save filename
#' @return Nothing
#' @export
#' @examples
#' filename<-system.file("extdata","sample.afs", package = "CCCA")
#' sampleAFS<-readAFS(filename)
#' writeAFS(sampleAFS,"sample.afs")
writeAFS<-function(data,fname){
    frame<-as.data.frame(do.call(cbind,data))
    chrs<-data$chr
    frame[,1]<-chrs
    frame
    write.table(frame,file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}

#' Order Bed
#'
#' sorts a data.frame based on the first and third columns
#' @param ret data frame
#' @return sorted data frame
#' @examples
#' orderBed(data.frame(c("a","a","c"),c(1,2,3),c(5,4,6)))
._orderBed<-function(ret){
    if(is.list(ret)){
        if(sum(!is.na(match(c("chr","start"),names(ret))))<2)
            stop("requires a chr and start member if list")
        or<-order(as.character(ret$chr),ret$start)
        name<-names(ret)
        ret<-lapply(ret,function(x) x[or])
        names(ret)<-name
        ret
    }
    else{
        ret[order(as.character(ret[,1]),ret[,3]),] # for matrix
    }    
}


#' Read AFS
#' Read an AFS file
#' @export
#' @param fname Save filename
#' @return An AFS object
#' @examples
#' filename<-system.file("extdata","sample.afs", package = "CCCA")
#' sampleAFS<-readAFS(filename)
readAFS<-function(fname){
    ret<-._orderBed(read.table(fname,header=TRUE))
    ret
}


#' sort data frame
#'
#' A small utility to functionalize the sorting of data frames
#' @examples
#' dd<-data.frame(initial1=LETTERS[runif(100,1,24)],
#'                initial2=LETTERS[runif(100,1,24)],
#'                age=floor(runif(100,21,35)))
#' sortDataFrame(dd,"initial2","age")
._sortDataFrame<-function(dd,...)    
    dd[do.call(order,lapply(list(...),function(x) dd[x])),]


#' Make UDM
#'
#' Generates a pile up matrix from a unified set of peaks and a list of
#' raw data sets.
#' @param peaks A preloaded bed data.frame which includes slots $chro
#' $start $end
#' @param rawdata a list of raw data files
#' @param n the number of nodes to use. If 0 then the parrallel package
#' is not used
#' @param verbose if not null then print each file that is found
#' @param clust pass in a cluster defined outside of makeUDM
#' @return A matrix of class UDM containing the pile up counts under each
#' peak in the AFS
#' @examples
#' rawdata<-sapply(c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' sampleAFS<-readAFS(system.file("extdata","sample.afs", package = "CCCA"))
#' sampleUDM<-makeUDM(sampleAFS,rawdata,n=0)
#' ## Parallel version
#' ## library(parallel)
#' ## cl<-makeForkCluster(3)
#' ## sampleUDM<-makeUDM(sampleAFS,rawdata,n=1,clust=cl)
#' ## stopCluster(cl)
#' @export
makeUDM<-function(peaks,rawdata,n=0,verbose=NULL,clust=NULL){
    for(file in rawdata){
        if(! file.exists(file)){
            stop("Can't find file ",file,".")
        }
        if(!is.null(verbose))
            print(paste("# Raw data file ",file," was found.",sep=""))
    }
    peakLength<-length(peaks$chr)
    ## Sort afs
    data<-._orderBed(peaks)    
    chroms<-as.character(data$chr)
    ## Parallell Implementation
    if(n>0){
        if(!requireNamespace("parallel"))
            stop("The parallel package must be loaded to run in parallel")
        ## if a cluster is passed in then make one now
        if(is.null(clust))
            cs<-makeForkCluster(n,renice=0)
        else
            cs<-clust
        ret<-matrix(unlist(parLapply(cs,rawdata,._getPileUp,data,chroms,
                                     peakLength)),nrow=peakLength)
        if(is.null(clust))
            stopCluster(cs)
    }
    ## Serial Implementation
    else{
        ret<-matrix(unlist(lapply(rawdata,._getPileUp,data,
                                  chroms,peakLength)),
                    nrow=peakLength)
    }
    ret
}

#' Write UDM
#'
#' Write the unified data matrix to a file.
#' @param data UDM object to be saved
#' @param fname Save filename
#' @return Nothing
#' @export
#' @examples
#' filename<-system.file("extdata","sample.udm", package = "CCCA")
#' sampleUDM<-readUDM(filename)
#' writeUDM(sampleUDM,"sample.udm")
writeUDM<-function(data,fname){
    write.table(data,file=fname,quote=FALSE,row.names=FALSE,sep="\t")
}

#' Read UDM
#'
#' Read the UDM from its tab format
#' @param fname Save filename
#' @export
#' @return A udm object
#' @examples
#' filename<-system.file("extdata","sample.udm", package = "CCCA")
#' sampleUDM<-readUDM(filename)
readUDM<-function(fname){
    ret<-read.table(fname,header=TRUE)
    ret
}


#' makeCCCA
#'
#' Generate a ChIP Component object from peaks, reads, and a list of
#' categories
#' @param dataSets list of files containing the raw reads in bed format
#' @param peakLists list of files containing the peaks reads in MACS'
#' xls format
#' @param categories the moniker to be attached to each of the data set files
#' @param pvalue the pvalue that the AFS peaks will be filtered by
#' @return a list of class ccca containing afs,udm, and prc
#' @examples 
#' dataSets<-sapply(c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' peakLists<-sapply(c("sample1.xls","sample2.xls","sample3.xls"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' categories<-c("s1","s2","s3")
#' ccca<-makeCCCA(dataSets,peakLists,categories,20)
#' ## Find dimensions that seperate data sets of interest
#' plot(ccca,1,2)
#' plot(ccca,1,3)
#' ccca<-local({
#'   nm<-CCCA::normalizePRC(ccca$prc)
#'   ccca<-addRegion(ccca,"s1",nm[,1]<(mean(nm[,1])-3*sd(nm[,1])))
#'   ccca<-addRegion(ccca,"s2",nm[,1]>(mean(nm[,1])+3*sd(nm[,1])))
#'   ccca<-addRegion(ccca,"s3",nm[,2]>(mean(nm[,2])+3*sd(nm[,2])))
#'   ccca<-addRegion(ccca,"s1.me",ccca$reg[,"s1"] & !(ccca$reg[,"s2"]
#'                | ccca$reg[,"s3"]))
#'   ccca<-addRegion(ccca,"s2.me",ccca$reg[,"s2"] & !(ccca$reg[,"s1"]
#'                | ccca$reg[,"s3"]))
#'   ccca<-addRegion(ccca,"s3.me",ccca$reg[,"s3"] & !(ccca$reg[,"s2"]
#'                | ccca$reg[,"s1"]))
#'   ccca
#' })
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome=BSgenome.Hsapiens.UCSC.hg19
#' addFasta(ccca,genome)
#' }
#' @export
makeCCCA<-function(dataSets,peakLists,categories,pvalue){
    afs<-makeAFS(peakLists,categories,pvalue=pvalue)
    udm<-makeUDM(afs,dataSets)
    prc<-makePRC(udm)
    ret<-list(afs=afs,udm=udm,prc=prc,fasta=NULL,reg=NULL,
              categories=categories)
    ret
}

#' Load CCCA
#'
#' Generate the ccca object from premade AFS and UDM files.
#' See \link[CCCA]{makeCCCA}
#' @param peaks the AFS tab file
#' @param heights The UDM tab file
#' @param categories an array of categories
#' @return A ccca object
#' @examples
#' afsfiles<-system.file("extdata", "sample.afs", package = "CCCA")
#' udmfiles<-system.file("extdata", "sample.udm", package = "CCCA")
#' ccca<-loadCCCA(afsfiles,udmfiles,c("s1","s2","s3"))
#' @export
loadCCCA<-function(peaks,heights,categories){
    afs<-readAFS(peaks)
    udm<-readUDM(heights)
    prc<-makePRC(udm)
    ret<-list(afs=afs,udm=udm,prc=prc,fasta=NULL,reg=NULL,
              categories=categories)
    ret

}

#' make PRC
#'
#' Normalizes a data matrix and then preformes quantile normalization
#' using \code{prcomp} on the transform of the normalized data. The
#' eigenvectors are returned in a list with the normalized data
#' @param data The matrix to be analyzed
#' @param norm The normalization method
#' norm may either be a user provided function or one of the
#' following
#' \describe{
#' \item{rowSumOne}{ all rows sum to 1}
#' \item{colSumOne}{ all columns sum to 1}
#' \item{normRow}{all rows have mean 0 and sd 1}
#' \item{normCol}{all cols have mean 0 and sd 1}
#' \item{rowsVarOne}{ all rows have sd 1}
#' \item{colsVarOne}{ all cols have sd 1}
#' \item{qn}{ quantile normalization}
#' \item{none}{ no normalization}
#' }
#' @return list($normData,$eigenVectors)
#' @examples
#' ## raw read file names
#' dataSets<-sapply(c("raw_sample1.bed","raw_sample2.bed","raw_sample3.bed"),
#'   function(file){
#'     system.file("extdata", file, package = "CCCA")
#'   })
#' ## load AFS
#' filename<-system.file("extdata","sample.afs", package = "CCCA")
#' sampleAFS<-readAFS(filename)
#' ## Generate UDM
#' udm<-makeUDM(sampleAFS,dataSets)
#' @export
makePRC<-function(data,norm="qn"){
    if( is.function(norm))
        normData<-norm(data)
    else
        normData<-switch(
            norm,
            rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
            colSumOne=apply(data,2, function(x) x/sum(x)),
            normRow=t(apply(data,1, function(x) (x-mean(x))/stats::var(x))),
            normCol=apply(data,2, function(x) (x-mean(x))/stats::var(x)),
            rowsVarOne=t(apply(data,1, function(x) x/stats::var(x))),
            colsVarOne=apply(data,2, function(x) x/stats::var(x)),
            qn=._qn(data),
            none=data,
            data)
    prc<-stats::prcomp(t(normData))$rotation
    ret<-list(normData=normData,eigenVectors=prc)
    ret
}

#' Normalize PRC
#'
#' Substract the mean and divide by the standard deviation.
#' @param x Vector or matrix to be normalized
#' @param ... additional arguments
#' @return an array wher len = length(x)
#' @examples
#' r<-runif(100,0,100)
#' normalize(r)
#' b=NULL
#' b$eigenVector<-matrix(r,10)
#' normalize(b)
#' @export
normalizePRC<-function(x,...){
    apply(x$eigenVector,2,._normalize)
}

#' map zip list
#'
#' Takes a list of lists and groups the nth member of each list
#' into a new one.
#' @param ... an arbitraraly long set of equal length lists
#' @examples
#' mapziplist(list(1,2,3),list("A","B","C"))
#' # [[1]]
#' # [1] "1" "A"
#' # 
#' # [[2]]
#' # [1] "2" "B"
#' # 
#' # [[3]]
#' # [1] "3" "C"
._mapziplist<-function(...){
    inlist<-list(...)
    l<-length(inlist[[1]])
    zip<-function(i)
        unlist(lapply(inlist,"[[",i))
    lapply(seq(l),zip)
}

#' Shift from zero
#' 
#' Checks if the result of subtracting width from value is zero
#'
#' @param value integer or double which cannot be less than 0
#' @param width integer or double
#' @return a range value+width value-width where value-width >0
._shiftFromZero<-function(summit){
    testBedSE<-cbind(summit-350,summit+350)
    x<-which(testBedSE[,1]<0)
    testBedSE[x,2]=testBedSE[x,2]-testBedSE[x,1]
    testBedSE[x,1]=0
    testBedSE
}

#' Add Region
#'
#' Add a region to the ccca object. These regions will be used as a mask to
#' identify seperate contexts.
#' @param x The object with a $reg value to be added
#' @param tag The moniker of the new region to be added
#' @param logic A logical list of equal length to other regions in x$reg
#' @param ... additional arguments
#' @return An array is appended to <ccca>$reg
#' @examples
#' ## See ?makeCCCA
#' @export
addRegion<-function(x, tag,logic,...){
    if(is.null(x$reg)){
        x$reg<-cbind(logic)
        colnames(x$reg)=tag
    }
    else{
        if(length(logic)!=dim(x$reg)[1])
            stop("Reg length not equal to region being added")
        logicp<-cbind(logic)
        colnames(logicp)<-tag
        x$reg<-cbind(x$reg,logicp)
    }
    x
}

#' Add Fasta
#'
#' Sets the fasta value of ccca to the neuclotide values of the supplied
#' genome between start and end. 
#' @param ccca The object used to generate the fasta strings
#' @param genome A BS.genome object 
#' @param width The width of fasta files, centered around the summit
#' of each peak in ccca
#' @param ... Extension aguments
#' @return ccca with fasta values
#' @examples
#' #' @examples
#' ## See ?makeCCCA
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' genome=BSgenome.Hsapiens.UCSC.hg19
#' addFasta(CCCA_prc,genome)
#' }
#' @export
addFasta<-function(ccca,genome,width=200,...){
    requireNamespace("Biostrings")
    if (is.null(ccca$afs$chr) | is.null(ccca$afs$start)) 
        stop("addFasta env list must contain afs$chr afs$start and afs$end")
    if (!requireNamespace("Biostrings")) 
        stop(paste0("Must install the Biostrings package from",
                    " Bioconductor.\n",
                    "source(\"https://bioconductor.org/biocLite.R\");",
                    " biocLite(\"Biostrings\")"))

    start<- (ccca$afs$start + ccca$afs$end)/2-floor(width/2)
    ccca$fasta <- getSeq(genome,
                         ccca$afs$chr,
                         start = start,
                         width = width)
    ccca
}

#' unify bed file
#'
#' unify bed file
._unifyBedFile<-function(dataset,overlapWidth, chr="chr",
                         summit="summit",name="name"){
    chrcats<-levels(dataset[,chr])
    l<-nrow(dataset)
    t1<-dataset[1:l-1,]
    t2<-dataset[2:l,]
    b<-(t1[,chr]==t2[,chr] & abs(as.numeric(t1[,summit])-
                                 as.numeric(t2[,summit]))<overlapWidth)
    prePeaks<-unlist(lapply(which(b==FALSE),function(x){c(x,x+1)}))
    peaks<-t(matrix(c(1,prePeaks[1:length(prePeaks)-1]),nrow=2))
    data<-._unityOutput(peaks,
                        as.integer(dataset[,chr]),
                        as.integer(dataset[,summit]),
                        as.integer(dataset[,name]))
    names(data)<-list(chr,summit,"matrix")
    data[[chr]]<-as.factor(data[[chr]])
    levels(data[[chr]])<-levels(dataset[,chr])
    colnames(data$matrix)<-levels(dataset[,name])
    od<-data.frame(chr=data[[chr]],summit=data[[summit]])
    od<-cbind(od,data$matrix)
    od
}

#' Quantile Normalization
#'
#' Applies quantile normalization to data in matrix form.
#' 
#' \itemize{
#' \item{1. Order each of the columns}
#' \item{2. Sum each orderd row and normalize by the width}
#' \item{3. Put the data back in order}
#' }
#' @param data matrix of data to be normalized
#' @return matrix of the same dimensions as data
#' @examples
#' matr<-cbind(rbind(5,2,3,4),rbind(4,1,4,2),rbind(3,4,6,8))
#' qn(matr)
._qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    matr<-matrix(unlist(lapply(seq(shape[2]),
                               function(i,x,y) x[y[,i],i],data,sequence)),
                 ncol=shape[2])
    ranks<-apply(matr,1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}

#' PCA tp Matrix Transformation
#'
#' Use before plotting the principal components
#' @param x the normalized input data and eigenvectors 
#' @param n normalizing function applied to the eigevectors
#' @return the dot product of the normalized data and eigenvectors
._pca2Matr<-function(x,n=function(x) x){
    if(is.null(x$normData) | is.null(x$eigenVectors))
        stop("x needs variables normData and eigenVectors")
    ## need to add tests for dimensions
    t(x$normData)%*%apply(x$eigenVectors,2,n)
}

#' Plot Principal Components
#'
#' Basic scatter-plot for visualizing similarities and differences
#' between cell types based on their principal components. To be used as
#' a template.
#' @param ccca class holding the UDM AFS and PRC
#' @param components two principal components
#' @export
plotPrincipalComponents<-function (ccca,components=c(1,2)){
    matrix<-ChipPComp:::._pca2Matr(ccca$prc)
    plot(matrix[,components],xlim=c(min(matrix[,components[1]]) -
                                    abs(min(matrix[,components[1]])) * 0.1,
                                    max(matrix[,components[1]]) +
                                    abs(max(matrix[,components[1]]) )* 0.3))
    text(matrix[,components], labels = ccca$categories, pos = 4)
}

#' Normalize
#'
#' privides unit normalization for numeric x
#' @param x numeric input
#' @return numeric of length x
#' @examples
#' x<-rnorm(100,10,2)
#' y<-normalize(x)
#' print(c(mean(y),sd(y)))
._normalize<-function(x){
    norm<-function(x){
        (x-mean(x))/(sqrt(stats::var(x)))}
    if(is.null(dim(x)))
        norm(x)
    else
        apply(x,2,norm)
}
