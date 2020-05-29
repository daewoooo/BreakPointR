#' Genotype the space between breakpoints
#'
#' Genotype the space between breakpoints with Fisher's exact test.
#'
#' @param breaks A \code{\link{GRanges-class}} object with breakpoint coordinates.
#' @param fragments A \code{\link{GRanges-class}} object with read fragments.
#' @param background The percent (e.g. 0.05 = 5\%) of background reads allowed for WW or CC genotype calls.
#' @param minReads The minimal number of reads between two breaks required for genotyping.
#' @return A \code{\link{GRanges-class}} object with genotyped breakpoint coordinates.
#' @author David Porubsky, Ashley Sanders, Aaron Taudt
#' @importFrom stats fisher.test
#' @importFrom S4Vectors endoapply
#' @export
#' @examples
#'## Get an example file 
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Load the file 
#'breakpoint.objects <- get(load(exampleFile))
#'## Genotype regions between breakpoints
#'gbreaks <- GenotypeBreaks(breakpoint.objects$breaks, breakpoint.objects$fragments)
#'
GenotypeBreaks <- function(breaks, fragments, background=0.05, minReads=10) {
    if (length(breaks)==0) {
        stop("argument 'breaks' is empty")
    }
  
    breakrange.list <- GenomicRanges::GRangesList()
    seqlevels(breakrange.list) <- seqlevels(breaks)
    for (chrom in unique(seqnames(breaks))) {
          
        # create ranges between the breakpoints -> start and stops in a dataframe, use this to genotype between
        breaks.strand <- breaks[seqnames(breaks) == chrom]
        strand(breaks.strand) <- '*'
        # Remove the non-used seqlevels
        breaks.strand <- keepSeqlevels(breaks.strand, value=chrom)
        breakrange <- GenomicRanges::gaps(breaks.strand)
        breakrange <- breakrange[strand(breakrange)=='*']
        
        ## pull out reads of each line, and genotype in the fragments
        strand(breakrange) <- '-'
        breakrange$Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
        strand(breakrange) <- '+'
        breakrange$Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
        strand(breakrange) <- '*'
        breakrange$readNo <- breakrange$Ws + breakrange$Cs
        
        ## bestFit genotype each region by Fisher Exact test
        fisher <- lapply(1:length(breakrange), function(x) genotype.fisher(cReads=breakrange$Cs[x], wReads=breakrange$Ws[x], roiReads=breakrange$readNo[x], background=background, minReads=minReads))
        fisher <- do.call(cbind, fisher)
      
        breakrange$genoT <- unlist(fisher[1,])
        breakrange$pVal <- unlist(fisher[2,])
        breakrange <- breakrange[!is.na(breakrange$genoT)]
      
        ## remove break if genotype is the same on either side of it
        equal.on.either.side <- c(breakrange$genoT[-length(breakrange)] == breakrange$genoT[-1], TRUE)
        breakrange.new <- breakrange[!equal.on.either.side]
        start(breakrange.new) <- end(breakrange[!equal.on.either.side])
        end(breakrange.new) <- start(breakrange[which(!equal.on.either.side)+1])
        breakrange.new$genoT <- paste(breakrange$genoT[!equal.on.either.side], breakrange$genoT[which(!equal.on.either.side)+1], sep='-')
      
        mcols(breakrange.new)[c('Ws','Cs','readNo','pVal')] <- NULL
  
        breakrange.list[[chrom]] <- breakrange.new
        
    }
    breakrange.new <- unlist(breakrange.list, use.names=FALSE)
    
    if (length(breakrange.new)) {
        ## Refine breakpoint regions to the highest deltaW in the given region
        hits <- GenomicRanges::findOverlaps(breakrange.new, breaks)
        ToRefine <- split(breaks[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))
        refined <- unlist(endoapply(ToRefine, RefineBreaks), use.names = FALSE)
        ranges(breakrange.new) <- ranges(refined)
        breakrange.new$deltaW <- refined$deltaW
        return(breakrange.new)
    } else {
        return(breakrange.new<-NULL)
    }  
}


#' Genotype read numbers
#'
#' Genotype read numbers with Fisher's exact test.
#' 
#' @param cReads Number of Crick reads.
#' @param wReads Number of Watson reads.
#' @param roiReads Total number of reads.
#' @param background Parameter for frequency of background reads.
#' @param minReads Minimal number of reads to perform the test.
#' @return A list with the $bestFit and $pval.
#' @author Ashley Sanders, David Porubsky, Aaron Taudt
#' 
genotype.fisher <- function(cReads, wReads, roiReads, background=0.02, minReads=10) {
    ## FISHER EXACT TEST
    result <- list(bestFit=NA, pval=NA)
    if (length(roiReads)==0) {
        return(result)
    }
    if (is.na(roiReads)) {
        return(result)
    }
    if ( roiReads >= minReads ) {
        m <- matrix(c(cReads, wReads, round(roiReads*(1-background)), round(roiReads*background)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
        CCpVal <- stats::fisher.test(m, alternative="greater")[[1]]
        m <- matrix(c(cReads, wReads, round(roiReads*0.5), round(roiReads*0.5)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Cs','Ws')))
        WCpVal <- 1 - stats::fisher.test(m, alternative="two.sided")[[1]]
        m <- matrix(c(wReads, cReads, round(roiReads*(1-background)), round(roiReads*background)), ncol=2, byrow=TRUE, dimnames=list(case=c('real', 'reference'), reads=c('Ws','Cs')))
        WWpVal <- stats::fisher.test(m, alternative="greater")[[1]]
        
        pVal <- c(wc=WCpVal, cc=CCpVal, ww=WWpVal)
        result <- list(bestFit=names(pVal)[which.min(pVal)], pval=min(pVal))
        return(result)
        
    } else { 
        return(result)
    }
}  

####################
# Helper functions #
####################

RefineBreaks <- function(gr) {
    maxDeltaW <- max(gr$deltaW)
    new.gr <- gr[gr$deltaW == maxDeltaW]
    new.gr <- mergeGR(new.gr)
    new.gr$deltaW <- maxDeltaW
    return(new.gr)
}

mergeGR <- function(gr) {
    gr <- sort(gr)
    #new.gr <- GRanges(seqnames=as.character(seqnames(gr))[1], ranges=IRanges(start=start(gr[1]), end=end(gr[length(gr)])))
    new.gr <- GRanges(seqnames=as.character(seqnames(gr))[1], ranges=IRanges(start=min(start(gr)), end=max(end(gr))))
    return(new.gr)
}
