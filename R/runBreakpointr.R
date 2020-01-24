#' Find breakpoints in Strand-seq data
#'
#' Find breakpoints in Strand-seq data. See section Details on how breakpoints are located.
#'
#' Breakpoints are located in the following way:
#' \enumerate{
#' \item calculate deltaWs chromosome-by-chromsome
#' \item localize breaks that pass zlim above the threshold
#' \item genotype both sides of breaks to confirm whether strand state changes
#' \item write a file of _reads, _deltaWs and _breaks in a chr fold -> can upload on to UCSC Genome browser
#' \item write a file for each index with all chromosomes included -> can upload on to UCSC Genome browser
#' }
#'
#' @param bamfile A file with aligned reads in BAM format.
#' @param ID A character string that will serve as identifier in downstream functions.
#' @param windowsize The window size used to calculate deltaWs, either number of reads or genomic size depending on \code{binMethod}.
#' @inheritParams readBamFileAsGRanges
#' @param binMethod Method used to calculate optimal number of reads in the window ("size", "reads"). By default \code{binMethod='size'}.
#' @inheritParams breakSeekr
#' @inheritParams GenotypeBreaks
#' @param maskRegions List of regions to be excluded from the analysis (tab-separated file: chromosomes start end).
#' @inheritParams confidenceInterval
#' @return A \code{\link{BreakPoint}} object.
#' @author David Porubsky, Ashley Sanders, Aaron Taudt
#' @importFrom utils flush.console
#' @importFrom S4Vectors subjectHits
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_bams", package="breakpointRdata")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Run breakpointR
#'brkpts <- runBreakpointr(exampleFile, chromosomes='chr22', pairedEndReads=FALSE)
#'
runBreakpointr <- function(bamfile, ID=basename(bamfile), pairedEndReads=TRUE, chromosomes=NULL, windowsize=1e6, binMethod="size", trim=10, peakTh=0.33, zlim=3.291, background=0.05, min.mapq=10, pair2frgm=FALSE, filtAlt=FALSE, genoT='fisher', minReads=20, maskRegions=NULL, conf=0.99) {

    ## check the class of the bamfile, make GRanges object of file
    if (!is(bamfile, 'GRanges')) {  
        ptm <- startTimedMessage("Reading file ", bamfile, " ...")
        suppressWarnings( fragments <- readBamFileAsGRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, min.mapq=min.mapq, pair2frgm=pair2frgm, filtAlt=filtAlt) )
        stopTimedMessage(ptm)
    } else {
        fragments <- bamfile
        bamfile <- 'CompositeFile'
    }
  
    filename <- basename(bamfile)
    message("Working on ", filename)
    
    ## remove reads from maskRegions
    if (!is.null(maskRegions)) {
        mask <- IRanges::findOverlaps(maskRegions, fragments)
        fragments <- fragments[-S4Vectors::subjectHits(mask)]
    } 
  
    reads.all.chroms <- fragments
  
    deltas.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(deltas.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(deltas.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    breaks.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(breaks.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(breaks.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    confint.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(confint.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(confint.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    counts.all.chroms <- GenomicRanges::GRangesList()
    GenomeInfoDb::seqlevels(counts.all.chroms) <- GenomeInfoDb::seqlevels(fragments)
    GenomeInfoDb::seqlengths(counts.all.chroms) <- GenomeInfoDb::seqlengths(fragments)
    for (chr in unique(seqnames(fragments))) {
        message("  Working on chromosome ",chr)
        fragments.chr <- fragments[seqnames(fragments)==chr]
        fragments.chr <- GenomeInfoDb::keepSeqlevels(fragments.chr, chr)
    
        ptm <- startTimedMessage("    calculating deltaWs ...")
        if (binMethod == "size") {
            ## normalize for size of chromosome one and for read counts of each chromosome
            if (windowsize > GenomeInfoDb::seqlengths(fragments)[chr]) {
                stopTimedMessage(ptm)
                next
            }
            tiles <- unlist(GenomicRanges::tileGenome(seqlengths(fragments)[chr], tilewidth = windowsize))
            counts <- GenomicRanges::countOverlaps(tiles, fragments.chr)
            reads.per.window <- max(10, round(mean(counts[counts>0], trim=0.05)))
          
            dw <- deltaWCalculator(frags=fragments.chr, reads.per.window=reads.per.window)
        } else if (binMethod == "reads") {
            ## normalize only for size of the chromosome 1
            #reads.per.window <- max(10, round(windowsize/(seqlengths(fragments)[1]/seqlengths(fragments)[chr]))) # scales the bin to chr size, anchored to chr1 (249250621 bp)
            ## do not normalize to the size of the chromosome 1
            reads.per.window <- windowsize 
            dw <- deltaWCalculator(frags=fragments.chr, reads.per.window=reads.per.window)
        }
        deltaWs <- dw[seqnames(dw)==chr]
        stopTimedMessage(ptm)
    
        ptm <- startTimedMessage("    finding breakpoints ...")
        breaks <- suppressMessages( breakSeekr(deltaWs, trim=trim, peakTh=peakTh, zlim=zlim) )
        stopTimedMessage(ptm)
    
        ## Genotyping
        if (length(breaks) > 0 ) {
            maxiter <- 10
            iter <- 1
            ptm <- startTimedMessage("    genotyping ",iter, " ...")
            utils::flush.console()
            newBreaks <- GenotypeBreaks(breaks=breaks, fragments=fragments, background=background, minReads=minReads, genoT=genoT)
            prev.breaks <- breaks  
            breaks <- newBreaks
            stopTimedMessage(ptm)
            
            while (length(prev.breaks) > length(newBreaks) && !is.null(breaks) ) {
                utils::flush.console()
                iter <- iter + 1
                ptm <- startTimedMessage("    genotyping ",iter, " ...")
                newBreaks <- GenotypeBreaks(breaks=breaks, fragments=fragments, background=background, minReads=minReads, genoT=genoT)
                prev.breaks <- breaks  
                breaks <- newBreaks
                stopTimedMessage(ptm)
                if (iter == maxiter) { break }
            }
            
        } else {
            newBreaks <- NULL # assigns something to newBreaks if no peaks found
            confint <- NULL
        }
    
        ### count reads between the breaks and write into GRanges
        if (is.null(newBreaks)) {
            Ws <- length(which(as.logical(strand(fragments.chr) == '-')))
            Cs <- length(which(as.logical(strand(fragments.chr) == '+')))
            chrRange <- GenomicRanges::GRanges(seqnames=chr, ranges=IRanges(start=1, end=GenomeInfoDb::seqlengths(fragments.chr)[chr]))
            counts <- cbind(Ws,Cs)
            mcols(chrRange) <- counts
            seqlengths(chrRange) <- GenomeInfoDb::seqlengths(fragments.chr)[chr]
      
            ## genotype entire chromosome
            WC.ratio <- (chrRange$Ws-chrRange$Cs)/sum(c(chrRange$Ws,chrRange$Cs))
            if (WC.ratio > 0.8) {
                state <- 'ww'
            } else if (WC.ratio < -0.8) {
                state <- 'cc'
            } else if (WC.ratio < 0.2 & WC.ratio > -0.2) {
                state <- 'wc'
            } else {
                state <- '?'
            }          
            chrRange$states <- state
            suppressWarnings( counts.all.chroms[[chr]] <- chrRange )
            
            confint <- NULL
      
        } else {
          
            breaks.strand <- newBreaks
            strand(breaks.strand) <- '*'
            breakrange <- GenomicRanges::gaps(breaks.strand)
            breakrange <- breakrange[strand(breakrange)=='*']
      
            ## pull out reads of each line, and genotype in the fragments
            strand(breakrange) <- '-'
            Ws <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- '+'
            Cs <- GenomicRanges::countOverlaps(breakrange, fragments)
            strand(breakrange) <- '*'
      
            ## pull out states for each region between the breakpoints
            concat.states <- paste(newBreaks$genoT, collapse="-")
            split.states <- unlist(strsplit(concat.states, "-"))
            states.idx <- seq(from=1,to=length(split.states), by=2)
            states.idx <- c(states.idx, length(split.states))
            states <- split.states[states.idx]
      
            counts <- cbind(Ws,Cs)
            mcols(breakrange) <- counts
            breakrange$states <- states
            suppressWarnings( counts.all.chroms[[chr]] <- breakrange )
            
            ## Confidence intervals
            ptm <- startTimedMessage("    confidence intervals ...")
            if (bamfile != 'CompositeFile') {
                # confint <- confidenceInterval(breaks = newBreaks, fragments = fragments, background = background, conf = conf)
                confint <- confidenceInterval.binomial(breaks = newBreaks, fragments = fragments, background = background, conf = conf)
            } else {
                # use nonbinomial test to calculate confidence intervals for composite file
                confint <- confidenceInterval.binomial(breaks = newBreaks, fragments = fragments, background = background, conf = conf)
            }
            stopTimedMessage(ptm)
        }
    
        ### write breaks and deltas into GRanges
        if (length(deltaWs) > 0) {
            deltas.all.chroms[[chr]] <- deltaWs[,'deltaW']  #select only deltaW metadata column to store
        }
        if (length(newBreaks) > 0) {
            breaks.all.chroms[[chr]] <- newBreaks
        }
        if (length(confint) > 0) {
            confint.all.chroms[[chr]] <- confint
        }
        
    }
    
    ## creating list of list where filename is first level list ID and deltas, breaks and counts are second list IDs
    deltas.all.chroms <- unlist(deltas.all.chroms, use.names=FALSE)
    breaks.all.chroms <- unlist(breaks.all.chroms, use.names=FALSE)
    confint.all.chroms <- unlist(confint.all.chroms, use.names=FALSE)
    counts.all.chroms <- unlist(counts.all.chroms, use.names=FALSE)
    
    #seqlevels(deltas.all.chroms) <- seqlevels(reads.all.chroms)
    #seqlevels(breaks.all.chroms) <- seqlevels(reads.all.chroms)
    #seqlevels(breaks.all.chroms) <- seqlevels(reads.all.chroms)
    #seqlevels(counts.all.chroms) <- seqlevels(reads.all.chroms)
    
    ## Estimate background reads
    counts.all.chroms[counts.all.chroms$states == 'ww'] -> ww
    counts.all.chroms[counts.all.chroms$states == 'cc'] -> cc
    if (length(ww) > 0) {
        #add pseudocounts to avoid working with zeros
        ww$Cs <- ww$Cs+1
        ww$Ws <- ww$Ws+1
        bg.estim.ww <-  sum(ww$Cs) / sum(ww$Ws)
    } else {
        bg.estim.ww <- 0
    }  
    
    if (length(cc) > 0) {
        #add pseudocounts to avoid working with zeros
        cc$Ws <- cc$Ws+1
        cc$Cs <- cc$Cs+1
        bg.estim.cc <- sum(cc$Ws) / sum(cc$Cs)
    } else {
        bg.estim.cc <- 0
    } 
    bg.estim <- mean(bg.estim.ww, bg.estim.cc) 
    
    ## Calculate reads per megabase
    tiles <- unlist(GenomicRanges::tileGenome(seqlengths(fragments), tilewidth = 1000000))
    counts <- GenomicRanges::countOverlaps(tiles, fragments)
    reads.MB <- round(median(counts))
    
    ## Calculate % of the genome covered
    red.frags <- GenomicRanges::reduce(fragments)
    perc.cov <- ( sum(as.numeric(width(red.frags))) / sum(as.numeric(seqlengths(red.frags))) )*100
    perc.cov <- round(perc.cov, digits = 2)
    
    library.metrics <- c(background.estimate=bg.estim, med.reads.per.MB=reads.MB, perc.coverage=perc.cov)

    ## save set parameters for future reference
    parameters <- c(windowsize=windowsize, binMethod=binMethod, trim=trim, peakTh=peakTh, zlim=zlim, background=background, minReads=minReads, confidence.int=conf)
    
    ## save brekpointR results into a single obejct
    fragments <- fragments[,'mapq'] #store only mapq values to save space
    #breaks.all.chroms$confint <- confint.all.chroms[,0] #store confidence intervals as a part of breaks object
    data.obj <- list(ID=ID, fragments=fragments, deltas=deltas.all.chroms, breaks=breaks.all.chroms, confint=confint.all.chroms, counts=counts.all.chroms, lib.metrics=library.metrics, params=parameters)
    class(data.obj) <- class.breakpoint
  
    return(data.obj)
}