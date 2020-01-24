#' Synchronize Strand-seq read directionality
#'
#' This function aims to synchronize strand directionality of reads that fall into WW and CC regions.
#' 
#' @param files2sync A list of files that contains \code{\link{BreakPoint}} objects.
#' @param collapseWidth A segment size to be collapsed with neighbouring segments.
#' @return A \code{\link{GRanges-class}} object that reads synchronized by directionality.
#' @importFrom S4Vectors endoapply
#' @author David Porubsky
#' @export
#' @examples
#'## Get some files that you want to load
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'files2sync <- list.files(exampleFolder, full.names=TRUE)[1]
#'synchronizeReadDir(files2sync=files2sync)
#'
synchronizeReadDir <- function(files2sync, collapseWidth=5000000) {
  
    ## helper functions
    #switch strand directionality for WW regions
    switchStrand <- function(gr) {
        max.strand <- names(which.max(BiocGenerics::table(strand(gr))))
        #switch strand if max directionality is Watson("-")
        if (max.strand == "-") {
            recode.gr <- gr
            strand(recode.gr[strand(gr) == "+"]) <- "-"
            strand(recode.gr[strand(gr) == "-"]) <- "+"
            return(recode.gr)
        } else {
            return(gr)
        }
    }
    
    ptm <- startTimedMessage("Synchronizing read directionality ...")

    ## load data of class BreakPoint
    data <- loadFromFiles(files2sync)

    allLibs.syncReads <- GenomicRanges::GRangesList()
    for (i in seq_along(data)) {
        lib.results <- data[[i]]
        region.counts <- lib.results$counts
        ## collapse regions smaller than collapseWidth
        #get index of each region smaller than collapseWidth
        indexes <- which(width(region.counts) < collapseWidth)
        #collapse regions flanked by the same state
        for (idx in indexes) {
            #skip if index correspond to a first or last range
            if (idx == 1 | idx == length(region.counts)) next
            #check if regions flanking selected regions are coming from the same chromosome
            if (as.character(seqnames(region.counts[idx-1])) == as.character(seqnames(region.counts[idx+1]))) {
                if (region.counts[idx-1]$states == region.counts[idx+1]$states) {
                    region.counts[idx]$states <- region.counts[idx-1]$states
                }
            }
        }
        
        ## For WW and CC regions separated by WC region take only the larger one
        regions.per.chr <- split(region.counts, as.character(seqnames(region.counts)))
        WWandCC.regions.grl <- endoapply(regions.per.chr, removeDoubleSCEs)
        #WWandCC.regions <- region.counts[region.counts$states != 'wc']
        #WWandCC.regions.grl <- GenomicRanges::split(WWandCC.regions, seqnames(WWandCC.regions))
        WWandCC.regions.collapsed <- endoapply(WWandCC.regions.grl, collapseBins)
        WWandCC.regions.collapsed <- unlist(WWandCC.regions.collapsed, use.names=FALSE)
        ## get fraction of reads with minority direction
        #WsandCs.df <- data.frame(Ws=WWandCC.regions$Ws, Cs=WWandCC.regions$Cs)
        #WsandCs.mat <- t(apply(df, 1, sort))
        #region.background <- WsandCs.mat[,1]/rowSums(WsandCs.mat)
        
        ## get fragments (reads) falling into WW and CC regions
        WWandCC.regions.frags <- IRanges::subsetByOverlaps(lib.results$fragments, WWandCC.regions.collapsed)
        WWandCC.regions.frags.grl <- GenomicRanges::split(WWandCC.regions.frags, seqnames(WWandCC.regions.frags))
        ## synchronize strand => max strand equals reference strand
        WWandCC.regions.frags.grl <- endoapply(WWandCC.regions.frags.grl, switchStrand)
        WWandCC.regions.frags.sync <- unlist(WWandCC.regions.frags.grl, use.names=FALSE)
        allLibs.syncReads[[i]] <- WWandCC.regions.frags.sync
    }
    syncReads <- unlist(allLibs.syncReads)
    syncReads <- GenomicRanges::sort(syncReads, ignore.strand=TRUE)
    
    stopTimedMessage(ptm)
    
    return(syncReads)
}
