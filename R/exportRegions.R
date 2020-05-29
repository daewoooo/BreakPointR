#' Function to print WC regions after breakpointR analysis
#'
#' @param datapath A path to that 
#' @param file A filename to print exported regions to.
#' @param collapseInversions Set to \code{TRUE} if you want to collapse putative inverted regions.
#' @param collapseRegionSize Upper range of what sized regions should be collapsed.
#' @param minRegionSize Minimal size of the region to be reported.
#' @param state A genotype of the regions to be exported ('ww', 'cc' or 'wc').
#' @return A \code{data.frame} object containing all regions with user defined 'state'.
#' @importFrom gtools mixedorder
#' @author David Porubsky
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'## To export regions genotyped as 'wc'
#'wc.regions <- exportRegions(datapath=exampleFolder, collapseInversions=FALSE, minRegionSize=5000000, state='wc') 

exportRegions <- function(datapath, file=NULL, collapseInversions=FALSE, collapseRegionSize=5000000, minRegionSize=5000000, state='wc') {

    ptm <- startTimedMessage("Exporting regions with ",state, " state ...")
  
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)

    ranges <- GenomicRanges::GRangesList()
    for (filename in files) {
        counts <- get(load(filename))$counts

        ## Collapse inversion !!!SO far tested only for wc region export!!!
        if (collapseInversions & collapseRegionSize > 0 & state=='wc') {
            grl <- GenomicRanges::split(counts, seqnames(counts)) #split GRanges into a list of ranges per chromosome
            gr.all <- GenomicRanges::GRanges()
            for (i in 1:length(grl)) {
                gr.chr <- grl[[i]]

                ## Get index of each region smaller than specified 'regionSize' with WW or CC inherirtance
                indexes <- which(width(gr.chr) < collapseRegionSize & (gr.chr$states == 'ww' | gr.chr$states == 'cc'))

                ## Those indexes (ww/cc regions) that are flanked by wc regions set to wc as well
                for (index in indexes) {
                    if (length(gr.chr[index-1]$states) == 0 | index == length(gr.chr)) next
                    if (gr.chr[index-1]$states == 'wc' & gr.chr[index+1]$states == 'wc') {
                        gr.chr[index]$states <- 'wc'
                    }
                }

                ## Collapse consecutive regions with WC regions
                cumsum(S4Vectors::runLength(S4Vectors::Rle(gr.chr$states))) -> ind.last  ##get indices of last range in a consecutive(RLE) run of the same haplotype 
                c(1,cumsum(S4Vectors::runLength(S4Vectors::Rle(gr.chr$states))) + 1) -> ind.first ##get indices of first range in a consecutive(RLE) run of the same haplotype
                ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices
                gr <- GenomicRanges::GRanges(seqnames=seqnames(gr.chr[ind.first]), ranges=IRanges(start=start(gr.chr[ind.first]), end=end(gr.chr[ind.last])), mcols=mcols(gr.chr[ind.first]))
                gr.all <- IRanges::append(gr.all, gr)
            }

            ## Filter only wc regios and store them GRangeslist per filename
            counts.filt <- gr.all[width(gr.all) >= minRegionSize & gr.all$mcols.states == state]
            if (length(counts.filt) > 0) {
                ranges[[filename]] <- counts.filt 
            }
        } else {
            counts.filt <- counts[width(counts) >= minRegionSize & counts$states == state]
            ranges[[filename]] <- counts.filt
        } 
    }

    ## Print selected wc regions into a text file
    ranges <- unlist(ranges)
    sort <- GenomicRanges::sort(ranges)
    filenames <- gsub(".RData", "", basename(names(ranges)))
    names(ranges) <- NULL
    ranges$filename <- filenames

    df <- as.data.frame(ranges)
    df2print <- df[,c('seqnames','start','end','filename')]
    
    if (!is.null(file)) {
        utils::write.table(df2print, file, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
        stopTimedMessage(ptm)
    } else {
        stopTimedMessage(ptm)
        return(df2print)
    }
}
