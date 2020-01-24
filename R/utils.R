#' Collapse consecutive bins with the same ID value
#'
#' Collapse consecutive bins with the same value defined in 'id.field'. 
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @param id.field A number of metadata column to use for region merging.
#' @return A \code{\link{GRanges-class}} object.
collapseBins <- function(gr, id.field=3) {
    ##get indices of last range in a consecutive(RLE) run of the same value
    ind.last <- cumsum(S4Vectors::runLength(S4Vectors::Rle(mcols(gr)[,id.field])))
    ##get indices of first range in a consecutive(RLE) run of the same value
    ind.first <- c(1,cumsum(S4Vectors::runLength(S4Vectors::Rle(mcols(gr)[,id.field]))) + 1)
    ind.first <- ind.first[-length(ind.first)]  ##erase last index from first range indices 
    collapsed.gr <- GenomicRanges::GRanges(seqnames=seqnames(gr[ind.first]), ranges=IRanges(start=start(gr[ind.first]), end=end(gr[ind.last])), mcols=mcols(gr[ind.first]))
    names(mcols(collapsed.gr)) <- names(mcols(gr[ind.first]))
    return(collapsed.gr)
}

#' Transform genomic coordinates
#'
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} with two additional metadata columns 'start.genome' and 'end.genome'.
transCoord <- function(gr) {
    cum.seqlengths <- cumsum(as.numeric(seqlengths(gr)))
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    names(cum.seqlengths.0) <- GenomeInfoDb::seqlevels(gr)
    gr$start.genome <- start(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    gr$end.genome <- end(gr) + cum.seqlengths.0[as.character(seqnames(gr))]
    return(gr)
}

#' Insert chromosome for in case it's missing
#' 
#' Add two columns with transformed genomic coordinates to the \code{\link{GRanges-class}} object. This is useful for making genomewide plots.
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} object with an additional metadata column containing chromosome name with 'chr'. 
insertchr <- function(gr) {
    mask <- which(!grepl('chr', seqnames(gr)))
    mcols(gr)$chromosome <- as.character(seqnames(gr))
    mcols(gr)$chromosome[mask] <- sub(pattern='^', replacement='chr', mcols(gr)$chromosome[mask])
    mcols(gr)$chromosome <- as.factor(mcols(gr)$chromosome)
    return(gr)
}

#' Process double SCE chromsomes: with internal WC region.
#' 
#' This function will take from a double SCE chromosome only  WW or CC region (Longer region is taken).
#'
#' @param gr A \code{\link{GRanges-class}} object.
#' @return The input \code{\link{GRanges-class}} object with only WW or CC region retained. 
removeDoubleSCEs <- function(gr) {
    if (any(gr$states == 'wc')) {
        wc.idx <- which(gr$states == 'wc')
        if (any(wc.idx > 1) & any(wc.idx < length(gr))) {
            ## Remove wc regions
            gr.new <- gr[-wc.idx ]
            ## Take strand state covering largest region
            gr.new.byState <- split(gr.new, gr.new$states)
            state.widths <- lapply(gr.new.byState, function(x) sum(width(x)))
            gr.new <- gr.new[gr.new$states == names(which.max(state.widths))]
            return(gr.new)
        } else {
            return(gr[-wc.idx ])
        }
    } else {
        return(gr)
    }
}