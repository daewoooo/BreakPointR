#' Estimate copy number (CN) in genomic regions
#'
#' This function calcualtes ratio between observed and expected number of reads in localized genomic regions as CN estimates.
#'
#' @inheritParams deltaWCalculator
#' @param regions A \code{\link{GRanges-class}} object containing regions defined by breakpoint boundaries.
#' @inheritParams breakSeekr
#' @return A \code{\link{GRanges-class}} object containing submitted regions with an extra metacolumn \code{CN.estim}.
#' @author David Porubsky

estimateCN <- function(frags, regions, zlim=3.291) {
    regions$CN.estim <- 0
    regions$zscore <- 0
    for (i in seq_along(regions)) {
        region <- regions[i]
        ## Skip regions larger than 1% of given chromosome size
        if (width(region) / seqlengths(region) <= 0.01) {
            observed.counts <- region$Ws + region$Cs
            ## Select reads covering genotyped region
            cov.region <- subsetByOverlaps(frags, region)
            ## Get number of covered bases by removing large gaps in coverage
            strand(cov.region) <- "*"
            coverage.gaps <- gaps(cov.region)
            ## Remove gaps from the start and the end of the selected region
            coverage.gaps <- coverage.gaps[strand(coverage.gaps) == "*"]
            coverage.gaps <- coverage.gaps[-c(1,length(coverage.gaps))]
            ## Calculate coverage of gaps in between the reads
            coverage.gaps$percCov <- ( width(coverage.gaps)/width(region) )*100
            filt <- quantile(coverage.gaps$percCov, 0.99) #remove 1% of the largest gaps
            coverage.gaps <- coverage.gaps[coverage.gaps$percCov > filt]
            ## Remove gaps from the total region size if they do not cover more 50% of the given region
            if (sum(width(coverage.gaps)) / width(region) <= 0.5) {
                cov.region <- width(region) - sum(width(coverage.gaps))
            } else {
                cov.region <- width(region)
            }
            #noCov.region <- sum(width(coverage.gaps[coverage.gaps$percCov > 1]))
            #cov.region <- width(region) - noCov.region #remove the largest gaps from the total region span
            ## Calculate expected number of reads (median) in equaly sized regions across given chromosome
            binned <- unlist(tileGenome(seqlengths=seqlengths(frags), tilewidth = cov.region))
            binned.counts <- countOverlaps(binned, frags)
            med.binned.counts <- median(binned.counts)
            ## Calculate zscore to detect regions we are confidently above or below average coverage
            zscore <- (observed.counts -  med.binned.counts) / sd(binned.counts) 
            ## Calcualte ratio of observed versus expected number of reads
            ratio <- observed.counts / med.binned.counts
            regions$CN.estim[i] <- round(ratio, digits = 2)
            regions$zscore[i] <-  zscore
        } else {
            ## For extremely large regions assign zscore 0
            regions$zscore[i] <-  0
        }
    }
    ## Mark regions that are above or below the set confidence level (zlim)
    regions$zscore[is.nan(regions$zscore)] <- 0 #Coerce zscore == NaN to 0
    regions$copy <- 'normal'
    regions$copy[regions$zscore > zlim] <- 'high'
    regions$copy[regions$zscore < -zlim] <- 'low'
    
    return(regions)
}   