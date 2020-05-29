#' Calculate deltaWs
#'
#' This function will calculate deltaWs from a \code{\link{GRanges-class}} object with read fragments.
#'
#' @param frags A \code{\link{GRanges-class}} with read fragments (see \code{\link{readBamFileAsGRanges}}).
#' @param reads.per.window Number of reads in each dynamic window.
#' @return The input \code{frags} with additional meta-data columns.
#' @import GenomicRanges
#' @importFrom BiocGenerics as.vector
#' @author Aaron Taudt
#' @seealso readBamFileAsGRanges
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_bams", package="breakpointRdata")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Load the file
#'fragments <- readBamFileAsGRanges(exampleFile, pairedEndReads=FALSE, chromosomes='chr22')
#'## Calculate deltaW values
#'dw <- deltaWCalculator(fragments)
#' 
deltaWCalculator <- function(frags, reads.per.window=100) {

    if (reads.per.window == 0) {
        stop("'reads.per.window' must be >= 1")
    }
    if (reads.per.window < 10) {
        warning("'reads.per.window' should at least be 10")
    }
    frags.split <- split(frags, seqnames(frags))
    reads.per.chrom <- vapply(frags.split, FUN=length, FUN.VALUE=numeric(1))
    chroms2parse <- names(reads.per.chrom)[reads.per.chrom>2*reads.per.window]
    chroms2skip <- setdiff(names(reads.per.chrom),chroms2parse)
    if (length(chroms2skip)>0) {
        warning(paste0("Not parsing chromosomes ",paste(chroms2skip, collapse=',')," because they do not have enough reads."))
    }
    if (length(chroms2parse)==0) {
        warning("None of the specified chromosomes has enough reads. Doing nothing.")
        return(GRanges())
    }

    frags.new <- GRangesList()
    for (chrom in chroms2parse) {
        f <- frags.split[[chrom]]
        f <- f[order(start(f))]
        f$pcsum <- cumsum(strand(f)=='+')
        f$mcsum <- cumsum(strand(f)=='-')
        f$preads <- c(rep(NA,reads.per.window),diff(BiocGenerics::as.vector(f$pcsum),lag=reads.per.window))
        f$mreads <- c(rep(NA,reads.per.window),diff(BiocGenerics::as.vector(f$mcsum),lag=reads.per.window))
        f$deltaW <- abs(c(diff(f$mreads,lag=reads.per.window),rep(NA,reads.per.window)))
        
        # Shift deltaWs to region between reads
        start.f <- end(f)
        end.f <- c(start(f)[-1],seqlengths(frags)[chrom])
        mask <- start.f < end.f
        f <- f[mask]
        start(f) <- start.f[mask]
        end(f) <- end.f[mask]
        frags.new[[chrom]] <- f
    }
    frags.new <- unlist(frags.new)
    names(frags.new) <- NULL
    # Replace NAs with 0 to avoid problems in downstream functions
    frags.new$deltaW[is.na(frags.new$deltaW)] <- 0
    frags.new$mreads[is.na(frags.new$mreads)] <- 0
    frags.new$preads[is.na(frags.new$preads)] <- 0

    return(frags.new)
}


#' Calculate deltaWs using various window sizes
#'
#' This function will calculate deltaWs from a \code{\link{GRanges-class}} object with read fragments.
#'
#' @param sizes User defined multiplications of the original window size.
#' @inheritParams deltaWCalculator
#' @return The input \code{frags} with additional meta-data columns.
#' @import GenomicRanges
#' @importFrom BiocGenerics as.vector
#' @author David Porubsky
#' @seealso deltaWCalculator

deltaWCalculatorVariousWindows <- function(frags, reads.per.window=100, sizes=c(2,4,6)) {
  
  dw.per.size <- list()
  ## Calculate deltaW for an original window size
  dw <- deltaWCalculator(frags=frags, reads.per.window=reads.per.window)
  dw.per.size[[1]] <- dw$deltaW
  for (i in seq_along(sizes)) {
    reads.per.window.rescaled <- reads.per.window * sizes[i]
    ## Make sure that the rescaled number of reads per window is not more than 5% of all available fragments
    if (reads.per.window.rescaled <= (length(frags) * 0.05)) {
      ## Calculate deltaW for a user defined window
      dw <- deltaWCalculator(frags=frags, reads.per.window=reads.per.window * sizes[i])
      ## Normalize calculated deltaW for a given window by the original window size
      dw.per.size[[length(dw.per.size) + 1]] <- dw$deltaW / sizes[i]
    }  
  }
  
  if (length(dw.per.size) > 1) {
    dw.per.size.matrix <- do.call(cbind, dw.per.size)
    ## Take max deltaW across all windows for any given window boundary
    #max.deltaW <- apply(dw.per.size.matrix, 1, max) #[EXPERIMENTAL]
    ## Take mean deltaW across all windows for any given window boundary
    max.deltaW <- apply(dw.per.size.matrix, 1, mean) #[EXPERIMENTAL]
  } else {
    max.deltaW <- unlist(dw.per.size)
  }  
  
  dw$deltaW <- max.deltaW
  return(dw)
}