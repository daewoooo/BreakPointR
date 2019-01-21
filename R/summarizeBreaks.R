#' Compile breakpoint summary table
#'
#' This function will calculate deltaWs from a \code{\link{GRanges-class}} object with read fragments.
#'
#' @param breakpoints A list containing breakpoints stored in \code{\link{GRanges-class}} object.
#' @return A \code{data.frame} of compiled breakpoints together with confidence intervals.
#' @author David Porubsky
#' @export
#' @examples 
#'## Get some files that you want to load
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'file <- list.files(exampleFolder, full.names=TRUE)[1]
#'breakpoints <- get(load(file))[c('breaks', 'confint')]
#'summarizeBreaks(breakpoints)

summarizeBreaks <- function(breakpoints) {
    breaks <- breakpoints$breaks
    confint <- breakpoints$confint
    if (length(breaks) > 0 & length(confint) == length(breaks)) {
        breaks.df <- as(breaks, 'data.frame')
        confint.df <- as(confint, 'data.frame')
        breaks.df <- breaks.df[,c('seqnames','start','end')]
        confint.df <- confint.df[,c('start','end','genoT')]
        breaksSummary <- cbind(breaks.df, confint.df)
        names(breaksSummary) <- c('seqnames','start','end','CI.start','CI.end','genoT')
        return(breaksSummary)
    } else if (length(breaks) > 0 & length(confint) == 0) {
        breaks.df <- as(breaks, 'data.frame')
        breaksSummary <- breaks.df[,c('seqnames','start','end','genoT')]
        return(breaksSummary)
    } else {
        return(NULL)
    }
}
