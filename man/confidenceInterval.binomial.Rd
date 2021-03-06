% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidenceInterval.R
\name{confidenceInterval.binomial}
\alias{confidenceInterval.binomial}
\title{Estimate confidence intervals for breakpoints}
\usage{
confidenceInterval.binomial(breaks, fragments, background = 0.02, conf = 0.99)
}
\arguments{
\item{breaks}{Genotyped breakpoints as outputted from function \code{\link{GenotypeBreaks}}.}

\item{fragments}{Read fragments from function \code{\link{readBamFileAsGRanges}}.}

\item{background}{The percent (e.g. 0.05 = 5\%) of background reads allowed for WW or CC genotype calls.}

\item{conf}{Desired confidence interval of localized breakpoints.}
}
\value{
A \code{\link{GRanges-class}} object of breakpoint ranges for a given confidence interval in \code{conf}.
}
\description{
Estimate confidence intervals for breakpoints by going outwards from the breakpoint read by read, and performing a binomial test of getting the observed or a more extreme outcome, given that the reads within the confidence interval belong to the other side of the breakpoint.
}
\examples{
\dontrun{ 
## Get an example file 
exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
## Load the file 
breakpoint.objects <- get(load(exampleFile))
## Calculate confidence intervals of genotyped breakpoints
confint <- confidenceInterval.binomial(breakpoint.objects$breaks, breakpoint.objects$fragments, background=0.02)}

}
\author{
Aaron Taudt, David Porubsky
}
