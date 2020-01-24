#' Create composite Strand-seq file
#'
#' This function will move through BAM files in a folder, read in each individual file and go through each chromosome, 
#' determine if the chromosome is WW or CC based on WCcutoff, reverse complement all reads in the WW file,
#' append to a new composite file for that chromosome, order the composite file of each chromosome based on position.
#'
#' @param file.list A list of BAM files to process.
#' @inheritParams readBamFileAsGRanges
#' @inheritParams GenotypeBreaks
#' @param WC.cutoff Percentage of WW or CC reads to consider chromosome being WW or CC
#' @return A \code{\link{GRanges-class}} object.
#' @importFrom gtools mixedsort
#' @importFrom BiocGenerics table
#' @import GenomicRanges
#' @author Ashley Sanders, David Porubsky
#'                    
createCompositeFile <- function(file.list, chromosomes=NULL, pairedEndReads=TRUE, pair2frgm=FALSE, min.mapq=10, filtAlt=FALSE, WC.cutoff=0.90, genoT='fisher', background=0.05) {

    message("Creating composite file from ", length(file.list), " bam files")

    composite.bam.grl <- GenomicRanges::GRangesList()
    for (bamfile in file.list) {
        message("Working on file ",bamfile)
        fragments <- suppressWarnings( readBamFileAsGRanges(bamfile, pairedEndReads=pairedEndReads, chromosomes=chromosomes, min.mapq=min.mapq, pair2frgm=pair2frgm, filtAlt=filtAlt) )
        if (length(fragments) > 0) { mcols(fragments)$lib <- bamfile } # appends file name to reads
        composite.bam <- GenomicRanges::GRangesList()
        for (chr in unique(seqnames(fragments))) {
            #message("Working on chromosome ",chr)
            chr.fragments <- fragments[seqnames(fragments)==chr]
            tab.counts <- BiocGenerics::table(strand(chr.fragments))
            plus.count <- tab.counts[1]
            minus.count <- tab.counts[2]

            # ## Genotyping based on W-C ratio
            # wcCall <- round( ( minus.count - plus.count ) / length(chr.fragments), digits=3 )
            # 
            # ## check if this is a pure (WW or CC) library
            # if( wcCall <= -WC.cutoff | wcCall >= WC.cutoff ) {
            #     ## if this is pure WW, reverse compliment all the reads
            #     if(wcCall <= -WC.cutoff) {
            #         chr.fragments.copy <- chr.fragments
            #         strand(chr.fragments)[strand(chr.fragments.copy) == '+'] <- '-'
            #         strand(chr.fragments)[strand(chr.fragments.copy) == '-'] <- '+'
            #     }
            #     composite.bam[[chr]] <- chr.fragments
            # }
            
            if (genoT == 'fisher') {
                ## Use Fisher Exact Test to genotype
                geno <- genotype.fisher(cReads=plus.count, wReads=minus.count, roiReads=plus.count+minus.count, background=background, minReads = 10)
            } else if (genoT == 'binom') {    
                ## Use binomial test to genotype
                geno <- genotype.binom(cReads=plus.count, wReads=minus.count, background=background, log = TRUE)
            } else {
                stop("Wrong argument 'genoT', genoT='fisher|binom'")
            } 
            
            ## check if this is a pure (WW or CC) library
            if (geno$bestFit == 'ww' | geno$bestFit == 'cc') {
                ## if this is pure WW, reverse compliment all the reads
                if (geno$bestFit == 'ww') {
                    chr.fragments.copy <- chr.fragments
                    strand(chr.fragments)[strand(chr.fragments.copy) == '+'] <- '-'
                    strand(chr.fragments)[strand(chr.fragments.copy) == '-'] <- '+'
                }
                composite.bam[[chr]] <- chr.fragments
            }
            
        }  
        composite.bam.grl[[bamfile]] <- unlist(composite.bam)
    }
    composite.data <- unlist(composite.bam.grl)
    names(composite.data) <- NULL
    seqlevels(composite.data) <- gtools::mixedsort( as.character(unique(seqnames(composite.data))) )
    composite.data <- GenomicRanges::sort(composite.data, ignore.strand=TRUE)
    return(composite.data)
}

