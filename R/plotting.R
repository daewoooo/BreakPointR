#' Plotting genome-wide ideograms \pkg{\link{breakpointR}}
#' 
#' This function will create genome-wide ideograms from a \code{\link{BreakPoint}} object.
#'
#' @param files2plot A list of files that contains \code{\link{BreakPoint}} objects or a single \code{\link{BreakPoint}} object.
#' @param file Name of the file to plot to.
#' @return A list with \code{\link[ggplot2:ggplot]{ggplot}} objects.
#'
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'exampleFile <- list.files(exampleFolder, full.names=TRUE)[1]
#'## Plot the file
#'plotBreakpoints(files2plot=exampleFile)

plotBreakpoints <- function(files2plot, file=NULL) {
    
    if (is(files2plot, class.breakpoint)) {
        numplots <- 1
    } else if (is.character(files2plot)) {
        numplots <- length(files2plot)
    } else {
        stop("Unsupported object class submitted!!!")
    }
  
    plots <- list()
    for (i in seq_len(numplots)) {
        if (is(files2plot, 'character')) {  
            data <- loadFromFiles(files2plot[i], check.class=class.breakpoint)[[1]]
        } else if (is(files2plot, class.breakpoint)) {  
            data <- files2plot
        } else {
            stop("Only 'BreakPoint' class object can be plotted")
        }
    
        filename <- data$ID
        ptm <- startTimedMessage("Plotting ", filename, " ...")
        
        bamfile <- data$ID
        reads <- data$fragments
        chroms2plot <- GenomeInfoDb::seqlevels(reads)
        breaks <- data$breaks
        counts <- data$counts
        lib.metrics <- data$lib.metrics
        lib.metrics <- round(lib.metrics, digits = 5)
        lib.metrics <- paste(names(lib.metrics), lib.metrics, sep = '=')
        lib.metrics <- paste(lib.metrics, collapse = "  |  ")
        
        #Skip chromosomes shorter then 5-times of the bin size 200kb
        if (any(seqlengths(reads) < 200000*5)) {
          message(" Skipping short chromosomes/contigs!")
          keep.chroms <- names(seqlengths(reads)[seqlengths(reads) >= 200000*5])
          reads <- GenomeInfoDb::keepSeqlevels(reads, keep.chroms, pruning.mode = 'coarse')
        }  
        
        binned.data <- unlist(GenomicRanges::tileGenome(seqlengths(reads), tilewidth = 200000))
        
        #counts overlaps between bins and our reads
        Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='-']) 
        Crickreads <- GenomicRanges::countOverlaps(binned.data, reads[strand(reads)=='+'])
        bothreads <- Watsonreads + Crickreads
        
        mcols(binned.data)$bothreads <- bothreads
        mcols(binned.data)$Watsonreads <- Watsonreads
        mcols(binned.data)$Crickreads <- Crickreads
      
        #transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
        cum.seqlengths <- cumsum(as.numeric(seqlengths(binned.data)))
        cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
        #get positions of ends of each chromosome to plot lones between the chromosomes
        if (length(cum.seqlengths) > 1) {
            chr.lines <- data.frame( y=cum.seqlengths[-length(cum.seqlengths)] )
        } else {
            chr.lines <- data.frame( y=0 )
        }    
        #get positions of each chromosomes names
        chr.label.pos <- round( cum.seqlengths.0 + (0.5 * seqlengths(binned.data) ) )
        names(chr.label.pos) <- gsub("chr", "", names(chr.label.pos)) #line to add to exclude chr
      
        #transform chromosome based coordinates into genomewide coordinates
        trans.reads <- transCoord(binned.data)
        trans.breaks <- transCoord(breaks)
        trans.counts <- transCoord(counts)
      
        dfplot.reads <- as.data.frame(trans.reads)
        dfplot.breaks <- as.data.frame(trans.breaks)
        dfplot.counts <- as.data.frame(trans.counts)
      
        my_theme <- theme(
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()
        )
      
        ### PLOT READS
      
        #get midpoint values for each genomic bin
        dfplot.reads$midpoint <- dfplot.reads$start.genome + ( (dfplot.reads$end.genome - dfplot.reads$start.genome) %/% 2 )
      
        #filter bins with extremely high amount of reads
        Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
        Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
        #set outlier bins to the limit
        dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
        dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
        
        #construct ggplot
        dfplot.reads$mCrickreads <- -dfplot.reads$Crickreads
        ggplt1 <- ggplot(dfplot.reads) +
            geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color="paleturquoise4", size=0.2) +
            geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color="sandybrown", size=0.2) +
            geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') + xlab(NULL) +
            ylab("Read counts") +
            scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) +
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            my_theme +
            ggtitle(bamfile, subtitle = lib.metrics)
        
        ### PLOT COUNTS
      
        ## Scaling size of the rectangle to amout of reads in a given region
        scale <- (dfplot.counts[,c('Ws','Cs')] / dfplot.counts$width) * 1000000
          
        ## filter regions small regions => hard to see on the plot
        outlier.W <- stats::quantile(scale$Ws, 0.9)
        outlier.C <- stats::quantile(scale$Cs, 0.9)
        
        #set.max.W <- round(max(scale$Ws[scale$Ws < outlier.W]))
        #set.max.C <- round(max(scale$Cs[scale$Cs < outlier.C]))  
        #scale$Ws[scale$Ws > outlier.W] <- set.max.W
        #scale$Cs[scale$Cs > outlier.C] <- set.max.C
        
        scale$Ws[scale$Ws >= outlier.W] <- outlier.W
        scale$Cs[scale$Cs >= outlier.C] <- outlier.C
            
        ## putting scaled and filtered regions into a data frame
        names(scale) <- c('W.scaled', 'C.scaled')
        df.W <- cbind(dfplot.counts, scaled=scale$W.scaled, fill.strand=rep('W', length(scale$W.scaled)) )
        df.C <- cbind(dfplot.counts, scaled=-scale$C.scaled, fill.strand=rep('C', length(scale$W.scaled)) )
        dfplot.counts <- rbind(df.W, df.C)
        
        dfplot.counts <- dfplot.counts[dfplot.counts$width > 20000,]
        
        #get midpoint values for each breakpoint
        dfplot.breaks$midpoint <- dfplot.breaks$start.genome + ( (dfplot.breaks$end.genome - dfplot.breaks$start.genome) %/% 2)
        
        ggplt2 <- ggplot(dfplot.counts) +
            geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax='scaled', fill='fill.strand')) +
            geom_linerange(data=chr.lines, aes_string(ymin=-Inf, ymax=Inf, x='y'), col='black') +
            scale_fill_manual(values=c("sandybrown","paleturquoise4"))
        if (nrow(dfplot.breaks) > 0) {
            ggplt2 <- ggplt2 + 
                geom_point(data=dfplot.breaks, aes_string(x='midpoint', y=0), size=5, color='red', shape=124, inherit.aes = FALSE) +
                xlab(NULL) +
                ylab("Breaks") +
                scale_x_continuous(expand = c(0,0)) +
                theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
                my_theme  
        } else {
            ggplt2 <- ggplt2 +
                xlab(NULL) +
                ylab("Breaks") +
                scale_x_continuous(expand = c(0,0)) +
                theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
                my_theme
        }
            
        ### PLOT STATES
        
        ggplt3 <- ggplot(dfplot.counts) +
            geom_rect(aes_string(xmin='start.genome', xmax='end.genome', ymin=0, ymax=10, fill='states')) +
            geom_linerange(data=chr.lines, aes_string(ymin=0, ymax=10, x='y'), col='black') +
            xlab("Chromosomes") +
            ylab("States") +
            scale_x_continuous(breaks=chr.label.pos, labels=names(chr.label.pos), expand = c(0,0)) +
            scale_fill_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) +
            theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
            my_theme
        
        ### PLOT ALL
      
        if (length(breaks)) {
            p <- suppressWarnings( cowplot::plot_grid(ggplt1, ggplt2, ggplt3, ncol=1, align="v", rel_heights = c(3,3,2)) )
        } else {
            p <- suppressWarnings( cowplot::plot_grid(ggplt1, ncol=1, align="v", rel_heights = 3) )
        }
        plots[[length(plots)+1]] <- p
        stopTimedMessage(ptm)
    }
  
    ### PRINT TO FILE
  
    ## printing to a file or returning plot object
    if (!is.null(file)) {
        message("Printing to PDF ",file)
    
        grDevices::pdf(file, width=max(10, length(chroms2plot)), height=5)
        bquiet = lapply(plots, print)
        d <- grDevices::dev.off()
    }
    return(plots)
}


#' Genome wide heatmap of template inheritance states
#'
#' Plot a genome-wide heatmap of template inheritance states from a \code{\link{BreakPoint}} object.
#'
#' @param files2plot A list of files that contains \code{\link{BreakPoint}} objects or a single \code{\link{BreakPoint}} object.
#' @param file Name of the file to plot to.
#' @param hotspots A \code{\link{GRanges-class}} object with locations of breakpoint hotspots.
#' @return A \code{\link[ggplot2:ggplot]{ggplot}} object.
#' 
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @importFrom grDevices pdf dev.off
#' @importFrom S4Vectors endoapply
#' @export
#' @examples 
#'## Get example BreakPoint objects to plot
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'exampleFiles <- list.files(exampleFolder, full.names=TRUE)
#'breakpoint.objects <- loadFromFiles(exampleFiles)
#'## Plot the heatmap
#'plotHeatmap(breakpoint.objects)
#'
plotHeatmap <- function(files2plot, file=NULL, hotspots=NULL) {

    if (is(files2plot, class.breakpoint)) {
        stop("Cannot make heatmap from only one object.")
    } else {
        numlibs2plot <- length(files2plot)
    } 
  
    ptm <- startTimedMessage("Preparing heatmap from ",numlibs2plot, " libraries ...")

    IDs <- list()
    grl <- GRangesList()
    breaks <- GRangesList()
    for (i in seq_len(numlibs2plot)) {  
        data <- loadFromFiles(files2plot[i], check.class=class.breakpoint)
        IDs[[i]] <- data[[1]]$ID
        grl[[i]] <- data[[1]]$counts
        suppressWarnings( breaks[[i]] <- data[[1]]$breaks )
    }
    
    #transform genomic ranges to genome-wide coordinates
    grl <- endoapply(grl, transCoord)
    
    GenomeInfoDb::seqlevels(breaks) <- GenomeInfoDb::seqlevels(grl)
    GenomeInfoDb::seqlengths(breaks) <- GenomeInfoDb::seqlengths(grl)

    # disjoin overlaping breaks
    breaks <- unlist(breaks, use.names=FALSE)
    ranges <- GenomicRanges::disjoin(breaks) # redefine ranges in df
    hits <- GenomicRanges::countOverlaps(ranges, breaks) # counts number of breaks overlapping at each range
    mcols(ranges)$hits <- hits

    disjoin.breaks <- transCoord(ranges)
    dfplot.disjoin.breaks <- as.data.frame(disjoin.breaks)

    #transform bin coordinates of each chromosome into genomewide coordinates (cumulative sum of bin coordintes)
    cum.seqlengths <- cumsum(as.numeric(seqlengths(grl[[1]])))
    cum.seqlengths.0 <- c(0,cum.seqlengths[-length(cum.seqlengths)])
    # Chromosome lines for heatmap
    label.pos <- round( cum.seqlengths.0 + 0.5 * seqlengths(grl[[1]]) )
    names(label.pos) <- gsub("chr", "", names(label.pos)) #line to add to exclude chr
    df.chroms <- data.frame(y=c(0,cum.seqlengths))
    
    my_theme <- theme(
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()  
    )

    # Plot breaks summary
    #ggplt1 <- ggplot(dfplot.disjoin.breaks)
    #ggplt1 <- ggplt1 + geom_rect(aes_string(ymin=0, ymax='hits', xmin='start.genome', xmax='end.genome'), fill="red", color="red")
    #ggplt1 <- ggplt1 + geom_vline(aes_string(xintercept='y'), data=df.chroms, col='black')
    #ggplt1 <- ggplt1 + 
    #    scale_y_continuous(expand = c(0,0)) +
    #    ylab("Break\ncounts") +
    #    my_theme
    
    # Data
    df <- list()
    for (i1 in seq_along(grl)) {
        df[[length(df)+1]] <- data.frame(start=grl[[i1]]$start.genome, end=grl[[i1]]$end.genome, seqnames=seqnames(grl[[i1]]), sample=IDs[[i1]], state=grl[[i1]]$states)
    }
    df <- do.call(rbind, df)

    ## PLOT
    if (is.null(hotspots)) {
        ggplt2 <- ggplot(df) +
            geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) +
            scale_y_continuous(breaks=label.pos, labels=names(label.pos), expand=c(0,0)) +
            coord_flip() +
            scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) +
            theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20)) +
            geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black')
    } else {
        hotspots.trans <- transCoord(hotspots)
        hotspot.midpoint <- hotspots.trans$start.genome + (hotspots.trans$end.genome - hotspots.trans$start.genome)/2
        hotspot.midpoint <- data.frame(y=hotspot.midpoint)
        ggplt2 <- ggplot(df) +
            geom_linerange(aes_string(ymin='start', ymax='end', x='sample', col='state'), size=5) +
            scale_y_continuous(breaks=label.pos, labels=names(label.pos)) +
            coord_flip() +
            scale_color_manual(values=c('cc'="paleturquoise4", 'wc'="olivedrab",'ww'="sandybrown",'?'="red")) +
            theme(panel.background=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=20)) +
            geom_hline(aes_string(yintercept='y'), data=df.chroms, col='black') +
            geom_hline(aes_string(yintercept='y'), data=hotspot.midpoint, col='red', alpha=0.5)
    }  
    stopTimedMessage(ptm)
    
    ## PRINT TO FILE
    ## printing to a file or returning plot object
    if (!is.null(file)) {
        message("Printing to PDF ",file)
        height.cm <- length(files2plot) * 0.25
        width.cm <- max(10, length(levels(df$seqnames))*5 )
        grDevices::pdf(file, width=width.cm, height=height.cm)
        print(ggplt2)
        d <- grDevices::dev.off()
        #p <- arrangeGrob(p)
        #ggsave(file=file, p, width=width.cm/5, height=height.cm/5, limitsize = F, device = "pdf")
    } 
    return(ggplt2)
}



#' Plotting chromosome specific ideograms \pkg{\link{breakpointR}}
#' 
#' This function will create chromsome specific enome-wide ideograms from a \code{\link{BreakPoint}} object.
#'
#' @param files2plot A list of files that contains \code{\link{BreakPoint}} objects or a single \code{\link{BreakPoint}} object.
#' @param plotspath Directory to store plots.
#' @param chromosomes Set specific chromosome(s) to be plotted.
#' @return A list with \code{\link[ggplot2:ggplot]{ggplot}} objects.
#'
#' @author David Porubsky
#' @import ggplot2
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
#'## Get an example file
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'exampleFiles <- list.files(exampleFolder, full.names=TRUE)
#'## Plot results
#'plotBreakpointsPerChr(exampleFiles, chromosomes='chr7')
#'
plotBreakpointsPerChr <- function(files2plot, plotspath=NULL, chromosomes=NULL) {

    if (is(files2plot, class.breakpoint)) {  
        numplots <- 1
        chroms.in.data <- GenomeInfoDb::seqlevels(files2plot$fragments)
        chrom.lengths <- GenomeInfoDb::seqlengths(files2plot$fragments)
        #Skip chromosomes shorter then 5-times of the bin size 200kb
        if (any(chrom.lengths < 200000*5)) {
            message("Will skip short chromosomes/contigs!")
            chroms.in.data <- names(chrom.lengths[chrom.lengths >= 200000*5])
        }
        #files2plot <- list(files2plot)
    } else if (is.character(files2plot)) {
        numplots <- length(files2plot)
        #get sequence levels and sequence lengths from a single BreakpointR obejct
        first.file <- loadFromFiles(files2plot[[1]], check.class=class.breakpoint)[[1]]
        chroms.in.data <- GenomeInfoDb::seqlevels(first.file$fragments)
        chrom.lengths <- GenomeInfoDb::seqlengths(first.file$fragments)
        #Skip chromosomes shorter then 5-times of the bin size 200kb
        if (any(chrom.lengths < 200000*5)) {
            message("Will skip short chromosomes/contigs!")
            chroms.in.data <- names(chrom.lengths[chrom.lengths >= 200000*5])
        }
    } else {
      stop("Unsupported object class submitted!!!")
    }
  
    if (is.null(chromosomes)) {
        chromosomes <- chroms.in.data
    }  
    chroms2plot <- intersect(chromosomes, chroms.in.data)
    if (length(chroms2plot)==0) {
        chrstring <- paste0(chromosomes, collapse=', ')
        stop('The specified chromosomes ', chrstring, ' do not exist in the data. Please try ', paste(paste0('chr',chromosomes), collapse=', '), ' instead.')
    }
  
    plots <- list()
    for (chr in chroms2plot) {
        ptm <- startTimedMessage("Plotting breakpoints for chromosome: ", chr, " ...")
    
        dfplot.reads.chr <- list()
        dfplot.breaks.chr <- list()
        for (i in seq_len(numplots)) {
        
            if (is(files2plot, 'character')) {  
                data <- loadFromFiles(files2plot[i], check.class=class.breakpoint)[[1]]
            } else if (is(files2plot, class.breakpoint)) {  
                data <- files2plot
            } else {
                stop("Only 'BreakPoint' class object can be plotted")
            }
              
            bamfile <- data$ID
            reads <- data$fragments
            breaks <- data$breaks
        
            #select reads and breaks for a single chromosome
            reads.chr <- GenomeInfoDb::keepSeqlevels(reads, chr, pruning.mode="coarse")
            breaks.chr <- GenomeInfoDb::keepSeqlevels(breaks, chr, pruning.mode="coarse")
        
            binned.data <- unlist(GenomicRanges::tileGenome(seqlengths(reads.chr), tilewidth = 200000))
        
            #add file IDs
            binned.data$ID <- bamfile
            if (length(breaks.chr) > 0) {
                breaks.chr$ID <- bamfile
                dfplot.breaks <- as.data.frame(breaks.chr)
                dfplot.breaks.chr[[i]] <- dfplot.breaks
            }
        
            #counts overlaps between bins and our reads
            Watsonreads <- GenomicRanges::countOverlaps(binned.data, reads.chr[strand(reads.chr)=='-']) 
            Crickreads <- GenomicRanges::countOverlaps(binned.data, reads.chr[strand(reads.chr)=='+'])
            bothreads <- Watsonreads + Crickreads
    
            mcols(binned.data)$bothreads <- bothreads
            mcols(binned.data)$Watsonreads <- Watsonreads
            mcols(binned.data)$Crickreads <- Crickreads
        
            dfplot.reads <- as.data.frame(binned.data)
        
            #filter bins with extremely high amount of reads
            Crickreads.outlier <- stats::quantile(dfplot.reads$Crickreads, 0.999)
            Watsonreads.outlier <- stats::quantile(dfplot.reads$Watsonreads, 0.999)
            #set outlier bins to the limit
            dfplot.reads$Crickreads[dfplot.reads$Crickreads >= Crickreads.outlier] <- Crickreads.outlier
            dfplot.reads$Watsonreads[dfplot.reads$Watsonreads >= Watsonreads.outlier] <- Watsonreads.outlier
        
            dfplot.reads.chr[[i]] <- dfplot.reads
        }
    
        dfplot.reads.chr <- do.call(rbind, dfplot.reads.chr)
        dfplot.breaks.chr <- do.call(rbind, dfplot.breaks.chr)
        
        my_theme <- theme(
        legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank(),
            strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm"))
        )
    
        ### PLOT READS
    
        #get midpoint values for each genomic bin
        dfplot.reads.chr$midpoint <- dfplot.reads.chr$start + ( (dfplot.reads.chr$end - dfplot.reads.chr$start) %/% 2 )
    
        #construct ggplot
        dfplot.reads.chr$mCrickreads <- -dfplot.reads.chr$Crickreads
        ggplt <- ggplot(dfplot.reads.chr) +
            geom_linerange(aes_string(ymin=0, ymax='mCrickreads', x='midpoint'), color="paleturquoise4", size=1) +
            facet_wrap(ID ~ seqnames, scales = "free", ncol=1) +
            geom_linerange(aes_string(ymin=0, ymax='Watsonreads', x='midpoint'), color="sandybrown", size=1) +
            xlab("Genomic position") +
            ylab("Read counts") +
            scale_x_continuous(expand = c(0,0))
        if (!is.null(dfplot.breaks.chr)) {
            #get midpoint values for each genomic bin
            dfplot.breaks.chr$midpoint <- dfplot.breaks.chr$start + ( (dfplot.breaks.chr$end - dfplot.breaks.chr$start) %/% 2)
            ggplt <- ggplt +
                geom_linerange(data=dfplot.breaks.chr, aes_string(ymin=-Inf, ymax=Inf, x='midpoint'), col='black') +
                facet_wrap(ID ~ seqnames, scales = "free", ncol=1) + my_theme
        } else {
            ggplt <- ggplt +
                facet_wrap(ID ~ seqnames, scales = "free", ncol=1) + my_theme
        }
        
        if (!is.null(plotspath)) {
            file.destination <- file.path(plotspath, paste0(chr, "_breakpoints.pdf"))
            grDevices::pdf(file.destination, width=10, height=numplots*2)
            print(ggplt)
            d <- grDevices::dev.off()   
        }
    
        plots[[chr]] <- ggplt
        stopTimedMessage(ptm)
    }
    return(plots)
}

