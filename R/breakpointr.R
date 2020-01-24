#' Main function for the \pkg{\link{breakpointR}} package
#'
#' This function is an easy-to-use wrapper to find breakpoints with \code{\link[breakpointR]{runBreakpointr}} in parallel, write the results to file, plot results and find hotspots.
#'
#' @param inputfolder Folder with BAM files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing files in \code{outputfolder} should be reused.
#' @param callHotSpots Search for regions of high abundance of breakpoints in single cells.

#' @inheritParams runBreakpointr
#' @inheritParams readBamFileAsGRanges
#' @inheritParams createCompositeFile
#' @inheritParams confidenceInterval

#' @return \code{NULL}
#' @author David Porubsky, Aaron Taudt, Ashley Sanders
#' @importFrom utils read.table
#' @import foreach
#' @import doParallel
#' @export
#' 
#' @examples
#'\dontrun{
#'## The following call produces plots and genome browser files for all BAM files in "my-data-folder"
#'breakpointr(inputfolder="my-data-folder", outputfolder="my-output-folder")}
#'
breakpointr <- function(inputfolder, outputfolder, configfile=NULL, numCPU=1, reuse.existing.files=FALSE, windowsize=1e6, binMethod="size", pairedEndReads=FALSE, pair2frgm=FALSE, chromosomes=NULL, min.mapq=10, filtAlt=FALSE, genoT='fisher', trim=10, peakTh=0.33, zlim=3.291, background=0.05, minReads=10, maskRegions=NULL, callHotSpots=FALSE, conf=0.99) {
 
#=======================
### Helper function ###
#=======================
  
## This parallel helper function performs breakpoint detection and export results into the pre-set locations.  
runBreakpointrANDexport <- function(file, datapath, browserpath, config) {
    savename <- file.path(datapath, paste0(basename(file),'.RData'))
    ## Find breakpoints
    if (!file.exists(savename)) {
        tC <- tryCatch({
            breakpoints <- runBreakpointr(bamfile=file, ID=basename(file), pairedEndReads=config[['pairedEndReads']], pair2frgm=config[['pair2frgm']], min.mapq=config[['min.mapq']], filtAlt=config[['filtAlt']], chromosomes=config[['chromosomes']], windowsize=config[['windowsize']], binMethod=config[['binMethod']], trim=config[['trim']], peakTh=config[['peakTh']], zlim=config[['zlim']], background=config[['background']], genoT=config[['genoT']], minReads=config[['minReads']], maskRegions=maskRegions, conf=config[['conf']])
        }, error = function(err) {
            stop(file,'\n',err)
        })  
        save(breakpoints, file=savename)
    } else {
        breakpoints <- get(load(savename))
    }
    ## Write BED file
    savename.breakpoints <- file.path(browserpath,paste0(basename(file), '_breakPoints.bed.gz'))
    if (!file.exists(savename.breakpoints)) {
        breakpointr2UCSC(index=basename(file), outputDirectory=browserpath, fragments=breakpoints$fragments, deltaWs=breakpoints$deltas, breakTrack=breakpoints$breaks, confidenceIntervals=breakpoints$confint)
    }
}  
  
#========================
### General variables ###
#========================
config <- NULL
if (is.character(configfile)) {
    ## Read config file ##
    errstring <- tryCatch({
        config <- readConfig(configfile)
        errstring <- ''
    }, error = function(err) {
        errstring <- paste0("Could not read configuration file ",configfile)
    })
    if (errstring!='') {
        stop(errstring)
    }
}

## Set createCompositeFile to FALSE [experimental]
createCompositeFile=FALSE

## Put options into list and merge with conf ##
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, 
               windowsize=windowsize, binMethod=binMethod, pairedEndReads=pairedEndReads, 
               pair2frgm=pair2frgm, chromosomes=chromosomes, min.mapq=min.mapq, 
               filtAlt=filtAlt, trim=trim, peakTh=peakTh, zlim=zlim, background=background, 
               genoT=genoT, minReads=minReads, createCompositeFile=createCompositeFile, 
               maskRegions=maskRegions, callHotSpots=callHotSpots, conf=conf)
config <- c(config, params[setdiff(names(params),names(config))])

## Checks user input ##
if (!config[['pairedEndReads']] & config[['pair2frgm']]) {
    message("Option pair2frgm=TRUE can be used only when pairedEndReads=TRUE\nSetting pair2frgm to FALSE!!!")
    config[['pair2frgm']] <- FALSE
}

if (config[['createCompositeFile']] & config[['callHotSpots']]) {
    message("If createCompositeFile=TRUE breakpoint hotspots can't be called\nSetting callHotSpots to FALSE!!!")
    config[['callHotSpots']] <- FALSE
}
    
## Helpers ##
windowsize <- config[['windowsize']]

## Set up the directory structure ##
datapath <- file.path(outputfolder,'data')
browserpath <- file.path(outputfolder,'browserfiles')
plotspath <- file.path(outputfolder,'plots')
breakspath <- file.path(outputfolder,'breakpoints')

## Delete old directory if desired ##
if (config[['reuse.existing.files']]==FALSE) {
    if (file.exists(outputfolder)) {
        message("Deleting old directory ",outputfolder)
        unlink(outputfolder, recursive=TRUE)
    }
}

## Creating required directories ##
if (!file.exists(outputfolder)) {
    dir.create(outputfolder)
}
if (!file.exists(datapath)) { dir.create(datapath) }
if (!file.exists(browserpath)) { dir.create(browserpath) }
if (!file.exists(plotspath)) { dir.create(plotspath) }
if (!file.exists(breakspath)) { dir.create(breakspath) }

## Write README file RR
savename <- file.path(outputfolder, 'README.txt')
cat("", file=savename)
cat("BreakpointR outputfolder contains the following folders.\n", file=savename, append=TRUE)
cat("========================================================\n", file=savename, append=TRUE)
cat("- breakpoints: UCSC browser compatible bedgraphs compiling all breakpoints across all single-cell libraries.\n", file=savename, append=TRUE)
cat("               List of all localized breakpoints in all single-cell libraries.\n", file=savename, append=TRUE)
cat("               Locations of breakpoint hotspots if 'callHotSpots=TRUE'.\n", file=savename, append=TRUE)
cat("- browserfiles: UCSC browser formated files with exported reads, deltaWs and breakPoints for every single-cell library.\n", file=savename, append=TRUE)
cat("                UCSC browser formated list of masked genomic regions if set option 'maskRegions'.\n", file=savename, append=TRUE)
cat("- data: RData files storing results of BreakpointR analysis for each single-cell library in an object of class 'BreakPoint'.\n", file=savename, append=TRUE)
cat("- plots: Genome-wide plots for selected chromsosome, genome-wide heatmap of strand states as well as chromosome specific read distribution together with localized breakpoints.\n", file=savename, append=TRUE)

## Make a copy of the config file ##
writeConfig(config = config, configfile=file.path(outputfolder, 'breakpointR.config'))

## Load user defined mask regions ##
if (!is.null(maskRegions)) {
    mask.df <- utils::read.table(maskRegions, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    if (is.null(chromosomes)) {
        maskRegions <- GenomicRanges::GRanges(seqnames=mask.df$V1, IRanges(start=mask.df$V2, end=mask.df$V3))   
    } else {
        if (any(mask.df$V1 %in% chromosomes)) {
            mask.df <- mask.df[mask.df$V1 %in% chromosomes,] #select only chromosomes submitted for analysis
            maskRegions <- GenomicRanges::GRanges(seqnames=mask.df$V1, IRanges(start=mask.df$V2, end=mask.df$V3))      
        } else {
            warning(paste0('Skipping maskRegions. Chromosomes names conflict: ', mask.df$V1[1], ' not equal ',chromosomes[1]))
            maskRegions <- NULL
        } 
    }  
}

#=====================================
# Find breakpoints and write BED files
#=====================================
files <- list.files(inputfolder, full.names=TRUE, recursive=FALSE, pattern='.bam$')
if (length(files) == 0) {
    stop("No BAM files present in an 'inputfolder'.")  
}

### Run breakpointR ###
if (createCompositeFile) {
    config[['numCPU']] <- 1 #always use only one CPU for composite file analysis
    fragments <- createCompositeFile(file.list=files, chromosomes=config[['chromosomes']], pairedEndReads=config[['pairedEndReads']], pair2frgm=config[['pair2frgm']], min.mapq=config[['min.mapq']], filtAlt=config[['filtAlt']], genoT=config[['genoT']], background=config[['background']])    
    runBreakpointrANDexport(file = fragments, datapath = datapath, browserpath = browserpath, config = config)
} else if (numCPU > 1) {
    ## Parallelization ##
    message("Using ",config[['numCPU']]," CPUs")
    cl <- parallel::makeCluster(config[['numCPU']])
    doParallel::registerDoParallel(cl)
  
    message("Finding breakpoints ...", appendLF=FALSE); ptm <- proc.time()
    temp <- foreach (file = files, .packages=c('breakpointR')) %dopar% {
        runBreakpointrANDexport(file = file, datapath = datapath, browserpath = browserpath, config = config)
    }
    
    parallel::stopCluster(cl)
    time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
} else {
    config[['numCPU']] <- 1 #if to use only one CPU or CPU argument not defined
    for (file in files) {
        runBreakpointrANDexport(file = file, datapath = datapath, browserpath = browserpath, config = config)
    }
}

## Write masked regions to BED file
if (!is.null(maskRegions)) {
    ranges2UCSC(gr=maskRegions, index='MaskedRegions', outputDirectory=browserpath, colorRGB='0,0,0')
}

## Compile all breaks using disjoin function
if (createCompositeFile==FALSE) {
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)
  
    #breaks.all.files <- GenomicRanges::GRangesList()
    #breaksConfInt.all.files <- GenomicRanges::GRangesList()
    breaks.all.files <- list()
    breaksConfInt.all.files <- list()
    summaryBreaks <- list()
    for (file in files) {
        data <- get(load(file))[c('breaks', 'confint')]
        summaryBreaks[[basename(file)]] <- summarizeBreaks(data)
        breakpoints <- data$breaks
        breaks.confint <- data$confint
        if (length(breakpoints)) {
            suppressWarnings( breaks.all.files[[file]] <- breakpoints ) #TODO check if this can be done without warnings
        }  
        if (length(breaks.confint)) {
            suppressWarnings( breaksConfInt.all.files[[file]] <- breaks.confint ) 
        }  
    }
    
    names(breaks.all.files) <- NULL
    breaks <- do.call(c, breaks.all.files)
    ranges.br <- GenomicRanges::disjoin(breaks) # redefine ranges in df
    hits <- GenomicRanges::countOverlaps(ranges.br, breaks) # counts number of breaks overlapping at each range
    mcols(ranges.br)$hits <- hits # appends hits as a metacolumn in ranges
    
    names(breaksConfInt.all.files) <- NULL
    breaks.CI <- do.call(c, breaksConfInt.all.files)
    ranges.CI <- GenomicRanges::disjoin(breaks.CI) # redefine ranges in df
    hits <- GenomicRanges::countOverlaps(ranges.CI, breaks) # counts number of breaks overlapping at each range
    mcols(ranges.CI)$hits <- hits # appends hits as a metacolumn in ranges
  
    ## write a bedgraph of data (overlapping breaks)
    breakpointr2UCSC(index='BreakpointSummary', outputDirectory=breakspath, breaksGraph=ranges.br)
    breakpointr2UCSC(index='BreakpointConfIntSummary', outputDirectory=breakspath, breaksGraph=ranges.CI)
} else {
    files <- list.files(datapath, pattern=".RData$", full.names=TRUE)
  
    summaryBreaks <- list()
    for (file in files) {
        data <- get(load(file))[c('breaks', 'confint')]
        summaryBreaks[[basename(file)]] <- summarizeBreaks(breakpoints=data)
    }
}

## Write all breakpoints to a file
summaryBreaks.df <- do.call(rbind, summaryBreaks)
summaryBreaks.df$filenames <- rownames(summaryBreaks.df)
write.table(summaryBreaks.df, file = file.path(breakspath, 'breakPointSummary.txt'), quote = FALSE, row.names = FALSE)

## Synchronize read directionality [EXPERIMENTAL] 
#files2sync <- list.files(datapath, pattern = ".RData", full.names = TRUE)
#syncReads <- synchronizeReadDir(files2sync)
#breakpointr2UCSC(index="syncReads", outputDirectory=browserpath, fragments=syncReads)

## Search for SCE hotspots ##
if (callHotSpots) {
    hotspots <- hotspotter(gr.list=breaks.all.files, bw = 1000000, pval = 1e-10)
    if (length(hotspots)) {
        ranges2UCSC(gr=hotspots, index='HotSpots', outputDirectory=breakspath, colorRGB='0,0,255')
    }
}

## Plotting
if (createCompositeFile==FALSE) {
    files2plot <- list.files(datapath, pattern = ".RData", full.names = TRUE)
    plotBreakpoints(files2plot=files2plot, file=file.path(plotspath, 'breaksPlot.pdf')) -> beQuiet
    plotBreakpointsPerChr(files2plot=files2plot, plotspath=plotspath) -> beQuiet
    if (callHotSpots) {
        plotHeatmap(files2plot=files2plot, file=file.path(plotspath, 'breaksHeatmapPlot.pdf'), hotspots=hotspots) -> beQuiet
    } else {
        plotHeatmap(files2plot=files2plot, file=file.path(plotspath, 'breaksHeatmapPlot.pdf')) -> beQuiet
    }  
}  

}
