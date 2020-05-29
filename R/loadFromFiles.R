#' Load \pkg{breakpointR} objects from file
#'
#' Wrapper to load \pkg{\link{breakpointR}} objects from file and check the class of the loaded objects.
#'
#' @param files A list of \code{\link{GRanges-class}} or \code{\link{BreakPoint}} objects or a vector of files that contain such objects.
#' @param check.class Any combination of \code{c('GRanges', 'BreakPoint')}. If any of the loaded objects does not belong to the specified class, an error is thrown.
#' @return A list of \code{\link{GRanges-class}} or \code{\link{BreakPoint}} objects.
#' @export
#' @examples
#'## Get some files that you want to load
#'exampleFolder <- system.file("extdata", "example_results", package="breakpointRdata")
#'exampleFiles <- list.files(exampleFolder, full.names=TRUE)
#'## Load the processed data                          
#'breakpoint.objects <- loadFromFiles(exampleFiles)
 
loadFromFiles <- function(files, check.class=c('GRanges', 'BreakPoint')) {

    # ptm <- startTimedMessage("Loading data from files ...")
    if (is.null(files)) {
        # stopTimedMessage(ptm)
        return(files)
    }
    if (any(! check.class %in% c('GRanges', class.breakpoint))) {
        stop("Argument 'check.class' must contain any combination of c('", paste0(c('GRanges', class.breakpoint), collapse="', '"), "').")
    }
    modellist <- list()
    if (is.character(files)) {
        for (file in files) {
            temp.env <- new.env()
            model <- get(load(file, envir=temp.env), envir=temp.env)
            if (! class(model) %in% check.class) {
                stop("File '", file, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[file]] <- model
        }
    } else if (class(files) %in% check.class) {
        modellist[[1]] <- files
    } else if (is.list(files)) {
        for (file in files) {
            model <- file
            if (! class(model) %in% check.class) {
                stop("List entry '", length(modellist)+1, "' does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
            }
            modellist[[length(modellist)+1]] <- model
        }
        names(modellist) <- names(files)
    } else if (! class(files) %in% check.class) {
        stop("Input does not contain an object of class ", paste0(check.class, collapse=' or '), ".")
    }
    # stopTimedMessage(ptm)
    return(modellist)
}
