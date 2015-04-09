# Unexported helper functions copied from AnnotationHub 1.6.0

cache.path <- function() AnnotationHub::hubCache()

cache.exists <- function(filePath) 
  !is.na(filePath) && isTRUE(file.info(filePath)$isdir)

cache.create <- function(dirPath) {
    if (!file.exists(dirPath))
        if (!dir.create(dirPath, showWarnings=FALSE, recursive=TRUE)) {
            warning(gettextf("unable to create %s", sQuote(dirPath)),
                    domain = NA)
        }
    else if (!file.info(dirPath)$isdir)
        stop("path exists but is not a directory: ", sQuote(dirPath))
    dirPath
}

# appends resources to cache path
cache.resources <- function(hubCache, path=character()) {
    url <- file.path(hubCache, "resources")
    if (length(path))
        url <- file.path(url, .pathToLocalPath(path))
    url
}