# Unexported helper functions copied from AnnotationHub

cache.path <- function() AnnotationHub::hubCache()

cache.exists <- function(path) 
  !is.na(path) && isTRUE(file.info(path)$isdir)

cache.create <- function(path) {
    if (!file.exists(path))
        if (!dir.create(path, showWarnings=FALSE, recursive=TRUE)) {
            warning(gettextf("unable to create %s", sQuote(path)), domain = NA)
        } else {
          AnnotationHub::setHubOption("CACHE", path)
        }
    else if (!file.info(path)$isdir)
        stop("path exists but is not a directory: ", sQuote(path))
    path
}

# appends resources to cache path
cache.resources <- function(hubCache, path=character()) {
  file.path(hubCache, "resources")
}