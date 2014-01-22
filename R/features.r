#' Load AnnotationHub features
#' 
#' Check for feature in local cache before AnnotationHub
#' 
#' @param name Name of AnnotationHub feature to load
#' @inheritParams feature.search

load.feature <- function(name, hub = NULL) {
  
  path <- AnnotationHub:::hubCache()
  cached <- dir(path, recursive = TRUE)
  
  cached.match <- grepl(name, make.names(cached))
  
  if (any(cached.match)) {
    load(file.path(path, cached[cached.match]))
    return(obj)
  } else {
    if (testBioCConnection()) {
      stop(name, " isn't cached and there is no internet connection.")
    }
    if (is.null(hub)) hub <- AnnotationHub()
    obj <- AnnotationHub:::.getResource(hub, name)
  }
  return(obj)
}


#' Search AnnotationHub for features
#' 
#' @param data.filter named list of vectors
#' @inheritParams hub.metadata
#' 
#' @export
#' 
#' @return A named list of DataFrames containing (at minimum) columns indicating
#' the Title and location (LocalPath) of features matching search query

hub.search <- function(data.filter, genome, md, online = FALSE) {
  
  if (missing(md)) md <- hub.metadata(online = online)
  if (is.atomic(data.filter)) data.filter <- list(data.filter)
  
  # Filter based on genome
  if (!missing(genome) & "genome" %in% colnames(md)) 
    md <- md[md$Genome %in% genome,]
  
  # Apply user specified filters
  mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
    hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
    which(rowSums(hits) == length(pattern))
  }
  
  filter.hits <- lapply(data.filter, mgrep, x = md$RDataPath)
  filter.hits <- lapply(filter.hits, function(x) DataFrame(md[x, ]))                          
  
  return(filter.hits)
}


#' Retrieve AnnotationHub metadata
#' 
#' @param online if TRUE search is conducted using latest AnnotationHub metadata, otherwise only the AnnotationHub cache directory is searched
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' @importFrom Biobase testBioCConnection

hub.metadata <- function(online = FALSE, cache.dir) {
  
  if (missing(cache.dir)) cache.dir <- AnnotationHub:::hubCache()
  
  # Catalog cached files
  cache.files <- dir(cache.dir, full.names = TRUE, recursive = TRUE)
  
  if (online) {
    # Retrieve latest features and identify which are already cached
    md <- metadata()
    local.path <- file.path(cache.dir, "resources", md$RDataPath)
    md$LocalPath <- ifelse(local.path %in% cache.files, local.path, NA)
  } else {
    md <- DataFrame(Title = feature.labels(cache.files), LocalPath = cache.files)
    rownames(md) <- NULL # DataFrame (1.20.6) doesn't respect row.names = NULL
  }
  
  return(md)
}

#' Create pretty feature labels from the full RDataPaths
#' 
#' @param x characer vector of full RDataPaths

feature.labels <- function(x) {
  out <- basename(x)
  out <- sapply(strsplit(out, "\\.", fixed = F), function(x) x[1])
  names(out) <- x
  return(out)  
}

