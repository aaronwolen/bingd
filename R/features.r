#' Load AnnotationHub features
#' 
#' Check for feature in local cache before AnnotationHub
#' 
#' @param path Path of feature to load
#' @inheritParams hub.search

load.feature <- function(path) {
  
  if (file.exists(path)) {
    # Load in temp environment so object name can be determined
    obj.env <- new.env()
    obj.name <- load(path, obj.env)
    obj <- get(obj.name, envir = obj.env)
  } else {
    stop("No such file:\n", path, call. = FALSE)
  }
  
  return(obj)
}



#' Search AnnotationHub for features
#' 
#' @inheritParams filter.features
#' @inheritParams hub.metadata
#' 
#' @export
#' 
#' @return A named list of DataFrames containing (at minimum) columns indicating
#' the Title and location (LocalPath) of features matching search query

hub.search <- function(query, genome, md, online = FALSE, cache.path = "default") {
  
  if (missing(md)) md <- hub.metadata(online = online, cache.path = cache.path)
  
  # Filter based on genome
  if (!missing(genome) & "genome" %in% colnames(md)) 
    md <- md[md$Genome %in% genome,]
  
  filter.features(query, md)
}



#' Retrieve information about AnnotationHub features
#' 
#' @param online if TRUE search is conducted using latest AnnotationHub metadata, otherwise only the AnnotationHub cache directory is searched
#' @param cache.path define custom location used for AnnotationHub's cached
#' directory, which should contain a 'resources' subdirectory. Normally this
#' shouldn't be changed.
#' @importFrom AnnotationHub AnnotationHub metadata
#' @importFrom Biobase testBioCConnection

hub.metadata <- function(online = FALSE, cache.path = "default") {
  
  if (cache.path == "default") cache.path <- AnnotationHub:::hubCache()
  
  # Catalog cached files
  if (!grepl("resources", cache.path)) cache.path <- file.path(cache.path, "resources")
  cache.files <- local.features(NULL, cache.path)
  
  if (online) {
    # Retrieve latest features and identify which are already cached
    md <- metadata()
    local.path <- file.path(cache.path, md$RDataPath)
    md$LocalPath <- ifelse(local.path %in% cache.files$LocalPath, local.path, NA)
  } else {
    md <- cache.files
  }
  
  return(md)
}



#' Retrieve information about local features
#' 
#' @param path Path of directory containing features
#' @inheritParams filter.features

local.features <- function(query = NULL, path) {
  
  files <- dir(path, full.names = TRUE, recursive = TRUE)
  
  flist <- DataFrame(Title = feature.labels(files), LocalPath = files)
  rownames(flist) <- NULL # DataFrame (1.20.6) doesn't respect row.names = NULL
  
  return(flist)
}

#' Filter list of features based on search terms
#'
#' @param query a vector of strings to filter features; may also be a named list
#' of vectors to perform multiple queries for different categories of features
#' @param file.list a \code{\link{DataFrame}} containing, at minimum,
#' \code{Title} and \code{LocalPath} columns

filter.features <- function(query, file.list) {
  
  if(!all(c("LocalPath", "Title") %in% names(file.list))) {
    stop("file.list must contain LocalPath and Title columns", call. = FALSE)
  }
  if (is.atomic(query)) query <- list(query)
  
  mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
    hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
    which(rowSums(hits) == length(pattern))
  }
  
  query.hits <- lapply(query, mgrep, x = file.list$LocalPath)
  query.hits <- lapply(query.hits, function(x) file.list[x, ])
  
  query.hits <- as.FeatureList(query.hits)
  return(query.hits)
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

