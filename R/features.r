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



#' Retrieve information about AnnotationHub features
#' 
#' @param online if TRUE search is conducted using latest AnnotationHub metadata, otherwise only the AnnotationHub cache directory is searched
#' @param path define custom location used for AnnotationHub's cached
#' directory, which should contain a 'resources' subdirectory. Normally this
#' shouldn't be changed.
#' 
#' @return A named list of DataFrames containing (at minimum) columns indicating
#' the Title and location (LocalPath) of features matching search query
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' 
#' @export

hub.features <- function(query = NULL, path, genome, online = FALSE) {

  if (missing(path)) {
    path <- AnnotationHub:::hubCache()  
    if (!grepl("resources", path)) path <- file.path(path, "resources")
  } 
  
  cached.files <- local.features(NULL, path)
  
  if (online) {
    # Retrieve latest features and identify which are already cached
    flist <- metadata()
    
    # Filter based on genome
    if (!missing(genome)) flist <- flist[flist$Genome == genome,]
    
    local.path <- file.path(path, flist$RDataPath)
    flist$LocalPath <- ifelse(local.path %in% cached.files$LocalPath, 
                              local.path, NA)
  } else {
    flist <- cached.files
  }
  
  if (is.null(query)) return(flist)
  filter.features(query, flist)
}



#' Retrieve information about local features
#' 
#' @param path Path of directory containing features
#' @inheritParams filter.features
#' 
#' @export

local.features <- function(query = NULL, path) {

  files <- dir(path, full.names = TRUE, recursive = TRUE)
  
  flist <- DataFrame(Title = feature.labels(files), LocalPath = files)
  rownames(flist) <- NULL # DataFrame (1.20.6) doesn't respect row.names = NULL
  
  if (is.null(query)) return(flist)
  filter.features(query, flist)
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

