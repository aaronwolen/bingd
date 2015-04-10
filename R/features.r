#' Load a local GRanges feature
#' 
#' \code{path} must point to an \code{R} object that contains a 
#' \code{GRanges} feature
#' 
#' @param path Path of feature to load

load.feature <- function(path) {
  
  if (file.exists(path)) {
    # Load in temp environment so object name can be determined
    obj.env <- new.env()
    obj.name <- load(path, obj.env)
    obj <- get(obj.name, envir = obj.env)
  } else {
    stop("No such file:\n", path, call. = FALSE)
  }
  
  if (class(obj) != "GRanges")
    stop("Features must be GRanges objects.")
  
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
#' the Title, location (LocalPath) of features matching search query and whether
#' the files are downloaded (Cached)
#' 
#' @param genome filter results based on genome version
#' @inheritParams filter.features
#' 
#' 
#' @export

hub.features <- function(query = NULL, path, genome, online = TRUE) {
  
  if (missing(path)) {
    path <- cache.path()
  } else {
    if (!cache.exists(path)) cache.create(path)
    options("AnnotationHub.Cache" = path)
  }
  
  cached.files <- local.features(NULL, cache.resources(path))
  if (!is.null(cached.files)) cached.files <- stack(cached.files)[-1]
  
  if (online) {
    # Retrieve latest feature
    ah <- AnnotationHub::AnnotationHub()
    f.files <- AnnotationHub::metadata(ah)
    
    # Filter based on genome
    if (!missing(genome)) f.files <- f.files[f.files$Genome == genome,]
    
    # Identify which features are already cached
    f.files$LocalPath <- file.path(cache.resources(path), f.files$RDataPath)
    f.files$Cached    <- file.exists(f.files$LocalPath)
  } else {
    if (is.null(cached.files)) stop("No files found in ", path)
    f.files <- cached.files
  }
  
  f.list <- FeatureList(f.files)
  if (is.null(query)) return(f.list)
  filter.features(f.list, query)
}



#' Retrieve information about local features
#' 
#' @param path Path of directory containing features
#' @inheritParams filter.features
#' 
#' @export

local.features <- function(query = NULL, path) {
  
  files <- dir(path, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(NULL)
  
  f.files <- DataFrame(Title = feature.labels(files), 
                       LocalPath = files, 
                       Cached = TRUE)
  rownames(f.files) <- NULL # DataFrame (1.20.6) doesn't respect row.names=NULL
  names(f.files$Title) <- NULL
  
  f.list <- FeatureList(f.files)
  if (is.null(query)) return(f.list)
  filter.features(f.list, query)
}



#' Filter FeatureList object based on search terms
#'
#' Search terms can be grouped together using a list of queries.
#'
#' @param query a vector of strings to filter features; may also be a named list
#' of vectors to perform multiple queries for different categories of features
#' @param object a \code{\link{FeatureList}} object comprising one or more
#' DataFrames that containing, at minimum, \code{Title} \code{LocalPath}, and
#' \code{Cached} columns

setGeneric("filter.features", 
  function(object, query) {
    standardGeneric("filter.features")
})
           
setMethod("filter.features", "FeatureList",
  function(object, query) {

    if (is.atomic(query)) query <- list(query)
    
    # FeatureList groupings are ignored in favor of query groupings
    object <- stack(object)[-1]
    
    # Feature group names
    if(is.null(names(query))) {
      f.names <- paste0("features", seq_along(query))
    } else {
      f.names <- names(query)
    }
      
    # multi-grep: pattern can be a vector of multiple character strings 
    mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
      hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
      hits <- .rowSums(hits, m = length(x), n = length(pattern))
      which(hits == length(pattern))
    }
    
    query.hits <- lapply(query, mgrep, x = object$LocalPath)
    object <- lapply(query.hits, function(x) object[x, ])
    
    return(FeatureList(object))
})



#' Create pretty feature labels from the full RDataPaths
#' 
#' Basically just strips off the extension from the basename
#' 
#' @param x characer vector of full RDataPaths

feature.labels <- function(x) {
  out <- basename(x)
  ext.pos <- sapply(gregexpr("\\.", out), function(x) tail(x, 1))
  out <- substr(out, 1, stop = ext.pos - 1)
  names(out) <- x
  return(out)  
}

