#' Load AnnotationHub features
#' 
#' Check for feature in local cache before AnnotationHub
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
#' @inheritParams filter.features
#' @importFrom AnnotationHub AnnotationHub metadata
#' 
#' @export

hub.features <- function(query = NULL, path, genome, online = FALSE) {
   
  if (missing(path)) path <- cache.path()
  if (!grepl("resources", path)) path <- file.path(path, "resources")
  if (!cache.exists(path)) cache.create(path)
  
  cached.files <- local.features(NULL, path)
  
  if (online) {
    # Retrieve latest feature
    f.files <- metadata()
    
    # Filter based on genome
    if (!missing(genome)) f.files <- f.files[f.files$Genome == genome,]
    
    # Identify which features are already cached
    f.files$LocalPath <- file.path(path, f.files$RDataPath)
    
    if (!is.null(cached.files)) {
      f.files$Cached <- ifelse(f.files$LocalPath %in% cached.files$LocalPath, T, F)  
    } else {
      f.files$Cached <- FALSE
    }
    
  } else {
    if (is.null(cached.files)) stop("No files found in ", path)
    f.files <- cached.files
  }
  
  if (is.null(query)) return(f.files)
  filter.features(query, f.files)
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
  
  f.list <- DataFrame(Title = feature.labels(files), 
                       LocalPath = files, 
                       Cached = TRUE)
  rownames(f.list) <- NULL # DataFrame (1.20.6) doesn't respect row.names=NULL
  names(f.list$Title) <- NULL
  f.list <- FeatureList(f.list)

  if (is.null(query)) return(f.list)
  filter.features(f.list, query)
}



#' Filter FeatureList object based on search terms
#'
#' Search terms can be grouped together using a list of queries.
#'
#' @param query a vector of strings to filter features; may also be a named list
#' of vectors to perform multiple queries for different categories of features
#' @param query a \code{\link{FeatureList}} object comprising one or more
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
    object <- subset(stack(object), select = -name)
    
    # Feature group names
    if(is.null(names(query))) {
      f.names <- paste0("features", seq_along(query))
    } else {
      f.names <- names(query)
    }
      
    # multi-grep: pattern can be a vector of multiple character strings 
    mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
      hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
      which(rowSums(hits) == length(pattern))
    }
    
    query.hits <- lapply(query, mgrep, x = object$LocalPath)
    object <- lapply(query.hits, function(x) object[x, ])
    
    return(FeatureList(object))
})



#' Create pretty feature labels from the full RDataPaths
#' 
#' @param x characer vector of full RDataPaths

feature.labels <- function(x) {
  out <- basename(x)
  out <- sapply(strsplit(out, "\\.", fixed = F), function(x) x[1])
  names(out) <- x
  return(out)  
}

