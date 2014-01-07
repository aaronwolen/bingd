#' Search AnnotationHub for features
#' 
#' @param filter named list of vectors
#' @inheritParams get.metadata
#' 
#' @export

feature.search <- function(data.filter, genome, md, online = NULL) {

  if (missing(md)) md <- get.metadata(online = online)
  
  # Filter based on genome
  if (!missing(genome) & "genome" %in% colnames(md)) 
    md <- md[md$Genome %in% genome,]
  
  # Apply user specified filters
  mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
    hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
    which(rowSums(hits) == length(pattern))
  }
  
  filter.hits <- lapply(data.filter, mgrep, x = md$RDataPath)
  filter.hits <- lapply(filter.hits, function(x) make.names(md$RDataPath[x]))
  
  return(filter.hits)
}


#' Retrieve AnnotationHub metadata
#' 
#' @param online logical, should metadata be retrieved online
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' @importFrom Biobase testBioCConnection

get.metadata <- function(online = NULL) {
  
  if (is.null(online)) online <- testBioCConnection()
  
  if (online) {
    md <- metadata()
  } else {
    path <- AnnotationHub:::hubCache()
    md <- dir(file.path(path, "resources"), recursive = TRUE)
    md <- DataFrame(RDataPath = md, Title = feature.labels(md))
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

