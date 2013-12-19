#' Search AnnotationHub for features
#' 
#' @param filter named list of vectors
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' @export

feature.search <- function(data.filter, genome, hub, md) {

  if (missing(md)) {
    if(missing(hub)) {
      hub <- AnnotationHub()
    }
    md <- metadata(hub)
  }
  
  # Filter based on genome
  if (!missing(genome)) md <- md[md$Genome %in% genome,]
  
  # Apply user specified filters
  mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
    hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
    which(rowSums(hits) == length(pattern))
  }
  
  filter.hits <- lapply(data.filter, mgrep, x = md$RDataPath)
  filter.hits <- lapply(filter.hits, function(x) make.names(md$RDataPath[x]))
  
  return(filter.hits)
}
