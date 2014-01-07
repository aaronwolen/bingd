#' Annotate GWAS object using AnnotationHub
#' 
#' @param gwas GWAS GRanges object
#' @param filter named list of vectors
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' @importFrom foreach getDoParWorkers
#' @importFrom parallel mclapply
#' @export
#' 
#' @examples
#' scz <- annotate.gwas(scz, data.filter = list(DNaseI = c("UwDnaseN", "broadPeak")))

annotate.gwas <- function(gwas, data.filter, hub = NULL, ...) {
  
  filter.hits <- feature.search(data.filter, genome(gwas)[1], ...)
  
  # Retrieve features and check for gwas overlaps
  message("Annotating GWAS markers...")
  overlaps <- mclapply(unlist(filter.hits), function(f)
                     gwas %over% load.feature(f, hub = hub),
                     mc.cores = getDoParWorkers())
  
  # Use metadata labels
  md <- get.metadata(online = FALSE)
  md.index <- sapply(unlist(filter.hits), grep, make.names(md$RDataPath))
  names(overlaps) <- md$Title[md.index]
  
  # It'd be nice if overlap results for each feature type could
  # be stored in different slots but for now all features are just
  # added as ordinary variables
  overlaps <- DataFrame(overlaps)
  
  # Features will be denoted by a .prefix
  filter.types <- rep(names(filter.hits), each = sapply(filter.hits, length))
  names(overlaps) <- paste0(".", filter.types, ".", names(overlaps))
  
  mcols(gwas) <- DataFrame(mcols(gwas), overlaps)
  
  return(gwas)
}