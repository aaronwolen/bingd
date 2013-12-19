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

annotate.gwas <- function(gwas, data.filter) {

  message("Creating AnnotationHub object...")
  hub <- AnnotationHub()
  md <- metadata(hub)
  
  filter.hits <- feature.search(data.filter, genome(gwas), md, hub)
  
  # Retrieve features and check for gwas overlaps
  overlaps <- mclapply(unlist(filter.hits), function(f)
                     gwas %over% AnnotationHub:::.getResource(hub, f),
                     mc.cores = getDoParWorkers())
  
  names(overlaps) <- md$Description[match(unlist(filter.hits), make.names(md$RDataPath))]
  
  # It'd be nice if overlap results for each feature type could
  # be stored in different slots but for now all features are just
  # added as ordinary variables
  overlaps <- DataFrame(overlaps)
  # Features will be denoted by a .prefix
  names(overlaps) <- paste0(".", names(overlaps))
  
  mcols(gwas) <- DataFrame(mcols(gwas), overlaps)
  
  return(gwas)
}