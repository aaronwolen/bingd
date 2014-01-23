#' Annotate GWAS object using AnnotationHub
#' 
#' @param gwas GWAS GRanges object
#' @param feature.list \code{FeatureList} object created with \code{\link{hub.search}}
#' 
#' @importFrom foreach getDoParWorkers
#' @importFrom parallel mclapply
#' @export
#' 
#' @examples
#' query <- list(DNaseI = c("Dnasen", "broadPeak", "rep1"))
#' feature.list <- hub.search(query, online = FALSE)
#' scz <- annotate.gwas(scz, feature.list)

annotate.gwas <- function(gwas, feature.list) {
  
  feature.list <- is.FeatureList(feature.list)
  
  overlaps <- lapply(feature.list, function(df) 
                     mclapply(df$LocalPath, function(f)
                              gwas %over% load.feature(f),
                              mc.cores = getDoParWorkers()))
  
  # Label with feature titles
  overlaps <- mapply(function(f, o) {
                       names(o) <- f$Title
                       DataFrame(o)
                     }, feature.list, overlaps, SIMPLIFY = FALSE)
  
  # It'd be nice if overlap results for each feature type could
  # be stored in different slots but for now all features are just
  # added as ordinary variables
  overlaps <- DataFrame(overlaps)
  
  # Features will be denoted by a .prefix
  names(overlaps) <- paste0(".", names(overlaps))
  
  mcols(gwas) <- DataFrame(mcols(gwas), overlaps)
    
  return(gwas)
}
