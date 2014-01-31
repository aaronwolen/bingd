#' Annotate GWAS object using AnnotationHub
#' 
#' @param object \code{\link{GWAS}} object
#' @param feature.list \code{FeatureList} object created with 
#' \code{\link{hub.features}} or \code{\link{local.features}}
#' 
#' @importFrom foreach getDoParWorkers
#' @importFrom parallel mclapply
#' @export
#' 
#' @examples
#' query <- list(DNaseI = c("Dnasen", "broadPeak", "rep1"))
#' feature.list <- hub.features(query)
#' scz <- annotate.gwas(scz, feature.list)


setGeneric("annotate.gwas", 
  function(object, feature.list) {
    standardGeneric("annotate.gwas")
})

setMethod("annotate.gwas", "GWAS", 
  function(object, feature.list) {
    
    feature.list <- is.FeatureList(feature.list)
    
    overlaps <- lapply(feature.list, function(df) 
                       mclapply(df$LocalPath, function(f)
                                object %over% load.feature(f),
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
    
    mcols(object) <- DataFrame(mcols(object), overlaps)
      
    return(object)
})
