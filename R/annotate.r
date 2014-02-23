#' Annotate GWAS makers with genomic features
#' 
#' @param object \code{\link{GWAS}} object
#' @param feature.list \code{FeatureList} object created with 
#' \code{\link{hub.features}} or \code{\link{local.features}}
#' 
#' @importFrom foreach getDoParWorkers
#' @importFrom parallel mclapply
#' @export
#' 
#' @return \code{\link{AnnotatedGWAS}} object
#' 
#' @examples
#' query <- list(DNaseI = c("Dnasen", "broadPeak", "rep1"))
#' feature.list <- hub.features(query)
#' scz <- annotate.gwas(scz, feature.list)


setGeneric("annotate.gwas", 
  function(object, feature.list) {
    standardGeneric("annotate.gwas")
})

setMethod("annotate.gwas", c(object = "GWAS", feature.list = "FeatureList"), 
  function(object, feature.list) {
    
    overlaps <- featureOverlaps(query = object, subject = feature.list)
    
    f.index <- lapply(overlaps, function(x) make.names(names(x)))
    overlaps <- do.call("DataFrame", as(overlaps, "list"))
    names(overlaps) <- unlist(f.index)
    
    mcols(object) <- DataFrame(mcols(object), overlaps)
    
    object <- new("AnnotatedGWAS", object, 
                  featureIndex = as(f.index, "SimpleList"))
      
    return(object)
})
