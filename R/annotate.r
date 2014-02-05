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

setMethod("annotate.gwas", c(object = "GWAS", feature.list = "FeatureList"), 
  function(object, feature.list) {

    overlaps <- overlapsAny(query = object, subject = feature.list)

    # It'd be nice if overlap results for each feature type could
    # be stored in different slots but for now all features are just
    # added as ordinary variables
    overlaps <- DataFrame(overlaps)
    
    # Features will be denoted by a .prefix
    names(overlaps) <- paste0(".", names(overlaps))
    
    mcols(object) <- DataFrame(mcols(object), overlaps)
      
    return(object)
})



#' Find features overlapping GWAS markers
#' 
#' Determine which GWAS markers overlap each \code{FeatureList} feature.
#' 
#'  @return A list of \code{DataFrame}s, the structure of which mirrors that of
#'  the supplied \code{FeatureList} object. Each element of the list is a
#'  \code{DataFrame} with one logical vector for each feature, indicating
#'  marker overlap. 
#' 
#' @param query \code{GWAS} object
#' @param subject \code{FeatureList} object
#' @param ... additional arguments passed to \code{\link[IRanges]{overlapsAny}}
#' 
#' @export
#' @aliases overlapsAny-method

setMethod("overlapsAny", c(query = "GWAS", subject = "FeatureList"),
  function(query, subject, ...) {
    
    result <- list()
    
    for (i in names(subject)) {
      
      paths  <- structure(subject[[i]]$LocalPath, names =  subject[[i]]$Title)
      
      result[[i]] <- mclapply(paths, function(p) query %over% load.feature(p),
                              mc.cores = getDoParWorkers())
    }
    
    return(result)
})
