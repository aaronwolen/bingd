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
#' @aliases overlapsAny-method
#' @export

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