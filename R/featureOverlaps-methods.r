#' Find features overlapping GWAS markers
#' 
#' Determine which GWAS markers overlap each \code{FeatureList} feature.
#'
#' @return 
#' A list of \code{DataFrame}s, the structure of which mirrors that of
#' the supplied \code{FeatureList} object. Each element of the list is a
#' \code{DataFrame} with one logical vector for each feature, indicating
#' marker overlap.
#'  
#' @param query \code{GWAS} object
#' @param subject \code{FeatureList} object
#' @inheritParams IRanges::findOverlaps
#' @param minoverlap see documentation for \code{maxgap} argument
#' 
#' @export

setGeneric("featureOverlaps", 
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"), ...) {
    standardGeneric("featureOverlaps")
})


setMethod("featureOverlaps", c("GWAS", "FeatureList"),
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"), ...) {

    type <- match.arg(type)
    
    result <- findOverlaps(query, subject, 
                           maxgap = maxgap, minoverlap = minoverlap, type = type,
                           select = "arbitrary", ...)
    
    result <- lapply(result, function(x) !is.na(x))
    result <- split(result, f = stack(LocalPath(subject))$ind)
    result <- lapply(result, DataFrame)
    
    result <- DataFrameList(result)
  
    return(result)
})
