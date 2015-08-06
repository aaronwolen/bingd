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
#' @exportMethod featureOverlaps

setGeneric("featureOverlaps", 
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"), ...) {
    standardGeneric("featureOverlaps")
})


#' @rdname featureOverlaps
setMethod("featureOverlaps", c("GWAS", "AnnotationHubList"),
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"), ...) {
    
    result <- lapply(subject, featureOverlaps, query = query, 
                     maxgap = maxgap, minoverlap = minoverlap, type = type, ...)

    DataFrameList(result)
})


setMethod("featureOverlaps", c("GWAS", "AnnotationHub"),
  function(query, subject, maxgap = 0L, minoverlap = 1L,
           type = c("any", "start", "end", "within"), ...) {
    
    result <- findOverlaps(query, subject, 
                           maxgap = maxgap, minoverlap = minoverlap, type = type,
                           select = "arbitrary", ...)
    
    result <- lapply(result, function(x) !is.na(x))
    DataFrame(result, row.names = names(query))
})
