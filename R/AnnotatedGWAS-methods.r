#' Access AnnotatedGWAS feature index
#' 
#' @param object \code{AnnotatedGWAS} object

setGeneric("featureIndex", function(object) {
  standardGeneric("featureIndex")
})

setMethod("featureIndex", "AnnotatedGWAS", function(object) object@featureIndex)


#' Extract list of features from AnnotatedGWAS object
#' 
#' @inheritParams featureIndex

setGeneric("features", function(object) {
  standardGeneric("features")
})

setMethod("features", "AnnotatedGWAS", function(object) {
  out <- lapply(featureIndex(object), function(x) 
                rename(mcols(object)[, x], structure(names(x), names = x)))
  return(DataFrameList(out))
})


#' Extract features from AnnotatedGWAS object
#' 
#' @inheritParams featureIndex

setGeneric("fcols", function(object) {
  standardGeneric("fcols")
})

setMethod("fcols", "AnnotatedGWAS", function(object) {
  
  out <- DataFrame(features(object))
  names(out) <- paste0(".", names(out))
  
  return(out)
})
