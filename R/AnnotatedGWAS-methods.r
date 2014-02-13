#' Access AnnotatedGWAS feature index
#' @param object \code{AnnotatedGWAS} object

setGeneric("featureIndex", function(object) {
  standardGeneric("featureIndex")
})

setMethod("featureIndex", "AnnotatedGWAS", function(object) object@featureIndex)


#' Extract list of features from AnnotatedGWAS object
#' @inheritParams featureIndex

setGeneric("features", function(object) {
  standardGeneric("features")
})

setMethod("features", "AnnotatedGWAS", function(object) {
  out <- lapply(featureIndex(object), function(x) mcols(object)[, x])
  return(DataFrameList(out))
})


#' Extract features from AnnotatedGWAS object
#' @inheritParams featureIndex

setGeneric("fcols", function(object) {
  standardGeneric("fcols")
})

setMethod("fcols", "AnnotatedGWAS", function(object) {
  out <- structure(as.list(features(object)), names = NULL)
  return(DataFrame(out))
})
