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
  out <- lapply(featureIndex(object), function(x) mcols(object)[x])
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


#' Consolidate each group of annotated features
#'
#' consolidate combines the features within each group to produce a new single
#' feature that comprises all of the original individual features
#' 
#' @param \code{AnnotatedGWAS} object
#' @return \code{AnnotatedGWAS} object
#' 
#' @export

setGeneric("consolidate", function(object) {
  standardGeneric("consolidate")
})

setMethod("consolidate", "AnnotatedGWAS", function(object) {
  
  # Consolidate
  new.features <- lapply(features(object), 
                         function(f) DataFrame(n = rowSums(as.matrix(f)) > 0))
  new.features <- DataFrame(new.features)
  names(new.features) <- sub("\\.n$", "", names(new.features))
  
  # Replace features
  old.features <- unlist(featureIndex(object))
  mcols(object) <- mcols(object)[, setdiff(names(mcols(object)), old.features)]
  mcols(object) <- DataFrame(mcols(object), new.features)
  
  # Update index
  new.index <- as.list(structure(names(new.features), names = names(new.features)))
  object@featureIndex <- as(new.index, "SimpleList")
  
  return(object)
})

#' Check if AnnotatedGWAS features have been consolidated
setGeneric("is.consolidated", function(object) standardGeneric("is.consolidated"))

setMethod("is.consolidated", "AnnotatedGWAS", function(object) {
  all(sapply(features(object), length) == 1)
})

