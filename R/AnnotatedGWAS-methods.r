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


#' Accessor for feature annotations stored in gwas object
#' 
#' Only needed until feature annotations are stored hierarchically by type
#' 
#' @param gwas a \code{\link{GRanges}} object
#' 
#' @export
#' 
#' @return 
#' List of data.frames in which each slot corresponds to a feature type
#' 

setGeneric("pull.features", function(gwas) {
  standardGeneric("pull.features")
})

setMethod("pull.features", "GRanges", function(gwas) {

  md.names <- names(mcols(gwas))
  feature.cols <- md.names[grep("^\\.", md.names)]
  
  if (length(feature.cols) == 0) return(feature.cols) 
  
  feature.vars <- do.call("rbind", strsplit(feature.cols, "\\."))
  feature.types <- feature.vars[, 2]
  feature.names <- split(feature.vars[, 3], feature.types)
  
  out <- split(feature.cols, feature.types)
  out <- lapply(out, function(x) mcols(gwas)[x])
  
  out <- mapply(function(x, y) {
                  names(x) <- y; x
                }, out, feature.names, SIMPLIFY = FALSE)
  return(out)
})


#' Check if a GWAS GRanges object is annotated
#' 
#' @inheritParams pull.features

setGeneric("is.annotated", function(gwas) {
  standardGeneric("is.annotated")
})

setMethod("is.annotated", "GRanges", function(gwas) {
  ifelse(length(pull.features(gwas)) > 0, TRUE, FALSE)
})

