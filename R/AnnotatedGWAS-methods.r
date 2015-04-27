#' Access AnnotatedGWAS feature index
#' @param object \code{AnnotatedGWAS} object

setGeneric("featureIndex", function(object) {
  standardGeneric("featureIndex")
})

setMethod("featureIndex", "AnnotatedGWAS", function(object) object@featureIndex)


#' Extract list of features from AnnotatedGWAS object
#' @inheritParams featureIndex
#' @exportMethod features

setGeneric("features", function(object) {
  standardGeneric("features")
})

#' @rdname features
setMethod("features", "AnnotatedGWAS", function(object) {  
  out <- lapply(featureIndex(object), function(x) mcols(object)[x])
  return(DataFrameList(out))
})


#' Extract features from AnnotatedGWAS object
#' @inheritParams featureIndex
#' @exportMethod fcols

setGeneric("fcols", function(object) {
  standardGeneric("fcols")
})

#' @rdname fcols
setMethod("fcols", "AnnotatedGWAS", function(object) {
  out <- structure(as.list(features(object)), names = NULL)
  return(DataFrame(out))
})


#' Consolidate each group of annotated features
#'
#' consolidate combines the features within each group to produce a new single
#' feature that comprises all of the original individual features
#' 
#' @param object \code{AnnotatedGWAS} object
#' @return \code{AnnotatedGWAS} object
#' 
#' @exportMethod consolidate

setGeneric("consolidate", function(object) {
  standardGeneric("consolidate")
})

#' @rdname consolidate
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

# Check if AnnotatedGWAS features have been consolidated
setGeneric("is.consolidated", function(object) standardGeneric("is.consolidated"))

setMethod("is.consolidated", "AnnotatedGWAS", function(object) {
  all(sapply(features(object), length) == 1)
})



#' @describeIn as.GWAS Drop annotations to coerce an \code{\link{AnnotatedGWAS}} object back to a \code{GWAS} object

setMethod("as.GWAS", signature = "AnnotatedGWAS", 
  function(object) {
    
    fvars <- unlist(object@featureIndex)
    mvars <- setdiff(names(mcols(object)), fvars)
    mcols(object) <- mcols(object)[mvars]
    
    or <- beta <- NULL
    if ("or"   %in% mvars) or   <- "or"
    if ("beta" %in% mvars) beta <- "beta"
    
    as.GWAS(as(object, "GRanges"),
            genome = GenomeInfoDb::genome(object)[1],
            marker = "marker",
            pvalue = "pvalue",
            zscore = "zscore",
            or     = or,
            beta   = beta)
})

