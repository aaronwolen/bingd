
# Define ------------------------------------------------------------------

#' AnnotatedGWAS-class
#' 
#' @aliases AnnotatedGWAS

setClass("AnnotatedGWAS", contains = "GRanges",  
         representation = representation(featureIndex = "SimpleList"))


# Validate ----------------------------------------------------------------

.validAnnotatedGWAS <- function(object) {
  
  # GWAS object is valid
  stopifnot(.validGWAS(object))
  
  # All indexed features are present
  index.features <- unlist(object@featureIndex)
  if (!all(index.features %in% names(mcols(object)))) {
    stop("GWAS object is missing features present in the featureIndex.\n",
         call. = FALSE)
  }
  
  # Features are logical
  class.features <- sapply(mcols(object)[, index.features], class)
  if (!all(class.features == "logical")) {
    stop("Annotated features must be logical.\n", call. = FALSE)
  }
  
  return(TRUE)
}

setValidity("AnnotatedGWAS", .validAnnotatedGWAS)
