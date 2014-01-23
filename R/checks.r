#' FeatureList object check
is.FeatureList <- function(x) {
   if (class(x)[1] != "FeatureList") {
    stop("Must provide a FeatureList object.")
  }
  return(x)
}
 