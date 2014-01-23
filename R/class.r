#' Create a \code{FeatureList} object
#' 
#' @param x a \code{\link{list}} of \code{\link{DataFrame}} objects

as.FeatureList <- function(x) {
  stopifnot(class(x) == "list")
  class(x) <- c("FeatureList", class(x))
  return(x)
}
