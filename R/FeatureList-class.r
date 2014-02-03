
# Define ------------------------------------------------------------------

setClass("FeatureList", contains="SimpleDataFrameList")


# Validate ----------------------------------------------------------------

.validFeatureList <- function(object) {
  .validColClass(object, "FeatureList", "Title",     "character")
  .validColClass(object, "FeatureList", "LocalPath", "character")
  .validColClass(object, "FeatureList", "Cached",    "logical")
  return(TRUE)
}

setValidity("FeatureList", .validFeatureList)

# Constructor -------------------------------------------------------------

#' @param ... either the names of \code{data.frames}, \code{DataFrames}, or a
#' list of either class, optionally named
#' @param names a character vector of names to assign to each element
#'
#' @return \code{FeatureList}

FeatureList <- function(..., names) {
  
  objs <- as.list(substitute(list(...))[-1])
    
  if (length(objs) == 1) {
    object <- eval(objs[[1]])
    
    if (is.list(object)) {
      out <- object
    } else if (class(object) %in% c("data.frame", "DataFrame")) {
      out <- list(object)
    } else {
      stop("Must provide data.frames, DataFrames, or a list of either class.")
    }
  } else {
    out <- lapply(objs, eval)
  }
  
  if (is.null(base::names(out))) {
    if (!missing(names)) {
      if (length(names) != length(out)) 
        stop("Length of names differs from length supplied objects.")
      base::names(out) <- names
    } else {
      base::names(out) <- paste0("features", seq_along(out))
    }
  }

 new("FeatureList", DataFrameList(out))
}

