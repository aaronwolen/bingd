
# Define ------------------------------------------------------------------

setClass("FeatureList", contains="SimpleDataFrameList")


# Validate ----------------------------------------------------------------

.validFeatureList <- function(object) {
  
  lapply(object, .validColClass, "FeatureList", "Title",     "character")
  lapply(object, .validColClass, "FeatureList", "LocalPath", "character")
  lapply(object, .validColClass, "FeatureList", "Cached",    "logical")

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
    object <- eval(objs[[1]], envir = parent.frame())
    
    if (is.list(object) & !is.data.frame(object)) {
      out <- object
    } else if (class(object) %in% c("data.frame", "DataFrame")) {
      out <- list(object)
    } else {
      stop("Must provide data.frames, DataFrames, or a list of either class.")
    }
  } else {
    out <- lapply(objs, eval, envir = parent.frame())
  }
  
  if (!missing(names)) {
    if (length(names) != length(out)) {
      stop("Length of names differs from length of supplied object(s)")
    } else {
     base::names(out) <- names 
    }
  } else {
    if (is.null(base::names(out))) {
      if (length(out) == 1) {
        base::names(out) <- "features"
      } else {
        base::names(out) <- paste0("features", seq_along(out))
      }
    }
  }

 new("FeatureList", DataFrameList(out))
}

