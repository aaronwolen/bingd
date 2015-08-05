setClass("AnnotationHubList",
         prototype = prototype(elementType = "AnnotationHub"),
         contains = "SimpleList")

AnnotationHubList <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]])) args <- args[[1]]
  if (!all(vapply(args, is, "AnnotationHub", FUN.VALUE = logical(1))))
            stop("all elements in '...' must be AnnotationHub objects")
  new("AnnotationHubList", S4Vectors::SimpleList(args))
}


setMethod("show", "AnnotationHubList", 
  function(object) {
    len <- length(object)
    cat("AnnotationHubList object with", len,
      ifelse(len == 1, "group", "groups"),
      "of features.\n\n")
    
    for (i in seq_len(len)) {
      cat("\nGroup ", i, " - ", names(object)[i], "\n", sep = "")
      cat(rep("-", getOption("width", 80) - 1), "\n", sep = "")
      AnnotationHub::show(object[[i]])
    }
})


setGeneric("cache", 
  function(x, ..., max.downloads = hubOption("MAX_DOWNLOADS")) {
    standardGeneric("cache")
})
   
setMethod("cache", "AnnotationHub",
  function(x, ..., max.downloads = hubOption("MAX_DOWNLOADS")) {
    AnnotationHub::cache(x, ..., max.downloads = max.downloads)
})

setMethod("cache", "AnnotationHubList",
  function(x, ..., max.downloads = hubOption("MAX_DOWNLOADS")) {
    if (is.null(names(x))) names(x) <- paste0("Group", seq_along(x))
    out <- setNames(vector(mode = "list", length(x)), names(x))
    
    for (g in names(x)) {
      cat("\nCaching", g, "\n")
      out[[g]] <- AnnotationHub::cache(x[[g]], ..., max.downloads = max.downloads)
    }
    out
})


setGeneric("hubCache", function(x) standardGeneric("hubCache"))
setMethod("hubCache", "AnnotationHub", function(x) AnnotationHub::hubCache(x))

setMethod("hubCache", "AnnotationHubList",
  function(x) lapply(x, AnnotationHub::hubCache))
