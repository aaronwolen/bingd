library(IRanges)
setClass("Features", 
  slots = list(Title = "character", 
               LocalPath = "character", 
               Cached = "logical",
               Metadata = "DataFrame"),
  prototype = prototype(DataFrame()),
  contains = "DataFrame")
    

ex.files <- dir("data/resources/", "RData",
                all.files=T, recursive=T, full.names=T)
x <- new("Features", 
    Title = basename(ex.files), 
    LocalPath = ex.files, 
    Cached = TRUE)

setMethod("show", "Features",
  function(object) {
    out <- DataFrame(Title = object@Title,
                     Cached = object@Cached)
    out <- rbind(out, object@Metadata)
    show(out)
})

x

setGeneric("filterFeatures", 
  function(file.list, query) {
    standardGeneric("filterFeatures")
})

setMethod("filterFeatures", signature = "Features",

  function(file.list, query) {
    browser()
    mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
      hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
      which(rowSums(hits) == length(pattern))
    }
  
    query.hits <- lapply(query, mgrep, x = file.list@LocalPath)
    query.hits <- lapply(query.hits, function(x) file.list[x, ])
  
    query.hits <- as.FeatureList(query.hits)
    return(query.hits)
})

