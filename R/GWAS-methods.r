# Accessors ---------------------------------------------------------------

#' pvalue accessor and setter
#' 
#' Gets and sets an object's pvalues
#' 
#' @param object an object to get or set the pvalue of
#' @export
setGeneric("pvalue", function(object) standardGeneric("pvalue"))

#' @describeIn pvalue Get and set the pvalues of a \code{GWAS} object
setMethod( "pvalue", "GWAS",          function(object) mcols(object)$pvalue)


#' marker accessor and setter
#' 
#' Gets and sets an object's markers
#' 
#' @param object an object to get or set the markers of
#' @export
setGeneric("marker", function(object) standardGeneric("marker"))

#' @describeIn marker Get and set a \code{GWAS} object's markers
setMethod( "marker", "GWAS",          function(object) mcols(object)$marker)


#' zscore accessor and setter
#' 
#' Gets and sets an object's zscores
#' 
#' @param object an object to get or set the zscore of
#' @export
setGeneric("zscore", function(object) standardGeneric("zscore"))

#' @describeIn zscore Get and set a \code{GWAS} object's zscores
setMethod( "zscore", "GWAS",          function(object) mcols(object)$zscore)
zvalue <- function(object) zscore(object)


# Summaries ---------------------------------------------------------------

# Summarize AnnotatedGWAS object
# @param object \code{AnnotatedGWAS} object
setMethod("summary", "AnnotatedGWAS", function(object) {
  
  f <- features(object)
  f.len <- sapply(f, ncol)
  
  cat("AnnotatedGWAS object contains", length(object), "markers and",
      sum(f.len), "genomic features.\n")

  cat("\nProportion of markers overlapping features:\n\n")
  out <- lapply(f, function(g) sapply(g, mean))
  out <- lapply(out, format, digits = 2)
  
  name.width <- max(nchar(unlist(lapply(f, names))))
  
  for (g in names(f)) {
    g.df <- structure(data.frame(out[[g]]), 
                      names = g,
                      row.names = format(names(out[[g]]), 
                                         width = name.width,
                                         justify = "right"))
                      
    print(g.df, quote = FALSE)
    cat("\n")
  }

})
