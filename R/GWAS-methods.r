# Accessors ---------------------------------------------------------------

#' Access GWAS pvalues
#' @param object \code{GWAS} or \code{AnnotatedGWAS} object
#' @export
setGeneric("pvalue", function(object) standardGeneric("pvalue"))
setMethod( "pvalue", "GWAS",          function(object) mcols(object)$pvalue)

#' Access GWAS markers
#' @param object \code{GWAS} or \code{AnnotatedGWAS} object
#' @export
setGeneric("marker", function(object) standardGeneric("marker"))
setMethod( "marker", "GWAS",          function(object) mcols(object)$marker)

#' Access GWAS z-scores
#' @param object \code{GWAS} or \code{AnnotatedGWAS} object
#' @export
setGeneric("zscore", function(object) standardGeneric("zscore"))
setMethod( "zscore", "GWAS",          function(object) mcols(object)$zscore)
zvalue <- function(object) zscore(object)


# Summaries ---------------------------------------------------------------

#' Summarize AnnotatedGWAS object
#' @param object \code{AnnotatedGWAS} object
#' @export

setGeneric("summary", function(object) standardGeneric("summary"))

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
