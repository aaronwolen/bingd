# Accessors ---------------------------------------------------------------

setGeneric("pvalue", function(object) standardGeneric("pvalue"))
setMethod( "pvalue", "GWAS",          function(object) mcols(object)$pvalue)
setMethod( "pvalue", "AnnotatedGWAS", function(object) mcols(object)$pvalue)
  
setGeneric("marker", function(object) standardGeneric("marker"))
setMethod( "marker", "GWAS",          function(object) mcols(object)$marker)
setMethod( "marker", "AnnotatedGWAS", function(object) mcols(object)$marker)

setGeneric("zscore", function(object) standardGeneric("zscore"))
setMethod( "zscore", "GWAS",          function(object) mcols(object)$zscore)
setMethod( "zscore", "AnnotatedGWAS", function(object) mcols(object)$zscore)
zvalue <- function(object) zscore(object)
