# Accessors ---------------------------------------------------------------

setGeneric("pvalue", function(object) standardGeneric("pvalue"))
setMethod( "pvalue", "GWAS",          function(object) mcols(object)$pvalue)
setMethod( "pvalue", "AnnotatedGWAS", function(object) mcols(object)$pvalue)
  
setGeneric("marker", function(object) standardGeneric("marker"))
setMethod( "marker", "GWAS",          function(object) mcols(object)$marker)
setMethod( "marker", "AnnotatedGWAS", function(object) mcols(object)$marker)
