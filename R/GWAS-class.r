
# Define ------------------------------------------------------------------

setClass("GWAS", contains="GRanges")


# Validate ----------------------------------------------------------------

.validGWAS <- function(object) {
  
  columns.req <- c(marker = "character", pvalue = "numeric")
  columns.opt <- c(or = "numeric", beta = "numeric")
  
  md <- mcols(object)
  
  # Required columns
  if (!all(names(columns.req) %in% names(md))) {
    return(cat("GWAS object metadata must contain all of the following columns:\n", 
         paste(names(columns.req), collapse = "\n")))
  }
  
  # Optional columns
  if (!any(names(columns.opt) %in% names(md))) {
    return(cat("GWAS object metadata must contain one of the following columns:\n", 
         paste(names(columns.opt), collapse = "\n")))
  }
  
  # Class information of required columns
  col.class <- c(columns.req, columns.opt)
  col.class <- col.class[which(names(col.class) %in% names(md))]
  
  md.class <- sapply(md[, names(col.class)], class)
  
  if (!identical(col.class[names(col.class)], md.class)) {
    class.info <- paste(names(col.class), col.class, sep = " = ")
    return(cat("Required GWAS metadata classes should be:\n",
         paste(class.info, collapse = "\n")))
  }
  
  return(TRUE)
}

setValidity("GWAS", .validGWAS)


# Constructors ------------------------------------------------------------

#' Create a GWAS object
#' 
#' @param snpid column label corresponding to SNP ID column
#' @param chr column label corresponding to chromosome column
#' @param pos
#' @param pval
#' @param or

setGeneric("as.GWAS", 
  function(object, marker, chr, bp, pvalue, or, beta) {
    standardGeneric("as.GWAS")
}) 

setMethod("as.GWAS", "data.frame", 
  function(object, marker, chr, bp, pvalue, or, beta) {
    
    object[[chr]] <- paste0("chr", gsub("[chr]", "", object[[chr]]))
    
    gr <- GRanges(object[[chr]], IRanges(object[[bp]], width = 1), strand = "*", 
                  marker = object[[marker]], pvalue = object[[pvalue]])
    
    if (!missing(or))   gr$or    <- object[[or]]
    if (!missing(beta)) gr$beta  <- object[[beta]]
    
    as.GWAS(gr)
})

setMethod("as.GWAS", "GRanges", 
  function(object, marker, chr, bp, pvalue, or, beta) {
      
    md <- mcols(object)
    if ("or"   %in% names(md)) stat <- "or"
    if ("beta" %in% names(md)) stat <- "beta"
    
    object$z <- ifelse(md[[stat]] > 1, qnorm(md$pvalue / 2), -qnorm(md$pvalue / 2))
    return(object)      
})
