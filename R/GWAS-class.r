
# Define ------------------------------------------------------------------

setClass("GWAS", contains="GRanges")


# Validate ----------------------------------------------------------------

.validGWAS <- function(object) {
  
  columns.req <- c(marker = "character", pvalue = "numeric")
  columns.opt <- c(or = "numeric", beta = "numeric")
  
  md <- mcols(object)
  
  # Required columns
  if (!all(names(columns.req) %in% names(md))) {
    stop("GWAS object metadata must contain all of the following columns:\n", 
         paste(names(columns.req), collapse = "\n"), call. = FALSE)
  }
  
  # Optional columns
  if (!any(names(columns.opt) %in% names(md))) {
    stop("GWAS object metadata must contain one of the following columns:\n", 
         paste(names(columns.opt), collapse = "\n"), call. = FALSE)
  }
  
  # Class information of required columns
  col.class <- c(columns.req, columns.opt)
  col.class <- col.class[which(names(col.class) %in% names(md))]
  
  md.class <- sapply(md[, names(col.class)], class)
  
  if (!identical(col.class[names(col.class)], md.class)) {
    class.info <- paste(names(col.class), col.class, sep = " = ")
    stop("Required GWAS metadata classes should be:\n",
         paste(class.info, collapse = "\n"), call. = FALSE)
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
    
    gr <- GRanges(object[[chr]], IRanges(object[[bp]], width = 1), strand = "*")
    
    # Attach all additional columns in on object as metadata
    mcols(gr) <- object[, -match(c(chr, bp), names(object))]
    
    as.GWAS(gr, marker = marker, pvalue = pvalue, or = or, beta = beta)
})

setMethod("as.GWAS", "GRanges", 
  function(object, marker, chr, bp, pvalue, or, beta) {
    
    # Rename required metadatda columns
    req.cols <- structure(c("marker", "pvalue"), names = c(marker, pvalue))
    if (!missing(or))   req.cols <- c(req.cols, structure("or",  names = or))
    if (!missing(beta)) req.cols <- c(req.cols, structure("beta", names = beta))
    if (missing(or) & missing(beta)) 
      stop("Must supply odds ratios (or) or beta values (beta).")
      
    mcols(object) <- rename(mcols(object), req.cols)
    
    if (is.factor(object$marker)) object$marker <- as.character(object$marker)
    
    # Calculate z-score
    object$z <- calc.z(object$pvalue, 
                       or   = ifelse(missing(or),   NULL, object$or),
                       beta = ifelse(missing(beta), NULL, object$beta))
    
    return(new("GWAS", object))      
})

