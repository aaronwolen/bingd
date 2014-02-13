
# Define ------------------------------------------------------------------

#' GWAS-class
#' 
#' The GWAS class is a container for representing the genetic markers used in a
#' genome-wide associate study (GWAS), along with their corresponding summary 
#' statistics.
#' 
#' The GWAS class is a \link[GenomicRanges]{GRanges} with the following required
#' metadata columns:
#' 
#' \describe{
#'  \item{\code{marker}}{Marker or SNP identifiers}
#'  \item{\code{pvalue}}{Marker association p-values}
#' }
#' 
#' A GWAS object must also contain a metadata column indicating the direction of
#' the minor allele's effect, which can be provided as one of the following:
#' 
#' \describe{
#'  \item{\code{or}}{Odds ratios}
#'  \item{\code{beta}}{Beta values (i.e., the regression coefficient)}
#' }
#' 
#' @aliases GWAS

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
#' Create a \code{\link{GWAS}} object from a \code{\link{data.frame}} or a
#' \code{\link[GenomicRanges]{GRanges}} object.
#' 
#' @param object a \code{\link{data.frame}} or 
#' \code{\link[GenomicRanges]{GRanges}} object containing the required fields
#' necessary to construct a \code{\link{GWAS}} object
#' @param genome genome abbreviation used by UCSC (e.g., hg19)
#' @param marker name of the column in \code{object} that contains the marker
#' (or SNP) identifiers
#' @param chr name of the column in \code{object} that contains the chromosome
#' names associated with each marker
#' @param bp name of the column in \code{object} that contains the genomic 
#' position (or start position) associated with each marker
#' @param pvalue name of the column in \code{object} that contains the GWAS 
#' association p-value for each marker
#' @param or name of the column in \code{object} that contains the GWAS odds
#' ratio for each marker
#' @param beta name of the column in \code{object} that contains the GWAS beta
#' values (or regression coefficients) for each marker
#'
#' @export

setGeneric("as.GWAS", 
  function(object, genome, marker, chr, bp, pvalue, or, beta) {
    standardGeneric("as.GWAS")
}) 


#' Create a GWAS object
#' 
#' Create a \code{\link{GWAS}} object from a \code{\link{data.frame}}
#' 
#' @inheritParams as.GWAS

setMethod("as.GWAS", "data.frame", 
  function(object, genome, marker, chr, bp, pvalue, or, beta) {
    
    object[[chr]] <- paste0("chr", gsub("[chr]", "", object[[chr]]))
    
    gr <- GRanges(object[[chr]], IRanges(object[[bp]], width = 1), strand = "*")
    
    # Attach all additional columns in on object as metadata
    mcols(gr) <- object[, -match(c(chr, bp), names(object))]
    
    as.GWAS(gr, genome = genome, marker = marker, 
            pvalue = pvalue, or = or, beta = beta)
})


#' Create a GWAS object
#' 
#' Create a \code{\link{GWAS}} object from a \link[GenomicRanges]{GRanges}
#' object
#' 
#' @inheritParams as.GWAS

setMethod("as.GWAS", "GRanges", 
  function(object, genome, marker, pvalue, or, beta) {
    
    # Rename required metadatda columns
    req.cols <- structure(c("marker", "pvalue"), names = c(marker, pvalue))
    if (!missing(or))   req.cols <- c(req.cols, structure("or",  names = or))
    if (!missing(beta)) req.cols <- c(req.cols, structure("beta", names = beta))
    if (missing(or) & missing(beta)) 
      stop("Must supply odds ratios (or) or beta values (beta).")
      
    mcols(object) <- rename(mcols(object), req.cols)
    
    if (is.factor(object$marker)) object$marker <- as.character(object$marker)
    
    # Genome information
    if (missing(genome)) {
      if (all(is.na(GenomicRanges::genome(object))))
        stop("Must supply UCSC genome abbreviation (e.g., hg19).")
    } else {
     GenomicRanges::genome(object) <- genome 
    }
     
    # Calculate z-score
    object$zscore <- calc.z(object$pvalue, 
                            or   = if (!missing(or))   object$or,
                            beta = if (!missing(beta)) object$beta)
    
    return(new("GWAS", object))      
})

