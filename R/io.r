#' Import GWAS results
#' 
#' @param file 
#' @param snpid column label corresponding to SNP ID column
#' @param chr column label corresponding to chromosome column
#' @param pos
#' @param pval
#' @param or
#' @param genome the genome identifier or assembly name
#' @param ... extra parameters passed to read.table
#' 
#' @return GRanges object
#' @export
#' 
#' @examples
#' scz <- import.gwas("data/pgc.scz.clump.txt.gz", snpid = "snpid",
#'                    chr = "hg18chr", pos = "bp", pval = "pval",
#'                    or = "or", genome = "hg18")

import.gwas <- function(file, snpid, chr, pos, pval, or, genome, ...) {
  
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, ...)
  
  vars <- list(required = c(snpid, chr, pos, pval, or, genome))
  vars$other <- setdiff(names(df), vars$required)
  
  df[[chr]] <- sub("chr", "", df[[chr]])
  
  gr <- GRanges(seqnames = paste0("chr", df[[chr]]), 
                ranges = IRanges(df[[pos]], width = 1),
                strand = "*",
                snpid = df[[snpid]], pval = df[[pval]], or = df[[or]])
  
  gr$z <- ifelse(gr$or > 1, qnorm(gr$pval / 2), -qnorm(gr$pval / 2))
  
  mcols(gr) <- df[vars$other]
  
  gr <- fix_chrs(gr, genome = genome, add.length = TRUE)
  
  return(gr)
}