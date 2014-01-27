library(GenomicRanges)

setClass("GWAS", 
         representation(markers = "GRanges",
                        snpid = "character",
                        pval = "numeric",
                        or = "numeric"),
         contains = "GRanges")
         

gwas.df <- read.table("data/pgc.scz.clump.txt.gz", 
                      header = T, stringsAsFactors = FALSE)

gwas.gr <- with(gwas.df, GRanges(paste0("chr", chr), 
                                 IRanges(pos, width = 1),
                                 strand = "*"))

gwas <- new("GWAS", 
            markers = gwas.gr, 
            snpid = gwas.df$snpid,
            pval = gwas.df$pval,
            or = gwas.df$or)

g

setGeneric("filterFeatures", 
  function(file.list, query) {
    standardGeneric("filterFeatures")
})


setMethod("makeGWAS", "GRanges",
  function(object, snpid, chr, pos, pval, or, genome) {
          
})

new.gwas <- function(data, snpid, chr, pos, pval, or, genome, ...) {
  
}

gwas.gr <- gwas
mcols(gwas) <- mcols(gwas)[, "z"]
new("GWAS", snpid = gwas$)


#' @param snpid column label corresponding to SNP ID column
#' @param chr column label corresponding to chromosome column
#' @param pos
#' @param pval
#' @param or