context("GWAS classes")

n.markers <- 100
test.snps <- paste0("rs", seq_len(n.markers))
test.chrs <- rep(paste0("chr", c(1:22, "X")), each = n.markers)
test.pos  <- round(runif(n.markers, 1, 100e6), 0)
test.pval <- runif(n.markers)
test.or   <- rnorm(n.markers)

test_that("GWAS object from data.frame", {
 
  test.df <- data.frame(test.snps, test.chrs, test.pos, test.pval, test.or)
  
  test.gwas <- as.GWAS(test.df,
                       marker = "test.snps", chr = "test.chrs", bp = "test.pos", 
                       pvalue = "test.pval", or = "test.or")
  
  expect_match(class(gwas), "GWAS")
})



test_that("GWAS object from GRanges", {
 
  test.gr <- GRanges(test.snps, IRanges(test.pos, width = 1), strand = "*",
                     test.snps, test.pval, test.or)
  
  gwas <- as.GWAS(test.gr,
                  marker = "test.snps", pvalue = "test.pval", or = "test.or")
  expect_match(class(gwas), "GWAS")
})
.