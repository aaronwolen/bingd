context("GWAS classes")

test.chrs <- rep(paste0("chr", c(1:22, "X")), each = 100)
n.markers <- length(test.chrs)

test.snps <- paste0("rs", seq_len(n.markers))
test.pos  <- round(runif(test.chrs, 1, 100e6), 0)
test.pval <- runif(n.markers)
test.or   <- rnorm(n.markers)

test_that("GWAS object from data.frame", {
 
  test.df <- data.frame(test.snps, test.chrs, test.pos, test.pval, test.or)
  
  gwas.df <- as.GWAS(test.df,
                     marker = "test.snps", chr = "test.chrs", bp = "test.pos", 
                     pvalue = "test.pval", or = "test.or")
  
  expect_match(class(gwas.df), "GWAS")
})



# Create GRanges object
test.gr <- GRanges(test.chrs, IRanges(test.pos, width = 1), strand = "*",
                   test.snps, test.pval, test.or)

test_that("GWAS creation fails without required information", {
  expect_error(as.GWAS(test.gr, marker = "test.snps", pvalue = "test.pval"))
  expect_error(as.GWAS(test.gr, marker = "test.snps", or = "test.or"))
  expect_error(as.GWAS(test.gr, pvalue = "test.pval", or = "test.or"))
})



# Create GWAS object
gwas.gr <- as.GWAS(test.gr, marker = "test.snps", 
                   pvalue = "test.pval", or = "test.or")

test_that("GWAS object from GRanges", {
  expect_match(class(gwas.gr), "GWAS")
})


# Create annotated GWAS object
features <- local.features(query, path = test.dir)
gwas.annot <- annotate.gwas(gwas.gr, feature.list = features)

test_that("GWAS object annotation", {
  
  expect_true(is.annotated(gwas.annot))
  
  # Pull annotated GWAS features
  f.df <- pull.features(gwas.annot)
  expect_match("DataFrame", sapply(f.df, class))
  
  # One column exists for every cached feature
  expect_equivalent(sapply(features, nrow), sapply(f.df, length))
})
