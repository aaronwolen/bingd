context("GWAS classes")

test.chrs <- rep(paste0("chr", c(1:22, "X")), each = 100)
n.markers <- length(test.chrs)

test.snps <- paste0("rs", seq_len(n.markers))
test.pos  <- round(runif(test.chrs, 1, 100e6), 0)
test.pval <- runif(n.markers)
test.or   <- rnorm(n.markers)

test.df <- data.frame(test.snps, test.chrs, test.pos, test.pval, test.or)

test.args <- c(genome = "hg19", 
               marker = "test.snps", chr = "test.chrs", bp = "test.pos", 
               pvalue = "test.pval", or = "test.or")

test_that("GWAS creation from GRanges fails without required information", {   
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-1])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-2])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-3])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-4])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-5])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-6])))
})

test_that("GWAS creation from data.frame", {
  gwas.df <- do.call("as.GWAS", c(list(test.df), test.args))  
  expect_match(class(gwas.df), "GWAS")
})


# Create GRanges object
test.gr <- GRanges(test.chrs, IRanges(test.pos, width = 1), strand = "*",
                   test.snps, test.pval, test.or)

test.args <- test.args[!names(test.args) %in% c("chr", "bp")]

test_that("GWAS creation from GRanges fails without required information", {
  expect_error(do.call("as.GWAS", c(list(test.gr), test.args[-1])))
  expect_error(do.call("as.GWAS", c(list(test.gr), test.args[-2])))
  expect_error(do.call("as.GWAS", c(list(test.gr), test.args[-3])))
  expect_error(do.call("as.GWAS", c(list(test.gr), test.args[-4])))
})


# Create GWAS object
gwas.gr <- do.call("as.GWAS", c(list(test.gr), test.args))

test_that("GWAS object from GRanges", {
  expect_match(class(gwas.gr), "GWAS")
})

test_that("GWAS object using annotated genome", {
  genome(test.gr) <- "hg19"
  gwas.gr <- do.call("as.GWAS", c(list(test.gr), test.args[-1]))
  expect_match(class(gwas.gr), "GWAS")
})



context("Create annotated GWAS object")

features <- local.features(query, path = test.dir)
gwas.annot <- annotate.gwas(gwas.gr, feature.list = features)

test_that("GWAS object annotation", {
  
  expect_match("AnnotatedGWAS", class(gwas.annot))
  
  # Access annotated GWAS features
  f.df <- features(gwas.annot)
  expect_match("DataFrame", sapply(f.df, class))
  
  # One column exists for every cached feature
  expect_equivalent(sapply(features, nrow), sapply(f.df, length))
})

test_that("Feature genome versions must match GWAS", {
  genome(gwas.gr) <- "hg18"
  expect_error(annotate.gwas(gwas.gr, feature.list = features))
})



context("Ennrichment analysis of GWAS objects")

log.pvals <- -log10(gwas.gr$pvalue)
thresh.levels <- seq(1, ceiling(max(log.pvals)), 1)

# Enrichment of unannotated GWAS object
enrich <- calc.enrich(gwas.gr, feature.list = features, 
                      stat = log.pvals, thresh.levels = thresh.levels)

test_that("Enrichment results are valid", {
  
  # Results contain expected number of rows
  enrich.rows <- sum(sapply(features, nrow)) * length(thresh.levels)
  expect_identical(nrow(enrich), enrich.rows)
  
  # feature column matches FeatureList names
  expect_match(unique(enrich$feature), names(features))
  
  # sample column matches FeatureList Titles
  expect_match(unique(enrich$sample), sapply(features, function(x) x$Title))
})


test_that("Same results for annotated and unannotated GWAS objects", {  
  enrich.annot <- calc.enrich(gwas.annot,
                            stat = log.pvals, thresh.levels = thresh.levels)
 expect_identical(enrich, enrich.annot)
})


test_that("Feature genome versions must match GWAS", {
  genome(gwas.gr) <- "hg18"
  expect_error(calc.enrich(gwas.gr, feature.list = features, 
                      stat = log.pvals, thresh.levels = thresh.levels))
})

