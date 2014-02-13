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



context("Create AnnotatedGWAS object")

f.list <- local.features(query, path = test.dir)
gwas.annot <- annotate.gwas(gwas.gr, feature.list = f.list)

test_that("GWAS object annotation", {  
  expect_match("AnnotatedGWAS", class(gwas.annot))
  
  # Annotating added the expected number of columns
  expect_equal(ncol(mcols(gwas.annot)),
               ncol(mcols(gwas.gr)) + sum(sapply(f.list, nrow)))
})

test_that("Features can be accessed with features()", {
  
  f <- features(gwas.annot)
  expect_match(class(f), "DataFrameList")
  
  expect_true(all(nrow(f) == length(gwas.annot)))
  expect_equal(length(f), length(f.list))
  # Group names match
  expect_equal(names(f), names(f.list))
  # Feature names match
  expect_match(unlist(colnames(f)),
               unlist(lapply(LocalPath(f.list), names)))
})

test_that("Features can be accessed with fcols()", {
  
  f <- fcols(gwas.annot)
  expect_match(class(f), "DataFrame")
  
  expect_equal(nrow(f), length(gwas.annot))
  expect_equal(ncol(f), nrow(unlist(f.list)))
  
  # Feature names match
  expect_match(names(f),
               unlist(lapply(LocalPath(f.list), names)))
})



test_that("Feature genome versions must match GWAS", {
  genome(gwas.gr) <- "hg18"
  expect_error(annotate.gwas(gwas.gr, feature.list = f.list))
})






context("Ennrichment analysis of GWAS objects")

log.pvals <- -log10(gwas.gr$pvalue)
thresh.levels <- seq(1, ceiling(max(log.pvals)), 1)

# Enrichment of unannotated GWAS object
enrich <- calc.enrich(gwas.gr, feature.list = f.list, 
                      stat = log.pvals, thresh.levels = thresh.levels)

test_that("Enrichment results are valid", {
  
  # Results contain expected number of rows
  enrich.rows <- sum(sapply(f.list, nrow)) * length(thresh.levels)
  expect_identical(nrow(enrich), enrich.rows)
  
  # feature column matches FeatureList names
  expect_match(unique(enrich$feature), names(f.list))
  
  # sample column matches FeatureList Titles
  expect_match(unique(enrich$sample), sapply(f.list, function(x) x$Title))
})


test_that("Same results for annotated and unannotated GWAS objects", {  
  enrich.annot <- calc.enrich(gwas.annot,
                            stat = log.pvals, thresh.levels = thresh.levels)
 expect_identical(enrich, enrich.annot)
})


test_that("Feature genome versions must match GWAS", {
  genome(gwas.gr) <- "hg18"
  expect_error(calc.enrich(gwas.gr, feature.list = f.list, 
                      stat = log.pvals, thresh.levels = thresh.levels))
})




context("Consolidation")

gwas.cons <- consolidate(gwas.annot)

test_that("Same number of feature groups after consolidation", {
  features(gwas.cons)
})

