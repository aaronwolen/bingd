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

test_that("GWAS creation from data.frame", {
  gwas.df <- do.call("as.GWAS", c(list(test.df), test.args))  
  expect_match(class(gwas.df), "GWAS")
})

test_that("GWAS creation with different signed statistics", {
  test.args <- test.args[-6]
  test.args <- c(test.args, zscore = "test.or")
  
  gwas.df <- do.call("as.GWAS", c(list(test.df), test.args))  
  expect_match(class(gwas.df), "GWAS")
  expect_identical(zscore(gwas.df), test.or)
})

test_that("GWAS creation from data.frame  fails without required information", {   
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-1])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-2])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-3])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-4])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-5])))
  expect_error(do.call("as.GWAS", c(list(test.df), test.args[-6])))
})


# Create GRanges object
test.gr <- GenomicRanges::GRanges(test.chrs, 
                                  IRanges::IRanges(test.pos, width = 1), 
                                  strand = "*",
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
  GenomeInfoDb::genome(test.gr) <- "hg19"
  gwas.gr <- do.call("as.GWAS", c(list(test.gr), test.args[-1]))
  expect_match(class(gwas.gr), "GWAS")
})

test_that("Feature genome versions must match GWAS", {
  GenomeInfoDb::genome(gwas.gr) <- "hg18"
  expect_error(annotate.gwas(gwas.gr, feature.list = f.list))
})


context("Overlap functions")

f.list <- local.features(query, path = test.dir)
f.ranges <- lapply(unlist(LocalPath(f.list)), load.feature)

test_that("findOverlaps results are accurate", {
  
  flist.overlaps <- findOverlaps(gwas.gr, f.list)
  expect_is(flist.overlaps, "HitsList")
  
  f.overlaps <- lapply(f.ranges, function(f) findOverlaps(test.gr, f))
  
  expect_true(all(mapply(identical, flist.overlaps, f.overlaps)))
})
  
test_that("featureOverlaps results are accurate", {
  
  flist.overlaps <- featureOverlaps(gwas.gr, f.list)
  expect_is(flist.overlaps, "DataFrameList")
  
  flist.overlaps <- DataFrame(as(flist.overlaps, "list"))
  
  f.overlaps <- lapply(f.ranges, function(f) overlapsAny(test.gr, f))
  expect_true(all(mapply(identical, flist.overlaps, f.overlaps)))
})

context("Create AnnotatedGWAS object")

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
  expect_equivalent(unlist(colnames(f)),
                    unlist(lapply(LocalPath(f.list), names)))
})

test_that("Features can be accessed with fcols()", {
  
  f <- fcols(gwas.annot)
  expect_match(class(f), "DataFrame")
  
  expect_equal(nrow(f), length(gwas.annot))
  expect_equal(ncol(f), nrow(unlist(f.list)))
  
  # Feature names match
  expect_equivalent(names(f),
                    unlist(lapply(LocalPath(f.list), names)))
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



context("Test feature consolidation")

test_that("Same number of feature groups after consolidation", {
  
  gwas.cons <- consolidate(gwas.annot)
  
  expect_true(is.consolidated(gwas.cons))
  expect_false(is.consolidated(gwas.annot))
  
  # All groups are still present
  expect_equal(length(features(gwas.cons)), length(features(gwas.annot)))
  expect_equal(names(features(gwas.cons)), names(features(gwas.annot)))
  
  # One feature per group
  expect_true(all(sapply(features(gwas.cons), length) == 1))

  # Group and feature names match
  expect_match(names(features(gwas.cons)),
               sapply(features(gwas.cons),  names))
})
