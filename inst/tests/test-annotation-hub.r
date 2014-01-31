context("AnnotationHub")

query <- list(DNase   = c("GM12878", "UwDnase", "narrowPeak"),
              Histone = c("GM12878", "H3k4me3", "narrowPeak"))

features <- hub.features(query, path = test.dir, online = TRUE)

test_that("Searching for uncached features", {
  
  # Specified cache directory was created
  expect_true(file.exists(test.dir))

  # Returns a FeatureList
  expect_identical(is.FeatureList(features), features)

  # Returns expected number of matches
  expect_equal(sapply(features, nrow), 
               c("DNase" = 2, "Histone" = 2))
})



test_that("Test features are uncached", {
  
  # Files aren't cached
  expect_equal(sum(sapply(features, function(x) x$Cached)), 0)
  feature.files <- unlist(lapply(features, function(x) x$LocalPath))
  expect_false(any(file.exists(feature.files)))
})
  


# Cache test features
features <- cache.features(features, test.dir)

test_that("Test features were cached", {
  
  # Returns a FeatureList
  expect_identical(is.FeatureList(features), features)
  
  # Files are now cached
  expect_true(all(file.exists(feature.files)))
  
  # FeatureList indicates files are cached
  expect_equal(sum(sapply(features, function(x) x$Cached)), 4)
})
