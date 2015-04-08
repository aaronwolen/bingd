context("AnnotationHub")

f.list <- hub.features(query, path = test.dir, online = TRUE)

test_that("Searching for uncached features", {
  
  # Specified cache directory was created
  expect_true(file.exists(test.dir))
  
  # Returns a FeatureList
  expect_match(class(f.list), "FeatureList")

  # Returns expected number of matches
  expect_equal(sapply(f.list, nrow), 
               c("DNase" = 2, "Histone" = 2))
})



test_that("Test features are uncached", {
  
  # FeatureList indicates files aren't cached
  expect_equal(sum(sapply(f.list, function(x) x$Cached)), 0)
  
  # Files don't exist
  feature.files <- unlist(lapply(f.list, function(x) x$LocalPath))
  expect_true(all(grepl(test.dir, feature.files)))
  expect_false(any(file.exists(feature.files)))
})
  


# Cache test features
suppressMessages(f.list <- cache.features(f.list, test.dir))

test_that("Test features were cached", {
  
  # Returns a FeatureList
  expect_true(.validFeatureList(f.list))
  
  # FeatureList indicates files are cached
  expect_equal(sum(IRanges::stack(f.list)$Cached), 4)
  
  # Files do exist
  feature.files <- IRanges::stack(f.list)$LocalPath
  expect_true(all(grepl(test.dir, feature.files)))
  expect_true(all(file.exists(feature.files)))
})
