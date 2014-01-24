context("Feature search")

query <- list(Tier1 = c("Gm12878", "Encode"))
cache.dir <- system.file("data/resources", package = "bingd")
required.cols <- c("Title", "LocalPath")

test_that("Offline AnnotationHub feature search", {
  
  list.hit <- hub.features(query, path = cache.dir, online = FALSE)
  atomic.hit <- hub.features(query[[1]], path = cache.dir, online = FALSE)
  
  expect_identical(list.hit[[1]], atomic.hit[[1]])
  
  expect_equal(names(list.hit[[1]]), required.cols)
  expect_equal(names(atomic.hit[[1]]), required.cols)  
})


test_that("Online AnnotationHub feature search", {
  
  online.hit <- hub.features(query, path = cache.dir, online = TRUE)
  offline.hit <- hub.features(query, path = cache.dir, online = FALSE)
  
  expect_true(all(required.cols %in% names(online.hit[[1]])))
  expect_equivalent(online.hit[[1]][, required.cols], offline.hit[[1]])
})
