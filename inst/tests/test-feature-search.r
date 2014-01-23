context("Feature search")

query <- list(Tier1 = c("Gm12878", "Encode"))
cache.dir <- system.file("data", package = "bingd")
required.cols <- c("Title", "LocalPath")

test_that("AnnotationHub feature search works offline", {
  
  list.hit <- hub.search(query, online = FALSE, cache.dir = cache.dir)
  atomic.hit <- hub.search(query[[1]], online = FALSE, cache.dir = cache.dir)
  
  expect_identical(list.hit[[1]], atomic.hit[[1]])
  
  expect_equal(names(list.hit[[1]]), required.cols)
  expect_equal(names(atomic.hit[[1]]), required.cols)  
})


test_that("AnnotationHub feature search works online", {
  
  online.hit <- hub.search(query, online = TRUE, cache.dir = cache.dir)
  
  expect_true(all(required.cols %in% names(online.hit[[1]])))
  expect_equal(nrow(online.hit[[1]]), 2)  
})
