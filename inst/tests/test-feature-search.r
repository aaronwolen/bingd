context("Feature search")

test_that("AnnotationHub feature search works offline", {
  query <- list(Tier1 = c("Gm12878", "Encode"))
  
  cache.dir <- file.path(system.file("data", package = "bingd"), "resources")
  
  list.hit <- hub.search(query, online = FALSE, cache.dir = cache.dir)
  atomic.hit <- hub.search(query[[1]], online = FALSE, cache.dir = cache.dir)
  
  expect_identical(list.hit[[1]], atomic.hit[[1]])
  
  expect_equal(names(list.hit[[1]]), c("Title", "LocalPath"))
  expect_equal(names(atomic.hit[[1]]), c("Title", "LocalPath"))  
})