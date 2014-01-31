context("Feature search")


test_that("Query format for local search", {
  
  # Query can be vector or a list
  loc.vec <- local.features(query["DNase"],   path = test.dir)
  loc.lst <- local.features(query[["DNase"]], path = test.dir)
  
  expect_identical(loc.vec[[1]], loc.lst[[1]])
})

test_that("Query format for online search", {
  
  # Query can be vector or a list
  hub.vec <- hub.features(query["DNase"],   path = test.dir, online = FALSE)
  hub.lst <- hub.features(query[["DNase"]], path = test.dir, online = FALSE)
  
  expect_identical(hub.vec[[1]], hub.lst[[1]])
})



test_that("Online AnnotationHub feature search", {
  
  # Modify query to be more general
  query <- lapply(query, Filter, f = function(x) x != "narrowPeak")
  
  # Online search identifies previously cached features
  offline <- hub.features(query, path = test.dir, online = FALSE)
  online  <- hub.features(query, path = test.dir, online = TRUE)
  
  online.cached <- lapply(online, function(x) x[x$Cached, names(offline[[1]])])
  expect_equivalent(online.cached, offline)
  
  # Offline features are a subset of online features
  offline.hits <- unlist(lapply(offline, function(x) x$LocalPath))
  online.hits  <- unlist(lapply(online,  function(x) x$LocalPath))
  
  expect_true(all(offline.hits %in% online.hits))
})
