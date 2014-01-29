context("Clean up")

test_that("Removed test cache directory", {
  unlink(test.dir, recursive = TRUE)
  expect_false(file.exists(test.dir))
})
