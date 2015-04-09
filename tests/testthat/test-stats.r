context("Test stats calculations")

p.null <- runif(1000)
p.alt  <- c(runif(800), runif(200, max = 0.01))

test_that("Warn about unreliable p-value threshold", {
  expect_warning(calc.threshold(p.null))
  expect_is(calc.threshold(p.alt), "numeric")
})
