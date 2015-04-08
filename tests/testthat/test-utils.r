context("Test group labeling")

df2 <- expand.grid(c(TRUE, FALSE), c(TRUE, FALSE))
df3 <- expand.grid(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE))

test_that("Combinations of two groups are labeled accurately", {
  expect_equal(label.groups(df2), 
               c("Var1.and.Var2", "Var2", "Var1", "neither"))
})

test_that("Combinations of three groups are labeled accurately", {
  expect_equal(label.groups(df3), 
               c("Var1.and.Var2.and.Var3",
                 "Var2.and.Var3", "Var1.and.Var3",
                 "Var3", "Var1.and.Var2", "Var2",
                 "Var1", "neither" ))
})
