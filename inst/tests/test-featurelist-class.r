context("FeatureList class")

df <- data.frame(LocalPath = rep(test.dir, 3),
                 Title = basename(test.dir),
                 Cached = TRUE,
                 stringsAsFactors = FALSE)

test_that("From data.frame", {

  fl1 <- FeatureList(df)
  expect_match(class(fl1), "FeatureList")
  expect_match(names(fl1), "features")
  
  fl2 <- FeatureList(df, df)
  expect_equivalent(class(fl2), "FeatureList")
  expect_equivalent(names(fl2), c("features1", "features2"))
  
  fl3 <- FeatureList(label1 = df, label2 = df)
  expect_equivalent(class(fl3), "FeatureList")
  expect_equivalent(names(fl3), c("label1", "label2"))
  
  expect_identical(FeatureList(df, df, names = c("label1", "label2")), fl3)
})

test_that("From DataFrame", {

  df <- IRanges::DataFrame(df)
  
  fl1 <- FeatureList(df)
  expect_match(class(fl1), "FeatureList")
  expect_match(names(fl1), "features")
  
  fl2 <- FeatureList(df, df)
  expect_equivalent(class(fl2), "FeatureList")
  expect_equivalent(names(fl2), c("features1", "features2"))
  
  fl3 <- FeatureList(label1 = df, label2 = df)
  expect_equivalent(class(fl3), "FeatureList")
  expect_equivalent(names(fl3), c("label1", "label2"))
  
  expect_identical(FeatureList(df, df, names = c("label1", "label2")), fl3)
})
