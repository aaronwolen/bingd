
# Define ------------------------------------------------------------------

setClass("FeatureList", contains="DataFrame")


# Validate ----------------------------------------------------------------

.validFeatureList <- function(object) {
  .validColClass(object, "FeatureList", "Title",     "character")
  .validColClass(object, "FeatureList", "LocalPath", "character")
  .validColClass(object, "FeatureList", "Cached",    "logical")
  return(TRUE)
}

setValidity("FeatureList", .validFeatureList)
