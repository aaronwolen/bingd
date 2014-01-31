# Create a test cache directory
make_test <- function(dir.name) {
  if (missing(dir.name)) dir.name <- paste(cache.path(), "test", sep = "-")
  if (file.exists(dir.name)) {
    warning(dir.name, " already exists.")
    rndm.ext <- paste(sample(letters, 3), collapse = "")
    dir.name <- paste(dir.name, rndm.ext, sep = "-")
  }
  return(dir.name)
}

test.dir <- make_test()

query <- list(DNase   = c("GM12878", "UwDnase", "narrowPeak"),
              Histone = c("GM12878", "H3k4me3", "narrowPeak"))

features <- hub.features(query, path = test.dir, online = TRUE)
