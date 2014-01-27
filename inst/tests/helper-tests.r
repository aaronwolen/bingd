# Create a test cache directory
test.dir <- paste(cache.path(), "test", sep = "-")

query <- list(DNase   = c("GM12878", "UwDnase", "narrowPeak"),
              Histone = c("GM12878", "H3k4me3", "narrowPeak"))

features <- hub.features(query, path = test.dir, online = TRUE)
