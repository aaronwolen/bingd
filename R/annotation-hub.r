# Unexported helper functions copied from AnnotationHub

cache.path <- function() AnnotationHub::hubCache()

cache.exists <- function(path) 
  !is.na(path) && isTRUE(file.info(path)$isdir)

cache.create <- function(path) {
    if (!file.exists(path))
        if (!dir.create(path, showWarnings=FALSE, recursive=TRUE)) {
            warning(gettextf("unable to create %s", sQuote(path)), domain = NA)
        } else {
          AnnotationHub::setHubOption("CACHE", path)
        }
    else if (!file.info(path)$isdir)
        stop("path exists but is not a directory: ", sQuote(path))
    path
}

# path of AnnotationHub database
ah_db <- function(path) {
  if (missing(path)) path <- cache.path()
  db <- dir(path, "annotationhub.sqlite3", full.names = TRUE)
  if (length(db) == 0) stop("AnnotationHub database not found in ", path)
  db
}
  
# list of feature files available in cache directory
cached_files <- function(path) {
  if (missing(path)) path <- cache.path()
  db <- ah_db(path)
  
  ids <- dplyr::src_sqlite(db) %>% 
    dplyr::tbl("resources") %>% 
    dplyr::select_(.dots = "id") %>%
    dplyr::collect()
  
  files <- setdiff(dir(path, full.names = TRUE), db)
  files <- files[basename(files) %in% as.character(ids$id)]
  setNames(files, paste0("AH", basename(files)))
}