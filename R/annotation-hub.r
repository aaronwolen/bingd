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
  
  ah.sql <- dplyr::src_sqlite(db)
  
  resources  <- dplyr::tbl(ah.sql, "resources") %>%
    dplyr::select_(.dots = c("id", "ah_id")) %>%
    dplyr::rename_(.dots = c(resource_id = "id"))
    
  rdatapaths <- dplyr::tbl(ah.sql, "rdatapaths") %>% 
    dplyr::select_(.dots = c("id", "resource_id", "rdatapath"))
  
  # filter rdatapaths for local files
  local.files <- dir(path, full.names = TRUE)
  rdatapaths  <- rdatapaths %>%
    dplyr::filter_(.dots = interp(~id %in% f, f = basename(local.files)))

  files <- resources %>%
    dplyr::inner_join(rdatapaths, by = "resource_id") %>% 
    dplyr::collect()
 
  # multiple resource_ids occur when files include an index file (*.tbi or *fai)
  # ignore these index files when listing cached files
  files <- files[!grepl("(tbi|fai)$", files$rdatapath), ]
  if (any(duplicated(files$ah_id))) stop("Bug! Duplicate ah_ids detected.")

  # match output of AnnotationHub:::.datapathIds()
  setNames(file.path(path, files$id), files$ah_id)
}

# retrieve all ids in database
db_ids <- function(db) {
  if (missing(db)) db <- ah_db()
  ids <- dplyr::src_sqlite(db) %>% 
    dplyr::tbl("resources") %>%
    dplyr::select_(.dots = "ah_id") %>%
    dplyr::collect()
  
  ids$ah_id
}

# retrieve metadata fields matching AnnotationHub::mcols
# TODO: Add rdataclass, sourceurl and sourcetype fields
db_metadata <- function(db, ids) {
  if (missing(db)) db <- ah_db()
  if (missing(ids)) ids <- db_ids(db)

  # ah_id %in% ids filtering fails if length(ids) == 1
  # convert ids to a data.frame and perform inner_joins as a workaround
  ids <- dplyr::data_frame(ah_id = ids)
  
  ah.sql <- dplyr::src_sqlite(db)
  
  mdata <- ah.sql %>% 
    dplyr::tbl("resources") %>% 
    dplyr::inner_join(ids, by = "ah_id", copy = TRUE) %>%
    dplyr::select_(.dots = c("id", .DB_RESOURCE_FIELDS)) %>%
    dplyr::collect()
  
  tags <- ah.sql %>%
    dplyr::tbl("tags") %>%
    dplyr::collect() %>%
    dplyr::filter_(.dots = interp(~ resource_id %in% x, x = mdata$id)) %>%
    dplyr::group_by_(.dots = "resource_id") %>%
    dplyr::summarize(tags = paste0(tag, collapse = ", "))
  
  mdata <- mdata %>%
    dplyr::left_join(tags, by = c("id" = "resource_id")) %>%
    dplyr::select_(.dots = "-id")
    
  DataFrame(mdata[-1], row.names = mdata[["ah_id"]])
}

