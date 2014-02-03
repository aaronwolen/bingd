setMethod("show", "FeatureList", 
  function(object) {
    len <- length(object)
    cat("FeatureList object with", len,
      ifelse(len == 1, "group", "groups"),
      "of features.\n\n")
    
    for (i in seq_len(len)) {
      cat("Group ", i, " - ", names(object)[i], 
          " (", ncol(object[[i]]), " columns):\n", sep = "")
      cat("  ..  ", nrow(object[[i]]), "features\n")
      cat("  ..  ", sum(object[[i]]$cached), "of which are cached\n")
    }
})


#' Download uncached features
#' 
#' @inheritParams annotate.gwas
#' @inheritParams local.features
#' 
#' @importFrom AnnotationHub AnnotationHub
#' @export
#' 
#' @return \code{FeatureList}

cache.features <- function(feature.list, path) {
  feature.list <- is.FeatureList(feature.list)
  if (missing(path)) path <- cache.path()
  
  is.uncached <- function(x) x[!x$Cached, "LocalPath"]
  uncached <- lapply(feature.list, is.uncached)
  
  uncached.files <- basename(unlist(uncached))
  
  hub <- AnnotationHub(hubCache = path)
  hub.files <- hub@snapshotPaths
    
  # Download uncached files
  cached.files <- character()
  for (file in uncached.files) {
    
    hub.file <- hub.files[match(file, basename(hub.files))]
    
    # Report unmatched files
    if (is.na(hub.file)) {
      warning(file, " not found on AnnotationHub.\b", call. = F)
    }
    
    dl <- AnnotationHub:::.downloadFile(hub, hub.file)
    
    if (dl == 0) {
      cached.files <- c(cached.files, file)
    } else {
      stop("Failed to download:\n\t", file, call. = FALSE)
    }
  }
  
  feature.list <- Map(function(f) {
    f$Cached[basename(f$LocalPath) %in% cached.files] <- TRUE; f
  }, feature.list)
  
  return(as.FeatureList(feature.list))
}
