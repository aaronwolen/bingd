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
      cat("  ..  ", sum(object[[i]]$Cached), "of which are cached\n")
    }
})


#' Download uncached features
#' 
#' @inheritParams annotate.gwas
#' @inheritParams local.features
#' 
#' @exportMethod cache.features
#' 
#' @return \code{FeatureList}

setGeneric("cache.features", 
  function(object, path) {
    standardGeneric("cache.features")
})
           
setMethod("cache.features", "FeatureList",
  function(object, path) {
    
    if (missing(path)) path <- cache.path()
    
    object <- stack(object)
    uncached.files <- basename(subset(object, !Cached)$LocalPath)
    
    hub <- AnnotationHub::AnnotationHub(hubCache = path)
    hub.files <- hub@snapshotPaths
      
    # Download uncached files
    cached.files <- character()
    for (file in uncached.files) {
      
      hub.file <- hub.files[match(file, basename(hub.files))]
      
      # Report unmatched files
      if (is.na(hub.file)) {
        warning(file, " not found on AnnotationHub.\b", call. = F)
      }
      
      # $ is the only exported method from AnnoHub 1.6 to download resources
      dl <- do.call("$", list(hub, names(hub.file)))
      
      if (length(dl) == 0) stop("Failed to download:\n\t", file, call. = FALSE)
      cached.files <- c(cached.files, file)
    }
    
    object$Cached[basename(object$LocalPath) %in% cached.files] <- TRUE

    FeatureList(split(object[-1], object$name))
})


# Accessors ---------------------------------------------------------------

#' Access LocalPath information from FeatureList object
#' @param object \code{FeatureList} object
#' @export 

setGeneric("LocalPath", function(object) standardGeneric("LocalPath"))
           
setMethod("LocalPath", "FeatureList",
  function(object) {
    Map(function(x) structure(x$LocalPath, names = x$Title), object)
})

#' Access Cached column from FeatureList object
#' @param object \code{FeatureList} object
#' @export 

setGeneric("Cached", function(object) standardGeneric("Cached"))
           
setMethod("Cached", "FeatureList",
  function(object) {
    Map(function(x) structure(x$Cached, names = x$Title), object)
})
