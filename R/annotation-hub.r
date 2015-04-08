# Helper functions for interfacing with AnnotationHub

cache.path <- function() AnnotationHub:::hubCache()

cache.exists <- function(path) AnnotationHub:::.hubCacheExists(path)

cache.create <- function(path) AnnotationHub:::.dirCreate(path)

# appends resources to cache path
cache.resources <- function(path) AnnotationHub:::.cacheResource(path)