#' Annotate GWAS object using AnnotationHub
#' 
#' @param gwas GWAS GRanges object
#' @param filter named list of vectors
#' 
#' @importFrom AnnotationHub AnnotationHub metadata
#' @export
#' 
#' @examples
#' scz <- annotate.gwas(scz, data.filter = list(DNaseI = c("UwDnaseN", "broadPeak")))

annotate.gwas <- function(gwas, data.filter) {

  message("Creating AnnotationHub object...")
  hub <- AnnotationHub()
  md <- metadata(hub)
  
  # Filter based on genome  
  md <- md[md$Genome %in% genome(gwas),]
  
  # Apply user specified filters
  mgrep <- function(pattern, x, ignore.case = TRUE, ...) {
    hits <- sapply(pattern, grepl, x = x, ignore.case = ignore.case, ...)
    which(rowSums(hits) == length(pattern))
  }
  
  filter.hits <- lapply(data.filter, mgrep, x = md$RDataPath)
  filter.hits <- lapply(filter.hits, function(x) make.names(md$RDataPath[x]))
  
  # Download/retrieve features
  get_features <- function(x, names) {
    out <- lapply(names, function(n) AnnotationHub:::.getResource(x, n))
    names(out) <- names
    return(out)
  }
  
  message("Retrieving features...")
  features <- lapply(filter.hits, get_features, x = hub)
  
  feature.labels <- lapply(features, function(x)
                          md$Description[match(names(x), 
                                               make.names(md$RDataPath))])

  features <- mapply(function(f, l) {
                       names(f) <- l; f
                     }, features, feature.labels, SIMPLIFY = FALSE)
  
  # Annotate GWAS data
  message("Annotating SNPs...")
  overlaps <- lapply(features, function(f)
                      sapply(f, overlapsAny, query = gwas))

  # It'd be nice if overlap results for each feature type could
  # be stored in different slots but for now all features are just
  # added as ordinary variables
  overlaps <- do.call("DataFrame", overlaps)
  
  # Features will be denoted by a .prefix
  names(overlaps) <- paste0(".", names(overlaps))
  
  mcols(gwas) <- DataFrame(mcols(gwas), overlaps)
  
  return(gwas)
}