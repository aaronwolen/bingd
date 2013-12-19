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
  
  # Retrieve features and check for gwas overlaps
  overlaps <- lapply(unlist(filter.hits), function(f)
                     gwas %over% AnnotationHub:::.getResource(hub, f))
  
  names(overlaps) <- md$Description[match(unlist(filter.hits), make.names(md$RDataPath))]
  
  # It'd be nice if overlap results for each feature type could
  # be stored in different slots but for now all features are just
  # added as ordinary variables
  overlaps <- DataFrame(overlaps)
  # Features will be denoted by a .prefix
  names(overlaps) <- paste0(".", names(overlaps))
  
  mcols(gwas) <- DataFrame(mcols(gwas), overlaps)
  
  return(gwas)
}