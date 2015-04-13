#' Calculate enrichment of features across thresholds
#' 
#' @param object \code{\link{GWAS}} object
#' @inheritParams annotate.gwas  
#' @inheritParams serial.enrich
#'
#' @description
#' \code{gwas} must be annotated with \code{\link{annotate.gwas}} OR a
#' \code{feature.list} must be provided.
#' 
#' @exportMethod calc.enrich       

setGeneric("calc.enrich", 
 function(object, feature.list, stat, thresh.levels) {
   standardGeneric("calc.enrich")
}) 

#' @rdname calc.enrich
setMethod("calc.enrich", c(object = "GWAS", feature.list = "FeatureList"), 
  function(object, feature.list, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  features <- LocalPath(feature.list)
  enrich <- list()
  
  for (i in names(features)) {
    
    enrich[[i]] <- parallel::mclapply(features[[i]], function(p) 
                               serial.enrich(object %over% load.feature(p),
                                             stat, thresh.levels),
                               mc.cores = foreach::getDoParWorkers())
  }
  
  return(format.enrich(enrich))
})


#' @rdname calc.enrich
setMethod("calc.enrich", c(object = "AnnotatedGWAS", feature.list = "missing"),
  function(object, feature.list, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  features <- features(object)
  enrich <- list()    
  
  for (i in names(features)) {
    
    enrich[[i]] <- parallel::mclapply(features[[i]], function(f) 
                               serial.enrich(f, stat, thresh.levels),
                               mc.cores = foreach::getDoParWorkers())
  }
  
  return(format.enrich(enrich))
})


#' Calculate the enrichment of a particular discrete feature across thresholds
#'
#' @param feature logical or character vector
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#' @param thresh.levels a numeric vector containing the various thresholds at
#' which \code{feature} enrichment is calculated
#' 
#' @return data.frame containing `enrichment`
#' @export

serial.enrich <- function(feature, stat, thresh.levels) {
  
  n <- length(stat)
  
  stat.hits <- vapply(thresh.levels, function(t) stat >= t, FUN.VALUE = logical(n))
  
  hits <- matrix(feature, nrow = 1) %*% stat.hits
  
  dplyr::data_frame(threshold = thresh.levels, count = hits[1,]) %>%
    dplyr::mutate_(prop = ~count / colSums(stat.hits)) %>%
    dplyr::mutate_(enrichment = ~prop / mean(feature))
}


# Convert list of enrichment results to tidy data.frame
format.enrich <- function(x) {
  tidyr::unnest(lapply(x, tidyr::unnest, col = "sample"), col = "feature")
} 
