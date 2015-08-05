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
setMethod("calc.enrich", c(object = "GWAS", feature.list = "AnnotationHub"), 
  function(object, feature.list, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  enrich <- parallel::mclapply(names(feature.list), function(p) 
                               serial.enrich(object %over% feature.list[[p]],
                                             stat, thresh.levels),
                               mc.cores = foreach::getDoParWorkers())
  
  setNames(enrich, names(feature.list)) %>% tidyr::unnest(col = "feature")
})



#' @rdname calc.enrich
setMethod("calc.enrich", c(object = "GWAS", feature.list = "AnnotationHubList"), 
  function(object, feature.list, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  if (is.null(names(feature.list))) 
    names(features.list) <- paste0("Group", seq_along(feature.list))
  
  enrich <- lapply(feature.list, calc.enrich, object = object, 
                   stat = stat, thresh.levels = thresh.levels)
  
  tidyr::unnest(enrich, col = "group")
})



#' @rdname calc.enrich
setMethod("calc.enrich", c(object = "AnnotatedGWAS", feature.list = "missing"),
  function(object, feature.list, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  fmat <- as.matrix(fcols(object))
  enrich <- serial.enrich(fmat, stat, thresh.levels)
  
  # add feature labels
  flist <- lapply(features(object), names)
  fnames <- rep(names(flist), elementLengths(flist))

  enrich <- dplyr::data_frame(feature = fnames, sample = unlist(flist)) %>%
    dplyr::right_join(enrich, by = "sample") 
  
  structure(enrich, class = "data.frame")
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
  
  if (is.vector(feature)) feature <- matrix(feature, ncol = 1)
  
  enrich <- list()
  enrich$count      <- t(crossprod(feature, stat.hits))
  enrich$prop       <- enrich$count / colSums(stat.hits)
  enrich$enrichment <- t(t(enrich$prop) / colMeans(feature))
  
  enrich <- lapply(enrich, as.data.frame.table) %>%
    tidyr::unnest(col = "variable") %>%
    tidyr::spread_("variable", "Freq") %>%
    dplyr::rename_("threshold" = "Var1", "sample" = "Var2") %>%
    dplyr::mutate_(threshold = ~thresh.levels[as.integer(threshold)]) %>%
    dplyr::select_("sample", "threshold", "count", "prop", "enrichment") %>%
    dplyr::mutate_(sample = ~as.character(sample)) %>%
    dplyr::arrange_("sample", "threshold")
  
  ns <- dplyr::n_distinct(enrich$sample)
  if (ns == 1) enrich <- dplyr::select_(enrich, ~-sample)
  
  return(enrich)
}


# Convert list of enrichment results to tidy data.frame
format.enrich <- function(x) {
  tidyr::unnest(lapply(x, tidyr::unnest, col = "sample"), col = "feature")
} 
