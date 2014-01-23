#' Calculate enrichment of all annotated features across thresholds
#' 
#' @param gwas GWAS \code{\link{GRanges}} object
#' @inheritParams annotate.gwas  
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#'
#' @description
#' \code{gwas} must be annotated with \code{\link{annotate.gwas}} OR a
#' \code{feature.list} must be provided.
#' 
#' @importFrom reshape2 melt
#' @export
#'         

calc.enrich <- function(gwas, feature.list, stat, thresh.levels) {

  annotated <- ifelse(is.annotated(gwas), TRUE, FALSE)
  
  if (!missing(feature.list)) {
    feature.list <- is.FeatureList(feature.list)
  } else {
    if (!annotated & missing(feature.list)) {
      stop("gwas is unannotated and no FeatureList was provided.", call. = FALSE)
    }
  }

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  if (annotated) {
    features <- pull.features(gwas.annot)
    labels <- lapply(features, names)
    overlap <- function(x) x
  } else {
    features <- lapply(feature.list, function(x) x$LocalPath)
    labels <- lapply(feature.list, function(x) x$Title)
    overlap <- function(x) gwas %over% load.feature(x)
  }
  
  enrich <- lapply(features, function(group)
                   mclapply(group, function(f) 
                              serial.enrich(overlap(f), stat, thresh.levels),
                            mc.cores = getDoParWorkers()))
      
  # Label with feature titles
  enrich <- mapply(function(e, l) {
                     names(e) <- l; e
                   }, enrich, labels, SIMPLIFY = FALSE)

  enrich <- melt(enrich, measure.vars = NULL)
  
  # Temp replacement for rename
  names(enrich)[which(names(enrich) == "L1")] <- "feature"
  names(enrich)[which(names(enrich) == "L2")] <- "sample"

  enrich$threshold <- factor(enrich$threshold)

  return(enrich)
}


#' Calculate the enrichment of a particular discrete feature across thresholds
#'
#' @param feature logical or character vector
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#' 
#' @return data.frame containing `enrichment`
#' @export

serial.enrich <- function(feature, stat, thresh.levels) {
  
  n <- length(stat)
  
  stat.hits <- vapply(thresh.levels, function(t) stat >= t, FUN.VALUE = logical(n))
  
  hits.mat <- feature & stat.hits
  
  out <- data.frame(threshold = thresh.levels, count = colSums(hits.mat))
  out$prop <- out$count / colSums(stat.hits)
  out$enrichment <- out$prop / mean(feature)
  
  return(out)
}
