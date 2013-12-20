#' Calculate enrichment of all annotated features across thresholds
#' 
#' @param gwas GWAS GRanges object
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#' 
#' @importFrom reshape2 melt
#' @importFrom plyr rename
#' @export
#'         

calc.enrich <- function(object, stat, thresh.levels) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  out <- lapply(features(object), function(f) 
                 lapply(f, serial.enrich, stat, thresh.levels))

  out <- melt(out, measure.vars = NULL)
  out <- rename(out, c(L1 = "feature", L2 = "sample"))
  out$threshold <- factor(out$threshold)
  return(out)
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
