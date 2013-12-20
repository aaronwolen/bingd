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
  thresh.levels <- unique(c(0, thresh.levels))
  thresh.levels <- matrix(rep(thresh.levels, each = n), nrow = n)

  hits.mat <- feature & (stat >= thresh.levels)
  
  out <- data.frame(threshold = thresh.levels[1,],
                        count = colSums(hits.mat),
                         prop = colMeans(hits.mat))
  out$enrichment <- out$prop / mean(feature)
  
  return(out)
}
