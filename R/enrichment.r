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
  out <- subset(out, feature == TRUE, select = -feature)
  out <- rename(out, c(L1 = "feature", L2 = "sample"))
  levels
  out$threshold <- factor(out$threshold, 
                          sort(as.numeric(levels(out$threshold))))
  return(out)
}


#' Calculate the enrichment of a particular discrete feature across thresholds
#'
#' @param feature logical or character vector
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#' 
#' @return data.frame containing `enrichment`
#' @importFrom plyr ddply
#' @export

serial.enrich <- function(feature, stat, thresh.levels) {
  
  # Perform counts across all thresholds
  base.counts <- data.frame(table(feature))
  names(base.counts) <- c("feature", "count")
  base.counts$prop <- base.counts$count / sum(base.counts$count)

  thresh.counts <- lapply(thresh.levels, function(t) table(feature[stat >= t]))
  names(thresh.counts) <- thresh.levels
  thresh.counts <- thresh.counts[sapply(thresh.counts, dim) > 0]
  
  thresh.counts <- lapply(thresh.counts, data.frame)

  # Re-format output
  out <- data.frame(threshold = rep(names(thresh.counts), 
                                    sapply(thresh.counts, nrow)))
  out <- data.frame(out, do.call("rbind", thresh.counts))
  names(out) <- c("threshold", "feature", "count")
  
  # Proportion of feature at each threshold
  out <- ddply(out, .(threshold), transform, 
               prop = count / sum(count, na.rm = TRUE))
  
  base.index <- match(out$feature, base.counts$feature)
  out$enrichment <- out$prop / base.counts$prop[base.index]
  
  return(out)
}
