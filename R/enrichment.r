#' Calculate enrichment of all annotated features across thresholds
#' 
#' @param gwas GWAS GRanges object
#' @param stat numeric vector containing statistic to which thresholds will 
#'        be applied for each enrichment calculation
#' @param data.filter optional, calculate enrichment among features that match data.filter to avoid annotating the gwas object, which can require lots of memory        
#' 
#' @importFrom reshape2 melt
#' @importFrom plyr rename
#' @export
#'         

calc.enrich <- function(object, stat, thresh.levels, data.filter, online = NULL) {

  if (missing(thresh.levels)) {
    thresh.levels <- quantile(stat, seq(0, 1, 0.1))
  }
  
  if (!missing(data.filter)) {
    
    features <- feature.search(data.filter, genome(object)[1], online = online)
    
    out <- lapply(features, function(f)
              mclapply(f, function(x)
                       serial.enrich(object %over% load.feature(x),
                       stat, thresh.levels),  
                       mc.cores = getDoParWorkers()))
    
    out <- mapply(function(x, y) {
                names(x) <- y; x
              }, out, features, SIMPLIFY = FALSE)
    
  } else {
    features <- features(object)
    
    out <- lapply(features, function(f) 
                  mclapply(f, serial.enrich, stat, thresh.levels,  
                           mc.cores = getDoParWorkers()))

  }
  
  out <- melt(out, measure.vars = NULL)
  out <- rename(out, c(L1 = "feature", L2 = "sample"))
  out$threshold <- factor(out$threshold)
  
  # Rename features if data.filter was used (this should be a separate function)
  if (!missing(data.filter)) {
    md <- get.metadata(online = FALSE)
    md.index <- sapply(out$sample, grep, make.names(md$RDataPath))
    out$sample <- md$Title[md.index]
  }
  
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
