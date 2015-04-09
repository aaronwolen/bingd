#' Convert p-values to z-scores
#' 
#' @param p a vector of p-values
#' @param or a vector of odds ratios
#' @param beta a vector of beta values
#' 
#' @description Must provide either odds ratios OR beta values.
#' @return numeric vector of z-scores

calc.z <- function(p, or = NULL, beta = NULL) {
  
  z <- qnorm(p / 2)
  
  if (!is.null(or))
    z[or < 1]   <- -z[or < 1]
  else if (!is.null(beta)) 
    z[beta < 0] <- -z[beta < 0]
  else 
    stop("Must provide odds ratios (or) or beta values (beta).\n")
  
  return(z)
}

#' Calculate p-value threshold for identifying risk SNPs
#' 
#' We select at least 100 SNPs with FDR < 0.05 or 1000 SNPs with FDR < 0.1. If 
#' this can't be done then we use 100 SNPs but warn the estimate is ureliable.
#' 
#' @inheritParams stats::p.adjust
#' @inheritParams calc.priors
#' 

calc.threshold <- function(p, verbose = FALSE, method = "BH") {

  q <- stats::p.adjust(p, method)
 
  if (sum(q < 0.05) >= 100) {
    thresh <- max(p[q < 0.05])
    warn <- FALSE
  }  

  if (sum(q < 0.1)  >= 1000) {
    thresh <- max(p[q < 0.1 ])
    warn <- TRUE
  }

  if (!exists("thresh")) {
    thresh <- sort(p)[100]  
    warn <- TRUE
  }
  
  if (verbose) report(format(thresh, scientific = TRUE), "Risk p-value threshold")
  if (warn) warning("p-value threshold may be unreliable.\n", call. = FALSE, immediate. = TRUE)
  
  return(thresh)
}

