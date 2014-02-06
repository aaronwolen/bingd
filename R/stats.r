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