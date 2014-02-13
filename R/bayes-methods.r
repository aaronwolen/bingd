#' Calculate prior probability of being an effective SNP
#'
#' @param object \code{AnnotatedGWAS} object
#' @param effect assumed effect size, used to shift the ideal z-score distribution
#' @param trace Non-negative integer. If positive, tracing information on the
#'  progress of the optimization is produced. See \code{\link[stats]{optim}} for
#'  details.
#'  
#'  @return A list with components:
#' \describe{
#'  \item{\code{p.e}}{P(E): prior probability of being an effective SNP}
#'  \item{\code{p.n}}{P(N): prior probability of a SNP having no effect}
#'  \item{\code{lambda}}{GWAS inflation factor (\eqn{\lambda})}
#' }

setGeneric("calc.pe", function(object, effect = 2, trace = 0) {
  standardGeneric("calc.pe")
})

setMethod("calc.pe", "AnnotatedGWAS", function(object, effect = 2, trace = 0) {
    
  # Ideal densities at the same points as observed z-scores
  z.ideal <- seq(-50, 50) / 10
  fz.n <- dnorm(z.ideal)
  fz.e <- (dnorm(z.ideal + effect) + dnorm(z.ideal - effect)) / 2
  
  # Density of observed z-scores
  z.obs <- density(zscore(object), 
                   bw = 0.1, from = -5, to = 5, n = length(z.ideal))
  
  # Calculate residual SS
  calc.resids.sq <- function(pars, z.vals , obs.z.density, fz.e) {
    results <- (1 - pars[2]) * dnorm(z.vals / pars[1]) / pars[1] + pars[2] * fz.e
    sum((results - obs.z.density)^2)
  }

  roots <- optim(par = c(1.05, .005), 
                 fn = calc.resids.sq, 
                 z.vals = z.ideal,
                 obs.z.density = z.obs$y,
                 fz.e = fz.e,
                 lower = c(1, 0), 
                 upper = c(1.2, .1), 
                 method = 'L-BFGS-B', 
                 control = list(trace = trace))

  return(list(p.e = roots$par[2], 
              p.n = 1 - roots$par[2], 
              lambda = roots$par[1]))  
})
