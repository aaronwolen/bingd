#' Calculate posterior probability of being an effective SNP
#' 
#' Integrate GWAS results with annotated genomic features
#' 
#' @param object \code{AnnotatedGWAS} object
#' @inheritParams calc.pr
#' @inheritParams calc.conditionals
#' 
#' @importFrom plyr ddply
#' @export

setGeneric("calc.bayes", 
  function(object, risk.thresh, adjust, effect = 2, verbose = FALSE, trace = 0) {
  standardGeneric("calc.bayes")
})

setMethod("calc.bayes", "AnnotatedGWAS", 
  function(object, risk.thresh, adjust, effect = 2, verbose = FALSE, trace = 0) {
  
  gwas.params <- calc.pr(object, trace = trace)
  if (verbose) report(gwas.params, "Prior probabilities")
  
  if (!missing(adjust)) {
    gwas.params$p.r <- gwas.params$p.r / (1 + adjust)
    if (verbose) report(gwas.params, "LD adjusted prior probabilities")
  } 
  
  cond.probs <- calc.conditionals(object, risk.thresh = risk.thresh, adjust = adjust)
  if (verbose) report(cond.probs, "Conditional probabilities")
  
  
  post.probs <- data.frame(label = label.groups(fcols(object)),
                          marker = marker(object),
                          zscore = zscore(object),
                             p.r = gwas.params$p.r,
                             p.n = gwas.params$p.n,
                           stringsAsFactors = FALSE)
  
  
  # density of zscores for risk SNPs and conditional on null density (inflated)
  fz_e <- function(z, l) (dnorm((2 - z) / l) + dnorm((-2 - z) / l)) / 2 / l
  fz_n <- function(z, l)  dnorm(z / l) / l
  
  post.probs <- transform(post.probs,
                          fz.e = fz_e(zscore, gwas.params$lambda), 
                          fz.n = fz_n(zscore, gwas.params$lambda))
  
  # Add conditional feature probabilities
  l.index <- match(post.probs$label, cond.probs$label)
  post.probs <- transform(post.probs,
                          p.f.r = cond.probs$prob[l.index],
                          p.f.n = cond.probs$prob.base[l.index])
  
  # Calculate posterior probabilities by feature combination groups                  
  post.probs <- ddply(post.probs, "label", transform,
                 post.prob.gwas = (fz.e * p.r) /
                                 ((fz.e * p.r) + fz.n * p.n),
                      post.prob = (p.f.r * fz.e * p.r) / 
                                 ((p.f.r * fz.e * p.r) + p.f.n * fz.n * p.n))  
  
  m.index <- match(marker(object), post.probs$marker)
  post.probs <- post.probs[m.index, c("label", "post.prob.gwas", "post.prob")]
  post.probs <- rename(post.probs, c(label = "feature.category"))
  mcols(object) <- DataFrame(mcols(object), post.probs)
  
  return(object)
})



#' Calculate prior probability of being a risk SNP
#'
#' @param object \code{AnnotatedGWAS} object
#' @param effect assumed effect size, used to shift the ideal z-score distribution
#' @param trace Non-negative integer. If positive, tracing information on the
#'  progress of the optimization is produced. See \code{\link[stats]{optim}} for
#'  details.
#'  
#'  @return A list with components:
#' \describe{
#'  \item{\code{p.r}}{P(R): prior probability of being a risk SNP}
#'  \item{\code{p.n}}{P(N): prior probability of a SNP having no effect}
#'  \item{\code{lambda}}{GWAS inflation factor (\eqn{\lambda})}
#' }

setGeneric("calc.pr", function(object, effect = 2, trace = 0) {
  standardGeneric("calc.pr")
})

setMethod("calc.pr", "AnnotatedGWAS", function(object, effect = 2, trace = 0) {
    
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

  out <- list(p.r = roots$par[2], 
              p.n = 1 - roots$par[2], 
           lambda = roots$par[1])
  out <- structure(out, class = c("gwas.priors", class(out)))
  return(out)  
})
