#' Calculate posterior probability of being an effective SNP
#' 
#' Integrate GWAS results with annotated genomic features
#' 
#' @param object \code{AnnotatedGWAS} object
#' @param adjust optional parameter to adjust P(R), as well as conditional
#' probabilities for linkage disequilibrium
#' @param prob.risk probability of being a risk SNP, this is normally calculated
#' automatically but can be supplied in situations where the result was 0
#' @inheritParams calc.priors
#' @inheritParams calc.conditionals
#' 
#' @importFrom plyr ddply
#' @export

setGeneric("calc.bayes", 
  function(object, risk.thresh = NULL, adjust = NULL, effect = 2, prob.risk = NULL,
           verbose = FALSE, trace = 0) {
  standardGeneric("calc.bayes")
})

setMethod("calc.bayes", "AnnotatedGWAS", 
  function(object, risk.thresh = NULL, adjust = NULL, effect = 2, prob.risk = NULL,  
           verbose = FALSE, trace = 0) {
  
  gwas.params <- calc.priors(object, effect, adjust, verbose, trace)
  
  if (!is.null(prob.risk)) {
    gwas.params$p.r <- prob.risk
    gwas.params$p.n <- 1 - prob.risk
  }
  
  if (gwas.params$p.r == 0) stop("P(R) is zero, the analysis cannot proceed.")
    
  cond.probs <- calc.conditionals(object, risk.thresh, adjust, verbose)
  
  post.probs <- data.frame(label = label.groups(fcols(object)),
                          marker = marker(object),
                          zscore = zscore(object),
                             p.r = gwas.params$p.r,
                             p.n = gwas.params$p.n,
                           stringsAsFactors = FALSE)
  
  # density of zscores for risk SNPs and conditional on null density (inflated)
  fz_r <- function(z, l) (dnorm((2 - z) / l) + dnorm((-2 - z) / l)) / 2 / l
  fz_n <- function(z, l)  dnorm(z / l) / l
  
  post.probs <- transform(post.probs,
                          fz.r = fz_r(zscore, gwas.params$lambda), 
                          fz.n = fz_n(zscore, gwas.params$lambda))
  
  # Add conditional feature probabilities
  l.index <- match(post.probs$label, cond.probs$label)
  post.probs <- transform(post.probs,
                          p.f.r = cond.probs$prob[l.index],
                          p.f.n = cond.probs$prob.base[l.index])
  
  # Calculate posterior probabilities by feature combination groups                  
  post.probs <- ddply(post.probs, "label", transform,
                 post.prob.gwas = (fz.r * p.r) /
                                 ((fz.r * p.r) + fz.n * p.n),
                      post.prob = (p.f.r * fz.r * p.r) / 
                                 ((p.f.r * fz.r * p.r) + p.f.n * fz.n * p.n))  
  
  m.index <- match(marker(object), post.probs$marker)
  post.probs <- post.probs[m.index, c("label", "post.prob.gwas", "post.prob")]
  post.probs <- rename(post.probs, c(label = "feature.category"))
  mcols(object) <- DataFrame(mcols(object), post.probs)
  
  return(object)
})
