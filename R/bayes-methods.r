#' Calculate posterior probability of being an effective SNP
#' 
#' Integrate GWAS results with annotated genomic features
#' 
#' @param object \code{AnnotatedGWAS} object
#' @param adjust optional parameter to adjust P(R), as well as conditional
#' probabilities for linkage disequilibrium
#' @inheritParams calc.priors
#' @inheritParams calc.conditionals
#' 
#' @importFrom plyr ddply
#' @export

setGeneric("calc.bayes", 
  function(object, risk.thresh = NULL, adjust = NULL, effect = 2, 
           verbose = FALSE, trace = 0) {
  standardGeneric("calc.bayes")
})

setMethod("calc.bayes", "AnnotatedGWAS", 
  function(object, risk.thresh = NULL, adjust = NULL, effect = 2, 
           verbose = FALSE, trace = 0) {
  
  gwas.params <- calc.priors(object, effect, adjust, verbose, trace)
    
  cond.probs <- calc.conditionals(object, risk.thresh, adjust, verbose)
  
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
