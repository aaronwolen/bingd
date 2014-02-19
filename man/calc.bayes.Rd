\name{calc.bayes}
\alias{calc.bayes}
\title{Calculate posterior probability of being an effective SNP}
\usage{
calc.bayes(object, risk.thresh = NULL, adjust = NULL, effect = 2,
  verbose = FALSE, trace = 0)
}
\arguments{
  \item{object}{\code{AnnotatedGWAS} object}

  \item{adjust}{optional parameter to adjust P(R), as well
  as conditional probabilities for linkage disequilibrium}

  \item{effect}{assumed effect size, used to shift the
  ideal z-score distribution}

  \item{trace}{Non-negative integer. If positive, tracing
  information on the progress of the optimization is
  produced. See \code{\link[stats]{optim}} for details.}

  \item{risk.thresh}{the p-value threshold applied to
  identify ``risky'' markers}
}
\description{
Integrate GWAS results with annotated genomic features
}

