% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/bayes-methods.r
\docType{methods}
\name{calc.bayes}
\alias{calc.bayes}
\alias{calc.bayes,AnnotatedGWAS-method}
\title{Calculate posterior probability of being an effective SNP}
\usage{
calc.bayes(object, risk.thresh = NULL, adjust = NULL, effect = 2,
  prob.risk = NULL, verbose = FALSE, trace = 0)

\S4method{calc.bayes}{AnnotatedGWAS}(object, risk.thresh = NULL,
  adjust = NULL, effect = 2, prob.risk = NULL, verbose = FALSE,
  trace = 0)
}
\arguments{
\item{object}{\code{AnnotatedGWAS} object}

\item{risk.thresh}{the p-value threshold applied to identify ``risky'' markers}

\item{adjust}{optional parameter to adjust P(R), as well as conditional
probabilities for linkage disequilibrium}

\item{effect}{assumed effect size, used to shift the ideal z-score distribution}

\item{prob.risk}{probability of being a risk SNP, this is normally calculated
automatically but can be supplied in situations where the result was 0}

\item{verbose}{print out intermediate statistics calculated during the analysis}

\item{trace}{Non-negative integer. If positive, tracing information on the
 progress of the optimization is produced. See \code{\link[stats]{optim}} for
 details.}
}
\description{
Integrate GWAS results with annotated genomic features
}

