% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/stats.r
\name{calc.threshold}
\alias{calc.threshold}
\title{Calculate p-value threshold for identifying risk SNPs}
\usage{
calc.threshold(p, verbose = FALSE, method = "BH")
}
\arguments{
\item{p}{numeric vector of p-values (possibly with \code{\link{NA}}s).
    Any other \R is coerced by \code{\link{as.numeric}}.}

\item{verbose}{print out intermediate statistics calculated during the analysis}

\item{method}{correction method}
}
\description{
We select at least 100 SNPs with FDR < 0.05 or 1000 SNPs with FDR < 0.1. If
this can't be done then we use 100 SNPs but warn the estimate is ureliable.
}

