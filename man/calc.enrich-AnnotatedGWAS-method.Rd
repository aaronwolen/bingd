\docType{methods}
\name{calc.enrich,AnnotatedGWAS-method}
\alias{calc.enrich,AnnotatedGWAS-method}
\title{Calculate enrichment of annotated features across thresholds}
\usage{
\S4method{calc.enrich}{AnnotatedGWAS}(object, stat, thresh.levels)
}
\arguments{
  \item{object}{\code{\link{GWAS}} object}

  \item{stat}{numeric vector containing statistic to which
  thresholds will be applied for each enrichment
  calculation}

  \item{thresh.levels}{a numeric vector containing the
  various thresholds at which \code{feature} enrichment is
  calculated}
}
\description{
Calculate enrichment of annotated features across
thresholds
}

