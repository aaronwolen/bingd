\docType{methods}
\name{calc.enrich,GWAS-method}
\alias{calc.enrich,GWAS-method}
\title{Calculate enrichment of FeatureList features across thresholds}
\usage{
\S4method{calc.enrich}{GWAS}(object, feature.list, stat, thresh.levels)
}
\arguments{
  \item{object}{\code{\link{GWAS}} object}

  \item{feature.list}{\code{FeatureList} object created
  with \code{\link{hub.features}} or
  \code{\link{local.features}}}

  \item{stat}{numeric vector containing statistic to which
  thresholds will be applied for each enrichment
  calculation}

  \item{thresh.levels}{a numeric vector containing the
  various thresholds at which \code{feature} enrichment is
  calculated}
}
\description{
Calculate enrichment of FeatureList features across
thresholds
}

