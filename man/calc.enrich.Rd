\name{calc.enrich}
\alias{calc.enrich}
\title{Calculate enrichment of features across thresholds}
\usage{
calc.enrich(object, feature.list = NULL, stat, thresh.levels)
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
\code{gwas} must be annotated with
\code{\link{annotate.gwas}} OR a \code{feature.list} must
be provided.
}

