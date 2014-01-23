\name{calc.enrich}
\alias{calc.enrich}
\title{Calculate enrichment of all annotated features across thresholds}
\usage{
calc.enrich(gwas, feature.list, stat, thresh.levels)
}
\arguments{
  \item{gwas}{GWAS \code{\link{GRanges}} object}

  \item{stat}{numeric vector containing statistic to which
  thresholds will be applied for each enrichment
  calculation}

  \item{feature.list}{\code{FeatureList} object created
  with \code{\link{hub.search}}}
}
\description{
\code{gwas} must be annotated with
\code{\link{annotate.gwas}} OR a \code{feature.list} must
be provided.
}

