\name{calc.enrich}
\alias{calc.enrich}
\title{Calculate enrichment of all annotated features across thresholds}
\usage{
calc.enrich(object, stat, thresh.levels, data.filter, online = NULL)
}
\arguments{
  \item{gwas}{GWAS GRanges object}

  \item{stat}{numeric vector containing statistic to which
  thresholds will be applied for each enrichment
  calculation}

  \item{data.filter}{optional, calculate enrichment among
  features that match data.filter to avoid annotating the
  gwas object, which can require lots of memory}
}
\description{
Calculate enrichment of all annotated features across
thresholds
}

