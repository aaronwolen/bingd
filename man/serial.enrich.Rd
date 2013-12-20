\name{serial.enrich}
\alias{serial.enrich}
\title{Calculate the enrichment of a particular discrete feature across thresholds}
\usage{
serial.enrich(feature, stat, thresh.levels)
}
\arguments{
  \item{feature}{logical or character vector}

  \item{stat}{numeric vector containing statistic to which
  thresholds will be applied for each enrichment
  calculation}
}
\value{
data.frame containing `enrichment`
}
\description{
Calculate the enrichment of a particular discrete feature
across thresholds
}

