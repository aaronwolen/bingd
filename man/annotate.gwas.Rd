\name{annotate.gwas}
\alias{annotate.gwas}
\title{Annotate GWAS object using AnnotationHub}
\usage{
annotate.gwas(object, feature.list)
}
\arguments{
  \item{object}{\code{\link{GWAS}} object}

  \item{feature.list}{\code{FeatureList} object created
  with \code{\link{hub.features}} or
  \code{\link{local.features}}}
}
\description{
Annotate GWAS object using AnnotationHub
}
\examples{
query <- list(DNaseI = c("Dnasen", "broadPeak", "rep1"))
feature.list <- hub.features(query)
scz <- annotate.gwas(scz, feature.list)
}

