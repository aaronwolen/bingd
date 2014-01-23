\name{annotate.gwas}
\alias{annotate.gwas}
\title{Annotate GWAS object using AnnotationHub}
\usage{
annotate.gwas(gwas, feature.list)
}
\arguments{
  \item{gwas}{GWAS GRanges object}

  \item{feature.list}{\code{FeatureList} object created
  with \code{\link{hub.search}}}
}
\description{
Annotate GWAS object using AnnotationHub
}
\examples{
query <- list(DNaseI = c("Dnasen", "broadPeak", "rep1"))
feature.list <- hub.search(query, online = FALSE)
scz <- annotate.gwas(scz, feature.list)
}

