\name{annotate.gwas}
\alias{annotate.gwas}
\title{Annotate GWAS object using AnnotationHub}
\usage{
annotate.gwas(gwas, data.filter)
}
\arguments{
  \item{gwas}{GWAS GRanges object}

  \item{filter}{named list of vectors}
}
\description{
Annotate GWAS object using AnnotationHub
}
\examples{
scz <- annotate.gwas(scz, data.filter = list(DNaseI = c("UwDnaseN", "broadPeak")))
}

