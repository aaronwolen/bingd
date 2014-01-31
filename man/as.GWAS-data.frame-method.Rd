\docType{methods}
\name{as.GWAS,data.frame-method}
\alias{as.GWAS,data.frame-method}
\title{Create a GWAS object}
\usage{
\S4method{as.GWAS}{data.frame}(object, marker, chr, bp, pvalue, or, beta)
}
\arguments{
  \item{object}{a \code{\link{data.frame}} or
  \code{\link[GenomicRanges]{GRanges}} object containing
  the required fields necessary to construct a
  \code{\link{GWAS}} object}

  \item{marker}{name of the column in \code{object} that
  contains the marker (or SNP) identifiers}

  \item{chr}{name of the column in \code{object} that
  contains the chromosome names associated with each
  marker}

  \item{bp}{name of the column in \code{object} that
  contains the genomic position (or start position)
  associated with each marker}

  \item{pvalue}{name of the column in \code{object} that
  contains the GWAS association p-value for each marker}

  \item{or}{name of the column in \code{object} that
  contains the GWAS odds ratio for each marker}

  \item{beta}{name of the column in \code{object} that
  contains the GWAS beta values (or regression
  coefficients) for each marker}
}
\description{
Create a \code{\link{GWAS}} object from a
\code{\link{data.frame}}
}

