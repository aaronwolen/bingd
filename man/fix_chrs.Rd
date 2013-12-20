\name{fix_chrs}
\alias{fix_chrs}
\title{Fix a GRanges object's chromosome information}
\usage{
fix_chrs(object, genome = "hg19", add.length = TRUE, db = NA)
}
\arguments{
  \item{object}{object of class \code{GRanges}}

  \item{add.length}{Add missing chromosome length
  information}

  \item{db}{Optional bioconductor organism package (e.g.,
  org.Hs.eg.db) to provide chromosome lengths so that an
  internet connection isn't necessary}
}
\description{
Edit a GRanges object's chromosome information by adding
chromosome lengths from UCSC and properly ordering the
seqnames factor levels. Adding lengths requires an internet
connection. If one is not available then \code{add.length}
should be set to \code{FALSE}.
}

