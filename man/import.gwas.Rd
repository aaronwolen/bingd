\name{import.gwas}
\alias{import.gwas}
\title{Import GWAS results}
\usage{
import.gwas(file, snpid, chr, pos, pval, or, genome, ...)
}
\arguments{
  \item{file}{}

  \item{snpid}{column label corresponding to SNP ID column}

  \item{chr}{column label corresponding to chromosome
  column}

  \item{pos}{}

  \item{pval}{}

  \item{or}{}

  \item{genome}{the genome identifier or assembly name}

  \item{...}{extra parameters passed to read.table}
}
\value{
GRanges object
}
\description{
Import GWAS results
}
\examples{
scz <- import.gwas("data/pgc.scz.clump.txt.gz", snpid = "snpid",
                   chr = "chr", pos = "pos", pval = "pval",
                   or = "or", genome = "hg19")
}

