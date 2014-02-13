\name{calc.pe}
\alias{calc.pe}
\title{Calculate posterior probability of being an effective SNP}
\usage{
calc.pe(object, effect = 2, trace = 0)
}
\arguments{
  \item{object}{\code{AnnotatedGWAS} object}

  \item{effect}{assumed effect size, used to shift the
  ideal z-score distribution}

  \item{trace}{Non-negative integer. If positive, tracing
  information on the progress of the optimization is
  produced. See \code{\link[stats]{optim}} for details.}
}
\value{
A list with components: \describe{ \item{\code{p.e}}{P(E):
posterior probability of being an effective SNP}
\item{\code{p.n}}{P(N): posterior probability of a SNP
having no effect} \item{\code{lambda}}{GWAS inflation
factor (\eqn{\lambda})} }
}
\description{
Calculate posterior probability of being an effective SNP
}

