\name{calc.pr}
\alias{calc.pr}
\title{Calculate prior probability of being a risk SNP}
\usage{
calc.pr(object, effect = 2, trace = 0)
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
A list with components: \describe{ \item{\code{p.r}}{P(R):
prior probability of being a risk SNP}
\item{\code{p.n}}{P(N): prior probability of a SNP having
no effect} \item{\code{lambda}}{GWAS inflation factor
(\eqn{\lambda})} }
}
\description{
Calculate prior probability of being a risk SNP
}

