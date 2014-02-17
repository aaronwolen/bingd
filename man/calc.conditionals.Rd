\name{calc.conditionals}
\alias{calc.conditionals}
\title{Calculate conditional probabilities for annotated features}
\usage{
calc.conditionals(object, risk.thresh, adjust)
}
\arguments{
  \item{object}{\code{AnnotatedGWAS} object}

  \item{risk.thresh}{the p-value threshold applied to
  identify ``risky'' markers}

  \item{adjust}{optional parameter to adjust baseline and
  conditional probabilities}
}
\value{
data.frame containing probabilities of observing each
combination of annotated features.
}
\description{
Calculate conditional probabilities for annotated features
}

