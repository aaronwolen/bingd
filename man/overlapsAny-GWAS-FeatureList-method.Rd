\docType{methods}
\name{overlapsAny,GWAS,FeatureList-method}
\alias{overlapsAny,GWAS,FeatureList-method}
\alias{overlapsAny-method}
\title{Find features overlapping GWAS markers}
\usage{
\S4method{overlapsAny}{GWAS,FeatureList}(query, subject, ...)
}
\arguments{
  \item{query}{\code{GWAS} object}

  \item{subject}{\code{FeatureList} object}

  \item{...}{additional arguments passed to
  \code{\link[IRanges]{overlapsAny}}}
}
\value{
A list of \code{DataFrame}s, the structure of which mirrors
that of the supplied \code{FeatureList} object. Each
element of the list is a \code{DataFrame} with one logical
vector for each feature, indicating marker overlap.
}
\description{
Determine which GWAS markers overlap each
\code{FeatureList} feature.
}

