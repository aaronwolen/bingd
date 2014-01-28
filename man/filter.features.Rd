\name{filter.features}
\alias{filter.features}
\title{Filter list of features based on search terms}
\usage{
filter.features(query, file.list)
}
\arguments{
  \item{query}{a vector of strings to filter features; may
  also be a named list of vectors to perform multiple
  queries for different categories of features}

  \item{file.list}{a \code{\link{DataFrame}} containing, at
  minimum, \code{Title} and \code{LocalPath} columns}
}
\description{
Filter list of features based on search terms
}

