\name{hub.search}
\alias{hub.search}
\title{Search AnnotationHub for features}
\usage{
hub.search(query, genome, md, online = FALSE, cache.dir = "default")
}
\arguments{
  \item{query}{a vector of strings to filter AnnotationHub
  features; may also be a named list of vectors to perform
  multiple queries}

  \item{online}{if TRUE search is conducted using latest
  AnnotationHub metadata, otherwise only the AnnotationHub
  cache directory is searched}

  \item{cache.dir}{define custom location used for
  AnnotationHub's cached directory, which should contain a
  'resources' subdirectory. Normally this shouldn't be
  changed.}
}
\value{
A named list of DataFrames containing (at minimum) columns
indicating the Title and location (LocalPath) of features
matching search query
}
\description{
Search AnnotationHub for features
}

