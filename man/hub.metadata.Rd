\name{hub.metadata}
\alias{hub.metadata}
\title{Retrieve AnnotationHub metadata}
\usage{
hub.metadata(online = FALSE, cache.dir = "default")
}
\arguments{
  \item{online}{if TRUE search is conducted using latest
  AnnotationHub metadata, otherwise only the AnnotationHub
  cache directory is searched}

  \item{cache.dir}{define custom location used for
  AnnotationHub's cached directory, which should contain a
  'resources' subdirectory. Normally this shouldn't be
  changed.}
}
\description{
Retrieve AnnotationHub metadata
}

