\name{label.groups}
\alias{label.groups}
\title{Generate labels based on column names of logical variables in a data.frame}
\usage{
label.groups(df, sep = ".and.", empty = "neither")
}
\description{
Generate labels based on column names of logical variables
in a data.frame
}
\examples{
label.groups(data.frame(a = c(T,F, T), b = c(F, T, T)))
}

