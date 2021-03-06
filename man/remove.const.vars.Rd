% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{remove.const.vars}
\alias{remove.const.vars}
\title{Remove constant variables}
\usage{
remove.const.vars(df)
}
\arguments{
\item{df}{data.frame from which constant variable will be removed.}
}
\value{
data.frame with removed constant variables.
}
\description{
remove constant variables from a data.frame.
}
\examples{
df <- data.frame(a = rnorm(5), b = 1:5, c = rep(1, 5))
df <- remove.const.vars(df)
}
