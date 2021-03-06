% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{write.excel}
\alias{write.excel}
\title{Copy data to clipboard suitable for Excel pasting.}
\usage{
write.excel(x, row.names = TRUE, col.names = TRUE, ...)
}
\arguments{
\item{x}{vector, data.frame or matrix to copy.}

\item{row.names}{copy with row names.}

\item{col.names}{copy with column names.}
}
\value{
nothing, just copy object to clipboard.
}
\description{
Copy data to clipboard suitable for Excel pasting.
}
\details{
column names in copied object will be shifted one column left.
}
\examples{
a <- 1:10
names(a) <- LETTERS[1:10]
\dontrun{write.excel(a)}
b <- data.frame(V1 = 1:10, V2 = rnorm(10), row.names = LETTERS[1:10])
\dontrun{write.excel(b)}
}
