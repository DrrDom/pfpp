% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ggplot.boxplot.stats}
\alias{ggplot.boxplot.stats}
\title{Calculate boxplot statistics for ggplot}
\usage{
ggplot.boxplot.stats(df, group_by, value.col, ncpu = 1)
}
\arguments{
\item{df}{input data.frame.}

\item{group_by}{column names or indices used for grouping data}

\item{value.col}{column name of index of a value used to calculate boxplot statistics}

\item{ncpu}{number of cores to calculate statistics. If more than 1 parallel::mclapply function will be used}
}
\value{
list of two data.frames, the first one with boxplot statistics and the second one with outliers
}
\description{
Calculate boxplot statistics for ggplot
}
\details{
uses boxplot.stats as backend. All columns not in group_by will be dropped in output 
data.frame with outliers
}
\examples{
d <- data.frame(a = sample(LETTERS[1:3], 100, replace = TRUE), 
                b = sample(LETTERS[1:3], 100, replace = TRUE), 
                c = sample(c(1:3, 10), 100, replace = TRUE),
                d = sample(c(1:3, 10), 100, replace = TRUE))
                
ggplot.boxplot.stats(d, c("a", "b"), 3L)
ggplot.boxplot.stats(d, c("a", "b"), "c")                 
}
