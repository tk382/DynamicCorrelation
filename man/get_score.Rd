% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_score}
\alias{get_score}
\title{Computes the score statistic}
\usage{
get_score(x, y1, y2)
}
\arguments{
\item{x}{covariate matrix to test}

\item{y1}{vector of first variable}

\item{y2}{vector of second variable}
}
\value{
double score statistic
}
\description{
Computes the score statistic
}
\examples{
x = rnorm(100)
y1 = rnorm(100)
y2 = rnorm(100)
q = get_score(x, y1, y2)
}
