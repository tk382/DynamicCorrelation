% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_cor_r}
\alias{get_cor_r}
\title{Helper function for get_est_H
Returns each element of matrix H}
\usage{
get_cor_r(rho12, rho23, rho13)
}
\arguments{
\item{rho12}{correlation between variable 1 and variable 2}

\item{rho23}{correlation between variable 2 and variable 3}

\item{rho13}{correlation between variable 1 and variable 3}
}
\value{
element (j,k) of H from elements (i,j), (j,k), (i,k) of Sigma
}
\description{
Helper function for get_est_H
Returns each element of matrix H
}
