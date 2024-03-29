% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DCARS.R
\name{prop_nonzero_matrix}
\alias{prop_nonzero_matrix}
\title{prop_nonzero_matrix}
\usage{
prop_nonzero_matrix(dat, PPI, offset = TRUE)
}
\arguments{
\item{dat}{matrix with rownames appearing in PPI}

\item{PPI}{two column matrix with entries that should appear in rownames of dat}

\item{offset}{should 1 be added when pospos is 0 or all}
}
\value{
\code{numeric} of weighted correlations for the sequence of weights given
}
\description{
The prop_nonzero_matrix function calculates the proportion of observations with nonzero values in both variables across a number of pairs of variables, gived in PPI
}
\examples{





}
