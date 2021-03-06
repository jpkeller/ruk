% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_aux.R
\name{varcov.eff}
\alias{varcov.eff}
\title{Creates exponential covariance matrix}
\usage{
varcov.eff(coords, pars, cov.model)
}
\arguments{
\item{coords}{coordinates.}

\item{pars}{A vector of three log covariance parameters: nugget, sill, and range.}

\item{cov.model}{character string}
}
\value{
Returns an \eqn{n x n} matrix.
}
\description{
Given covariance parameters, creates a covariance matrix
	according to an expoential model.
}
\details{
If the sill and range are zero
}
\author{
Paul Sampson
}

