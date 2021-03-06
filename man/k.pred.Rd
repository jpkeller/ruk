% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_aux.R
\name{k.pred}
\alias{k.pred}
\title{Kriging predictions}
\usage{
k.pred(reg.dat, dim.X, beta, log.cov.pars, cov.model)
}
\arguments{
\item{reg.dat}{matrix of regression data}

\item{dim.X}{TBA}

\item{beta}{TBA}

\item{log.cov.pars}{TBA}
}
\description{
Internal function for computing kriging predictions (i.e. conditional expectations).
}
\details{
Not intended to be called directly. See 
 \code{\link{predict.rlikfit}} for
	making kriging predictions from an \code{rlikfit} object.
}
\keyword{internal}

