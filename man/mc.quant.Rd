\name{mc.quant}
\alias{mc.quant}

\title{Compute GEV Quantiles from Markov Chains}
\description{
  Compute gev quantiles from samples stored within a Markov chain,
  corresponding to specified probabilities in the upper tail.
}

\usage{
mc.quant(post, p, lh = c("gev", "gpd"))
}

\arguments{
  \item{post}{A Markov chain generated using \code{posterior},
    containing samples of gev parameters.}
  \item{p}{A numeric vector of upper tail probabilities.}
  \item{lh}{Specify ``gev'' or ``gpd'' likelihood.}
}

\value{
  A matrix with \eqn{n} rows and \eqn{m} columns, where \eqn{n} is the
  number of samples stored within the chain, and \eqn{m} is the
  length of the vector \eqn{p}.
  If \eqn{m = 1} the dimension is dropped (i.e. a vector of length
  \eqn{n} is returned).
  The \code{(i,j)}th entry contains the gev quantile coresponding to the
  upper tail probability \code{p[j]}, evaluated at the parameters
  within sample \code{i}.
  
  If a linear trend on the location has been implemented, the quantiles
  correspond to the distribution obtained when the trend parameter is
  zero.
}
\details{
  See the user's guide. 
}

\seealso{\code{\link{posterior}}}

\keyword{manip}

