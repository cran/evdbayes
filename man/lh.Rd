\name{lh}
\alias{gevlik}
\alias{gpdlik}
\alias{pplik}
\alias{oslik}

\title{Calculate Log-likelihoods}
\description{
  Calculate log-likelihoods for the gev, order statistics or
  point process models.
}
\usage{
pplik(par, data, thresh, noy, trend, exact = FALSE)
gevlik(par, data, trend)
gpdlik(par, data, trend)
oslik(par, data, trend) 
}
\arguments{
  \item{par}{If \code{trend} is missing, should be a numeric vector of
    length three, containing the location, scale and shape parameters.
    If \code{trend} is not missing, should be a numeric vector of length
    four, containing the location intercept, scale, shape and location
    trend parameters, in that order.}
  \item{data}{For \code{pplik}, \code{gevlik} and \code{gpdlik}; a
    non-empty numeric vector containing the data at which the likelihood
    is evaluated, possibly containing missing values. For \code{oslik};
    a numeric matrix (see the user's guide).}
  \item{thresh}{Threshold. Typically a single number or a vector of
    the same length as \code{data}.}
  \item{noy}{Number of years/periods of observations, excluding
    any missing values.}
  \item{trend}{Trend vector (optional). If given, should be the same
    length as \code{data} for \code{pplik} and \code{gevlik}.
    For \code{oslik}, should contain one value for each row
    of \code{data}.}
  \item{exact}{In general, the point process likelihood includes an
    approximation to an integral. If \code{exact} is \code{TRUE}, every
    value in \code{trend} and \code{thresh} is used for the
    approximation.}
}
\details{
  See the user's guide. 
}
\note{
  These functions are essentially internal, and need not be called
  by the user. They are documented only because their arguments
  (excluding \code{par}) can be passed to
  \code{\link{posterior}}.
}

\seealso{\code{\link{posterior}}, \code{\link{prior.prob}}}

\keyword{models}

