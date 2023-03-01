\name{mposterior}
\alias{mposterior}

\title{Maximizing Posterior Distributions}
\description{
  Maximizing prior and posterior distibutions for the
  location (with optional trend), scale and shape parameters under
  the gev, order statistics or point process models.
}
    
\usage{
mposterior(init, prior, lh = c("none", "gev", "gpd", "pp", "os"),
    method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
    lower = -Inf, upper = Inf, control = list(), hessian = FALSE, ...)
}

\arguments{
  \item{init}{Numeric vector of length three/four, giving the initial
    values for the optimization.}
  \item{prior}{An object of class \code{"evprior"}, constructed
    using \code{prior.prob}, \code{prior.quant} or
    \code{prior.norm}.}
  \item{lh}{A character string specifying the likelihood; either
    \code{"gev"} for gev, \code{"gpd"} for gpd, \code{"os"} for order
    statistics, \code{"pp"} for Poisson process or \code{"none"} for none
    (the default). The latter can be used to maximize the
    prior distribution.}
  \item{method}{The method to be used. See \code{\link{optim}}.}
  \item{lower, upper}{Bounds on the variables for the \code{"L-BFGS-B"}
    method. See \code{\link{optim}}.}
  \item{control}{A list of control parameters. See \code{\link{optim}}.}
  \item{hessian}{Logical. See \code{\link{optim}}.}
  \item{...}{Arguments to the likelihood. Should include \code{data}
    unless \code{lh} is \code{"none"}. Should also include \code{thresh}
    and \code{noy} if \code{lh} is \code{"pp"}.
    Should include the vector \code{trend} if a linear trend on
    the location is implemented.
    See \code{\link{pplik}} for details.}
}


\value{
  A list. See \code{\link{optim}}.
}

\seealso{\code{\link{pplik}}, \code{\link{posterior}}, \code{\link{prior.prob}}}

\keyword{optimize}

