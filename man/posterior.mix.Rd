\name{posterior.mix}
\alias{posterior.mix}

\title{MCMC Sampling of Posterior Distributions with mixture}
\description{
  Constructing MCMC samples of prior and posterior distibutions for
  the location (with optional trend), scale and shape parameters under
  the gev, order statistics or point process models. Note this function
  is only effective for a mixture defined in the prior distribution.
}
\usage{
posterior.mix(n, init, prior, lh = c("none", "gev", "gpd"), \dots, psd,
    pMassProb, normPi0, xitilde, pMass = 0, cv, burn = 0, thin = 1)
}

\arguments{
  \item{n}{The run-length; the number of sampled vectors
    (excluding \code{init}).}
  \item{init}{Numeric vector of length three/four, giving the initial
    values for the chain, taken to be iteration zero.}
  \item{prior}{An object of class \code{"evprior"}, constructed
    using \code{prior.prob}, \code{prior.quant} or
    \code{prior.norm}.}
  \item{lh}{A character string specifying the likelihood; either
    \code{"gev"} for gev, \code{"os"} for order statistics,
    \code{"pp"} for Poisson process or \code{"none"} for none
    (the default). The latter can be used to sample from the
    prior distribution.}
  \item{\dots}{Arguments to the likelihood. Should include \code{data}
    unless \code{lh} is \code{"none"}. Should also include \code{thresh}
    and \code{noy} if \code{lh} is \code{"pp"}.
    Should include the vector \code{trend} if a linear trend on
    the location is implemented.
    See \code{\link{pplik}} for details.}
  \item{psd}{A vector of length three/four containing standard
    deviations for proposal distributions.}
  \item{pMassProb}{The probability for the point Mass.}
  \item{normPi0}{The constant such as the point mass prior distribution
    is a distribution.}
  \item{xitilde}{The mode of the marginal posterior distribution for the
  shape parameter - from a conventional Bayesian analysis i.e. without a
  point Mass.}
  \item{pMass}{The point mass for the shape parameter.}
  \item{cv}{A vector of length two containing probabilities for
    quantiles matching.}
  \item{burn}{The burn-in period (an integer); the first \code{burn}
    iterations (including \code{init}) are excluded from the chain.}
  \item{thin}{The thinning interval (an integer); iteration \eqn{k}
    is stored only if \eqn{k} mod \code{thin} is zero (and if \eqn{k}
    greater than or equal to \code{burn}).}
}
\details{
  See the user's guide. 
}

\value{
  A matrix with \code{1+floor(n/thin)-burn} rows.
  Row labels give the iteration numbers.
  Column labels give parameter names.
  
  An attribute \code{ar} is also returned.
  This is a matrix containing acceptence rates in the first row
  (the number of proposals accepted divided by the number of
  iterations) and ``external rates'' in the second (the number of
  proposals that resulted in a zero likelihood, divided by the
  number of iterations).
}

\seealso{\code{\link{pplik}}, \code{\link{prior.prob}}}

\keyword{models}

