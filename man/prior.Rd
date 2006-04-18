\name{prior}
\alias{prior.prob}
\alias{prior.quant}
\alias{prior.norm}
\alias{prior.loglognorm}

\title{Construction of Prior Distributions}
\description{
  Constructing prior distibutions for the location, scale and shape
  parameters using normal, beta or gamma distributions.
  A linear trend for the location can also be specified, using a
  prior normal distribution centered at zero for the trend
  parameter.
  
}
\usage{
prior.prob(quant, alpha, trendsd = 0)
prior.quant(prob = 10^-(1:3), shape, scale, trendsd = 0)
prior.norm(mean, cov, trendsd = 0)
prior.loglognorm(mean, cov, trendsd = 0)
}

\arguments{
  \item{quant, alpha}{Numeric vectors of length three and four
    respectively.
    Beta prior distibutions are placed on probability ratios
    corresponding to the quantiles given in \code{quant}.}
  \item{prob, shape, scale}{Numeric vectors of length three.
    Gamma prior distibutions, with parameters \code{shape} and
    \code{scale}, are placed on quantile differences corresponding
    to the probabilities given in \code{prob}.}
  \item{mean, cov}{The prior distibution for the location, log(scale)
    and shape is taken to be trivariate normal, with mean \code{mean}
    (a numeric vector of length three) and covariance matrix \code{cov}
    (a symmetric positive definite three by three matrix).}
  \item{trendsd}{The standard deviation for the marginal normal prior
    distribution (with mean zero) placed on the linear trend parameter
    for the location. If this is zero (the default) a linear trend is
    not implemented.}
}
\details{
  See the user's guide. 
}

\value{
  Returns an object of class \code{"evprior"}, which is essentially
  just a list of the arguments passed.
}

\seealso{\code{\link{posterior}}, \code{\link{pplik}}}

\examples{
mat <- diag(c(10000, 10000, 100))
prior.norm(mean = c(0,0,0), cov = mat, trendsd = 10)
prior.quant(shape = c(38.9,7.1,47), scale = c(1.5,6.3,2.6))
prior.prob(quant = c(85,88,95), alpha = c(4,2.5,2.25,0.25))
}
\keyword{models}

