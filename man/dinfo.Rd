\name{dinfo}
\alias{ibeta}
\alias{igamma}

\title{Information for Beta and Gamma Distributions}
\description{
  Show means, variances and modes for beta and gamma distributions.
}

\usage{
ibeta(mean, var, shape1, shape2)
igamma(mean, var, shape, scale)
}

\arguments{
  \item{mean, var}{Numeric vectors giving means and variances.}
  \item{shape1, shape2}{Numeric vectors. See \code{\link{dbeta}}.}
  \item{shape, scale}{Numeric vectors. See \code{\link{dgamma}}.}  
}
\details{
  For \code{ibeta}, either both of \code{mean} and \code{var} or both of
  \code{shape1} and \code{shape2} must be specified.
  For \code{igamma}, either both of \code{mean} and \code{var} or both of
  \code{shape} and \code{scale} must be specified.
  The pair of vectors that are passed to each function define a set
  of beta/gamma distributions.
  If one vector is shorter than the other, the shorter vector is
  replicated.
}

\value{
  A matrix with five columns and \eqn{n} rows, where \eqn{n} is the
  length of the longest argument.
  If \eqn{n = 1} the dimension is dropped (i.e. a vector of length
  five is returned).
  The columns contain the means, variances, modes, and the
  shape/scale parameters of the specified distributions.
  If a mode is \code{NA}, it does not exist, or it is not unique, or
  it does not occur in the interior of the support.
  If an entire row is \code{NA}, the corresponding arguments do not
  lead to a valid distribution.
}

\examples{
ibeta(shape1 = 5, shape2 = 4)
ibeta(mean = seq(0.1,0.9,0.2), var = 0.03)
igamma(shape=c(38.9,7.1,47), scale=c(1.5,6.3,2.6))
}

\seealso{\code{\link{dbeta}}, \code{\link{dgamma}}}

\keyword{distribution}

