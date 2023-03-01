\name{rl.pred}
\alias{rl.pred}

\title{Return Level Plots for GEV Predictive Distributions}
\description{
  Produce return level plots depicting prior and posterior predictive
  gev distributions.
}

\usage{
rl.pred(post, qlim, npy, lh = c("gev", "gpd"), period = 1, lty = 1, col = 1,
xlab = "return period", ylab = "return level", \dots)
}

\arguments{
  \item{post}{A Markov chain generated using \code{posterior},
    containing samples from the corresponding prior/posterior
    distribution.}
  \item{qlim}{A vector of length two, giving the limits for
    the quantiles at which the predictive probabilities are
    calculated.}
  \item{npy}{The Number of observation Per Year (in average). If ``gev''
  likelihood, ``npy'' is supposed to be equal to 1 i.e. annual maxima.} 
  \item{lh}{The likelihood.} 
  \item{period}{A vector of integers. One curve is plotted for
    each element of \code{period}. The \eqn{i}th curve depicts
    the probabilities that that quantiles will be exceeded over
    the next \code{period[i]} periods.}
  \item{lty}{Passed to \code{matplot}.}
  \item{col}{Passed to \code{matplot}.}
  \item{xlab, ylab}{Labels for the x and y axes.}
  \item{\dots}{Other arguments passed to \code{matplot}.}
}

\value{
  The first two arguments to \code{matplot} are returned invisibly
  as a list.  
  
  If a linear trend on the location has been implemented, the plot
  corresponds to the distribution obtained when the trend parameter is
  zero.
}
\details{
  See the user's guide. 
}

\seealso{\code{\link{matplot}}, \code{\link{posterior}}}

\keyword{hplot}

