\name{rl.pst}
\alias{rl.pst}

\title{Return Level Plots Depicting Distributions of GEV Quantiles}
\description{
  Produce return level plots depicting prior and posterior
  distributions of gev quantiles.
}

\usage{
rl.pst(post, npy, lh = c("gev", "gpd"), ci = 0.9, lty = c(2,1), col = c(2,1),
xlab = "return period", ylab = "return level",  \dots)
}

\arguments{
  \item{post}{A Markov chain generated using \code{posterior},
    containing samples from the corresponding prior/posterior
    distribution.}
  \item{npy}{The Number of observation Per Year (in average). If ``gev''
  likelihood, ``npy'' is supposed to be equal to 1 i.e. annual maxima.}
  \item{lh}{The likelihood.}
  \item{ci}{The confidence coefficient for the plotted
    prior/posterior probability interval.}
  \item{lty}{Passed to \code{matplot}. The first and
    second values specify the line type of the probability
    interval and the median line respectively.}
  \item{col}{Passed to \code{matplot}. The first and
    second values specify the colour of the probability
    interval and the median line respectively.}
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

