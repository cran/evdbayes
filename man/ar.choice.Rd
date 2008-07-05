\name{Accept Rate}
\alias{ar.choice}
\title{Compute Suited Proposal Standard Deviations}
\description{
  Compute suited proposal standard deviations for the MCMC algorithm. 
}
\usage{
ar.choice(init, prior, lh = c("none","gev","gpd","pp","os"), ..., psd,
ar = rep(.4, npar), n = 1000, tol = rep(.05, npar))
}
\arguments{
  \item{init}{a numeric vector for the starting value of the MCMC
    algorithm.}
  \item{prior}{A prior model. See function \code{\link{prior.prob}},
    \code{\link{prior.quant}}, \code{\link{prior.norm}} and
    \code{\link{prior.loglognorm}}.}
  \item{lh}{The likelihood function. Should be one of ``none'', ``gev'',
    ``gpd'', ``pp'' and ``os''.}
  \item{\dots}{Optional arguments to be passed to the
    \code{\link{posterior}} function.}
  \item{psd}{The initials proposal standard deviations.}
  \item{ar}{Optional. The objective accept rates - default is
    \code{rep(.4, npar)}.}
  \item{n}{Optional. The length of the simulated Markov Chains.}
  \item{tol}{Optional. The tolerance for the convergence test.}
}
\value{
  Return a list with two arguments. ``psd'': the suited proposal
  standard deviations and ``ar'': the accept rates related to these
  proposal standard deviations.
}
\details{
  The suited proposal standard deviations (psd) are computed through
  trial and error processes. Proposal standard deviations are
  fundamental to ensure good mixing properties for the Markov Chains.
  
  For this purpose, there exits a thumb rule:
  \emph{``In small dimensions, aim at an average acceptance rate of
    50. In large dimensions, at an average acceptance rate of 25. (Gelman
    et al., 1995)''}.

  For numerical conveniences, the trial and error process is more
  accurate with small initial starting psd.
}
\references{
  Gelman, A. and Roberts, G. and Gilks, W. (1995) \emph{Efficient
  Metropolis Jumping Rules.} Oxford University Press. 
}
\examples{
data(rainfall)
prrain <- prior.quant(shape = c(38.9, 7.1, 47), scale = c(1.5, 6.3,
  2.6))
n <- 10000; t0 <- c(43.2, 7.64, 0.32);
s <- ar.choice(init = t0, prior = prrain, lh = "pp", data = rainfall,
  thresh = 40, noy = 54, psd = rep(0.01, 3))
##Alec Stephenson choose suited psd = c(2, .2, .07) which is really
##close to values computed by ``ar.choice'' 
}
\author{Mathieu Ribatet}
\keyword{htest}
