\name{Jumps Definition}
\alias{jacFun}
\alias{movTyp1}
\alias{movTyp2}
\title{Functions to define reversible jumps}
\description{
  These are function usefull to define reversible jumps for the
  algorithm. There is currently three functions: jacFun which evaluates
  the jacobian, movType1 which returns a special proposal of type 1 and
  movTyp2 which returns a special proposal of type 2.
}
\usage{
jacFun(xi, pMass, cv)
movTyp1(prow, pMass, cv)
movTyp2(prow, propShape, cv)
}
\arguments{
  \item{xi}{The shape parameter at which the jacobian must be
    evaluated.}
  \item{pMass}{The point Mass i.e. the fixed value for the shape
    parameter. See reference guide.}
  \item{cv}{A numeric value for jumping rule. See reference guide.}
  \item{prow}{The current state of the Markov Chain.}
  \item{propShape}{The proposal shape candidate.}
}
\value{
  \code{movTyp1} and \code{movTyp2} returns a vector of lenght two. It
  corresponds to location and scale proposals.

  \code{jacFun} return the value of the Jacobian at value \code{xi}.
}
\details{
  These are special fucntions but can be overwritten to define new
  jumping rules.

  See user's guide.
}
\references{
  Green, P. (1995) Reversible jump Markov chain Monte Carlo computation and
  Bayesian model determination. Biometrica \bold{82}:711--732.
}
\author{Mathieu Ribatet}
\keyword{htest}
