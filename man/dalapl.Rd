\name{dalapl}
\alias{dalapl}
\alias{palapl}
\alias{ralapl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density and random draws from the asymmetric Laplace distribution }
\description{
  \code{dalapl} evaluates the probability density function,
  \code{palapl} the cumulative probability function
  and \code{ralapl} generates random draws.
}
\usage{
dalapl(x, th=0, scale=1, alpha=0, logscale=FALSE)

palapl(x, th=0, scale=1, alpha=0)

ralapl(n, th=0, scale=1, alpha=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of values at which to evaluate the pdf/cdf}
  \item{n}{ Number of random draws}
  \item{th}{ Location parameter (mode) }
  \item{scale}{ Scale parameter (proportional to variance)}
  \item{alpha}{ Asymmetry parameter, must be between -1 and 1}
  \item{logscale}{If TRUE the log-pdf is returned}
}
\details{
  For x<=th the asymmetric Laplace pdf is

  0.5*exp(-abs(th-x)/(sqrt(scale)*(1+alpha)))/sqrt(scale)

  and for x>th it is

  0.5*exp(-abs(th-x)/(sqrt(scale)*(1-alpha)))/sqrt(scale)
  }
\value{
  \code{dalapl} returns the density function,
  \code{palapl} the cumulative probability,
  \code{ralapl} random draws.
}
\author{ David Rossell }
\examples{
library(mombf)
e <- ralapl(n=10^4, th=1, scale=2, alpha=0.5)
thseq <- seq(min(e),max(e),length=1000)
hist(e, main='', breaks=30, prob=TRUE)
lines(thseq, dalapl(thseq, th=1, scale=2, alpha=0.5), col=2)
}
\keyword{ distribution }

