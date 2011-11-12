\name{bbPrior}
\alias{bbPrior}
\alias{unifPrior}
\alias{binomPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Priors on model space for variable selection problems }
\description{
  \code{unifPrior}
  implements a uniform prior (equal a priori probability for all
  models). \code{binomPrior} implements a Binomial prior.
  \code{bbPrior} implements a Beta-Binomial prior. 
}
\usage{
unifPrior(sel, logscale = TRUE)
binomPrior(sel, prob = 0.5, logscale = TRUE)
bbPrior(sel, alpha = 1, beta = 1, logscale = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sel}{ Logical vector indicating which variables are included in
    the model}
  \item{logscale}{ Set to \code{TRUE} to return the log-prior
    probability.}
  \item{prob}{ Success probability for the Binomial prior}
  \item{alpha}{ First parameter of the Beta-Binomial prior, which is equivalent
    to specifying a Beta(alpha,beta) prior on \code{prob}.}
  \item{beta}{ First parameter of the Beta-Binomial prior, which is equivalent
    to specifying a Beta(alpha,beta) prior on \code{prob}.}
}
\value{
  Prior probability of the specified model
}
\author{ David Rossell }
\examples{
library(mombf)
sel <- c(TRUE,TRUE,FALSE,FALSE)
unifPrior(sel,logscale=FALSE)
binomPrior(sel,prob=.5,logscale=FALSE)
bbPrior(sel,alpha=1,beta=1,logscale=FALSE)
}
\keyword{ distribution }
