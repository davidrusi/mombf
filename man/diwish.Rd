\name{diwish}
\alias{diwish}
\title{Density for Inverse Wishart distribution}
\description{
  \code{diwish} returns the density for the inverse Wishart(nu,S)
  evaluated at Sigma.
}
\usage{
diwish(Sigma, nu, S, logscale=FALSE)
}
\arguments{
  \item{Sigma}{Positive-definite matrix}
  \item{nu}{ Degrees of freedom of the inverse Wishart}
  \item{S}{ Scale matrix of the inverse Wishart}
  \item{logscale}{ If \code{logscale==TRUE} the log-density is returned}
}
\value{
  Inverse Wishart(nu,S) density evaluated at Sigma
}
\author{ David Rossell }
\seealso{ \code{\link{dpostNIW}} for the Normal-IW posterior density}
\examples{
library(mombf)
Sigma= matrix(c(2,1,1,2),nrow=2)
diwish(Sigma,nu=4,S=diag(2))
}
\keyword{ distribution }

