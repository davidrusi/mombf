\name{dpostNIW}
\alias{dpostNIW}
\alias{rpostNIW}
\title{Posterior Normal-IWishart density}
\description{
  dpostNIW evalutes the posterior Normal-IWishart density at (mu,Sigma).
  rpostNIW draws independent samples.
  This posterior corresponds to a Normal model for the data

  x[i,] ~ N(mu,Sigma) iid i=1,...,n

  under conjugate priors

  mu | Sigma ~ N(mu0, g Sigma)
  Sigma      ~ IW(nu0, S0)
}
\usage{
dpostNIW(mu, Sigma, x, g=1, mu0=rep(0,length(mu)), nu0=nrow(Sigma)+1, S0,
  logscale=FALSE)

rpostNIW(n, x, g=1, mu0=0, nu0, S0, precision=FALSE)
}
\arguments{
  \item{mu}{Vector of length p}
  \item{Sigma}{p x p positive-definite covariance matrix}
  \item{x}{n x p data matrix (individuals in rows, variables in columns)}
  \item{g}{Prior dispersion parameter for mu}
  \item{mu0}{Prior mean for mu}
  \item{nu0}{Prior degrees of freedom for Sigma}
  \item{S0}{Prior scale matrix for Sigma, by default set to I/nu0}
  \item{logscale}{set to TRUE to get the log-posterior density}
  \item{n}{Number of samples to draw}
  \item{precision}{If set to \code{TRUE}, samples from the precision
    matrix (inverse of Sigma) are returned instead}
}
\value{
\code{dpostNIW} returns the Normal-IW posterior density evaluated at
(mu,Sigma).

\code{rpostNIW} returns a list with two elements. The first element are
posterior draws for the mean. The second element are posterior draws for
the covariance (or its inverse if \code{precision==TRUE}). Only
lower-diagonal elements are returned (\code{Sigma[lower.tri(Sigma,diag=TRUE)]}).
}
\author{ David Rossell }
\seealso{ \code{\link{diwish}} for the inverse Wishart prior density,
\code{\link{marginalNIW}} for the integrated likelihood under a
  Normal-IW prior}
\examples{
#Simulate data
x= matrix(rnorm(100),ncol=2)
#Evaluate posterior at data-generating truth
mu= c(0,0)
Sigma= diag(2)
dpostNIW(mu,Sigma,x=x,g=1,nu0=4,log=FALSE)
}
\keyword{ distribution }

