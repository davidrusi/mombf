\name{dpostNIW}
\alias{dpostNIW}
\title{Posterior Normal-IWishart density}
\description{
  Posterior Normal-IWishart density evaluated at (mu,Sigma).
  This posterior corresponds to a Normal model for the data

  x[i] ~ N(mu,Sigma) iid i=1,...,n

  under conjugate priors

  mu | Sigma ~ N(mu0, g Sigma)
  Sigma^{-1} ~ W(nu0, S0)
}
\usage{
dpostNIW(mu, Sigma, x, g=1, mu0=rep(0,length(mu)), nu0=nrow(Sigma)+1, S0,
  logscale=FALSE)
}
\arguments{
  \item{mu}{Vector of length p}
  \item{Sigma}{p x p positive-definite covariance matrix}
  \item{x}{n x p data matrix (individuals in rows, variables in columns)}
  \item{g}{Prior dispersion parameter for mu}
  \item{mu0}{Prior mean for mu}
  \item{nu0}{Prior degrees of freedom for Sigma^{-1}}
  \item{S0}{Prior scale matrix for Sigma^{-1}, by default set to I/nu0}
  \item{logscale}{set to TRUE to get the log-posterior density}
}
\value{
Normal-IW posterior density evaluated at (mu,Sigma)
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

