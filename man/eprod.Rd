\name{eprod}
\alias{eprod}
\title{ Expectation of a product of powers of Normal or T random
  variables }
\description{
  Compute the mean of prod(x)^power when x follows T_dof(mu,sigma)
  distribution (dof= -1 for multivariate Normal).

}
\usage{
eprod(m, S, power = 1, dof = -1) 
}
\arguments{
  \item{m}{Location parameter}
  \item{S}{Scale matrix. For multivariate T with dof>2 the covariance is
    S*dof/(dof-2). For the multivariate Normal the covariance is S.}
  \item{power}{Power that the product is raised to}
  \item{dof}{Degrees of freedom of the multivariate T. Set to -1 for the
  multivariate Normal.}
}
\details{
  The calculation is based on the computationally efficient approach by Kan (2008).
}
\value{
  Expectation of the above-mentioned product
}
\references{
Kan R. From moments of sum to moments of product. Journal of
  Multivariate Analysis 99 (2008), 542-554.
}
\author{ John Cook }
\examples{
#Check easy independence case
m <- c(0,3); S <- matrix(c(2,0,0,1),ncol=2)

eprod(m, S, power=2)

(m[1]^2+S[1][1])*(m[2]^2+S[2][2])
}
\keyword{ models }
\keyword{ distrib }

