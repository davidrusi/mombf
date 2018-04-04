\name{marginalNIW}
\alias{marginalNIW}
\alias{marginalNIW-methods}
\alias{marginalNIW,missing,ANY,matrix,numeric,missing-method}
\alias{marginalNIW,matrix,missing,missing,missing,missing-method}
\alias{marginalNIW,matrix,missing,missing,missing,vector-method}
\alias{marginalNIW,missing,list,list,numeric,missing-method}

\title{
Marginal likelihood under a multivariate Normal likelihood and a conjugate
Normal-inverse Wishart prior.
}
\description{
The argument \code{z} can be used to specify cluster allocations. If
left missing then the usual marginal likelihood is computed, else it is
computed conditional on the clusters (this is equivalent to the product
of marginal likelihoods across clusters)
}
\usage{
marginalNIW(x, xbar, samplecov, n, z, g,  mu0=rep(0,ncol(x)),
nu0=ncol(x)+4, S0, logscale=TRUE)
}
\arguments{
  \item{x}{Data matrix (individuals in rows, variables in
    columns). Alternatively you can leave missing and specify
    \code{xbar}, \code{samplecov} and \code{n} instead}
  \item{xbar}{Either a vector with column means of \code{x} or a list
    where each element corresponds to the column means for each cluster}
  \item{samplecov}{Either the sample covariance matrix \code{cov(x)} or
    a list where each element contains the covariance for each clsuter}
  \item{n}{Either an integer indicating the sample size \code{nrow(x)}
    or a vector indicating the cluster counts \code{table(z)}}
  \item{z}{Optional argument specifying cluster allocations}
  \item{g}{Prior dispersion parameter for mu}
  \item{mu0}{Prior mean for mu}
  \item{nu0}{Prior degrees of freedom for Sigma}
  \item{S0}{Prior scale matrix for Sigma, by default set to I/nu0}
  \item{logscale}{set to TRUE to get the log-posterior density}
}
\value{
If \code{z} is missing the integrated likelihood under a Normal-IW
prior. If \code{z} was specified then the product of integrated
likelihoods across clusters
}
\details{
  The function computes

  p(x)= int p(x | mu,Sigma) p(mu,Sigma) dmu dSigma

  where p(x[i])= N(x[i]; mu,Sigma) iid i=1,...,n

  p(mu | Sigma)= N(mu; mu0, g Sigma)
  p(Sigma)=     IW(Sigma; nu0, S0)
}
\author{
David Rossell
}
\seealso{
\code{\link{dpostNIW}} for the posterior Normal-IW density.
}
\examples{
#Simulate data
x= matrix(rnorm(100),ncol=2)

#Integrated likelihood under correct model
marginalNIW(x,g=1,nu0=4,log=FALSE)

#Integrated likelihood under random cluster allocations
z= rep(1:2,each=25)
marginalNIW(x,z=z,g=1,nu0=4,log=FALSE)
}
\keyword{ models }
\keyword{ htest }
