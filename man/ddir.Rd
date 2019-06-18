\name{ddir}
\alias{ddir}
\title{ Dirichlet density }
\description{
  Evaluate the density of a Dirichlet distribution
}
\usage{
ddir(x, q, logscale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Vector or matrix containing the value at which to evaluate
    the density. If a matrix, the density is evaluated for each
    row. Rows are renormalized to ensure they add up to 1}
  \item{q}{ Dirichlet parameters. Must have the same length as
    \code{ncol(x)}, or length 1 (in which case a symmetric Dirichlet
    density is valuated)}
  \item{logscale}{ For \code{logscale==TRUE}, \code{dimom} returns the
    natural log of the prior density}
}
\value{
Density of a Dirichlet(q) distribution evaluated at each row of \code{x}
}
\author{ David Rossell }
\examples{
library(mombf)
x= matrix(c(1/3,2/3,.5,.5),nrow=2,byrow=TRUE)
ddir(x,q=2)
}
\keyword{ distribution }

