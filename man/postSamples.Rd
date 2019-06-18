\name{postSamples}
\alias{postSamples}
\alias{postSamples-methods}
\alias{postSamples,mixturebf-method}

\title{
Extract posterior samples from an object
}
\description{
  Obtain posterior model probabilities after running Bayesian model selection
}
\usage{
  postSamples(object)
}
\arguments{
  \item{object}{Object containing posterior samples, e.g. of class
    mixture bf as returned by \code{bfnormmix}}
}
\value{
  For objects of class \code{mixturebf}, a list with one element for
  each considered number of mixture components.

  Each element in the list
  contains posterior samples on the mixture weights (eta) and other
  component-specific parameters such as means (mu) and Cholesky
  decomposition of the inverse covariance matrix (cholSigmainv)
}
\author{
David Rossell
}
\examples{
#See help(bfnormmix)
}
\keyword{ models }
\keyword{ distribution }
