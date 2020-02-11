\name{postProb}
\alias{postProb}
\alias{postProb-methods}
\alias{postProb,msfit-method}
\alias{postProb,mixturebf-method}

\title{
Obtain posterior model probabilities
}
\description{
  Obtain posterior model probabilities after running Bayesian model selection
}
\usage{
  postProb(object, nmax, method='norm')
}
\arguments{
  \item{object}{Object of class msfit returned by \code{modelSelection}
    or class mixturebf, e.g. returned by \code{bfnormmix}}
  \item{nmax}{Maximum number of models to report (defaults to no max)}
  \item{method}{Only when \code{class(object)} is msfit.
    For 'norm' probabilities are obtained by renormalizing the
  stored integrated likelihoods, for 'exact' they are given by the proportion
  of MCMC visits to each model. 'norm' has less variability but can be biased
  if the chain has not converged.}
}
\value{
  A \code{data.frame} with posterior model probabilities in column pp.
  Column modelid indicates the indexes of the selected covariates (empty
  for the null model with no covariates).
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelection}} to perform model selection
}
\examples{
#See help(modelSelection)
}
\keyword{ models }
\keyword{ distribution }
