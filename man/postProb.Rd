\name{postProb}
\alias{postProb}
\alias{postProb-methods}
\alias{postProb,msfit-method}

\title{
Obtain posterior model probabilities
}
\description{
  Obtain posterior model probabilities after running Bayesian model selection
}
\usage{
  postProb(object, nmax)
}
\arguments{
  \item{object}{Object of class msfit, e.g. as returned by \code{modelSelection}}
  \item{nmax}{Maximum number of models to report (defaults to no max)}
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
