\name{msfit-class}
\Rdversion{1.1}
\docType{class}
\alias{msfit-class}
\alias{msfit}
\alias{show,msfit-method}

\title{Class "msfit"}
\description{
Stores the output of Bayesian variable selection, as produced by
function \code{modelSelection}.
The class extends a list, so all usual methods for lists also work for
\code{msfit} objects, e.g. accessing elements, retrieving names etc.

Some additional methods are provided for printing information on screen,
computing posterior probabilities or sampling from the posterior of
regression coefficients, as indicated below.
}
\section{Objects from the Class}{
  Typically objects are automatically created by a call to \code{modelSelection}.
  Alternatively, objects can be created by calls of the form
  \code{new("msfit",x)} where \code{x} is a list with the adequate
  elements (see slots).
}
\section{Slots}{
  The class extends a list with elements:
  \describe{

    \item{postSample}{\code{matrix} with posterior samples for the model
      indicator. \code{postSample[i,j]==1}
      indicates that variable j was included in the model in the MCMC
      iteration i}
    \item{postOther}{\code{postOther}
      returns posterior samples for parameters other than the model
      indicator, i.e. basically hyper-parameters.
      If hyper-parameters were fixed in the model specification, \code{postOther} will be empty.}
    \item{margpp}{Marginal posterior probability for inclusion of each
      covariate. This is computed by averaging marginal post prob for
      inclusion in each Gibbs iteration, which is much more accurate than
      simply taking \code{colMeans(postSample)}}.
    \item{postMode}{Model with highest posterior probability amongst all those visited}
    \item{postModeProb}{Unnormalized posterior prob of posterior mode (log scale)}
    \item{postProb}{Unnormalized posterior prob of each visited model (log
      scale)}
    \item{coef}{Estimated coefficients (via posterior mode) for highest
      posterior probability model}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "msfit")}: Displays general information about the object. }
    \item{postProb}{\code{signature(object = "msfit")}: Extracts
      posterior model probabilities. }
    \item{rnlp}{\code{signature(object = "msfit")}: Obtain posterior
      samples for regression coefficients. }
  }
}
\references{
Johnson VE, Rossell D. Non-Local Prior Densities for Default Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B, 2010, 72, 143-170

Johnson VE, Rossell D. Bayesian model selection in high-dimensional
settings. Journal of the American Statistical Association, 107, 498:649-660.
}
\author{
David Rossell
}
\seealso{
See also \code{\link{modelSelection}} and \code{\link{rnlp}}.
}
\examples{
showClass("msfit")
}
\keyword{classes}
