\name{msfit-class}
\Rdversion{1.1}
\docType{class}
\alias{msfit-class}
\alias{msfit}
\alias{show,msfit-method}
\alias{msfit.coef}
\alias{msfit.plot}
\alias{msfit.predict}
\alias{coefByModel}
\alias{coefByModel-methods}
\alias{coefByModel,msfit-method}


\title{Class "msfit"}
\description{
Stores the output of Bayesian variable selection, as produced by
function \code{modelSelection}.
The class extends a list, so all usual methods for lists also work for
\code{msfit} objects, e.g. accessing elements, retrieving names etc.

Methods are provided to compute posterior probabilities,
obtaining regression coefficient estimates and posterior intervals
(both via Bayesian model averaging and for individual models),
and sampling from their posterior distribution, as indicated below.
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
    \item{family}{Residual distribution, i.e. argument \code{family}
      when calling \code{modelSelection}}
    \item{p}{Number of variables}
    \item{priors}{Priors specified when calling \code{modelSelection}}
    \item{ystd}{For internal use. Stores the response variable,
      standardized if \code{center} or \code{scale} were set to
      \code{TRUE}}
    \item{xstd}{For internal use. Stores the covariates,
      standardized if \code{center} or \code{scale} were set to
      \code{TRUE}}
    \item{stdconstants}{For internal use. If \code{center} or
      \code{scale} were set to \code{TRUE}, stores the sample mean
      and standard deviation of the outcome and covariates}
    \item{call}{Stores info about the call, the formula used (if any),
      splines used etc}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{Obtains posterior means and intervals via Bayesian model
      averaging}
    \item{coefByModel}{Obtains posterior means and intervals for
      individual models}
    \item{plot}{Shows estimated posterior inclusion probability for each
      parameter vs. number of MCMC iterations}
    \item{predict}{Obtains posterior means and intervals for given
      covariate values. These are posterior intervals for the mean, not
      posterior predictive intervals for the outcome}
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
