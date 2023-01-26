\name{icfit-class}
\Rdversion{1.1}
\docType{class}
\alias{icfit-class}
\alias{icfit}
\alias{show,icfit-method}
\alias{icfit.coef}
\alias{icfit.predict}
\alias{icfit.summary}



\title{Class "icfit"}
\description{
Stores the output of the search for the model with best information
criterion value, e.g. produced by \code{bestBIC}, \code{bestBIC},
\code{bestAIC} or \code{bestIC}.
The class extends a list, so all usual methods for lists also work for
\code{icfit} objects, e.g. accessing elements, retrieving names etc.

Methods are provided to extract coefficients, predictions, confidence
intervals and summary information about the best model.
}
\section{Objects from the Class}{
  icfit objects are automatically created by a call to
  \code{bestBIC} or similar.
}
\section{Slots}{
  The class extends a list with elements:
  \describe{

  \item{topmodel}{names of the variables in the top model}

  \item{topmodel.fit}{top model as fitted by glm}

  \item{models}{data frame with the information criterion for all models
  (when enumeration is feasible) or those visited by an MCMC model
  search in modelSelection (when enumeration is not feasible)}

  \item{varnames}{the names of all variables in the design matrix}

  \item{msfit}{Output of modelSelection (used to search the top model)}

  }
}
\section{Methods}{
  \describe{
    \item{coef}{glm fit for the top model}
    \item{confint}{Confidence intervals under the top model}
    \item{predict}{Predictions for the top model (predict.glm)}
    \item{show}{\code{signature(object = "icfit")}: Displays general information about the object. }
  }
}
\author{
David Rossell
}
\seealso{
See also \code{\link{bestBIC}}.
}
\examples{
showClass("icfit")
}
\keyword{classes}
