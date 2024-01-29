\name{msfit_ggm-class}
\Rdversion{1.1}
\docType{class}
\alias{msfit_ggm-class}
\alias{msfit_ggm}
\alias{show,msfit_ggm-method}
\alias{msfit_ggm.coef}


\title{Class "msfit_ggm"}
\description{
Stores the output of Bayesian Gaussian graphical model selection and
averaging, as produced by function \code{modelSelectionGGM}.
The class extends a list, so all usual methods for lists also work for
\code{msfit_ggm} objects, e.g. accessing elements, retrieving names etc.

Methods are provided to obtain parameter estimates, posterior intervals
(Bayesian model averaging), and posterior probabilities of parameters
being non-zero

}
\section{Objects from the Class}{
  Objects are created by a call to \code{modelSelectionGGM}.
}
\section{Slots}{
  The class extends a list with elements:
  \describe{

    \item{postSample}{Sparse matrix (\code{dgCMatrix}) with posterior
      samples for the Gaussian precision (inverse covariance)
      parameters. Each row is a posterior sample. Within each row, only
      the upper-diagonal of the precision matrix is stored in a flat
      manner. The row and column indexes are stored in indexes}

    \item{indexes}{For each column in postSample, it indicates the row
    and column of the precision matrix}

    \item{p}{Number of variables}
    
    \item{priors}{Priors specified when calling \code{modelSelection}}

  }
}
\section{Methods}{
  \describe{
    \item{coef}{Obtain BMA posterior means, intervals and posterior
      probability of non-zeroes}

    \item{plot}{Shows estimated posterior inclusion probability for each
      parameter vs. number of MCMC iterations. Only up to the first 5000
      parameters are shown}

    \item{show}{\code{signature(object = "msfit_ggm")}:
      Displays general information about the object. }
  }
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelectionGGM}}
}
\examples{
showClass("msfit_ggm")
}
\keyword{classes}
