\name{icov}
\alias{icov}
\title{ Extract estimated inverse covariance }
\description{
  
  Extract the estimated inverse covariance from an \code{msfit_ggm}
  object.

  Bayesian model averaging is used, optionally entries with posterior
  probability of being non-zero below a threshold are set to 0.

}
\usage{

icov(fit, threshold) 

}
\arguments{
  
  \item{fit}{Object of class \code{msfit_ggm}, returned by
    \code{modelSelectionGGM}}
  
  \item{threshold}{Entries with posterior probability of being non-zero
  below threshold are set to 0. If missing this argument is ignored and 
  no entries are set to exact zeroes. When the goal is to identify
  zeroes, a sensible default is \code{threshold=0.95}}

  
}
\details{

  The inverse covariance is obtained via Bayesian model averaging, using
  posterior samples of Omega. When \code{threshold} is specified,
  entries in the BMA estimate are set to zero, which may result in a non
  positive-definite matrix.
  
}
\value{

  Estimated inverse covariance matrix.

}
\author{ David Rossell }
\seealso{
  \code{\link{modelSelectionGGM}}, 
  \code{coef.msfit_ggm} for Bayesian model averaging estimates and
  intervals. 
}
\examples{
## See modelSelectionGGM
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
