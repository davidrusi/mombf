\name{plotprior}
\alias{plotprior}
\alias{plotprior-methods}
\alias{plotprior,cilfit-method}

\title{
Plot estimated marginal prior inclusion probabilities
}
\description{
  Plot marginal prior inclusion probabilities as estimated by cil versus
  regression coefficients for the treatment(s) equation(s)
}
\usage{
  plotprior(object, xlab, ylab, ylim=c(0,1), ...)
}
\arguments{
  \item{object}{Object of class cilfit returned by \code{cil}}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{ylim}{y-axis limits}
  \item{...}{Other arguments passed on to \code{plot}}
}
\value{
  A plot of prior inclusion probabilities vs treatment regression
  coefficients (dots). The line shows the (empirical Bayes) fit
}
\author{
David Rossell
}
\seealso{
\code{\link{cil}} 
}
\examples{
#See help(cil)
}
\keyword{ models }
\keyword{ hplot }
