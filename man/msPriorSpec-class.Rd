\name{msPriorSpec-class}
\Rdversion{1.1}
\docType{class}
\alias{msPriorSpec-class}

\title{Class "msPriorSpec"}
\description{
Stores the prior distributions to be used for Bayesian variable selection
in normal regression models.
This class can be used to specify the prior on non-zero regression coefficients,
the model indicator or the nuisance parameters.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("msPriorSpec", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{priorType}:}{Object of class \code{"character"}. \code{"coefficients"} indicates
that the prior is for the non-zero regression coefficients.
\code{"modelIndicator"} that it is for the model indicator,
and \code{"nuisancePars"} that it is for the nuisance parameteres.
Several prior distributions are available for each choice of \code{priorType},
and these can be speicified in the slot \code{priorDist}. }

    \item{\code{priorDistr}:}{Object of class \code{"character"}.
If \code{priorType=="coefficients"}, \code{priorDistr} can be equal to "pMOM", "piMOM", or "peMOM"
(product moment, product inverse moment or product exponential moment prior, respectively).
If \code{priorType=="modelIndicator"}, \code{priorDistr} can be equal to "uniform" or "binomial"
to specify a uniform prior (all models equaly likely a priori) or a binomial prior. For a binomial prior,
the prior inclusion probability for any single variable must be specified in slot \code{priorPars}.
If \code{priorType=="nuisancePars"}, \code{priorDistr} must be equal to "invgamma". This corresponds to an
inverse gamma distribution for the residual variance, with parameters
specified in the slot \code{priorPars}.}

    \item{\code{priorPars}:}{Object of class \code{"vector"}, where each element must be named.
For \code{priorDistr=='pMOM'}, there must be an element "r" (MOM power
is 2r).
For any \code{priorDistr} there must be either an element "tau" indicating
the prior dispersion or elements "a.tau" and "b.tau" specifying an
inverse gamma hyper-prior for "tau".
Optionally, there may be an element "tau.adj" indicating the prior
dispersion for the adjustment variables (i.e. not undergoing variable
selection). If not defined, "tau.adj" is set to 0.001 by default.
For \code{priorDistr=='binomial'}, there must be either an element "p" specifying the prior inclusion probability
for any single covariate, or a vector with elements "alpha.p" and
"beta.p" specifying a Beta(alpha.p,beta.p) hyper-prior on p.
For \code{priorDistr=='invgamma'} there must be elements "alpha" and "lambda". The prior for the residual variance
is an inverse gamma with parameteres \code{.5*alpha} and \code{.5*lambda}.
  }
  }
}
\section{Methods}{
No methods defined with class "msPriorSpec" in the signature.
}
\references{
Prior Densities for Default Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B, 2010, 72, 143-170
}
\author{
David Rossell.
}
\note{
When new instances of the class are created a series of check are performed to ensure that a valid prior
specification is produced.
}
\seealso{
See also \code{\link{modelSelection}} for an example of defining an instance of the class
and perform Bayesian model selection.
}
\examples{
showClass("msPriorSpec")
}
\keyword{classes}
