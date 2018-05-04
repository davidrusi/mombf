\name{mixturebf-class}
\Rdversion{1.1}
\docType{class}
\alias{mixturebf-class}
\alias{mixturebf}
\alias{show,mixturebf-method}
\alias{coef.mixturebf}


\title{Class "mixturebf"}
\description{
Stores the output of Bayesian model selection for mixture models,
e.g. as produced by function \code{bfnormmix}.

Methods are provided for retrieving the posterior probability of a given
number of mixture components, posterior means and posterior samples of
the mixture model parameters.
}
\section{Objects from the Class}{
  Typically objects are automatically created by a call to \code{bfnormmix}.
}
\section{Slots}{
  The class has the following slots:
  \describe{

    \item{postprob}{data.frame containing posterior probabilities for
      different numbers of components (k) and log-posterior probability of
      a component being empty (contain no individuals)}

    \item{p}{Number of variables in the data to which the model was fit}

    \item{n}{Number of observations in the data to which the model was fit}

    \item{priorpars}{Prior parameters used when fitting the model}
    
    \item{postpars}{Posterior parameters for a 1-component mixture,
    e.g. for a Normal mixture the posterior is N(mu1,Sigma/prec)
    IW(nu1,S1)}
  
    \item{mcmc}{For each considered value of k, posterior samples for
    the parameters of the k-component model are stored}
  }
}
%  \S3method{coef}{mixturebf}(object, ...)
\section{Methods}{
  \describe{
    \item{coef}{Computes posterior means for all parameters}
    
    \item{show}{\code{signature(object = "mixturebf")}: Displays general
      information about the object. }
    
    \item{postProb}{\code{signature(object = "mixturebf")}: Extracts
      posterior model probabilities, Bayes factors and posterior
      probability of a cluster being empty}

    \item{postSamples}{\code{signature(object = "mixturebf")}: Extracts
      posterior samples}
  }
}
\references{
Fuquene J., Steel M.F.J., Rossell D. On choosing mixture components via
  non-local priors. 2018. arXiv
}
\author{
David Rossell
}
\seealso{
See also \code{\link{bfnormmix}}
}
\examples{
showClass("mixturebf")
}
\keyword{classes}
