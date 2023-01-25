\name{getBIC}
\alias{getAIC}
\alias{getAIC-methods}
\alias{getAIC,msfit-method}
\alias{getBIC}
\alias{getBIC-methods}
\alias{getBIC,msfit-method}
\alias{getEBIC}
\alias{getEBIC-methods}
\alias{getEBIC,msfit-method}
\alias{getIC}
\alias{getIC-methods}
\alias{getIC,msfit-method}

\title{
Obtain AIC, BIC, EBIC or other general information criteria (getIC)
}
\description{
  Extract information criteria from an msfit object.
}
\usage{
  getAIC(object)

  getBIC(object)

  getEBIC(object)

  getIC(object)
}
\arguments{
  \item{object}{Object of class msfit returned by \code{modelSelection}}
}
\value{
  BIC or EBIC values for all models enumerated / visited by modelSelection
}
\details{
  Let p be the total number of parameters and n the sample size. The BIC of a model k with p_k
  parameters is

  - 2 L_k + p_k log(n)

  the AIC is
  
  - 2 L_k + p_k 2
  
  the EBIC is

  - 2 L_k + p_k log(n) + 2 log(p choose p_k)

  and a general information criterion with a given model size penalty

  - 2 L_k + p_k penalty


  For getBIC() to work, object must be the result returned by
  modelSelection setting priorCoef=bic() and priorDelta=modelunifprior()

  For getEBIC() it is priorCoef=bic() and priorDelta=modelbbprior()

  For getAIC() it is priorCoef=aic() and priorDelta=modelunifprior()

  For getIC() it is priorCoef=ic() and priorDelta=modelunifprior()
  
  Function modelSelection returns the log posterior probability of a
  model, postProb = log(m_k) + log(prior k), where m_k is the marginal
  likelihood of the model and prior k its prior probability.
  
  When running function modelSelection with priorCoef=bicprior()
  and priorDelta=modelunifprior(), the BIC approximation is used for
  m_k, that is
  
  log(m_k) = L_k - 0.5 * p_k log(n)

  and all models are equally likely a priori, log(prior k)= p
  log(1/2). Then the BIC can be easily recovered

  BIC_k= -2 * [postProb + p log(2)]

  When using priorCoef=bicprior() and priorDelta=modelbbprior(),
  log(prior k)= - log(p+1) - log(p choose p_k), hence

  EBIC_k= -2 * [postProb + log(p+1)].
 
}
\author{
David Rossell
}
\seealso{
\code{\link{modelSelection}} to perform model selection
}
\examples{
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x \%*\% theta + rnorm(100)
ybin <- y>0

#Obtain BIC
ms= modelSelection(ybin, x=x, priorCoef=bicprior(),
priorDelta=modelunifprior(), family='binomial')
getBIC(ms)

#Obtain EBIC
ms2= modelSelection(ybin, x=x, priorCoef=bicprior(),
priorDelta=modelbbprior(), family='binomial')
getEBIC(ms2)

}
\keyword{ models }
\keyword{ distribution }
