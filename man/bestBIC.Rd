\name{bestBIC}
\alias{bestAIC}
\alias{bestBIC}
\alias{bestEBIC}
\alias{bestIC}

\title{
Model with best AIC, BIC, EBIC or other general information criteria (getIC)
}
\description{
  Search for the regression model attaining the best value of the
  specified information criterion
}
\usage{
  bestAIC(...)

  bestBIC(...)

  bestEBIC(...)

  bestIC(..., penalty)
}
\arguments{
  \item{...}{Arguments passed on to \code{modelSelection}. The first and
    main argument is a model formula, see the examples}
  \item{penalty}{General information penalty. For example, since the AIC
    penalty is 2, bestIC(...,penalty=2) is the same as bestAIC(...)}
}
\value{Object of class \code{icfit}. Use (coef, summary,
  confint, predict) to get inference for the top model,
  and \code{help(icfit-class)} for more details on the returned object.
}
\details{
  For details on the information criteria see help(getBIC).

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

#BIC for all models (the intercept is also selected in/out)
fit= bestBIC(y ~ x[,1] + x[,2])
fit

#Same, but setting the BIC's log(n) penalty manually
#change the penalty for other General Info Criteria
#n= nrow(x)
#fit= bestIC(y ~ x[,1] + x[,2], penalty=log(n))

summary(fit) #usual GLM summary

coef(fit) #MLE under top model

#confint(fit) #conf int under top model (requires MASS package)



}
\keyword{ models }
\keyword{ distribution }
