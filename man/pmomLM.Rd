\name{pmomLM}
\alias{pmomLM}
\alias{pmomPM}
\alias{emomLM}
\alias{emomPM}
\alias{pplPM}
\alias{ppmodel}
\title{ Bayesian variable selection and model averaging for linear
  and probit models via non-local priors. }
\description{
  Variable selection for linear and probit models, providing a sample from the
  joint posterior of the model and regression coefficients.
  \code{pmomLM} and \code{pmomPM} implement product Normal MOM and 
  heavy-tailed product MOM as prior distribution for linear and probit
  model coefficients (respectively).
  \code{emomLM} and \code{emomPM} set an eMOM prior.

  \code{pplPM} finds the value of the prior dispersion parameter \code{tau} 
  minimizing posterior expected predictive loss (Gelfand and Ghosh, 1998) for the Probit model, 
  i.e. can be used to automatically set up \code{tau}.

  \code{ppmodel} returns the proportion of visits to each model.
}
\usage{
pmomLM(y, x, xadj, center=FALSE, scale=FALSE, niter=10^4, thinning=1, 
burnin=round(niter/10), priorCoef, priorDelta, priorVar, initSearch='greedy', 
verbose=TRUE)
pmomPM(y, x, xadj, niter=10^4, thinning=1, burnin=round(niter/10),
priorCoef, priorDelta, initSearch='greedy', verbose=TRUE)

emomLM(y, x, xadj, center=FALSE, scale=FALSE, niter=10^4, thinning=1, 
burnin=round(niter/10), priorCoef, priorDelta, priorVar, initSearch='greedy', 
verbose=TRUE)
emomPM(y, x, xadj, niter=10^4, thinning=1, burnin =round(niter/10),
priorCoef, priorDelta, initSearch='greedy', verbose=TRUE)

pplPM(tauseq=exp(seq(log(.01),log(2),length=20)), kPen=1, y, x, xadj, niter=10^4,
thinning=1, burnin=round(niter/10), priorCoef, priorDelta, priorVar,
initSearch='greedy', mc.cores=1)

ppmodel(nlpfit)
}
\arguments{
\item{y}{Vector with observed responses. For \code{pmomLM} this must be
  a numeric vector. For \code{pmomPM} it can either be a logical vector,
a factor with 2 levels or a numeric vector taking only two distinct values.}
\item{x}{Design matrix with all potential predictors which are to
  undergo variable selection.}
\item{xadj}{Design matrix for adjustment covariates, i.e. variables
  which are included in the model with probability 1. For instance,
  \code{xadj} can be used to force the inclusion of an intercept in the model.}
\item{center}{If \code{center==TRUE}, \code{y} and \code{x} are centered to have zero mean, therefore eliminating the need to include an intercept term in x}
\item{scale}{If \code{scale==TRUE}, \code{y} and columns in \code{x} are scaled to have standard deviation 1}
\item{niter}{Number of MCMC sampling iterations}
\item{thinning}{MCMC thinning factor, i.e. only one out of each \code{thinning} iterations are reported. Defaults to thinning=1, i.e. no thinning}
\item{burnin}{Number of burn-in MCMC iterations. Defaults to
  \code{.1*niter}. Set to 0 for no burn-in}
\item{priorCoef}{Prior distribution for the coefficients. Must be object of class \code{msPriorSpec} with slot \code{priorType} set to 'coefficients'. Possible values for slot \code{priorDistr} are 'pMOM', 'piMOM' and 'peMOM'}
\item{priorDelta}{Prior on model indicator space. Must be object of
  class \code{msPriorSpec} with slot \code{priorType} set to
  'modelIndicator'. Possible values for slot \code{priorDistr} are
  'uniform' and 'binomial'. For 'binomial', you can either set the prior
probability 'p' or specify a Beta-binomial prior by specifying the
parameters 'alpha.p','beta.p'.}
\item{priorVar}{Prior on residual variance. Must be object of class \code{msPriorSpec} with slot \code{priorType} set to 'nuisancePars'. Slot \code{priorDistr} must be equal to 'invgamma'}
\item{initSearch}{Algorithm to refine
  \code{deltaini}. \code{initSearch=='greedy'} uses a greedy Gibbs
  sampling search. \code{initSearch=='SCAD'} sets \code{deltaini} to the
  non-zero elements in a SCAD fit with cross-validated regularization
  parameter. \code{initSearch=='none'} initializes to the null model
  with no variables in \code{x} included.}
\item{verbose}{Set \code{verbose==TRUE} to print iteration progress}
\item{tauseq}{Grid of \code{tau} values for which the posterior predictive loss
  should be evaluated.}
\item{kPen}{Penalty term specifying the relative importance of
  deviations from the observed data vs deviation from the posterior
  predictive. \code{kPen} can be set either to a
  numeric value or to \code{'msize'} to set penalty equal to the average
  model size. Loss is Dev(yp,yhat) +
  kPen*Dev(yhat,yobs), where yp: draw from post predictive, yobs:
  observed data and yhat is E(yp|yobs).
}
\item{mc.cores}{Allows for parallel computing. \code{mc.cores} is the number of processors to use. 
  Setting \code{mc.cores>1} requires the \code{parallel} package.}
\item{nlpfit}{Non-local prior model fit, as returned by \code{pmomLM},
  \code{pmomPM}, \code{emomLM} or \code{emomPM}.}
}
\details{
  The implemented MCMC scheme makes proposals from the joint posterior of
  (delta[i],theta[i]) given all other parameters and the data, where
  delta[i] is the indicator for inclusion/exclusion of covariate i and
  theta[i] is the coefficient value.
  In contrast with some model fitting options implemented in
  \code{modelSelection}, here the scheme is exact. However, sampling the
  coefficients can adversely affect the mixing when covariates are
  very highly correlated.
  In practice, the mixing seems to be reasonably good for correlations
  up to 0.9.

  \code{pmomPM} uses the scheme of Albert & Chib (1993) for probit models.
}
\value{
  \code{pmomLM} and \code{pmomPM} returns a list with elements
  
\item{postModel}{\code{matrix} with posterior samples for the model
  indicator. \code{postModel[i,j]==1}
  indicates that variable j was included in the model in the MCMC
  iteration i}
\item{postCoef1}{\code{matrix} with posterior samples for coefficients
  associated to \code{x}}
\item{postCoef2}{\code{matrix} with posterior samples for coefficients
  associated to \code{xadj}}
\item{postPhi}{\code{vector} with posterior samples for residual variance}
\item{postOther}{\code{postOther}
  returns posterior samples for other parameters, i.e. basically hyper-parameters.
  Currently the prior precision parameter \code{tau}}
\item{margpp}{Marginal posterior probability for inclusion of each
  covariate. This is computed by averaging marginal post prob for
  inclusion in each MCMC iteration, which is much more accurate than
  simply taking \code{colMeans(postModel)}}.

\code{pplPM} returns a list with elements

\item{optfit}{Probit model fit using \code{tauopt}. It is the result of a call to \code{pmomPM}.}
\item{PPL}{\code{data.frame} indicating for each value in \code{tauseq} the posterior predictive loss (\code{PPL}=\code{G+P}),
the goodness-of-fit (\code{G}) and penalty terms (\code{P})}, the
average number of covariates in the model (\code{msize}) including \code{xadj}
and the smoothed \code{sPPL} obtained via a \code{gam} fit.
\item{tauopt}{Value of \code{tau} minimizing the PPL}
}
\references{
  Johnson V.E., Rossell D. Non-Local Prior Densities for Default
  Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B,
  2010, 72, 143-170.
  
  Johnson V.E., Rossell D. Bayesian model selection in high-dimensional
  settings. Technical report. 2011
  See http://rosselldavid.googlepages.com for technical reports.
  
  Albert, J. and Chib, S. (1993) Bayesian analysis of binary and
  polychotomous response data. Journal of the American Statistical
  Association, 88, p669-679

  Gelfand, A. and Ghosh, S. (1998) Model choice: A minimum posterior predictive loss approach.
  Biometrika, 85, p1-11.
}
\author{ David Rossell, Donatello Telesca }
\seealso{
 For more details on the prior specification see \code{\link{msPriorSpec-class}}
 To compute marginal densities for a given model see
  \code{\link{pmomMarginalK}}, \code{\link{pmomMarginalU}},
  \code{\link{pimomMarginalK}}, \code{\link{pimomMarginalU}}.
}
\examples{
#Simulate data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
xadj <- rep(1,nrow(x))
theta <- matrix(c(1,1,0),ncol=1)
y <- 10*xadj + x \%*\% theta + rnorm(100)

#Beta-binomial prior on model space
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)

#Non-informative prior on residual variance
priorVar <- igprior(alpha=.01,lambda=.01)

#Product MOM prior with tau=0.3 on x coefficients
#Non-informative prior on xadj coefficients
priorCoef <- momprior(tau=0.3, tau.adj=10^6)

mom0 <- pmomLM(y=y,x=x,xadj=xadj,center=FALSE,scale=FALSE,niter=1000,
priorCoef=priorCoef,priorDelta=priorDelta,priorVar=priorVar)
round(colMeans(mom0$postModel),2)
round(colMeans(mom0$postCoef1),2)
round(colMeans(mom0$postCoef2),2)

#Alternative prior: hyper-prior on tau
priorCoef <- new("msPriorSpec",priorType='coefficients',priorDistr='pMOM',
priorPars=c(a.tau=1,b.tau=.135,tau.adj=10^6,r=1)) #hyper-prior
mom1 <-
pmomLM(y=y,x=x,xadj=xadj,center=FALSE,scale=FALSE,niter=1000,
priorCoef=priorCoef,priorDelta=priorDelta,priorVar=priorVar)
mean(mom1$postOther)  #posterior mean for tau

#Probit model
n <- 500; rho <- .25; niter <- 1000
theta <- c(.4,.6,0); theta.adj <- 0
V <- diag(length(theta)); V[upper.tri(V)] <- V[lower.tri(V)] <- rho
x <- rmvnorm(n,rep(0,length(theta)),V); xadj <- matrix(1,nrow=nrow(x),ncol=1)
lpred <- as.vector(x \%*\% matrix(theta,ncol=1) + xadj \%*\% matrix(theta.adj,ncol=1))
p <- pnorm(lpred)
y <- runif(n)<p

mom2 <- pmomPM(y=y,x=x,xadj=xadj,niter=1000,priorCoef=priorCoef,
priorDelta=priorDelta,initSearch='greedy')
colMeans(mom2$postCoef1)
coef(glm(y ~ x + xadj -1, family=binomial(link='probit')))
}
\keyword{ models }
\keyword{ htest }

