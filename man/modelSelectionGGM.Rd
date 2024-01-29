\name{modelSelectionGGM}
\alias{modelSelectionGGM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian variable selection for linear models via non-local priors. }
\description{
  Bayesian model selection for linear, asymmetric linear,
  median and quantile regression under
  non-local or Zellner priors. p>>n can be handled.

  modelSelection enumerates all models when feasible
  and uses a Gibbs scheme otherwise.
  See \code{coef} and \code{coefByModel} for estimates and posterior
  intervals of regression coefficients, and \code{rnlp} for posterior samples.

  modelsearchBlockDiag seeks the highest posterior
  probability model using an iterative block search.
}
\usage{


modelSelectionGGM(y, priorCoef=normalidprior(tau=1), 
priorModel=modelbinomprior(1/ncol(y)), 
priorDiag=exponentialprior(lambda=1), center=TRUE, scale=TRUE, 
almost_parallel= FALSE, sampler='Gibbs', niter=10^3, 
burnin= round(niter/10), pbirth=0.5, nbirth, 
Omegaini='glasso-ebic', verbose=TRUE)

}
\arguments{
  
  \item{y}{Data matrix}
  
  \item{priorCoef}{Prior on off-diagonal entries of the precision
    matrix, conditional on their not being zero (slab)}
  
  \item{priorModel}{Prior probabilities on having non-zero diagonal
    entries}
  
  \item{priorDiag}{Prior on diagonal entries of the precision matrix}
  
  \item{center}{If \code{TRUE}, the columns of \code{y} will be centered
    to zero mean}
  
  \item{scale}{If \code{TRUE}, the columns of \code{y} will be scaled to
    unit sample variance}

  \item{almost_parallel}{Use almost parallel algorithm sampling from each 
column independently and using an MH step}
  
  \item{sampler}{Posterior sampler. Options are "Gibbs", "birthdeath"
  and "zigzag"}
  
  \item{niter}{Number of posterior samples to be obtained}

  \item{pbirth}{Probability of a birth move. Ignored unless
    \code{sampler=="birthdeath"}}

  \item{nbirth}{Number of birth/death updates to perform for each row of
  the precision matrix. Defaults to \code{ncol(y)}}
  
  \item{burnin}{The first burnin samples will be discarded}
  
  \item{Omegaini}{Initial value of the precision matrix Omega. "null"
  sets all off-diagonal entries to 0. "glasso-bic" and "glasso-ebic" use
  GLASSO with regularization parameter set via BIC/EBIC,
  respectively. Alternatively, \code{Omegaini} can be a matrix}

  \item{verbose}{Set \code{verbose==TRUE} to print iteration progress}

  
}
\details{

  Let Omega be the inverse covariance matrix.
  A spike-and-slab prior is used. Specifically,
  independent priors are set on all Omega[j,k], and then a
  positive-definiteness truncation is added.
  
  The prior on diagonal entries Omega[j,j] is given by \code{priorDiag}.
  Off-diagonal Omega[j,k] are equal to zero with probability given by
  \code{modelPrior} and, when non-zero, they are
  
  Independent spike-and-slab priors are set on the off-diagonal entries of Omega,
  i.e. Omega[j,k]=0 with positive probability (spike) and otherwise
  arises from the distribution indicated in \code{priorCoef} (slab).
  
}
\value{
  Posterior inference on the inverse covariance of \code{y}.
  Object of class \code{msfit_ggm}, which extends a list with elements
  
\item{postSample}{Posterior samples for the upper-diagonal entries of
  the precision matrix. Stored as a sparse matrix, see package Matrix
  to utilities to work with such matrices}

\item{p}{Number of columns in \code{y}}

\item{priors}{List storing the priors specified when calling
  \code{modelSelectionGGM}}

}
\author{ David Rossell }
\seealso{
  \code{\link{msfit_ggm-class}} for further details on the output.
  \code{icov} for the estimated precision (inverse covariance) matrix.
  \code{coef.msfit_ggm} for Bayesian model averaging estimates and
  intervals. 
}
\examples{

#Simulate data with p=3
Th= diag(3); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)

z= matrix(rnorm(1000*3), ncol=3)
y= z %*% chol(sigma)

#Obtain posterior samples
fit= modelSelectionGGM(y, scale=FALSE)

#Parameter estimates, intervals, prob of non-zero
coef(fit)

#Estimated inverse covariance
icov(fit)

#Estimated inverse covariance, entries set to 0
icov(fit, threshold=0.95)

#Shows first posterior samples
head(fit$postSample)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
