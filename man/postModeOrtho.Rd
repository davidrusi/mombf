\name{postModeOrtho}
\alias{postModeOrtho}
\alias{postModeBlockDiag}
\title{Bayesian model selection and averaging under block-diagonal X'X for linear models.}
\description{
  postModeOrtho is for diagonal X'X,
  postModeBlockDiag for the more general block-diagonal X'X,
  where X is the matrix with predictors.

  Both functions return the model of highest posterior probability of any given
  size using an efficient search algorithm. This sequence of models includes
  the highest posterior probability model (HPM).
  Posterior model probabilities, marginal variable inclusion probabilities
  and Bayesian model averaging estimates are also provided.
  The unknown residual variance is integrated out using an exact deterministic
  algorithm of low computational cost (see details in reference).
}
\usage{
postModeOrtho(y, x, priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(1,1),
priorVar=igprior(0.01,0.01), bma=FALSE, includeModels, maxvars=100)

postModeBlockDiag(y, x, blocks, priorCoef=zellnerprior(tau=nrow(x)),
priorDelta=modelbinomprior(p=1/ncol(x)),priorVar=igprior(0.01,0.01), bma=FALSE,
maxvars=100, momcoef)

}
\arguments{
   \item{y}{Outcome}

   \item{x}{Matrix with predictors. If an intercept is desired x should include
   a column of 1's.}

   \item{blocks}{Factor or integer vector of length ncol(x) indicating the block
   that each column in x belongs to.}

   \item{priorCoef}{Prior distribution for the coefficients. Object created
    with \code{momprior}, \code{imomprior}, \code{emomprior} or \code{zellnerprior}.}

   \item{priorDelta}{Prior on model space. Use \code{modelbbprior()} for Beta-Binomial prior,
   \code{modelbinomprior(p)} for Binomial prior with prior inclusion
   probability \code{p},
   \code{modelcomplexprior} for Complexity prior,
   or \code{modelunifprior()} for Uniform prior}

   \item{priorVar}{Inverse gamma prior on residual variance, created with \code{igprior()}}

   \item{bma}{Set to TRUE to obtain marginal inclusion probabilities and
     Bayesian model averaging parameter estimates for each column of x.}

   \item{includeModels}{Models that should always be included when
   computing posterior model probabilities. It must be a list, each
   element in the list corresponds to a model and must be a logical or
   numeric vector indicating the variables in that model}

   \item{maxvars}{The search for the HPM is restricted to models with up to maxvars variables
   (note: posterior model probabilities and BMA are valid regardless of maxvars)}

   \item{momcoef}{optional argument containing pre-computed coefficients needed to obtain
   the marginal likelihood under the pMOM prior.
   A first call to postModeBlockDiag returns these coefficients,
   thus this argument is useful to speed up successive calls.}
}
\details{
  The first step is to list a sequence of models with 0,...,maxvars variables which,
  under fairly general conditions listed in Papaspiliopoulos & Rossell (2016),
  is guaranteed to include the HPM.
  Then posterior model probabilities are computed for all these models to determine
  the HPM, evaluate the marginal posterior of the residual variance on a grid,
  and subsequently compute the marginal density p(y) via adaptive quadrature.
  Finally this adaptive grid is used to compute marginal inclusion probabilities
  and Bayesian model averaging estimates.
  For more details see Papaspiliopoulos & Rossell (2016).
}
\value{
List with elemants
\item{models}{data.frame indicating the variables included in the sequence of models found
during the search of the HPM, and their posterior probabilities. The model with highest
posterior probability in this list is guaranteed to be the HPM.}

\item{phi}{data.frame containing an adaptive grid of phi (residual variance) values and
their marginal posterior density p(phi|y).}

\item{logpy}{log-marginal density p(y), i.e. normalization constant of p(phi|y).}

\item{bma}{Marginal posterior inclusion probabilities and Bayesian model averaging estimates
for each column in x.}

\item{postmean.model}{Coefficient estimates conditional on each of the models in \code{models}}

\item{momcoef}{If a MOM prior was specified in priorCoef, momcoef stores some coefficients needed
to compute its marginal likelihood}

}
\references{
  Papaspiliopoulos O., Rossell D. Scalable Bayesian variable selection
  and model averaging under block-orthogonal design. 2016
}
\author{ David Rossell }
\examples{
#Simulate data
set.seed(1)
p <- 400; n <- 410
x <- scale(matrix(rnorm(n*p),nrow=n,ncol=p),center=TRUE,scale=TRUE)
S <- cov(x)
e <- eigen(cov(x))
x <- t(t(x \%*\% e$vectors)/sqrt(e$values))
th <- c(rep(0,p-3),c(.5,.75,1)); phi <- 1
y <- x \%*\% matrix(th,ncol=1) + rnorm(n,sd=sqrt(phi))

#Fit
priorCoef=zellnerprior(tau=n); priorDelta=modelbinomprior(p=1/p); priorVar=igprior(0.01,0.01)
pm.zell <- postModeOrtho(y,x=x,priorCoef=priorCoef,priorDelta=priorDelta,priorVar=priorVar,
bma=TRUE)

#Best models
head(pm.zell$models)

#Posterior probabilities for sequence of models
nvars <- sapply(strsplit(as.character(pm.zell$models$modelid),split=','),length)
plot(nvars,pm.zell$models$pp,ylab='post prob',xlab='number of vars',ylim=0:1,xlim=c(0,50))

#Marginal posterior of phi
plot(pm.zell$phi,type='l',xlab='phi',ylab='p(phi|y)')

#Marginal inclusion prob & BMA estimates
plot(pm.zell$bma$margpp,ylab='Marginal inclusion prob')
plot(pm.zell$bma$coef,ylab='BMA estimate')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

