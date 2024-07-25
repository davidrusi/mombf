\name{localnulltest}
\alias{localnulltest}
\alias{localnulltest_fda}
\alias{localnulltest_givenknots}
\alias{localnulltest_fda_givenknots}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Local variable selection }
\description{
  Learn whether covariate effects are zero at given coordinates using
  Bayesian model selection or information criteria.

  Use \code{coef} to extract estimates and posterior
  probabilities for local effects.
}
\usage{

localnulltest(y, x, z, x.adjust, localgridsize, localgrid,
nbaseknots=20, nlocalknots=c(5,10,15), localknots,
basedegree=3, cutdegree=0,
usecutbasis=TRUE, priorCoef=normalidprior(),
priorGroup=priorCoef, priorDelta=modelbbprior(),
mc.cores=min(4,length(nlocalknots)), return.mcmc=FALSE, verbose=FALSE,
...)

localnulltest_fda(y, x, z, x.adjust, function_id,
Sigma='AR/MA', localgridsize, localgrid, nbaseknots=20,
nlocalknots=c(5,10,15), localknots,
basedegree=3, cutdegree=0, usecutbasis=TRUE,
priorCoef=momprior(), priorGroup=groupmomprior(),
priorDelta=modelbbprior(), mc.cores=min(4,length(nlocalknots)),
return.mcmc=FALSE, verbose=FALSE, ...)

localnulltest_givenknots(y, x, z, x.adjust, localgridsize,
localgrid, nbaseknots=20, nlocalknots=10, localknots,
basedegree=3, cutdegree=0,
usecutbasis=TRUE, priorCoef=normalidprior(),
priorGroup=priorCoef, priorDelta=modelbbprior(),
verbose=FALSE, ...)

localnulltest_fda_givenknots(y, x, z, x.adjust, function_id,
Sigma='AR/MA', localgridsize, localgrid, nbaseknots=20,
nlocalknots=10, localknots,
basedegree=3, cutdegree=0, usecutbasis=TRUE,
priorCoef=momprior(), priorGroup=groupmomprior(),
priorDelta=modelbbprior(), verbose=FALSE, ...)

}
\arguments{
\item{y}{Vector with the outcome variable}
\item{x}{Numerical matrix with covariate values}
\item{z}{Matrix with d-dimensional coordinates (d>=1$ for each entry in \code{y}, and d columns)}
\item{x.adjust}{Optionally, further adjustment covariates to be included
  in the model with no testing being performed}
\item{function_id}{Function identifier. It is assumed that one observes
  multiple functions over z, this is the identifier of each individual
  function}
\item{Sigma}{Error covariance. By default 'identity', other options are
  'MA', 'AR' or 'AR/MA' (meaning that BIC is used to choose between MA
  and AR). Alternatively the user can supply a function such that
  \code{Sigma(z[i,],z[j,])} returns the within-function \code{cov(y[i,], y[j,])}}
\item{localgridsize}{Local test probabilities will be returned for a
  grid of \code{z} values of size \code{localgridsize} for each
  dimension. If \code{ncol(z)==1} then \code{localgridsize} defaults to
  100, else to 10}
\item{localgrid}{Regions at which tests will be performed. Defaults to
  dividing each \code{[min(z[,i]),  max(z[,i])]} into 10 equal
  intervals. If provided, \code{localgrid} must be a list with one entry
  for each \code{z[,i]}, containing a vector with the desired grid for that \code{z[,i]}}
\item{nbaseknots}{Number of knots for the spline approximation to the
  baseline effect of \code{x} on \code{y}}
\item{nlocalknots}{Number of knots for the basis capturing the local effects. 
Ignored if \code{localknots} is specified}
\item{localknots}{Knots to be used for the local tests. The same knots are used for all 
columns in \code{z}, so make sure that all columns span the same range of values.
For a multi-resolution analysis, \code{localknots} should be a list where each entry
specifyes the knots for one resolution}
\item{basedegree}{Degree of the spline approximation to the baseline}
\item{cutdegree}{Degree of the cut spline basis used for testing}
\item{usecutbasis}{If \code{FALSE}, then the basis is not cut and a
  standard spline basis is returned (not recommended unless you know
  what you're doing)}
\item{priorCoef}{Prior on the coefficients, passed on to
  \code{modelSelection}}
\item{priorGroup}{Prior on grouped coefficients, passed on to
  \code{modelSelection}}
\item{priorDelta}{Prior on the models, passed on to
  \code{modelSelection}}
\item{mc.cores}{If package parallel is available on your system and
  \code{nlocalknots} has several entries defining several resolution
  levels, they will be run in parallel on \code{mc.cores}}
\item{return.mcmc}{Set to \code{TRUE} to return the MCMC output from \code{modelSelection}}
\item{verbose}{If \code{TRUE} some progress information is printed}
\item{...}{Other arguments to be passed on to \code{modelSelection},
  e.g. \code{family='binomial'} for logistic regression}
}
\details{
  Local variable selection considers the model

  \deqn{y_i= \beta_0(z_i) + sum_{j=1}^p \beta_j(z_i, x_i) + e_i}

  \eqn{\beta_0(z_i)} is the baseline mean
  
  \eqn{\beta_j(z_i,x_i)} is local effect of covariate j at coordinate z_i
  
  \eqn{e_i} a Gaussian error term assumed either independent or with a
  covariance structure given by Sigma. If assuming independence it is
  possible to consider alternatives to Gaussianity,
  e.g. set \code{family='binomial'} for logistic regression
  or \code{family='poisson'} for Poisson regression

  Note: a sum-to-zero type constraint is set on \eqn{\beta_1(z_i,x_i)} so
  that it defines a deviation from the baseline mean \eqn{\beta_0(z_i)}

  We model \eqn{\beta_0} using B-splines of degree \code{basedegree} with
  \code{nbaseknots} knots.
  We model \eqn{\beta_j} using B-splines of degree \code{cutdegree} with
  \code{nlocalknots}. Using \code{cutdegree=0} runs fastest is usually
  gives similar inference than higher degrees, and is hence recommended
  by default.
  
}
\value{
  Object of class \code{localtest}, which extends a list with elements

\item{covareffects}{Estimated local covariate effects at different
  \code{z} values, 0.95 posterior intervals and posterior probability
  for the existence of an effect}
\item{pplocalgrid}{Posterior probabilities for the existence of an
  effect for regions of \code{z} values. Do not use these unless you
  know what you're doing}
\item{covareffects.mcmc}{MCMC output used to build covareffects. Only
  returned if \code{return.mcmc=TRUE}}
\item{ms}{Objects of class \code{msfit} returned by \code{modelSelection}}
\item{pp_localknots}{Posterior probability for each resolution level
  (value of \code{nlocalknots})}
\item{Sigma}{Input parameter}
\item{nlocalknots}{Input parameter}
\item{basedegree}{Input parameter}
\item{cutdegree}{Input parameter}
\item{knots}{Input parameters}
\item{regionbounds}{List with region bounds defined by the local testing
  knots at each resolution level}

}
\author{ David Rossell }
\examples{

#Simulate outcome and 2 covariates
#Covariate 1 has local effect for z>0
#Covariate 2 has no effect for any z

truemean= function(x,z) {
    ans= double(nrow(x))
    group1= (x[,1]==1)
    ans[group1]= ifelse(z[group1] <=0, cos(z[group1]), 1)
    ans[!group1]= ifelse(z[!group1]<=0, cos(z[!group1]), 1/(z[!group1]+1)^2)
    return(ans)
}

n= 1000
x1= rep(0:1,c(n/2,n/2))
x2= x1 + rnorm(n)
x= cbind(x1,x2)
z= runif(n,-3,3)
m= truemean(x,z)
y= truemean(x,z) + rnorm(n, 0, .5)

#Run localnulltest with 10 knots
fit0= localnulltest(y, x=x, z=z, nlocalknots=10, niter=1000)

#Estimated covariate effects and posterior probabilities
b= coef(fit0)
b

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

