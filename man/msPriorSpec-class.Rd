\name{msPriorSpec-class}
\Rdversion{1.1}
\docType{class}
\alias{msPriorSpec}
\alias{msPriorSpec-class}
\alias{aic}
\alias{bic}
\alias{bicprior}
\alias{ic}
\alias{emomprior}
\alias{exponentialprior}
\alias{icarplusprior}
\alias{imomprior}
\alias{momprior}
\alias{zellnerprior}
\alias{normalidprior}
\alias{groupmomprior}
\alias{groupemomprior}
\alias{groupimomprior}
\alias{groupzellnerprior}
\alias{modelunifprior}
\alias{modelbinomprior}
\alias{modelbbprior}
\alias{modelcomplexprior}
\alias{igprior}

\title{Class "msPriorSpec"}
\description{
Stores the prior distributions to be used for Bayesian variable selection
in normal regression models.
This class can be used to specify the prior on non-zero regression coefficients,
the model indicator or the nuisance parameters.
}
\usage{
aic()
bic()
bicprior()
ic(penalty)

momprior(taustd=1, tau, tau.adj=10^6, r=1)
imomprior(tau, tau.adj=10^6)
emomprior(tau, tau.adj=10^6)
zellnerprior(taustd=1, tau, tau.adj=10^6)
normalidprior(taustd=1, tau, tau.adj=10^6)
icarplusprior(a=0.5, taustd=1, tau.adj=10^6)

exponentialprior(lambda = 1)

groupmomprior(taustd=1, tau, tau.adj=10^6)
groupimomprior(tau, tau.adj=10^6)
groupemomprior(tau, tau.adj=10^6)
groupzellnerprior(taustd=1, tau, tau.adj=10^6)

modelunifprior()
modelbinomprior(p=0.5)
modelbbprior(alpha.p=1, beta.p=1)
modelcomplexprior(c=1)

igprior(alpha=.01, lambda=.01)
}
\arguments{
  \item{penalty}{Penalty on model dimension, i.e. for the AIC penalty=2}
  \item{tau}{Prior dispersion parameter for covariates undergoing selection}
  \item{taustd}{Prior dispersion parameter for covariates undergoing selection.
    It is calibrated so that `taustd=1` equals the unit information prior.}
  \item{tau.adj}{Prior variance in Normal prior for covariates not undergoing selection}
  \item{r}{MOM prior parameter is \code{2*r}}
  \item{a}{The icarplus prior precision matrix is a P + (1-a) tau I, where P is
  an ICAR precision matrix and tau I a normalidprior precision matrix}
  \item{p}{Prior inclusion probability for binomial prior on model space}
  \item{alpha.p}{Beta-binomial prior on model space has parameters alpha.p, beta.p}
  \item{beta.p}{Beta-binomial prior on model space has parameters
    alpha.p, beta.p}
  \item{c}{Under the Complexity prior the prior
    probability of having k variables in the model is proportional to \code{1/p^(ck)}}
  \item{alpha}{Inverse gamma prior has parameters alpha/2, lambda/2}
  \item{lambda}{\code{igprior} defines an inverse gamma prior with
    parameters alpha/2, lambda/2. \code{exponentialprior} defines an
    exponential prior with rate parameter \code{lambda}}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("msPriorSpec",
  ...)}, but it is easier to use creator functions.

For priors on regression coefficients use \code{momprior},
\code{imomprior} or \code{emomprior}.
For prior on model space \code{modelunifprior}, \code{modelbinomprior}
\code{modelbbprior}, or \code{modelcomplexprior}.
For prior on residual variance use \code{igprior}.
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
If \code{priorType=="coefficients"}, \code{priorDistr} can be equal to
"pMOM", "piMOM", "peMOM", "zellner", "normalid", "groupMOM" or "groupzellner"
(product moment, product inverse moment, product exponential moment, Zellner prior, normal prior with \eqn{\Sigma=\mathbf{I}}{V=I}, respectively).
If \code{priorType=="modelIndicator"}, \code{priorDistr} can be equal to "uniform" or "binomial"
to specify a uniform prior (all models equaly likely a priori) or a
binomial prior, or to "complexity" for the Complexity prior of Castillo
et al 2015. For a binomial prior,
the prior inclusion probability for any single variable must be
specified in slot \code{priorPars['p']}. For a beta-binomial prior, the
Beta hyper-prior parameters must be in \code{priorPars['alpha.p']} and
\code{priorPars['beta.p']}.
For the Complexity prior, the prior parameter must be in the slot
\code{priorPars['c']}.
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
\details{
  DISCUSSION OF PRIOR ON PARAMETERS

  Let beta=(beta_1,...,beta_p) be the regression coefficients for
  individual variables and delta=(delta_1,...,delta_q) those for grouped
  variables (e.g. factors or smooth terms in modelSelection).

  momprior, emomprior, imomprior, zellnerprior, normalid and icarplus
  can be priors on both beta or delta. For further information see the vignette.

  groupzellnerprior is the prior density on delta

  \deqn{p_z(\delta; \tau)= \prod_j N(\delta_j; 0, (\tau/p_j)) (X_j'X_j)^{-1}}

  where \eqn{X_j} are the design matrix columns associated to \eqn{delta_j} and p_j=ncol(X_j)
  is the number of covariates in the group (for groupmomprior, the term in the
  denominator is (p_j +2) instead of p_j). A default
  tau=n=nrow(X_j) mimics the unit information prior and implies that the
  ratio of variance explained by X_j / residual variance is expected to be
  1 a priori. To set the dispersion in terms of unit information prior, taustd
  is also available.

  groupmomprior adds a quadratic MOM penalty

  p_m(delta; tau)= p_z(delta; tau * n) prod_j delta_j'X_j'X_jdelta_j ncol(X_j)/(tau * n * p_j / (p_j + 2))

  and analogously for eMOM and iMOM. Note that unlike groupzellnerprior, the nrow(X_j)
  factor is already included in the code. This is done to give user introduced tau values
  a roughly similar meaning between momprior and groupmomprior.


  DISCUSSION OF PRIOR ON MODELS

Under the uniform prior, the prior probability of any model is
1 / number of models.

Under the Binomial, Beta-Binomial and Complexity priors a model with k
out of K active variables has prior probability
P(Z=k) / (K choose k), where
where Z ~ Binomial(K,p),
Z ~ BetaBinomial(K,alpha.p,beta.p)
or for the Complexity prior P(Z=k) proportional to \code{1/K^(c*k)}.
}
\references{
Johnson VE, Rossell D. Non-Local Prior Densities for Default Bayesian Hypothesis Tests. Journal of the Royal Statistical Society B, 2010, 72, 143-170

Johnson VE, Rossell D. Bayesian model selection in high-dimensional
settings. Journal of the American Statistical Association, 107,
498:649-660.

Rossell D, Abril O, Bhattacharya A. Approximate Laplace approximations for scalable model selection (2021). Journal of the Royal Statistical Society B, 83, 4, 853-879
}
\author{
David Rossell
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
