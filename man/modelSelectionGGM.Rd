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
global_proposal= 'regression', prob_global=0.5,
tempering= 0.5, truncratio= 100,
save_proposal= FALSE, niter=10^3, burnin= round(niter/10), 
updates_per_iter= ceiling(sqrt(ncol(y))), updates_per_column= 10,
sampler='Gibbs', pbirth=0.75, pdeath=0.5*(1-pbirth),
bounds_LIT, Omegaini='glasso-ebic', verbose=TRUE)

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

  \item{global_proposal}{Either \code{'none'}, \code{'regression'} or 
  \code{'in-sample'}. 
  If \code{'none'}, serial MCMC is used as specified by \code{sampler}. 
  If \code{'regression'}, 
  MCMC uses Metropolis-Hastings where models are proposed for each column
  using regression posterior model probabilities for that column.
  If \code{'in-sample'}, each column's proposal use posterior model
  probabilities given a precision matrix estimate for the other columns}

  \item{prob_global}{Probability of proposing a sample from a
  global proposal, see details. This argument is ignored if
  \code{global_proposal == "none"}.
  }

  \item{tempering}{If \code{global_proposal != 'none'}, the posterior model
  probabilities of the proposal distribution are raised to the power
  indicated by \code{tempering} (set to 1 for no tempering)}

  \item{truncratio}{In the global proposal, any model's proposal probability 
  is >= prob(top model) / truncratio. This ensures bounded weight ratios
  in the MH step, to improve poor mixing when the current state has low
  proposal probability, often at the cost of decreasing the acceptance rate.
  If \code{truncratio <= 0}, no truncation is done}

  \item{save_proposal}{If \code{TRUE}, the global proposals are saved
  in \code{proposal} (a list with p entries, one per column) and the
  corresponding proposal densities in \code{proposaldensity}. Neither
  are typically needed, as they were already used to produce the
  posterior samples in \code{postSample}}
  
  \item{sampler}{Posterior sampler used when \code{global=="none"},
  and also to run the parallel proposals when \code{global!="none"}. 
  Options are \code{"Gibbs"} for Gibbs sampling,\code{"birthdeath"} for
  birth-death-swap, and \code{"LIT"} for the locally-informed truncated 
  algorithm of Zhou et al (2022)}
  
  \item{niter}{Number of posterior samples to be obtained. Each iteration
  consists of selecting a column of the precision matrix at random and 
  making \code{updates_per_column} updates to its entries}

  \item{burnin}{The first burnin samples will be discarded}

  \item{updates_per_iter}{An iteration consists of selecting
  \code{updates_per_iter} columns at random, and proposing
  \code{updates_per_column} edge updates within each column}

  \item{updates_per_column}{See \code{updates_per_iter}}

  \item{pbirth}{Probability of a birth move. The probability of a swap move 
  is \code{1-pbirth-pdeath}. Ignored unless \code{sampler=="birthdeath"}}

  \item{pdeath}{Probability of a death move. Ignored unless
    \code{sampler=="birthdeath"}}  

  \item{bounds_LIT}{log-proposal density bound for the locally-informed LIT 
  algorithm of Zhou et al. These bound the proposal density for a death move
  and for a birth move. By  default, \code{bounds_LIT} is
  \code{log(c("lbound_death"=1/p, "ubound_death"= 1, "lbound_birth"=1/p, "ubound_birth"=p))}
  }
 
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
  
  Independent spike-and-slab priors are set on the off-diagonal entries of 
  Omega, i.e. Omega[j,k]=0 with positive probability (spike) and otherwise
  arises from the distribution indicated in \code{priorCoef} (slab).

  Inference is based on MCMC posterior sampling. All sampling algorithms
  proceed by updating Omega[,k] given y and Omega[,-k]
  (of course, Omega[k,] is also set to Omega[,k]).
  Omega[,k] is updated by first updating the set of non-zero entries
  (i.e. edges in the graphical model) using either Gibbs sampling
  or a proposal distribution (see below), and then the non-zero
  entries of Omega[,k] are updated from their exact posterior given
  the current set of edges.

  An MCMC iteration consists of iterating over \code{updates_per_iter} columns 
  chosen at random and, for each column, do \code{updates_per_column} 
  proposals.

  If \code{global_proposal == "none"}, a local MCMC proposal is used to 
  update what entries in Omega[,k] are zero. 
  Local means that the proposed model is a neighbour of the current model,
  e.g. only one edge is added / killed.
  Available local kernels are Gibbs, birth-death-swap and LIT (Zhou et al 2022).

  If \code{global_proposal == "regression"}, a global proposal is used.
  A new model (set of non-zero entries in Omega[,k]) is proposed based 
  on the posterior distribution of a linear regression of y[,k] on the other
  covariates. 
  \code{prob_global} indicates the probability of using the global proposal,
  otherwise a local proposal is used.
  Proposal probabilities are tempered by raising them to the power
  \code{tempering}. Further, any model with probability 
  below prob(top model) / truncratio is assigned proposal probability
  prob(top model) / truncratio.
 
}
\value{
  Posterior inference on the inverse covariance of \code{y}.
  Object of class \code{msfit_ggm}, which extends a list with elements
  
\item{postSample}{Posterior samples for the upper-diagonal entries of
  the precision matrix. Stored as a sparse matrix, see package Matrix
  to utilities to work with such matrices}

\item{prop_accept}{If \code{almost_parallel} is \code{TRUE}, a vector
  with the proportion of accepted edge proposals. Note that Omega[,k]
  is always updated from its exact conditional posterior, regardless
  of the edge proposal being accepted or rejected.}

\item{proposal}{If \code{almost_parallel} and \code{save_proposal} are
  \code{TRUE}, this is a list with one entry per column of Omega,
  containing the proposed values of each column}

\item{proposaldensity}{log-proposal density for the samples in
  \code{proposal}. Entry (i,j) stores the log-proposal density for
  proposed sample i of column j}

\item{margpp}{Rao-Blackwellized estimates of posterior marginal
  inclusion probabilities. Only valid when using the Gibbs algorithm}

\item{priors}{List storing the priors specified when calling
  \code{modelSelectionGGM}}

\item{p}{Number of columns in \code{y}}

\item{indexes}{Indicates what row/column of Omega is stored in each
  column of \code{postSample}}

\item{samplerPars}{MCMC sampling parameters}

\item{almost_parallel}{Stores the input argument \code{almost_parallel}}

}
\author{ David Rossell }
\seealso{
  \code{\link{msfit_ggm-class}} for further details on the output.
  \code{icov} for the estimated precision (inverse covariance) matrix.
  \code{coef.msfit_ggm} for Bayesian model averaging estimates and
  intervals. 
}
\references{
  Zhou, Quan, et al. Dimension-free mixing for high-dimensional Bayesian 
  variable selection. JRSS-B 2022, 84, 1751-1784
}
\examples{

#Simulate data with p=3
Th= diag(3); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)

z= matrix(rnorm(1000*3), ncol=3)
y= z %*% chol(sigma)

#Obtain posterior samples
#y has few columns, initialize the MCMC at the sample precision matrix
#Else, leave the argument Omegaini in modelSelectionGGM empty
Omegaini= solve(cov(y))
fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini)

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
