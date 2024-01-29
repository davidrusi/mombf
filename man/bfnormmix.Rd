\name{bfnormmix}
\alias{bfnormmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of Normal mixture components under Normal-IW and
  Non-local priors}
\description{
  Posterior sampling and Bayesian model selection to choose the number
  of components k in multivariate Normal mixtures.

  \code{bfnormmix} computes posterior probabilities under non-local
  MOM-IW-Dir(q) priors, and also for local Normal-IW-Dir(q.niw) priors.
  It also computes posterior probabilities on cluster occupancy
  and posterior samples on the model parameters for several k.
}
\usage{
bfnormmix(x, k=1:2, mu0=rep(0,ncol(x)), g, nu0, S0, q=3, q.niw=1,
B=10^4, burnin= round(B/10), logscale=TRUE, returndraws=TRUE, verbose=TRUE)
}
\arguments{
\item{x}{n x p input data matrix}
\item{k}{Number of components}
\item{mu0}{Prior on mu[j] is N(mu0,g Sigma[j])}
\item{g}{Prior on mu[j] is N(mu0,g Sigma[j]). This is a critical MOM-IW prior
  parameter that specifies the separation between
  components deemed practically relevant. It defaults to assigning 0.95
  prior probability to any pair of mu's giving a bimodal mixture, see details}
\item{S0}{Prior on Sigma[j] is IW(Sigma_j; nu0, S0)}
\item{nu0}{Prior on Sigma[j] is IW(Sigma_j; nu0, S0)}
\item{q}{Prior parameter in MOM-IW-Dir(q) prior}
\item{q.niw}{Prior parameter in Normal-IW-Dir(q.niw) prior}
\item{B}{Number of MCMC iterations}
\item{burnin}{Number of burn-in iterations}
\item{logscale}{If set to TRUE then log-Bayes factors are returned}
\item{returndraws}{If set to \code{TRUE} the MCMC posterior draws under
  the Normal-IW-Dir prior are returned for all \code{k}}
\item{verbose}{Set to \code{TRUE} to print iteration progress}
}
\details{
  The likelihood is

  p(x[i,] | mu,Sigma,eta)= sum_j eta_j N(x[i,]; mu_j,Sigma_j)

  The Normal-IW-Dir prior is

  Dir(eta; q.niw) prod_j N(mu_j; mu0, g Sigma) IW(Sigma_j; nu0, S0)

  The MOM-IW-Dir prior is

  \deqn{d(\mu,A) Dir(\eta; q) \prod_j N(\mu_j; \mu0, g \Sigma_j) IW(\Sigma_j; \nu_0, S0)}

  where

  \deqn{d(\mu,A)= [\prod_{j<l} (\mu_j-\mu_l)' A (\mu_j-\mu_l)]}

  and A is the average of \eqn{\Sigma_1^{-1},...,\Sigma_k^{-1}}. Note that
  one must have q>1 for the MOM-IW-Dir to define a non-local prior.

  By default the prior parameter g is set such that

  P( (mu[j]-mu[l])' A (mu[j]-mu[l]) < 4)= 0.05.

  The reasonale when Sigma[j]=Sigma[l] and eta[j]=eta[l]
  then (mu[j]-mu[l])' A (mu[j]-mu[l])>4 corresponds to a bimodal
  density. That is, the default g focuses 0.95 prior prob on a degree of
  separation between components giving rise to a bimodal mixture density.

  \code{bfnormmix} computes posterior model probabilities under the
  MOM-IW-Dir and Normal-IW-Dir priors using MCMC output. As described in
  Fuquene, Steel and Rossell (2018) the estimate is based on the
  posterior probability that one cluster is empty under each possible k.
}
\value{
  A list with elements
\item{k}{Number of components}
\item{pp.momiw}{Posterior probability of k components under a
  MOM-IW-Dir(q) prior}
\item{pp.niw}{Posterior probability of k components under a
  Normal-IW-Dir(q.niw) prior}
\item{probempty}{Posterior probability that any one cluster is empty under a
  MOM-IW-Dir(q.niw) prior}
\item{bf.momiw}{Bayes factor comparing 1 vs k components under a
  MOM-IW-Dir(q) prior}
\item{logpen}{log of the posterior mean of the MOM-IW-Dir(q) penalty term}
\item{logbf.niw}{Bayes factor comparing 1 vs k components under a
  Normal-IW-Dir(q.niw) prior}
}
\references{
Fuquene J., Steel M.F.J., Rossell D. On choosing mixture components via
  non-local priors. 2018. arXiv
}
\author{ David Rossell }
\examples{
x <- matrix(rnorm(100*2),ncol=2)

bfnormmix(x=x,k=1:3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

