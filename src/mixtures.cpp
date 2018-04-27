#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "cstat.h"
//#include "do_mombf.h"
#include "mixtures.h"


SEXP normalmixGibbsCI(SEXP Sx, SEXP Sn, SEXP Sp, SEXP Sncomp, SEXP Sz, SEXP Smu0, SEXP Sg, SEXP Snu0, SEXP SS0, SEXP Sq, SEXP SB, SEXP Sburnin, SEXP Sverbose) {
  //Posterior sampling for Normal mixture models under a Normal-IWishart-Dir prior. Also posterior probability of one empty cluster, required for Bayes factor calculation
  //
  //Likelihood p(x[i,] | mu,Sigma,eta)= sum_j eta_j N(x[i,]; mu_j,Sigma_j)
  //Prior: p(mu_j, Sigma_j)= N(mu_j; mu0, g Sigma) IW(Sigma_j; nu0, S0) indep j=1,...,k
  //       p(eta)= Dir(eta; q)
  //Input
  // - x: n x p data matrix, individuals in rows and variables in columns
  // - ncomp: number of components
  // - z: initial cluster allocations (integer vector of length n, each z[i] must be in [1,ncomp])
  // - mu0, g, S0, nu0, q: prior parameters
  // - B: number of MCMC iterations
  // - burnin: number of burn-in iterations
  //
  //Output
  // - pponeempty: average posterior probability that one cluster is empty. Let n_j be the number of individuals in cluster j, then mean P(n_j=0 |y) across j=1,...,k
  // - eta: MCMC draws for mixing probabilities. A matrix with B-burnin rows and k columns.
  // - mu: MCMC draws for means, B-burnin rows and p*k columns. Each column is ordered mu_1, mu_2,...,mu_ncomp
  // - cholSigmainv: MCMC draws for Cholesky decomposition of inverse covariances. B-burnin rows, k*p*(p+1)/2 columns. Ordered Sigma_1,...,Sigma_ncomp
  int niter, nelemcov;
  double *pponeempty, *eta, *mu, *cholSigmainv;
  SEXP ans;

  niter= INTEGER(SB)[0] - INTEGER(Sburnin)[0];
  nelemcov= (INTEGER(Sp)[0])*(INTEGER(Sp)[0]+1)/2;

  PROTECT(ans= allocVector(VECSXP, 4));
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, 1));
  pponeempty= REAL(VECTOR_ELT(ans,0));

  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, INTEGER(Sncomp)[0] * niter)); //posterior draws for mixing probabilities
  eta= REAL(VECTOR_ELT(ans,1));

  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, INTEGER(Sncomp)[0] * INTEGER(Sp)[0] * niter)); //means
  mu= REAL(VECTOR_ELT(ans,2));

  SET_VECTOR_ELT(ans, 3, allocVector(REALSXP, INTEGER(Sncomp)[0] * nelemcov * niter)); //covariances
  cholSigmainv= REAL(VECTOR_ELT(ans,3));

  normalmixGibbsC(pponeempty, eta, mu, cholSigmainv, REAL(Sx), INTEGER(Sn), INTEGER(Sp), INTEGER(Sncomp), INTEGER(Sz), REAL(Smu0), REAL(Sg), INTEGER(Snu0), REAL(SS0), REAL(Sq), INTEGER(SB), INTEGER(Sburnin), INTEGER(Sverbose));

  UNPROTECT(1);
  return ans;
}


void normalmixGibbsC(double *pponeempty, double *eta, double *mu, double *cholSigmainv, double *x, int *n, int *p, int *ncomp, int *z, double *mu0, double *g, int *nu0, double *S0, double *q, int *B, int *burnin, int *verbose) {
  bool posdef;
  int b, i, j, l, idxeta=0, idxmu=0, idxSigma=0, idxctr, nelemmu, nelemSigma;
  double shrin, w, ginv, nplusginv, *dm, **xsum, *muz, *zcount, *qpost, **cholSigma, **Sinv, **cholSinv, **cholSigmainvcur, ***S;

  //Initialize
  (*pponeempty)= 0;
  nelemmu= (*ncomp)*(*p); nelemSigma= (*ncomp)*(*p)*(*p +1)/2;
  ginv= 1/(*g);

  dm= dvector(1,*p); zcount= dvector(1,*ncomp); qpost= dvector(1,*ncomp); muz= dvector(1,*p);
  xsum= dmatrix(1,*p,1,*ncomp); cholSigma= dmatrix(1,*p,1,*p); Sinv= dmatrix(1,*p,1,*p); cholSinv= dmatrix(1,*p,1,*p); cholSigmainvcur= dmatrix(1,*p,1,*p); S= darray3(1,*ncomp,1,*p,1,*p);

  for (l=1; l<=(*ncomp); l++) { for (j=1; j<=(*p); j++) { xsum[j][l]= 0; } }
  for (i=0; i<(*n); i++) {
    zcount[z[i]] +=1;
    for (j=1; j<=(*p); j++) xsum[j][z[i]] += x[i + (*n)*(j-1)];  //xsum[j,l] contains sum of variable j in component l
  }
  for (l=1; l<=(*ncomp); l++) { qpost[l]= zcount[l] + (*q); }

  //sum of squares and cross-products within cluster
  sumsqbyclus(x, *n, *p, z, *ncomp, false, S);
  for (i=1; i<=(*p); i++) { for (j=i; j<=(*p); j++) { for (l=1; l<=(*ncomp); l++) { S[l][i][j]+= S0[(i-1)*(*p)+j-1]; } } } //add prior scale matrix
  for (l=1; l<=(*ncomp); l++) { //add term (xbar - mu0) (xbar-mu0)' * n/(1+n*g)
    shrin= ((double) zcount[l]) / (1.0 + (*g) * ((double) zcount[l]));
    for (j=1; j<=(*p); j++) { dm[j]= (xsum[j][l]/((double) zcount[l]) - mu0[j]); }
    for (i=1; i<=(*p); i++) { for (j=i; j<=(*p); j++) { S[l][i][j] += shrin * dm[i] * dm[j]; } }
  }

  for (b=0; b<(*B); b++) {
    //Sample means (mu) and covariances (Sigma)
    for (l=1; l<=(*ncomp); l++) {
      //Sample Sigma
      //choldc_inv(S[l], *p, cholSinv, &posdef);  //cholSinv is the inverse of chol(S[l]). Then S[l]^{-1}= t(cholSinv) %*% cholSinv
      inv_posdef_upper(S[l], *p, Sinv, &posdef);
      choldc(Sinv, *p, cholSinv, &posdef);
      rwishartC(cholSigmainvcur, *nu0 + (int) zcount[l], cholSinv, *p, true);
      if (b>=(*burnin)) { //store sampled value
	idxctr= idxSigma + (l-1)*(*p)*(*p +1)/2;
	for (i=1; i<=(*p); i++) { for (j=i; j<=(*p); j++) { cholSigmainv[idxctr]= cholSigmainvcur[j][i]; idxctr++; } }
      }
      //Sample mu= rmvnorm(1, colMeans(x)*w + mu0*(1-w), Sigma/(n+1/g))
      w= ((double) zcount[l])/((double) zcount[l] + ginv);
      cholS_inv(cholSigmainvcur, *p, cholSigma);
      for (j=1; j<=(*p); j++) muz[j]= rnormC(0,1);
      Atx(cholSigma, muz, mu-1 + idxmu + (l-1)*(*p), 1, *p, 1, *p); //mu= t(cholSigma) %*% muz
      nplusginv= sqrt(((double) zcount[l]) + ginv);
      for (j=1; j<=(*p); j++) {
	mu[idxmu +j-1 + (l-1)*(*p)] /= nplusginv; //mu= (t(cholSigma) %*% muz) / sqrt(n+1/g)
	mu[idxmu +j-1 + (l-1)*(*p)] += (xsum[j][l] / ((double) zcount[l])) * w + mu0[j-1] * (1-w); //mu= [(t(cholSigma) %*% muz) / sqrt(n+1/g)] + [colMeans(x)*w + mu0*(1-w)]
      }
    }
    //Sample latent cluster indicators (z)
    //cluster prob: careful, requires new version of dmvnormC where the input is t(cholsinv) instead of cholsinv. Try to program overloaded function.
    //if cluster changes, update xsum
    //if cluster changes, update S
    //Sample mixing weights (eta)
    rdirichlet(eta+idxeta, zcount+1, ncomp);
    if (b>=(*burnin)) {
      idxeta+= (*ncomp);
      idxmu+= nelemmu;
      idxSigma+= nelemSigma;
    }

  }
  free_dvector(dm,1,*p); free_dvector(zcount,1,*ncomp); free_dvector(qpost,1,*ncomp); free_dvector(muz,1,*p);
  free_dmatrix(xsum,1,*p,1,*ncomp);
  free_dmatrix(cholSigma,1,*p,1,*p); free_dmatrix(Sinv,1,*p,1,*p); free_dmatrix(cholSinv,1,*p,1,*p); free_dmatrix(cholSigmainvcur,1,*p,1,*p); free_darray3(S,1,*ncomp,1,*p,1,*p);
}
