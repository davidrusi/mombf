#include "cstat.h"
#include "mixtures.h"


// [[Rcpp::export]]
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
  // - logpen: log pMOM penalty term at each MCMC iteration, a vector with B-burnin elements. This is the sum of log-Mahalanobis distances

  // sum_{j<l} log [ (mu_j - mu_l)' A (mu_j - mu_l) ] + sum_l dmvnorm(mu_l; mu0, A) - dmvnorm(mu_l; mu0; g Sigma_l^{-1})

  // where A= sum_l Sigma^{-1}/(g * ncomp).

  // - eta: MCMC draws for mixing probabilities. A matrix with B-burnin rows and k columns.
  // - mu: MCMC draws for means, B-burnin rows and p*k columns. Each column is ordered mu_1, mu_2,...,mu_ncomp
  // - cholSigmainv: MCMC draws for Cholesky decomposition of inverse covariances. B-burnin rows, k*p*(p+1)/2 columns. Ordered Sigma_1,...,Sigma_ncomp
  int niter, nelemcov;
  double *logpen, *pponeempty, *eta, *mu, *cholSigmainv;
  SEXP ans;

  niter= INTEGER(SB)[0] - INTEGER(Sburnin)[0];
  nelemcov= (INTEGER(Sp)[0])*(INTEGER(Sp)[0]+1)/2;

  PROTECT(ans= Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(ans, 0, Rf_allocVector(REALSXP, 1));
  pponeempty= REAL(VECTOR_ELT(ans,0));

  SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, niter));
  logpen= REAL(VECTOR_ELT(ans,1));

  SET_VECTOR_ELT(ans, 2, Rf_allocVector(REALSXP, INTEGER(Sncomp)[0] * niter)); //posterior draws for mixing probabilities
  eta= REAL(VECTOR_ELT(ans,2));

  SET_VECTOR_ELT(ans, 3, Rf_allocVector(REALSXP, INTEGER(Sncomp)[0] * INTEGER(Sp)[0] * niter)); //means
  mu= REAL(VECTOR_ELT(ans,3));

  SET_VECTOR_ELT(ans, 4, Rf_allocVector(REALSXP, INTEGER(Sncomp)[0] * nelemcov * niter)); //covariances
  cholSigmainv= REAL(VECTOR_ELT(ans,4));

  normalmixGibbsC(pponeempty, logpen, eta, mu, cholSigmainv, REAL(Sx), INTEGER(Sn), INTEGER(Sp), INTEGER(Sncomp), INTEGER(Sz), REAL(Smu0), REAL(Sg), INTEGER(Snu0), REAL(SS0), REAL(Sq), INTEGER(SB), INTEGER(Sburnin), INTEGER(Sverbose));

  UNPROTECT(1);
  return ans;
}


void normalmixGibbsC(double *pponeempty, double *logpen, double *eta, double *mu, double *cholSigmainv, double *x, int *n, int *p, int *ncomp, int *z, double *mu0, double *g, int *nu0, double *S0, double *q, int *B, int *burnin, int *verbose) {
  bool posdef;
  int b, bidx, i, j, jj, l, idxj, idxeta=0, idxmu=0, idxSigma=0, idxctr, nelemmu, nelemSigma, znew, *zcount;
  double det, detA, detprior, ginv, gsqrt, logprior, maxi, maxpponeempty, nplusginv, shrin, sumi, tmp, w;
  double *pponeemptybyclus, *maxlogpempty, *dd, *dm, **A, **logpempty, **xbar, **xsum, *muz, *zprob, *qpost;
  double **cholA, **cholSigma, **Sinv, **cholSinv, **cholSigmainvcur, **cholSigmainvprior, **logpclus, ***S, ***crossprodx;

  //Allocate memory
  zcount= ivector(1,*ncomp);
  dd= dvector(1,(*ncomp)*(*ncomp -1)/2);
  pponeemptybyclus= dvector(1,*ncomp); maxlogpempty= dvector(1,*ncomp); dm= dvector(1,*p); zprob= dvector(1,*ncomp); qpost= dvector(1,*ncomp); muz= dvector(1,*p);
  A= dmatrix(1,*p,1,*p); logpempty= dmatrix(1,*B - *burnin,1,*ncomp); xbar= dmatrix(1,*ncomp,1,*p); xsum= dmatrix(1,*ncomp,1,*p); logpclus= dmatrix(1,*ncomp,1,*n);
  cholA= dmatrix(1,*p,1,*p); cholSigma= dmatrix(1,*p,1,*p); Sinv= dmatrix(1,*p,1,*p); cholSinv= dmatrix(1,*p,1,*p); cholSigmainvcur= dmatrix(1,*p,1,*p); cholSigmainvprior= dmatrix(1,*p,1,*p);
  S= darray3(1,*ncomp,1,*p,1,*p); crossprodx= darray3(1,*ncomp,1,*p,1,*p);

  //Initialize
  (*pponeempty)= 0;
  nelemmu= (*ncomp)*(*p); nelemSigma= (*ncomp)*(*p)*(*p +1)/2;
  ginv= 1/(*g); gsqrt= sqrt(*g);
  for (l=1; l<=(*ncomp); l++) { maxlogpempty[l]= R_NegInf; }

  //sum of squares and cross-products within cluster
  for (l=1; l<=(*ncomp); l++) {
    for (j=1; j<=(*p); j++) {
      xsum[l][j]= crossprodx[l][j][j]= 0;
      for (jj=j+1; jj<=(*p); jj++) { crossprodx[l][j][jj]= 0; }
    }
  }
  for (i=0; i<(*n); i++) {
    (zcount[z[i]])++;
    for (j=1; j<=(*p); j++) {
      idxj= i + (*n)*(j-1);
      xsum[z[i]][j] += x[idxj];
      crossprodx[z[i]][j][j] += pow(x[idxj],2.0);
      for (jj=j+1; jj<=(*p); jj++) { crossprodx[z[i]][j][jj] += x[idxj] * x[i + (*n)*(jj-1)]; }
    }
  }
  crossprod2sumsq_byclus(crossprodx,xsum,zcount,ncomp,p,S,xbar); //convert column sums and cross-prod into means and sum of squares

  //Add prior scale matrix S0 and term (xbar - mu0) (xbar-mu0)' * n/(1+n*g) to S
  for (l=1; l<=(*ncomp); l++) {
    shrin= ((double) zcount[l]) / (1.0 + (*g) * ((double) zcount[l]));
    for (j=1; j<=(*p); j++) { dm[j]= xbar[l][j] - mu0[j-1]; }
    for (i=1; i<=(*p); i++) { for (j=i; j<=(*p); j++) { S[l][i][j] += shrin * dm[i] * dm[j] + S0[(i-1)*(*p)+j-1]; } }
  }

  for (b=0; b<(*B); b++) {
    //Sample mixing weights (eta)
    for (l=1; l<=(*ncomp); l++) { qpost[l]= ((double) zcount[l]) + (*q); }
    rdirichlet(eta+idxeta, qpost+1, ncomp);

    //For each cluster, sample means (mu) and covariances (Sigma)
    logprior= 0;
    for (j=1; j<=(*p); j++) { for (jj=j; jj<=(*p); jj++) { A[j][jj]= 0; } }  //initialize A= sum_l Sigma_l^{-1}
    for (l=1; l<=(*ncomp); l++) {
      //Sample Sigma
      inv_posdef_upper(S[l], *p, Sinv, &posdef);
      choldc(Sinv, *p, cholSinv, &posdef);
      rwishartC(cholSigmainvcur, *nu0 + zcount[l], cholSinv, *p, true);
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
	mu[idxmu +j-1 + (l-1)*(*p)] += xbar[l][j] * w + mu0[j-1] * (1-w); //mu= [(t(cholSigma) %*% muz) / sqrt(n+1/g)] + [colMeans(x)*w + mu0*(1-w)]
      }
      //Compute log p(y | cluster=l)
      det= choldc_det(cholSigmainvcur,*p);
      dmvnormmatC(logpclus[l]+1, x, *n, *p, mu-1 + idxmu + (l-1)*(*p), cholSigmainvcur, det, false, 1); //Recall: Sigma^{-1}= cholSigmainv %*% t(cholSigmainv)
      //Add Sigma^{-1}/(g*ncomp) to A
      if (b>=(*burnin)) {
	addcholStcholS(cholSigmainvcur, *p, (*g) * ((double) (*ncomp)), A);
	detprior= det / pow(*g, *p);
	for (j=1; j<=(*p); j++) { for (jj=1; jj<=j; jj++) { cholSigmainvprior[j][jj]= cholSigmainvcur[j][jj] / gsqrt; } }
        logprior -= dmvnormC(mu-1 + idxmu + (l-1)*(*p), *p, mu0-1, cholSigmainvprior, detprior, false, true, false); //log dmvnorm(mu_l, mu0, g*Sigma_l)
      }
    } //end for each cluster

    //pMOM penalty term
    bidx= b - *burnin;
    if (b>=(*burnin)) {
      //Sum of log-Mahalanobis distances between mu's
      choldc(A, *p, cholA, &posdef);
      detA= choldc_det(cholA, *p);
      mahaldist(mu + idxmu, *ncomp, *p, cholA, true, dd); //dd= squared Mahalanobis distances (mu_j - mu_l)' A (mu_j - mu_l)

      logpen[bidx]= 0;
      for (j=1; j<= (*ncomp)*(*ncomp -1)/2; j++) logpen[bidx] += log(dd[j]);
      //Contribution from the prior (sum_l log dmvnorm(mu_l, mu0, A^{-1}) - log dmvnorm(mu_l, mu0, g*Sigma_l)
      for (l=1; l<=(*ncomp); l++) {
	logprior += dmvnormC(mu-1 + idxmu + (l-1)*(*p), *p, mu0-1, cholA, detA, false, true, false); //log dmvnorm(mu_l, mu0, A^{-1})
      }
      logpen[bidx] += logprior;
    }

    //Sample latent clusters (z)
    if (b>=(*burnin)) { for (l=1; l<=(*ncomp); l++) { logpempty[bidx+1][l]= 0; } }
    for (i=1; i<= *n; i++) {
      logpclus[1][i]+= log(eta[idxeta]);
      maxi= logpclus[1][i];
      for (l=2; l<= *ncomp; l++) {
	logpclus[l][i]+= log(eta[idxeta+l-1]);
	if (logpclus[l][i]> maxi) { maxi= logpclus[l][i]; }
      }
      sumi= 0;

      for (l=1; l<= *ncomp; l++) { sumi += exp(logpclus[l][i] - maxi); }
      for (l=1; l<= *ncomp; l++) {
	zprob[l]= exp(logpclus[l][i] - maxi)/sumi;
	if (b>=(*burnin)) { logpempty[bidx +1][l] += max_xy(log(1-zprob[l]), -230.25850929940458); }  //prevent values <log(1e-100)= -230.258
      }
      znew= rdisc(zprob+1, *ncomp) + 1;
      if (znew != z[i-1]) {
        //Update xsum and cross-products in S
        for (j=1; j<=(*p); j++) {
	  idxj= i-1 + (*n)*(j-1);
	  xsum[z[i-1]][j] -= x[idxj]; xsum[znew][j] += x[idxj];
	  tmp= pow(x[idxj],2); crossprodx[z[i-1]][j][j] -= tmp; crossprodx[znew][j][j] += tmp;
	  for (jj=j+1; jj<=(*p); jj++) {
	    tmp= x[idxj] * x[i-1 + (*n)*(jj-1)];
 	    crossprodx[z[i-1]][j][jj] -= tmp;
 	    crossprodx[znew][j][jj] += tmp;
	  }
	}
	//Update cluster counts
	zcount[znew]++;
	zcount[z[i-1]]--;
	z[i-1]= znew;
      }
    } //end sample latent clusters

    crossprod2sumsq_byclus(crossprodx,xsum,zcount,ncomp,p,S,xbar); //convert column sums and cross-prod into means and sum of squares
    //Add S0 + (xbar - mu0) (xbar-mu0)' * n/(1+n*g) to S
    for (l=1; l<=(*ncomp); l++) {
      shrin= ((double) zcount[l]) / (1.0 + (*g) * ((double) zcount[l]));
      for (j=1; j<=(*p); j++) { dm[j]= xbar[l][j] - mu0[j-1]; }
      for (j=1; j<=(*p); j++) { for (jj=j; jj<=(*p); jj++) { S[l][j][jj] += shrin * dm[j] * dm[jj] + S0[(j-1)*(*p)+jj-1]; } }
    }

    if (b>=(*burnin)) {
      for (l=1; l<=(*ncomp); l++) { maxlogpempty[l]= max_xy(maxlogpempty[l], logpempty[bidx +1][l]); }
      idxeta+= (*ncomp);
      idxmu+= nelemmu;
      idxSigma+= nelemSigma;
    }

  } //free for b

  //Post-process probability of one empty cluster
  maxpponeempty= R_NegInf;
  for (l=1; l<=(*ncomp); l++) { pponeemptybyclus[l]= 0; }
  for (l=1; l<=(*ncomp); l++) {
    for (b=1; b<=(*B - *burnin); b++) { pponeemptybyclus[l] += exp(logpempty[b][l]-maxlogpempty[l]); }
    pponeemptybyclus[l]= log(pponeemptybyclus[l]) + maxlogpempty[l] - log((double) (*B - *burnin));
    if (pponeemptybyclus[l]> maxpponeempty) { maxpponeempty= pponeemptybyclus[l]; }
  }
  (*pponeempty)= 0;
  for (l=1; l<=(*ncomp); l++) { (*pponeempty)+= exp(pponeemptybyclus[l] - maxpponeempty); }
  (*pponeempty)= maxpponeempty + log(*pponeempty) - log((double) *ncomp);


  free_ivector(zcount,1,*ncomp);
  free_dvector(dd,1,(*ncomp)*(*ncomp -1)/2);
  free_dvector(pponeemptybyclus,1,*ncomp); free_dvector(maxlogpempty,1,*ncomp); free_dvector(dm,1,*p); free_dvector(zprob,1,*ncomp); free_dvector(qpost,1,*ncomp); free_dvector(muz,1,*p);
  free_dmatrix(A,1,*p,1,*p); free_dmatrix(logpempty,1,*B - *burnin,1,*ncomp); free_dmatrix(xbar,1,*ncomp,1,*p); free_dmatrix(xsum,1,*ncomp,1,*p); free_dmatrix(logpclus, 1,*ncomp,1,*n);
  free_dmatrix(cholA,1,*p,1,*p); free_dmatrix(cholSigma,1,*p,1,*p); free_dmatrix(Sinv,1,*p,1,*p); free_dmatrix(cholSinv,1,*p,1,*p); free_dmatrix(cholSigmainvcur,1,*p,1,*p); free_dmatrix(cholSigmainvprior,1,*p,1,*p);
  free_darray3(S,1,*ncomp,1,*p,1,*p); free_darray3(crossprodx,1,*ncomp,1,*p,1,*p);
}


//For all clusters l, convert column sums xsum[l] into means xbar[l] and cross-prod crossprodx[l] into sum of squares S[l]
//Input
// - crossprodx[l] contains cross-products X'X in cluster l=1,...,nclus. crossprodx is an [1..nclus][1..p][1..p] array
// - xsum[l] contains colSums(X) in cluster l=1,...,nclus. xsum is an [1..nclus][1..p] matrix
// - zcount[l] contains number of individuals in cluster l=1,...,nclus
// - nclus: number of clusters
// - p: number of columns in X
//Output
// - S[l] contains sums of squares (X-colMeans(X))' (X-colMeans(X)) for cluster l=1,...,nclus. S is an [1..nclus][1..p][1..p] array
// - xbar[l] contains colMeans(X) in cluster l=1,...,nclus. xbar is an [1..nclus][1..p] matrix
void crossprod2sumsq_byclus(double ***crossprodx, double **xsum, int *zcount, int *nclus, int *p, double ***S, double **xbar) {
  int l;
  for (l=1; l<=(*nclus); l++) {
    crossprod2sumsq(crossprodx[l],xsum[l],zcount[l],*p,S[l],xbar[l],false);
  }
}


//Add (cholS %*% t(cholS))/divideby to matrix A. All matrices must be p x p. Only upper diagonal elements are computed
void addcholStcholS(double **cholS, int p, double divideby, double **A) {
  int i,j,l;
  double tmp;
  for (i=1; i<=p; i++) {
    for (j=i; j<=p; j++) {
      tmp= 0;
      for (l=1; l<=i; l++) { tmp += cholS[i][l] * cholS[j][l]; }
      A[i][j] += tmp / divideby;
    }
  }
}
