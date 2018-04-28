#ifndef MIXTURES_H
#define MIXTURES_H 1

#include <R.h>
#include <Rinternals.h>

extern "C" {

  //Mixture models
  SEXP normalmixGibbsCI(SEXP Sx, SEXP Sn, SEXP Sp, SEXP Sncomp, SEXP Sz, SEXP Smu0, SEXP Sg, SEXP Snu0, SEXP SS0, SEXP Sq, SEXP SB, SEXP Sburnin, SEXP Sverbose);

}

void normalmixGibbsC(double *pponeempty, double *eta, double *mu, double *cholSigmainv, double *x, int *n, int *p, int *ncomp, int *z, double *mu0, double *g, int *nu0, double *S0, double *q, int *B, int *burnin, int *verbose);

void crossprod2sumsq_byclus(double ***crossprodx, double **xsum, int *zcount, int *nclus, int *p, double ***S, double **xbar);

#endif /* MIXTURES_H */
