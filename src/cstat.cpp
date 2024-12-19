/***********************************************************
 Basic, input-output and matrix manipulation

 Authors: Peter Mueller, Stephen Morris, David Rossell
          (some routines obtained from other sources)
 Edits: P. Roebuck
***********************************************************/

#include "cstat.h"

//static const char interface_c_sccs_id[] = "%W%";
//static const char mess_c_sccs_id[] = "%W%";
//static const char nrutil_c_sccs_id[] = "%W%";
//static const char vector_c_sccs_id[] = "%W%";
//static const char rand_c_sccs_id[] = "@(#)$Workfile: rand.c$ $Revision: 5$";
//static const char cstat_c_sccs_id[] = "@(#)$Workfile: cstat.c$ $Revision: 2011-12-10$";


/*
 * Much of this code has undocumented assumptions, the least of which
 * is IEC-559 / IEEE-754 standard compliance.
 */


/*
 * Globals
 */
static long is1 = 123456789, is2 = 981963;
static int cstat_set = 0;

FILE *ifile, *ofile;
int nv = 0;

long Xm1, Xm2;
long Xa1, Xa2;
long Xcg1[32], Xcg2[32];
long Xa1w, Xa2w;
long Xig1[32], Xig2[32];
long Xlg1[32], Xlg2[32];
long Xa1vw, Xa2vw;
long Xqanti[32];


/*
 * Private Routines
 */

/* Same as calling stop() method in R */
static void _cstatfatal(void)
{
    Rf_error("internal error occurred in R package 'mombf'");
    /*NOTREACHED*/
}


/* Writes error message and aborts */
static void _cstaterror(const char *proc,
                        const char *act,
                        const char *what)
{
    //assert(proc != NULL);
    //assert(act != NULL);
    //assert(what != NULL);

    REprintf("\n ** Error ");
    if (proc[0] != '\0') {
        REprintf("in function '%s', ", proc);
    }
    if (act[0] != '\0') {
        REprintf("trying to %s ", act);
    }
    if (what[0] != '\0') {
        REprintf("'%s'", what);
    }
    REprintf("\n ** (from a function in 'cstat.c').\n");
    _cstatfatal();
    /*NOTREACHED*/
}


/*
 * Public Routines
 */

/************************************************************************
                         STATISTICAL FUNCTIONS
************************************************************************/

/* Sample mean of elements 0 through lim of vector x */
double meani(const int *x,
             int lim)
{
    int i;
    double value = 0.0;

    for (i = 0; i <= lim; i++) { value += x[i]; }
    value *= 1.0 / (lim+1.0);
    return value;
}


double wmeani(const int *x,
              int lim,
              const double *w)
{
    int i;
    double value = 0.0;
    double wtot = 0.0;

    for (i = 0; i <= lim; i++) {
        value += w[i] * x[i];
        wtot += w[i];
    }
    value *= 1.0 / wtot;
    return value;
}


/* Sample mean of elements 0 through lim (both included) of vector x */
double meanx(const double *x,
             int lim)
{
    int i;
    double value = 0.0;

    for (i = 0; i <= lim; i++) { value += x[i]; }
    value *= 1.0 / (lim+1.0);
    return value;
}


double wmeanx(const double *x,
              int lim,
              const double *w)
{
    int i;
    double value = 0.0;
    double wtot = 0.0;

    for (i = 0; i <= lim; i++) {
        value += w[i] * x[i];
        wtot += w[i];
    }
    value *= 1.0 / wtot;
    return value;
}


/* Sample variance of elements 0 through lim of vector x */
double vari(const int *x,
            int lim,
            bool unbiased)
{
    int i;
    double value = 0.0;

    for (i = 0; i <= lim; i++) { value += pow((double) x[i], 2.0) / (1.0+lim);  }
    value -= pow(meani(x, lim), 2);
    if (unbiased && lim>0) {
        value *= (1.0+lim) / (0.0+lim);
    }
    return value;
}


double wvari(const int *x,
             int lim,
             const double *w)
{
    int i;
    double value = 0.0;
    double wtot = 0.0;

    for (i = 0; i <= lim; i++) {
      value += w[i] * pow((double) x[i], 2.0);
      wtot += w[i];
    }
    value = value/wtot - pow(wmeani(x, lim, w), 2);
    return value;
}


/* Sample variance of elements 0 through lim of vector x */
double varx(const double *x,
            int lim,
            bool unbiased)
{
    int i;
    double value = 0.0;

    for (i = 0; i <= lim; i++) { value += pow(x[i], 2) / (1.0+lim); }
    value -= pow(meanx(x, lim), 2);
    if (unbiased && lim>0) {
        value *= (1.0+lim) / (lim+0.0);
    }
    return value;
}


double wvarx(const double *x,
             int lim,
             const double *w)
{
    int i;
    double value = 0.0;
    double wtot = 0.0;

    for (i = 0; i <= lim; i++) {
        value += w[i] * pow(x[i], 2);
        wtot += w[i];
    }
    value = value/wtot - pow(wmeanx(x, lim, w), 2);
    return value;
}


/* Compute coefficient of variation of x[ini..fi] */
double cv(const double *x,
          int ini,
          int fi)
{
    int i;
    double m = 0.0;
    double s = 0.0;
    double ans;

    for (i = ini; i <= fi; i++) {
        double value;

        value = x[i];
        m += value;
        s += value * value;
    }
    m = m / (1.0+fi-ini);
    s = s / (0.0+fi-ini) - m*m*(1.0+fi-ini)/(0.0+fi-ini);
    ans = sqrt(s) / m;
    return(ans);
}


/* Compute inverse coefficient of variation of x[ini..fi] */
double cvinv(const double *x,
             int ini,
             int fi)
{
    int i;
    double m = 0.0;
    double s = 0.0;
    double ans;

    for (i = ini; i <= fi; i++) {
        double value;

        value = x[i];
        m += 1.0 / value;
        s += 1.0 / (value * value);
    }
    m = m / (1.0+fi-ini);
    s = s / (0.0+fi-ini) - m*m*(1.0+fi-ini)/(0.0+fi-ini);
    ans = sqrt(s) / m;
    return(ans);
}


/* Column means */
/* x is assumed to be in row order (x[0], x[1] ... x[ncol-1] are elem in 1st row) */
void colMeans(double *m,
              const double *x,
              int nrow,
              int ncol)
{
    int i, j;

    for (j = 0; j < ncol; j++) {
        m[j] = 0.0;
    }
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            m[j] += x[i*ncol+j];
        }
    }
    for (j = 0; j < ncol; j++) {
        m[j] = m[j] / (nrow+0.0);
    }
}


/* x is assumed to be in row order (x[0], x[1] ... x[ncol-1] are elem in 1st row) */
void colVar(double *v,
            const double *x,
            int nrow,
            int ncol)
{
    int i, j;
    double *m;
    double *m2;

    m  = dvector(0, ncol-1);
    m2 = dvector(0, ncol-1);

    for (j = 0; j < ncol; j++) {
        m[j] = m2[j] = 0.0;
    }
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            double value;

            value = x[i*ncol+j];
            m[j]  += value;
            m2[j] += value * value;
        }
    }
    for (j = 0; j < ncol; j++) {
        m[j] = m[j] / (0.0+nrow);    /* :TBD: unnecessary addition - change to cast? */
        v[j] = m2[j] / (nrow-1.0) - m[j]*m[j]*(nrow+0.0)/(nrow-1.0);
    }

    free_dvector(m,  0, ncol-1);
    free_dvector(m2, 0, ncol-1);
}


/* x is assumed to be in row order (x[0], x[1] ... x[ncol-1] are elem in 1st row) */
void colCV(double *cv,
           const double *x,
           int nrow,
           int ncol)
{
    int i;
    int j;
    double *m;
    double *s;

    m = dvector(0, ncol);
    s = dvector(0, ncol);

    for (j = 0; j < ncol; j++) {
        m[j] = s[j] = 0.0;
    }
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            double value;

            value = x[i*ncol+j];
            m[j] += value;
            s[j] += value * value;
        }
    }
    for (j = 0; j < ncol; j++) {
        m[j] = m[j] / (nrow+0.0);
        s[j] = s[j] / (nrow-1.0) - m[j]*m[j]*(nrow+0.0)/(nrow-1.0);
        cv[j] = sqrt(s[j]) / m[j];
    }

    free_dvector(m, 0, ncol);
    free_dvector(s, 0, ncol);
}


/* x is assumed to be in row order (x[0], x[1] ... x[ncol-1] are elem in 1st row) */
void colCVinv(double *cv,
              const double *x,
              int nrow,
              int ncol)
{
    int i;
    int j;
    double *m;
    double *s;

    m = dvector(0, ncol);
    s = dvector(0, ncol);

    for (j = 0; j < ncol; j++) {
        m[j] = s[j] = 0.0;
    }
    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            double value;

            value = x[i*ncol+j];
            m[j] += 1.0 / value;
            s[j] += 1.0 / (value * value);
        }
    }
    for (j = 0; j < ncol; j++) {
        m[j] = m[j] / (nrow+0.0);
        s[j] = s[j] / (nrow-1.0) - m[j]*m[j]*(nrow+0.0)/(nrow-1.0);
        cv[j] = sqrt(s[j]) / m[j];
    }

    free_dvector(m, 0, ncol);
    free_dvector(s, 0, ncol);
}


//Sample covariance matrix (divided by n).
//Input:
// - x is an n * p matrix arranged by columns (0-indexed).
// - n,p: dimensions of x
// - lowertri: if false only diagonal and upper-diagonal elements of S are returned
//Output: S[1][1],...,S[p][p]
void covxvec(double *x, int n, int p, bool lowertri, double **S) {
  int i, j, l, idxj, idxl;
  double *xbar, *m2;
  xbar= dvector(1,p); m2= dvector(1,p);
  //Variances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    xbar[j]= 0;
    for (i=0; i<n; i++) {
      xbar[j]+= x[i+idxj];
      m2[j]+= x[i+idxj]*x[i+idxj];
    }
    xbar[j]= xbar[j]/((double) n);
    S[j][j]= m2[j]/((double) n) - xbar[j]*xbar[j];
  }
  //Covariances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    for (l=j+1; l<=p; l++) {
      idxl= (l-1)*n;
      S[j][l]= 0;
      for (i=0; i<n; i++) { S[j][l]+= x[i+idxj]*x[i+idxl]; }
      S[j][l]= S[j][l]/((double) n) - xbar[j]*xbar[l];
    }
  }
  if (lowertri) {
    for (j=1; j<=p; j++) { for (l=1; l<j; l++) { S[l][j]= S[j][l]; } }
  }
  free_dvector(xbar,1,p); free_dvector(m2,1,p);
}



//Matrix with sum of squares and cross-product, that is n * cov(x)
//Input: x is an n * p matrix arranged by columns (0-indexed).
//Output: stored in S[1][1],...,S[p][p]
void sumsq(double *x, int n, int p, bool lowertri, double **S) {
  int i, j, l, idxj, idxl;
  double *xbar, *m2;
  xbar= dvector(1,p); m2= dvector(1,p);
  //Variances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    xbar[j]= 0;
    for (i=0; i<n; i++) {
      xbar[j]+= x[i+idxj];
      m2[j]+= x[i+idxj]*x[i+idxj];
    }
    xbar[j]= xbar[j]/((double) n);
    S[j][j]= m2[j] - ((double) n)*xbar[j]*xbar[j];
  }
  //Covariances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    for (l=j+1; l<=p; l++) {
      idxl= (l-1)*n;
      S[j][l]= 0;
      for (i=0; i<n; i++) { S[j][l]+= x[i+idxj]*x[i+idxl]; }
      S[j][l]= S[j][l] - ((double) n)*xbar[j]*xbar[l];
    }
  }
  if (lowertri) {
    for (j=1; j<=p; j++) { for (l=1; l<j; l++) { S[l][j]= S[j][l]; } }
  }
  free_dvector(xbar,1,p); free_dvector(m2,1,p);
}

//Sum of squares and cross-products matrix by clusters indicated in z[0],...,z[n-1]. Each z[i] must be in [1,nclus]
//Input:
// - x is an n * p matrix arranged by columns (0-indexed).
// - n,p: dimensions of x
// - z: cluster allocations (0-indexed)
// - nclus: number of clusters
// - lowertri: if false only diagonal and upper-diagonal elements of S are returned
// Output: S[l] is a p x p matrix containing sum of squares and cross-products for individuals with z[i]=l
void sumsqbyclus(double *x, int n, int p, int *z, int nclus, bool lowertri, double ***S) {
  int i, j, l, kk, idxj, idxl, *zcount;
  double **xbar, **m2;
  zcount= ivector(1,nclus);
  xbar= dmatrix(1,nclus,1,p); m2= dmatrix(1,nclus,1,p);
  for (kk=1; kk<=nclus; kk++) { zcount[kk]= 0; }
  for (i=0; i<n; i++) { (zcount[z[i]])++; }
  //Variances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    for (kk=1; kk<=nclus; kk++) { xbar[kk][j]= 0; }
    for (i=0; i<n; i++) {
      xbar[z[i]][j]+= x[i+idxj];
      m2[z[i]][j]+= x[i+idxj]*x[i+idxj];
    }
    for (kk=1; kk<=nclus; kk++) {
      xbar[kk][j]= xbar[kk][j]/((double) zcount[kk]);
      S[kk][j][j]= m2[kk][j] - ((double) zcount[kk])*xbar[kk][j]*xbar[kk][j];
    }
  }
  //Covariances
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    for (l=j+1; l<=p; l++) {
      idxl= (l-1)*n;
      for (kk=1; kk<=nclus; kk++) { S[kk][j][l]= 0; }
      for (i=0; i<n; i++) { S[z[i]][j][l]+= x[i+idxj]*x[i+idxl]; }
      for (kk=1; kk<=nclus; kk++) { S[kk][j][l]= S[kk][j][l] - ((double) zcount[kk])*xbar[kk][j]*xbar[kk][l]; }
    }
  }
  if (lowertri) {
    for (j=1; j<=p; j++) { for (l=1; l<j; l++) { for (kk=1; kk<=nclus; kk++) { S[kk][l][j]= S[kk][j][l]; } } }
  }
  free_ivector(zcount,1,nclus);  free_dmatrix(xbar,1,nclus,1,p); free_dmatrix(m2,1,nclus,1,p);
}


//Column sums colSums(x) and cross-products X'X
//Input:
// - x is an n * p matrix arranged by columns (0-indexed)
// - lowertri: if false only diagonal and upper-diagonal elements of S are returned
//
//Output:
// - Column sums stored in sums[1],...,sums[n]
// - Cross-products stored in S[1][1],...,S[p][p]
void sum_crossprod(double *x, int n, int p, bool lowertri, double *sumx, double **crossprodx) {
  int i, j, l, idxj, idxl;
  for (j=1; j<=p; j++) {
    idxj= (j-1)*n;
    sumx[j]= crossprodx[j][j]= 0;
    for (i=0; i<n; i++) {
      sumx[j]+= x[i+idxj];
      crossprodx[j][j]+= x[i+idxj]*x[i+idxj];
    }
    for (l=j+1; l<=p; l++) {
      idxl= (l-1)*n;
      crossprodx[j][l]= 0;
      for (i=0; i<n; i++) { crossprodx[j][l]+= x[i+idxj]*x[i+idxl]; }
    }
  }
  if (lowertri) {
    for (j=1; j<=p; j++) { for (l=1; l<j; l++) { crossprodx[l][j]= crossprodx[j][l]; } }
  }
}


//Convert column sums xsum into means xbar and cross-prod crossprodx into sum of squares S
//Input
// - crossprodx: matrix indexed [1..p][1..p] containing cross-products X'X.
// - xsum: vector indexed [1..p] containing colSums(X)
// - n: number of rows in X
// - p: number of columns in X
// - lowertri: if false only diagonal and upper-diagonal elements of S are returned
//Output
// - S: a p x p matrix containing sums of squares (X-colMeans(X))' (X-colMeans(X))
// - xbar: a p-vector containing colMeans(X)
//
//Note: if n==0, then xbar and S are set to 0
void crossprod2sumsq(double **crossprodx, double *xsum, int n, int p, double **S, double *xbar, bool lowertri) {
  int j, jj;
  if (n>0) {
    for (j=1; j<=p; j++) {
      xbar[j]= xsum[j] / ((double) n);
      S[j][j]= crossprodx[j][j] - xbar[j] * xsum[j];
      for (jj=j+1; jj<=p; jj++) {
        S[j][jj]= crossprodx[j][jj] - xbar[j] * xsum[jj];
      }
    }
  } else {
    for (j=1; j<=p; j++) {
      xbar[j]= S[j][j]= 0;
      for (jj=j+1; jj<=p; jj++) { S[j][jj]= 0; }
    }
  }
  if (lowertri) {
    for (j=1; j<p; j++) { for (jj=1; jj<j; jj++) { S[j][jj]= S[jj][j]; } }
  }
}




/************************************************************************
                         BASIC BAYESIAN MODELS
************************************************************************/

/********************************************
 *         normal_normal
 ********************************************/
/*
 * prior: N(x; mpr, r1*Spr)
 * likl:  N(y; x, r2*Slik)
 * p: dimensionality
 * returns:  post N(x; mpo,Spo)
 * Spo = (1/r1*Spr_inv + 1/r2*Slik_inv)^-1
 * mpo = Spo*(1/r1*Spr*mpr + 1/r2*Slik*y)
 * NOTE: input vectors and matrices must start at position 1, not 0
 */
void nn_bayes(double *mpo,
              double **Spo,
              double **Spo_inv,
              int p,
              double r1,
              double *mpr,
              double **Spr_inv,
              double r2,
              double *y,
              double **Slik_inv)
{
  bool posdef;
  double *z;

    z = dvector(1, p);

    rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, Spo_inv, 1, p, 1, p);
    inv_posdef(Spo_inv, p, Spo, &posdef);
    rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z, 1, p, 1, p);
    Ax(Spo, z, mpo, 1, p, 1, p);

    free_dvector(z, 1, p);
}


/*
 * Same as nn_bayes, but returns a draw only
 * returns:
 * theta: draw from the posterior N(theta; mpo, Spo):
 */
void nn_bayes_rand(double *theta,
                   int p,
                   double r1,
                   double **Spr_inv,
                   double *mpr,
                   double r2,
                   double **Slik_inv,
                   double *y)
{
    bool posdef;
    double *z;
    double *m;
    double **S;
    double **S_inv;
    double **cholS;

    /* Allocate memory */
    z = dvector(0, p-1);
    m = dvector(0, p-1);
    S = dmatrix(0, p-1, 0, p-1);
    S_inv = dmatrix(0, p-1, 0, p-1);
    cholS = dmatrix(0, p-1, 0, p-1);

    rA_plus_sB(1.0/r1, Spr_inv, 1.0/r2, Slik_inv, S_inv, 1, p, 1, p);
    inv_posdef(S_inv, p, S, &posdef);
    rAx_plus_sBy(1.0/r1, Spr_inv, mpr, 1.0/r2, Slik_inv, y, z, 1, p, 1, p);
    Ax(S, z, m, 1, p, 1, p);

    choldc(S, p, cholS, &posdef);
    rmvnormC(theta, p, m, cholS);

    free_dvector(z, 0, p-1);
    free_dvector(m, 0, p-1);
    free_dmatrix(S, 0, p-1, 0, p-1);
    free_dmatrix(S_inv, 0, p-1, 0, p-1);
    free_dmatrix(cholS, 0, p-1, 0, p-1);
}


/*
 * Compute normal-normal integral
 *   N(x;beta,rx*Vx) * N(beta;mpr,rpr*Vpr) with respect to beta
 *
 *    Vxinv, detVx   - inverse and determinant of Vx
 *    Vprinv, detVpr - inverse and determinant of Vpr
 *    p              - dimensionality
 *    logscale       - if not 0, result is returned in log scale
 */
double nn_integral(const double *x,
                   const double *rx,
                   double **Vxinv,
                   const double *detVx,
                   const double *mpr,
                   const double *rpr,
                   double **Vprinv,
                   const double *detVpr,
                   const int *p,
                   const int *logscale)
{
    bool posdef;
    double detsum;
    double **Vsum;
    double **Vsuminv;
    double **cholVsum;
    double *m;
    double ans;

    m = dvector(1, *p);
    Vsum = dmatrix(1, *p, 1, *p);
    Vsuminv = dmatrix(1, *p, 1, *p);
    cholVsum = dmatrix(1, *p, 1, *p);

    rA_plus_sB(1.0/(*rx), Vxinv, 1.0/(*rpr), Vprinv, Vsuminv, 1, *p, 1, *p);
    choldc_inv(Vsuminv, *p, cholVsum, &posdef);
    detsum = choldc_det(cholVsum, *p);
    inv_posdef_chol(cholVsum, *p, Vsum);
    rAx_plus_sBy(1.0/(*rx), Vxinv, x, 1.0/(*rpr), Vprinv, mpr, m, 1, *p, 1, *p);
    ans = xtAy(m, Vsum, m, 1, *p) -
          xtAy(x, Vxinv, x, 1, *p) -
          xtAy(mpr, Vprinv, mpr, 1, *p);

    ans = 0.5 * ans -
          0.5 * ((*p+0.0)*LOG_M_2PI + log(*detVx) + log(*detVpr) - log(detsum));
    /* :BUG?: Returns exponential when logscale flag enabled */
    if (*logscale != 0) {
        ans = exp(ans);
    }

    free_dvector(m, 1, *p);
    free_dmatrix(Vsum, 1, *p, 1, *p);
    free_dmatrix(Vsuminv, 1, *p, 1, *p);
    free_dmatrix(cholVsum, 1, *p, 1, *p);

    return(ans);
}


/*
 * Fits classical multiple linear regression
 *
 * Input
 *    - y: response variable y[1..n]
 *    - X: design matrix X[1..n][1..p]
 *    - n: number of observations
 *    - p: number of covariates
 *    - useXtX: if set to 0 the inverse of X'X is computed,
 *              otherwise the supplied value is used
 *  Ouput
 *    - b: least-squares estimate for regression coefficients
 *    - XtX, invXtX: X'X and its inverse (if useXtX==0 they're ouput param,
 *                                        otherwise they're input)
 *    - Xty: vector X'y (if useXtX==0 it's an output param, otherwise it's input)
 *    - s: residual variance (dividing by n-p)
 *    - ypred: predicted values i.e. X'b
 */
void lm(double *b,
        double **XtX,
        double **invXtX,
        double *Xty,
        double *s,
        double *ypred,
        const double *y,
        double **X,
        const int *n,
        const int *p,
        const int *useXtX)
{
  bool posdef;
  int i;

    //assert(b != NULL);
    //assert(XtX != NULL);
    //assert(invXtX != NULL);
    //assert(Xty != NULL);
    //assert(s != NULL);
    //assert(ypred != NULL);
    //assert(y != NULL);
    //assert(X != NULL);
    //assert(n != NULL);
    //assert(p != NULL);
    //assert(useXtX != NULL);

    if (*n < *p) {
        errorC("lm", "linear model with more variables than observations", 0);
        /*NOTREACHED*/
    }

    if (*useXtX == 0) {
        AtB(X, 1, *n, 1, *p, X, 1, *n, 1, *p, XtX);
        inv_posdef(XtX, *p, invXtX,&posdef);
        Atx(X, y, Xty, 1, *n, 1, *p); //X'y
    }

    Ax(invXtX, Xty, b, 1, *p, 1, *p); /* least squares estimate */
    Ax(X, b, ypred, 1, *n, 1, *p);    /* predicted values */

    (*s) = 0.0;
    for (i = 1; i <= (*n); i++) {
        double value;

        value = y[i] - ypred[i];
        (*s) += value * value;
    }
    (*s) = (*s) / (*n - *p);
}


/*
 * Bayesian conjugate multiple linear regression
 *    y ~ N(X'beta,sigma^2)
 *    beta ~ N(mpr,sigma^2*Spr) (if tauprior<=0)
 *    beta ~ N(mpr,tauprior*sigma^2*(X'X)^{-1}) (if tauprior>0)
 *           e.g. tauprior==n gives unit information prior
 *    sigma^2 ~ IG(0.5*nu0,0.5*s0)  (nu0: prior sample size; s0: prior sum of squares)
 * Input
 *      B -  number of posterior samples to draw (can be 0)
 *      y -  response variable y[1..n]
 *      X -  design matrix X[1..n][1..p]
 *      n -  number of observations
 *      p -  number of covariates
 *      useXtX -  if set to 0 the inverse of X'X is computed,
 *                otherwise the supplied value is used
 *      mpr, Spr_inv, tauprior -  prior parameters for beta
 *      nu0, s0 -  prior for sigma^2 is IG(0.5*nu0,0.5*s0)
 *  Ouput
 *      bpost -  matrix (B rows, p cols) with samples from the posterior
 *               of beta. Starts at bpost[1].
 *      spost -  vector (B rows) with samples from the posterior of sigma^2.
 *               Starts at spost[1].
 *      b, Vb -  posterior for regression coef is N(b,sigma^2*Vb)
 *      a_s, b_s -  posterior for sigma^2 is IG(a_s,b_s)
 *      XtX, invXtX -  X'X and its inverse (if useXtX==0 they're ouput param,
 *                                          otherwise they're input)
 *      Xty -  vector X'y (if useXtX==0 it's an output param, otherwise it's input)
 */
void lmbayes(double *bpost,
             double *spost,
             double *b,
             double **Vb,
             double *a_s,
             double *b_s,
             double **XtX,
             double **invXtX,
             double *Xty,
             int *B,
             double *y,
             double **X,
             int *n,
             int *p,
             int *useXtX,
             double *mpr,
             double **Spr_inv,
             double *tauprior,
             double *nu0,
             double *s0)
{
  bool posdef;
  double *b_ls;
  double s_ls;
  double **Vb_inv;

    //assert(bpost != NULL);
    //assert(spost != NULL);
    //assert(b != NULL);
    //assert(Vb != NULL);
    //assert(a_s != NULL);
    //assert(b_s != NULL);
    //assert(XtX != NULL);
    //assert(invXtX != NULL);
    //assert(Xty != NULL);
    //assert(B != NULL);
    //assert(y != NULL);
    //assert(X != NULL);
    //assert(n != NULL);
    //assert(p != NULL);
    //assert(useXtX != NULL);
    //assert(mpr != NULL);
    //assert(Spr_inv != NULL);
    //assert(tauprior != NULL);
    //assert(nu0 != NULL);
    //assert(s0 != NULL);

    if (*useXtX == 0) {
        AtB(X, 1, *n, 1, *p, X, 1, *n, 1, *p, XtX);
        inv_posdef(XtX, *p, invXtX, &posdef);
        Atx(X, y, Xty, 1, *n, 1, *p); //X'y
    }

    b_ls = dvector(1, *p);

    {
        /*const*/ int one = 1;
        double *ypred;

        ypred = dvector(1, *n);
        lm(b_ls, XtX, invXtX, Xty, &s_ls, ypred, y, X, n, p, &one);  //least-squares fit
        free_dvector(ypred, 1, *n);
    }

    /* Posterior for sigma^2 */
    *a_s = 0.5 * (*nu0 + *n);
    *b_s = 0.5 * (*s0 + (*n - *p)*s_ls);

    /* Posterior for beta */
    Vb_inv = dmatrix(1, *p, 1, *p);
    if (*tauprior > 0) {
        nn_bayes(b, Vb, Vb_inv, *p, *tauprior, mpr, XtX, 1.0, b_ls, XtX);
    } else {
        nn_bayes(b, Vb, Vb_inv, *p, 1.0, mpr, Spr_inv, 1.0, b_ls, XtX);
    }

    if (*B > 0) {             //posterior samples
        bool posdef;
        double *zeroes, **cholVb;
        int i, j;

        cholVb = dmatrix(1, *p, 1, *p);
        choldc(Vb, *p, cholVb,&posdef); //cholesky decomp of posterior covar for beta
        zeroes = dvector(1, *p);
        for (i = 1; i <= (*p); i++) {
            zeroes[i] = 0.0;
        }
        for (i = 1; i <= (*B); i++) {
            spost[i] = 1.0 / rgammaC(*a_s, *b_s);
            rmvnormC(bpost+(i-1)*(*p), *p, zeroes, cholVb);
            for (j = 1; j <= (*p); j++) {
                bpost[(i-1)*(*p)+j] = bpost[(i-1)*(*p)+j]*sqrt(spost[i])+b[j];
            }
        }
        free_dvector(zeroes, 1, *p);
        free_dmatrix(cholVb, 1, *p, 1, *p);
    }

    free_dvector(b_ls, 1, *p);
    free_dmatrix(Vb_inv, 1, *p, 1, *p);
}


/*
 * Bayesian conjugate multiple linear regression with known var
 *  y ~ N(X'beta,sigma^2)
 *  beta ~ N(mpr,sigma^2*Spr) (if tauprior<=0)
 *  beta ~ N(mpr,tauprior*sigma^2*(X'X)^{-1}) (if tauprior>0)
 *         e.g. tauprior==n gives unit information prior
 * Input
 *    sigma: residual standard deviation
 *    B: number of posterior samples to draw (can be 0)
 *    y: response variable y[1..n]
 *    X: design matrix X[1..n][1..p]
 *    n: number of observations
 *    p: number of covariates
 *    useXtX: if set to 0 the inverse of X'X is computed,
 *            otherwise the supplied value is used
 *    mpr, Spr_inv, tauprior: prior parameters for beta
 *  Ouput
 *    bpost: matrix (B rows, p cols) with samples from the posterior
 *           of beta. Starts at bpost[1].
 *    b, Vb: posterior for regression coef is N(b,sigma^2*Vb)
 *    XtX, invXtX: X'X and its inverse (if useXtX==0 they're ouput param,
 *                                      otherwise they're input)
 *    Xty: vector X'y (if useXtX==0 it's an output param, otherwise it's input)
 */
void lmbayes_knownvar(double *bpost,
                      double *b,
                      double **Vb,
                      double **XtX,
                      double **invXtX,
                      double *Xty,
                      double *sigma,
                      int *B,
                      double *y,
                      double **X,
                      int *n,
                      int *p,
                      int *useXtX,
                      double *mpr,
                      double **Spr_inv,
                      double *tauprior)
{
  bool posdef;
    double *b_ls;
    double s_ls;
    double **Vb_inv;

    //assert(bpost != NULL);
    //assert(b != NULL);
    //assert(Vb != NULL);
    //assert(XtX != NULL);
    //assert(invXtX != NULL);
    //assert(Xty != NULL);
    //assert(sigma != NULL);
    //assert(B != NULL);
    //assert(y != NULL);
    //assert(X != NULL);
    //assert(n != NULL);
    //assert(p != NULL);
    //assert(useXtX != NULL);
    //assert(mpr != NULL);
    //assert(Spr_inv != NULL);
    //assert(tauprior != NULL);

    if (*useXtX == 0) {
        AtB(X, 1, *n, 1, *p, X, 1, *n, 1, *p, XtX);
        inv_posdef(XtX, *p, invXtX, &posdef);
        Atx(X, y, Xty, 1, *n, 1, *p); //X'y
    }

    b_ls = dvector(1, *p);

    {
        double *ypred;
        /*const*/ int one = 1;

        ypred = dvector(1, *n);
        lm(b_ls, XtX, invXtX, Xty, &s_ls, ypred, y, X, n, p, &one);  //least-squares fit
        free_dvector(ypred, 1, *n);
    }

    Vb_inv = dmatrix(1, *p, 1, *p);   //posterior for beta
    if (*tauprior > 0) {
        nn_bayes(b, Vb, Vb_inv, *p, *tauprior, mpr, XtX, 1.0, b_ls, XtX);
    } else {
        nn_bayes(b, Vb, Vb_inv, *p, 1.0, mpr, Spr_inv, 1.0, b_ls, XtX);
    }

    if (*B > 0) {             //posterior samples
        bool posdef;
        double *zeroes, **cholVb;
        int i;

        cholVb = dmatrix(1, *p, 1, *p);
        choldc(Vb, *p, cholVb, &posdef); //cholesky decomp of posterior covar for beta
        zeroes = dvector(1, *p);
        for (i = 1; i <= (*p); i++) {
            zeroes[i] = 0.0;
        }
        for (i = 1; i <= (*B); i++) {
            int j;

            rmvnormC(bpost+(i-1)*(*p), *p, zeroes, cholVb);
            for (j = 1; j <= (*p); j++) {
                bpost[(i-1)*(*p)+j] = bpost[(i-1)*(*p)+j]*(*sigma)+b[j];
            }
        }
        free_dvector(zeroes, 1, *p);
        free_dmatrix(cholVb, 1, *p, 1, *p);
    }

    free_dvector(b_ls, 1, *p);
    free_dmatrix(Vb_inv, 1, *p, 1, *p);
}



/******************************************************************************
                            MESSAGES AND ERROR HANDLING
******************************************************************************/

void print_iterprogress(int *iter, int *niter, int *everyiter) {
  if ((*iter % *everyiter) == 0) {
    Rcout << "\r";  //Carriage return (deletes current line)
    Rcout << 100 * (*iter) / (*niter) << "%"; //% iteration progress
  }
}




void nrerror(const char *proc,
             const char *act,
             const char *what)
{
    _cstaterror(proc, act, what);
    /*NOTREACHED*/
}


void errorC(const char *module,
            const char *msg,
            int nr)
{
    REprintf("\n *** ERROR # %d in %s***\n %s\n", nr, module, msg);
    _cstatfatal();
    /*NOTREACHED*/
}





/******************************************************************************
                          MEMORY ALLOCATION
******************************************************************************/

/* Allocate int vector with subscript range v[nl..nh] */
int *ivector(int nl,
             int nh)
{
    int *v;
    size_t count = nh-nl+1;

    //assert(count >= 0);

    nv += count;
    v = (int *) calloc(count, sizeof(int));
    if (v == NULL) {
        nrerror("ivector", "allocate an int vector", "");
        /*NOTREACHED*/
    }
    return v-nl;
}


/* Allocate char vector with subscript range v[nl..nh] */
char *charvector(int nl,
                int nh)
{
    char  *v;
    size_t count = nh-nl+1;

    //assert(count >= 0);

    //v = malloc( sizeof(char) * ( count ) );
    v = (char *) calloc(count, sizeof(char));
    if (v == NULL) {
        nrerror("charvector", "allocate a character vector", "");
        /*NOTREACHED*/
    }
    return v-nl;
}


float *vector(int nl,
              int nh)
{
    float *v;
    size_t count = nh-nl+1;

    //assert(count >= 0);

    v = (float *) calloc(count, sizeof(float));
    if (v == NULL) {
        nrerror("vector", "allocate a float vector", "");
        /*NOTREACHED*/
    }
    return v-nl;
}


/* Allocate double vector with subscript range v[nl..nh] */
double *dvector(int nl,
                int nh)
{
    double  *v;
    size_t count = nh-nl+1;

    //assert(count >= 0);

    nv += count;
    v = (double *) calloc(count, sizeof(double));
    if (v == NULL) {
        nrerror("dvector", "allocate a double vector", "");
        /*NOTREACHED*/
    }
    return v-nl;
}


/* Allocate int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **imatrix(int nrl,
              int nrh,
              int ncl,
              int nch)
{
    int **m;
    size_t nrow = nrh-nrl+1;
    size_t ncol = nch-ncl+1;
    int i;

    //assert(nrow >= 0);
    //assert(ncol >= 0);

    nv += nrow * ncol;

    /* Allocate pointers to rows */
    m = (int **) calloc(nrow, sizeof(int *));
    if (m == NULL) {
        nrerror("imatrix", "allocate an int matrix (1st dim)", "");
        /*NOTREACHED*/
    }
    m -= nrl;

    /* For each row pointer... */
    for (i = nrl; i <= nrh; i++) {
        /* Allocate columns for individual row */
        m[i] = (int *) calloc(ncol, sizeof(int));
        if (m[i] == NULL) {
            nrerror("imatrix", "allocate an int matrix (2nd dim)", "");
            /*NOTREACHED*/
        }
        m[i] -= ncl;
    }
    return m;
}


/* Allocate double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl, int nrh, int ncl, int nch) {
    double **m;
    size_t nrow = nrh-nrl+1;
    size_t ncol = nch-ncl+1;
    int i;

    nv += nrow * ncol;

    /* Allocate pointers to rows */
    m = (double **) calloc(nrow, sizeof(double *));
    if (m == NULL) {
        nrerror("dmatrix", "allocate a double matrix (1st dim)", "");
        /*NOTREACHED*/
    }
    m -= nrl;

    /* For each row pointer... */
    for (i = nrl; i <= nrh; i++) {
        /* Allocate columns for individual row */
        m[i] = (double *) calloc(ncol, sizeof(double));
        if (m[i] == NULL) {
            nrerror("dmatrix", "allocate a double matrix (2nd dim)", "");
            /*NOTREACHED*/
        }
        m[i] -= ncl;
    }
    return m;
}


/* Allocates int array with subscript range a[0..n1-1][0..n2-1][0..n3-1] */
int ***iarray3(int n1, int n2, int n3) {
  int ***a, i, j;

    /* Allocate pointers to first dimension */
    a = (int ***) calloc(n1, sizeof(int **));
    if (a == NULL) {
        nrerror("iarray3", "allocate a 3dim int array (1st dim)", "");
        /*NOTREACHED*/
    }

    /* Allocate pointers to second dimension and set pointers */
    a[0] = (int **) calloc(n1 * n2, sizeof(int *));
    if (a[0] == NULL) {
        nrerror("iarray3", "allocate a 3dim int array (2nd dim)", "");
        /*NOTREACHED*/
    }
    for (i = 1; i < n1; i++) {
        a[i] = a[i-1] + n2;
    }

    /* Allocate pointers to third dimension and set pointers */
    a[0][0] = (int *) calloc(n1 * n2 * n3, sizeof(int));
    if (a[0][0] == NULL) {
        nrerror("iarray3", "allocate a 3dim int array (3rd dim)", "");
        /*NOTREACHED*/
    }
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
            a[i][j] = a[0][0] + n2*n3*i + j*n3;
        }
    }
    return a;
}

/* Allocates 3-way double array with subscript range [ini1..fi1][ini2..fi2][ini3..fi3] */
double ***darray3(int ini1, int fi1, int ini2, int fi2, int ini3, int fi3) {
    double ***a;
    int i;
    size_t ndim1 = fi1-ini1+1;

    /* Allocate pointers to first dimension */
    a = (double ***) calloc(ndim1, sizeof(double **));
    if (a == NULL) {
        nrerror("darray3", "allocate a 3dim double array (1st dim)", "");
        /*NOTREACHED*/
    }
    a -= ini1;

    for (i=ini1; i<=fi1; i++) a[i]= dmatrix(ini2,fi2,ini3,fi3);
    return a;
}

/* OLD VERSION OF darray3.
// Allocates double array with subscript range a[0..n1-1][0..n2-1][0..n3-1] */
// double ***darray3(int n1, int n2, int n3) {
//    double ***a;
//    int i, j;
//
//    /* Allocate pointers to first dimension */
//    a = (double ***) calloc(n1, sizeof(double **));
//    if (a == NULL) {
//        nrerror("darray3", "allocate a 3dim double array (1st dim)", "");
//        /*NOTREACHED*/
//    }
//
//    /* Allocate pointers to second dimension and set pointers */
//    a[0] = (double **) calloc(n1 * n2, sizeof(double *));
//    if (a[0] == NULL) {
//        nrerror("darray3", "allocate a 3dim double array (2nd dim)", "");
//        /*NOTREACHED*/
//    }
//    for (i = 1; i < n1; i++) {
//        a[i] = a[i-1] + n2;
//    }
//
//    /* Allocate pointers to third dimension and set pointers */
//    a[0][0] = (double *) calloc(n1 * n2 * n3, sizeof(double));
//    if (a[0][0] == NULL) {
//        nrerror("darray3", "allocate a 3dim double array (3rd dim)", "");
//        /*NOTREACHED*/
//    }
//    for (i = 0; i < n1; i++) {
//        for (j = 0; j < n2; j++) {
//            a[i][j] = a[0][0] + n2*n3*i + j*n3;
//        }
//    }
//    return a;
//}



/* Free char vector allocated with charvector() */
void free_charvector(char *v,
                  int nl,
                  int nh)
{

    //if ((v+nl) != NULL) {
        free((char  *) (v+nl));
    //}
}

/* Free int vector allocated with ivector() */
void free_ivector(int *v,
                  int nl,
                  int nh)
{

    //if ((v+nl) != NULL) {
        free((char  *) (v+nl));
    //}
    nv -= (nh-nl+1);
}


/* Free float vector allocated with vector() */
void free_vector(float *v,
                 int nl,
                 int nh)
{

    //if ((v+nl) != NULL) {
        free((char  *) (v+nl));
    //}
    nv -= (nh-nl+1);
}


/* Free double vector allocated with dvector() */
void free_dvector(double *v, int nl, int nh) {

    //if ((v+nl) != NULL) {
        free((char  *) (v+nl));
    //}
    nv -= (nh-nl+1);
}


/* Free int matrix allocated by imatrix() */
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch) {
    int i;
    size_t nrow = nrh-nrl+1;
    size_t ncol = nch-ncl+1;

    for (i = nrh; i >= nrl; i--) {
        //if ((m[i]+ncl) != NULL) {
            free((char *) (m[i]+ncl));
        //}
    }
    //if ((m+nrl) != NULL) {
        free((char *) (m+nrl));
    //}
    nv -= ncol * nrow;
}


/* Free double matrix allocated by dmatrix() */
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch) {
    int i;
    size_t nrow = nrh-nrl+1;
    size_t ncol = nch-ncl+1;

    for (i = nrh; i >= nrl; i--) {
        //if ((m[i]+ncl) != NULL) {
            free((char *) (m[i]+ncl));
        //}
    }
    //if ((m+nrl) != NULL) {
        free((char *) (m+nrl));
    //}
    nv -= ncol * nrow;
}


/* Free int array allocated by iarray3() */
void free_iarray3(int ***a, int n1, int n2, int n3) {
    free((char*) (a[0][0]));
    free((char*) (a[0]));
    free((char*) (a));
}

/* Free double array allocated by old version of darray3(), with indexes [ini1..fi1][ini2..fi2][ini3..fi3] */
void free_darray3(double ***a, int ini1, int fi1, int ini2, int fi2, int ini3, int fi3) {
  int i;
  for (i= fi1; i>=ini1; i--) {
    free_dmatrix(a[i],ini2,fi2,ini3,fi3);
  }
  //if ((a+ini1) != NULL) { 
    free((char*) (a+ini1)); 
  //}
}

//OLD VERSION OF free_darray3.
///* Free double array allocated by old version of darray3(), with indexes [0..n1][0..n2][0..n3] */
//void free_darray3(double ***a, int n1, int n2, int n3) {
//    free((char*) (a[0][0]));
//    free((char*) (a[0]));
//    free((char*) (a));
//}


/************************************************************************
                          MATHEMATICAL FUNCTIONS
************************************************************************/



double lfact(int n)
{

    if (n > 254)
    {
        double x = n + 1.0;
        return (x - 0.5)*log(x) - x + 0.5*log(2.0*M_PI) + 1.0/(12.0*x);
    }
    else
    {
        double lf[] =
        {
            0.000000000000000,
            0.000000000000000,
            0.693147180559945,
            1.791759469228055,
            3.178053830347946,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.604602902745251,
            12.801827480081469,
            15.104412573075516,
            17.502307845873887,
            19.987214495661885,
            22.552163853123421,
            25.191221182738683,
            27.899271383840894,
            30.671860106080675,
            33.505073450136891,
            36.395445208033053,
            39.339884187199495,
            42.335616460753485,
            45.380138898476908,
            48.471181351835227,
            51.606675567764377,
            54.784729398112319,
            58.003605222980518,
            61.261701761002001,
            64.557538627006323,
            67.889743137181526,
            71.257038967168000,
            74.658236348830158,
            78.092223553315307,
            81.557959456115029,
            85.054467017581516,
            88.580827542197682,
            92.136175603687079,
            95.719694542143202,
            99.330612454787428,
            102.968198614513810,
            106.631760260643450,
            110.320639714757390,
            114.034211781461690,
            117.771881399745060,
            121.533081515438640,
            125.317271149356880,
            129.123933639127240,
            132.952575035616290,
            136.802722637326350,
            140.673923648234250,
            144.565743946344900,
            148.477766951773020,
            152.409592584497350,
            156.360836303078800,
            160.331128216630930,
            164.320112263195170,
            168.327445448427650,
            172.352797139162820,
            176.395848406997370,
            180.456291417543780,
            184.533828861449510,
            188.628173423671600,
            192.739047287844900,
            196.866181672889980,
            201.009316399281570,
            205.168199482641200,
            209.342586752536820,
            213.532241494563270,
            217.736934113954250,
            221.956441819130360,
            226.190548323727570,
            230.439043565776930,
            234.701723442818260,
            238.978389561834350,
            243.268849002982730,
            247.572914096186910,
            251.890402209723190,
            256.221135550009480,
            260.564940971863220,
            264.921649798552780,
            269.291097651019810,
            273.673124285693690,
            278.067573440366120,
            282.474292687630400,
            286.893133295426990,
            291.323950094270290,
            295.766601350760600,
            300.220948647014100,
            304.686856765668720,
            309.164193580146900,
            313.652829949878990,
            318.152639620209300,
            322.663499126726210,
            327.185287703775200,
            331.717887196928470,
            336.261181979198450,
            340.815058870798960,
            345.379407062266860,
            349.954118040770250,
            354.539085519440790,
            359.134205369575340,
            363.739375555563470,
            368.354496072404690,
            372.979468885689020,
            377.614197873918670,
            382.258588773060010,
            386.912549123217560,
            391.575988217329610,
            396.248817051791490,
            400.930948278915760,
            405.622296161144900,
            410.322776526937280,
            415.032306728249580,
            419.750805599544780,
            424.478193418257090,
            429.214391866651570,
            433.959323995014870,
            438.712914186121170,
            443.475088120918940,
            448.245772745384610,
            453.024896238496130,
            457.812387981278110,
            462.608178526874890,
            467.412199571608080,
            472.224383926980520,
            477.044665492585580,
            481.872979229887900,
            486.709261136839360,
            491.553448223298010,
            496.405478487217580,
            501.265290891579240,
            506.132825342034830,
            511.008022665236070,
            515.890824587822520,
            520.781173716044240,
            525.679013515995050,
            530.584288294433580,
            535.496943180169520,
            540.416924105997740,
            545.344177791154950,
            550.278651724285620,
            555.220294146894960,
            560.169054037273100,
            565.124881094874350,
            570.087725725134190,
            575.057539024710200,
            580.034272767130800,
            585.017879388839220,
            590.008311975617860,
            595.005524249382010,
            600.009470555327430,
            605.020105849423770,
            610.037385686238740,
            615.061266207084940,
            620.091704128477430,
            625.128656730891070,
            630.172081847810200,
            635.221937855059760,
            640.278183660408100,
            645.340778693435030,
            650.409682895655240,
            655.484856710889060,
            660.566261075873510,
            665.653857411105950,
            670.747607611912710,
            675.847474039736880,
            680.953419513637530,
            686.065407301994010,
            691.183401114410800,
            696.307365093814040,
            701.437263808737160,
            706.573062245787470,
            711.714725802289990,
            716.862220279103440,
            722.015511873601330,
            727.174567172815840,
            732.339353146739310,
            737.509837141777440,
            742.685986874351220,
            747.867770424643370,
            753.055156230484160,
            758.248113081374300,
            763.446610112640200,
            768.650616799717000,
            773.860102952558460,
            779.075038710167410,
            784.295394535245690,
            789.521141208958970,
            794.752249825813460,
            799.988691788643450,
            805.230438803703120,
            810.477462875863580,
            815.729736303910160,
            820.987231675937890,
            826.249921864842800,
            831.517780023906310,
            836.790779582469900,
            842.068894241700490,
            847.352097970438420,
            852.640365001133090,
            857.933669825857460,
            863.231987192405430,
            868.535292100464630,
            873.843559797865740,
            879.156765776907600,
            884.474885770751830,
            889.797895749890240,
            895.125771918679900,
            900.458490711945270,
            905.796028791646340,
            911.138363043611210,
            916.485470574328820,
            921.837328707804890,
            927.193914982476710,
            932.555207148186240,
            937.921183163208070,
            943.291821191335660,
            948.667099599019820,
            954.046996952560450,
            959.431492015349480,
            964.820563745165940,
            970.214191291518320,
            975.612353993036210,
            981.015031374908400,
            986.422203146368590,
            991.833849198223450,
            997.249949600427840,
            1002.670484599700300,
            1008.095434617181700,
            1013.524780246136200,
            1018.958502249690200,
            1024.396581558613400,
            1029.838999269135500,
            1035.285736640801600,
            1040.736775094367400,
            1046.192096209724900,
            1051.651681723869200,
            1057.115513528895000,
            1062.583573670030100,
            1068.055844343701400,
            1073.532307895632800,
            1079.012946818975000,
            1084.497743752465600,
            1089.986681478622400,
            1095.479742921962700,
            1100.976911147256000,
            1106.478169357800900,
            1111.983500893733000,
            1117.492889230361000,
            1123.006317976526100,
            1128.523770872990800,
            1134.045231790853000,
            1139.570684729984800,
            1145.100113817496100,
            1150.633503306223700,
            1156.170837573242400,
        };
        return lf[n];
    }
}


/*
 * log-double factorial
 * Returns x*(x-2)*...*k, where k is the first term <2.
 * Corresponds to double factorial for integer x
 */
double ldoublefact(double x)
{
    int i;
    double ans = 0.0;

    /* :TBD: converts from double to integer without noting rounding intent */
    for (i = ((int) x); i >= 2; i -= 2) {
      ans += log((double) i);
    }
    return(ans);
}


/*
-----------------------------------------------------------------------
            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
-----------------------------------------------------------------------
     WRITTEN BY ALFRED H. MORRIS
          NAVAL SURFACE WARFARE CENTER
          DAHLGREN, VIRGINIA
--------------------------
     D = 0.5*(LN(2*PI) - 1)
--------------------------
*/
double gamln(double *a)
{
static double c0 = .833333333333333e-01;
static double c1 = -.277777777760991e-02;
static double c2 = .793650666825390e-03;
static double c3 = -.595202931351870e-03;
static double c4 = .837308034031215e-03;
static double c5 = -.165322962780713e-02;
static double d = .418938533204673e0;
static double gamln,t,w;
static int i,n;
static double T1;
/*
     ..
     .. Executable Statements ..
*/
    if(*a > 0.8e0) goto S10;
    gamln = gamln1(a)-log(*a);
    return gamln;
S10:
    if(*a > 2.25e0) goto S20;
    t = *a-0.5e0-0.5e0;
    gamln = gamln1(&t);
    return gamln;
S20:
    if(*a >= 10.0e0) goto S40;
    n = (int)(*a - 1.25e0);
    //n = (long)(*a - 1.25e0);
    t = *a;
    w = 1.0e0;
    for(i=1; i<=n; i++) {
        t -= 1.0e0;
        w = t*w;
    }
    T1 = t-1.0e0;
    gamln = gamln1(&T1)+log(w);
    return gamln;
S40:
    t = pow(1.0e0/ *a,2.0);
    w = (((((c5*t+c4)*t+c3)*t+c2)*t+c1)*t+c0)/ *a;
    gamln = d+w+(*a-0.5e0)*(log(*a)-1.0e0);
    return gamln;
}


/*
-----------------------------------------------------------------------
     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
-----------------------------------------------------------------------
*/
double gamln1(double *a)
{
static double p0 = .577215664901533e+00;
static double p1 = .844203922187225e+00;
static double p2 = -.168860593646662e+00;
static double p3 = -.780427615533591e+00;
static double p4 = -.402055799310489e+00;
static double p5 = -.673562214325671e-01;
static double p6 = -.271935708322958e-02;
static double q1 = .288743195473681e+01;
static double q2 = .312755088914843e+01;
static double q3 = .156875193295039e+01;
static double q4 = .361951990101499e+00;
static double q5 = .325038868253937e-01;
static double q6 = .667465618796164e-03;
static double r0 = .422784335098467e+00;
static double r1 = .848044614534529e+00;
static double r2 = .565221050691933e+00;
static double r3 = .156513060486551e+00;
static double r4 = .170502484022650e-01;
static double r5 = .497958207639485e-03;
static double s1 = .124313399877507e+01;
static double s2 = .548042109832463e+00;
static double s3 = .101552187439830e+00;
static double s4 = .713309612391000e-02;
static double s5 = .116165475989616e-03;
static double gamln1,w,x;
/*
     ..
     .. Executable Statements ..
*/
    if(*a >= 0.6e0) goto S10;
    w = ((((((p6**a+p5)**a+p4)**a+p3)**a+p2)**a+p1)**a+p0)/((((((q6**a+q5)**a+
      q4)**a+q3)**a+q2)**a+q1)**a+1.0e0);
    gamln1 = -(*a*w);
    return gamln1;
S10:
    x = *a-0.5e0-0.5e0;
    w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0)/(((((s5*x+s4)*x+s3)*x+s2)*x+s1)*x
      +1.0e0);
    gamln1 = x*w;
    return gamln1;
}


/*
 * Returns the derivative of the natural log of the gamma function,
 *    e.g. gamma'(x)/gamma(x)
 * - x must be positive and finite
 */
double digamma(double x)
{
    double stirling[] = {
        -8.333333333333333e-02, 8.333333333333333e-03,
        -3.968253968253968e-03, 4.166666666666667e-03,
        -7.575757575757576e-03, 2.109279609279609e-02,
        -8.333333333333334e-02, 4.432598039215686e-01,
        -3.053954330270120e+00, 2.645621212121212e+01,
        -2.814601449275362e+02, 3.607510546398047e+03
    };

    const double lower = 1.0e-8;
    const double upper = 19.5;
    const double euler_one = 0.422784335098467139393488;
    double x_inv;
    double x_pow;
    double ans;
    int i;

    if (x <= 0) {
        errorC("digamma", "argument must be positive", 1);
        /*NOTREACHED*/
    }

    if (x < lower) {
        ans = -1.0 / x - 1.0/(1.0+x) + euler_one;
        return(ans);
    }

    ans = 0.0;
    while (x < upper) {
        ans = ans - 1.0/x;
        x = x + 1.0;
    }

    x_inv = 1.0 / x;
    ans = ans + log(x) - 0.5*x_inv;

    x_inv = x_inv * x_inv;
    x_pow = x_inv;

    for (i = 0; i < 12; i++) {
        ans = ans + stirling[i] * x_pow;
        x_pow = x_pow * x_inv;
    }

    return(ans);
}


/* Bernoulli numbers of even order from 2 to 60 */
static double bernou[30] = {
    1.0/6.0,
   -1.0/30.0,
    1.0/42.0,
   -1.0/30.0,
    5.0/66.0,
   -691.0/2730.0,
    7.0/6.0,
   -3617.0/510.0,
    43867.0/798.0,
   -174611.0/330.0,
   854513.0/138.0,
   -236364091.0/2730.0,
   8553103.0/6.0,
   -23749461029.0/870.0,
   8615841276005.0/14322.0,
   -7709321041217.0/510.0,
   2577687858367.0/6.0,
   -1.371165521e13,
    4.883323190e14,
   -1.929657934e16,
    8.416930476e17,
   -4.033807185e19,
    2.115074864e21,
   -1.208662652e23,
    7.500866746e24,
   -5.038778101e26,
    3.652877648e28,
   -2.849876930e30,
    2.386542750e32,
   -2.139994926e34
};


double trigamma(double x)
{
    if (x > 1.0e-5) {
        /* fast approx to trigamma */
        return(1 / (x * x) +
               1 / ((x+1) * (x+1)) +
               1 / ((x+2) * (x+2)) +
               1 / (x+3) +
               0.5 / ((x+3) * (x+3)) +
               1 / (6.0 * pow(x+3.0, 3.0)));
    } else {
        /* slower computation */
        return(polygamma(x, 1, 0.0001, 100, 5, 1));
    }
}


/*
 * polygamma function of a real positive x of order n (n==1 is trigamma etc.).
 * setting low=0.0001, high=100 and terms=5 usually gives good results
 * nfact is n! e.g. 1 for trigamma etc.
 * no checks are made here on the suitability of arguments
 */
double polygamma(double x,
                 long n,
                 double low,
                 double high,
                 long terms,
                 double nfact)
{
    double asign;
    double ans = 0.0;
    double nd = (double) n;
    double nexp;
    double ser = 0.0;
    double t0;
    double x2_inv;
    long i;

    asign = (n % 2) ? 1.0 : -1.0;
    if (x < low) {
      return(asign * nfact / nd * pow(x, - nd) * (1.0 + nd * 0.5 / x));
    }
    nexp = - nd - 1.0;
    while (x < high) {
        ans = ans + asign * nfact * pow(x, nexp);
        x = x + 1.0;
    }
    t0 = nfact / nd * pow(x, - nd);
    ser = t0 * ( 1.0 + nd * 0.5 / x);
    x2_inv = pow(x, -2.0);
    for (i = 0; i < terms; i++) {
        if (n == 1) {
            t0 = t0 * x2_inv;
        } else {
            t0 = (2.0 * i + nd + 3.0) / (2.0 * i + 4.0) * (2.0 * i + nd + 2.0) / (2.0 * i + 3.0) * t0 * x2_inv;
        }
        ser = ser + bernou[i] * t0;
    }
    ans = ans + asign * ser;
    return(ans);
}


/* log of Beta function */
double lnbeta(double a,
              double b)
{
    double c = a + b;
    return(gamln(&a) + gamln(&b) - gamln(&c));
}


/*
 * Used by pbetaC: Evaluates continued fraction for incomplete beta function
 * by modified Lentz's method (x5.2).
 */
double betacf(double a,
              double b,
              double x)
{
    double aa;
    double c;
    double d;
    double del;
    double h;
    double qab;
    double qam;
    double qap;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;
    const int MAXIT = 100;
    int m;
    int m2;

    /* These q's will be used in factors that occur in coe cients (6.4.6). */
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;

    /* First step of Lentz's method */
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) {
        d = FPMIN;
    }
    d = 1.0 / d;
    h = d;
    for (m = 1; m <= MAXIT; m++) {
        m2 = 2 * m;
        aa = m * (b-m) * x / ((qam+m2) * (a+m2));
        /* One step (the even one) of the recurrence */
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
        /* Next step of the recurrence (the odd one) */
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) {
            c = FPMIN;
        }
        d = 1.0 / d;
        del = d * c;
        h *= del;
        /* Are we done? */
        if (fabs(del-1.0) < EPS) {
            break;
        }
    }
    if (m > MAXIT) {
        nrerror("betacf", "", "a or b too big, or MAXIT too small");
        /*NOTREACHED*/
    }
    return(h);
}


double lnchoose(double n,
                double k)
{
  return -lnbeta(1.0+n-k,1.0+k) - log(1.0+n);
}


double choose(double n,
              double k)
{
    return exp(lnchoose(n, k));
}


double logit(double x)
{
    return(log(x / (1-x)));
}


double ilogit(double x)
{
    return(1.0 / (1.0 + exp(-x)));
}


/* Returns 1.0 if x>=0, -1.0 if x<0 */
double dsign(double x)
{
    return (x >= 0) ? 1.0 : -1.0;
}


/* Returns 1.0 if x>0, 0 if x==0, -1.0 if x<0 */
double isign(int x)
{
    if (x == 0) {
        return 0.0;
    }
    return (x > 0) ? 1.0 : -1.0;
}

/**************************************************************/
/* MATRIX MANIPULATION                                        */
/**************************************************************/

//Copy A[model,model] into Aout
void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model) {
  int i, j;
  arma::SpMat<short>::iterator it, it2;

  for (it= model->begin(), i=0; it != model->end(); ++it, i++) {
    for (it2= model->begin(), j=0; it2 != model->end(); ++it2, j++) {
        Aout->at(i,j)= A->at(it.row(), it2.row());
    }
  }
}


/************************************************************************
                            VECTOR ALGEBRA
************************************************************************/

//Squared Mahalanobis distances (x[i,] - x[j,])' [ cholA %*% t(cholA) ] (x[i,] - x[j,])
// Input
// - x: an n x p matrix, formatted as a 0-indexed vector containing the rows in x, that is x[1,],x[2,],...,x[n,]
// - n: number of rows in x
// - p: number of columns in x
// - cholA: Cholesky decomposition of the scale matrix (p x p matrix, upper diagonal elements ignored)
//
// Output
// - ans: ans[1],...,ans[n*(n-1)/2] contain squared Mahalanobis distances between (1,2), (1,3), etc.
void mahaldist(double *x, int n, int p, double **cholA, bool squared, double *ans) {
  int i, j, l, idx;
  double **z;

  z= dmatrix(1,n,1,p);

  //Compute z= t(cholA) %*% x
  for (i=1; i<=n; i++) {
    idx= (i-1)*p;
    for (j=1; j<=p; j++) {
      z[i][j]= 0;
      for (l=j; l<=p; l++) { z[i][j] += cholA[l][j] * x[idx +l-1]; }
    }
  }

  //Euclidean distances on z
  idx= 1;
  for (i=1; i<=n; i++) {
    for (j=i+1; j<=n; j++) {
      ans[idx]= 0;
      for (l=1; l<=p; l++) { ans[idx] += pow(z[i][l] - z[j][l],2.0); }
      if (!squared) ans[idx]= sqrt(ans[idx]);
      idx++;
    }
  }

  free_dmatrix(z,1,n,1,p);
}



void grid(double x0,
          double xn,
          int n,
          double *x)
{
    int i;
    double dx;
    double xi;

    //assert(x != NULL);

    dx = (xn - x0) / (n - 1.0);
    xi = x0;
    for (i = 0; i < n; i++) {
        x[i] = xi;
        xi += dx;
    }
}


/*
 * Multiply matrix A[1..p][1..q] by scalar r.
 * Store results in matrix B.
 */
void rA(double r,
        double **A,
        double **B,
        int rowini,
        int rowfi,
        int colini,
        int colfi)
{
    int i;
    int j;

    //assert(A != NULL);
    //assert(B != NULL);

    for (i = rowini; i <= rowfi; i++) {
        for (j = colini; j <= colfi; j++) {
            B[i][j] = r * A[i][j];
        }
    }
}


/*
 * Sum matrix A[rowini..rowfi][colini..colfi] + B[rowini..rowfi][colini..colfi].
 * Store results in C.
 */
void A_plus_B(double **A,
              double **B,
              double **C,
              int rowini,
              int rowfi,
              int colini,
              int colfi)
{
    int i;
    int j;

    //assert(A != NULL);
    //assert(B != NULL);
    //assert(C != NULL);

    for (i = rowini; i <= rowfi; i++) {
        for (j = colini; j <= colfi; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}


/*
 * Add matrix A[rowini..rowfi][colini..colfi] (multiplied by scalar r) and
 * matrix B[rowini..rowfi][colini..colfi] (multiplied by scalar s).
 * Store results in C.
 */
void rA_plus_sB(double r,
                double **A,
                double s,
                double **B,
                double **C,
                int rowini,
                int rowfi,
                int colini,
                int colfi)
{
    int i, j;

    //assert(A != NULL);
    //assert(B != NULL);
    //assert(C != NULL);

    for (i = rowini; i <= rowfi; i++) {
        for (j = colini; j <= colfi; j++) {
            C[i][j] = (r * A[i][j]) + (s * B[i][j]);
        }
    }
}


/*
 * Sum the product of scalar r, matrix A[rowini..rowfi][colini..colfi],
 * and vector x[colini..colfi] with the product of scalar s,
 * matrix B[rowini..rowfi][colini..colfi], and vector y[colini..colfi].
 * Store results in vector z[rowini..rowfi].
 */
void rAx_plus_sBy(double r,
                  double **A,
                  const double *x,
                  double s,
                  double **B,
                  const double *y,
                  double *z,
                  int rowini,
                  int rowfi,
                  int colini,
                  int colfi)
{
    int i, j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(B != NULL);
    //assert(y != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += (r * A[i][j] * x[j]) + (s * B[i][j] * y[j]);
        }
    }
}


/*
 * Multiply matrix A[ini..fi][ini..fi] by vector x[ini..fi]
 * and add vector y[ini..fi].
 * Store result in vector z.
 */
void Ax_plus_y(double **A,
               const double *x,
               const double *y,
               double *z,
               int ini,
               int fi)
{
    int i, j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(y != NULL);
    //assert(z != NULL);

    for (i = ini; i <= fi; i++) {
        z[i] = y[i];
        for (j = ini; j <= fi; j++) {
            z[i] += A[i][j] * x[j];
        }
    }
}


/*
 * Multiply matrix A[ini..fi][ini..fi] by vector x[ini..fi].
 * Store result in vector z.
 */
void xA(const double *x,
        double **A,
        double *z,
        int ini,
        int fi)
{
    int i, j;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(z != NULL);

    for (i = ini; i <= fi; i++) {
        z[i] = 0.0;
        for (j = ini; j <= fi; j++) {
            z[i] += A[j][i] * x[j];
        }
    }
}


/*
 * Multiply matrix A[rowini..rowfi][colini..colfi] by vector x[colini..colfi].
 * Store result in vector z.
 */
void Ax(double **A,
        const double *x,
        double *z,
        int rowini,
        int rowfi,
        int colini,
        int colfi)
{
    int i;
    int j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[i][j] * x[j];
        }
    }
}


/*
 * Multiply (implicit matrix) vector A by vector x[colini..colfi].
 * Store result in vector z.
 */
void Avecx(const double *A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi) {
    int i;
    int j;
    int nrow = rowfi - rowini + 1;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[i + j*nrow] * x[j];
        }
    }
}


/*
 * Multiply z= A[,sel] %*% x, i.e. choose elements indicated by sel[0],...,sel[*nsel-1]
 */
void Aselvecx(const double *A, const double *x, double *z, int rowini, int rowfi, int *sel, int *nsel) {
  int i;
  int j;
  int nrow = rowfi - rowini + 1;

  for (i = rowini; i <= rowfi; i++) {
    z[i] = 0.0;
    for (j = 0; j < (*nsel); j++) {
      z[i] += A[i + sel[j]*nrow] * x[j];
    }
  }
}


/*
 * Multiply z= t(A[,sel]) %*% x, i.e. choose columns in A indicated by sel[0],...,sel[*nsel-1]
 */
void Atselvecx(const double *A, const double *x, double *z, int rowini, int rowfi, int *sel, int *nsel) {
  int i;
  int j;
  int nrow = rowfi - rowini + 1, idx;

  for (j=0; j<(*nsel); j++) {
    z[j]= 0.0;
    idx= sel[j]*nrow;
    for (i=rowini; i<=rowfi; i++) {
      z[j] += A[i + idx] * x[i];
    }
  }

}


/*
 * Multiply transposed (implicit matrix) vector A by vector x[colini..colfi].
 * Store result in vector z.
 */
void Atvecx(const double *A,
            const double *x,
            double *z,
            int rowini,
            int rowfi,
            int colini,
            int colfi)
{
    int i;
    int j;
    int ncol = colfi - colini + 1;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[j + i*ncol] * x[j];
        }
    }
}


/*
 * Returns sum of multiplying matrix A[ini..fi][ini..fi]
 * by transposed vector x[ini..fi] by vector y[ini..fi].
 */
double xtAy(const double *x,
            double **A,
            const double *y,
            int ini,
            int fi)
{
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(y != NULL);

    for (i = ini; i <= fi; i++) {
        for (j = ini; j <= fi; j++) {
            z += A[i][j] * x[j] * y[i];
        }
    }
    return(z);
}


/*
 * Returns sum of multiplying symmetric matrix A[ini..fi][ini..fi]
 * by transposed vector x[ini..fi] by vector x[ini..fi].
 *
 * Note: Faster than xtAy() for symmetric A (saves 25%-50% operations).
 */
double quadratic_xtAx(const double *x,
                      double **A,
                      int ini,
                      int fi)
{
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);

    for (i = ini; i <= fi; i++) {
        z += A[i][i] * x[i] * x[i];
        for (j = i+1; j <= fi; j++) {
            z += 2 * A[i][j] * x[i] *x[j];
        }
    }
    return(z);
}


/*
 * Returns sum of multiplying symmetric (implicit matrix) vector A[sel][sel]
 * by transposed vector x[sel] by vector x[sel] for quadratic forms.
 *     ncolA: number of columns in A.
 *     nsel : length of vector sel.
 *     sel  : vector with indexes for positions in x and (rows,cols) in A
 *            to be used in the operation.
 *
 * Note: Faster than xtAy() for symmetric A (saves 25%-50% operations).
 */
double quadratic_xseltAselxsel(const double *x,
                               crossprodmat *A,
                               const int *ncolA,
                               const int *nsel,
                               const int *sel)
{
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(ncolA != NULL);
    //assert(nsel != NULL);
    //assert(sel != NULL);

    for (i = 0; i <= (*nsel)-1; i++) {
        int i_sel;

        i_sel = sel[i];
        z += (A->at(i_sel * (*ncolA) + i_sel)) * x[i_sel] * x[i_sel];
        for (j = i+1; j <= (*nsel)-1; j++) {
            int j_sel;

            j_sel = sel[j];
            z += 2 * (A->at(i_sel * (*ncolA) + j_sel)) * x[i_sel] * x[j_sel];
        }
    }
    return(z);
}


/*
 * Returns sum of multiplying symmetric (implicit matrix) vector A[sel][sel]
 * by transposed vector x[0..nsel] by vector x[0..nsel] for quadratic forms.
 *     ncolA: number of columns in A.
 *     nsel : length of vector sel.
 *     sel  : vector with indexes for (rows,cols) in A to be used in the operation.
 *
 * Same as above but subset is only for A.
 * Note: Faster than xtAy() for symmetric A (saves 25%-50% operations).
 */
double quadratic_xtAselx(const double *x,
                         crossprodmat *A,
                         const int *ncolA,
                         const int *nsel,
                         const int *sel)
{
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(ncolA != NULL);
    //assert(nsel != NULL);
    //assert(sel != NULL);

    for (i = 0; i <= (*nsel)-1; i++) {
        int i_sel;

        i_sel = sel[i];
        z += (A->at(i_sel * (*ncolA) + i_sel)) * x[i] * x[i];
        for (j = i+1; j <= (*nsel)-1; j++) {
	  z += 2 * (A->at(i_sel * (*ncolA) + sel[j])) * x[i] * x[j];
        }
    }
    return(z);
}


/*
 * Returns sum of multiplying symmetric matrix A[sel][sel]
 * by transposed vector x[sel] by vector x[sel] for quadratic forms.
 *     ini  : element in A are indexed from ini to *nsel-1
 *     nsel : length of vector sel.
 *     sel  : vector with indexes for positions in x and (rows,cols) in A
 *            to be used in the operation.
 *
 * Same as above but subset is only for x.
 * Note: Faster than xtAy() for symmetric A (saves 25%-50% operations).
 */
double quadratic_xseltAxsel(const double *x,
                            double **A,
                            int ini,
                            const int *nsel,
                            const int *sel)
{
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(nsel != NULL);
    //assert(sel != NULL);

    for (i = 0; i <= (*nsel)-1; i++) {
        int rowA;
        int i_sel;

        rowA = ini + i;
        i_sel = sel[i];
        z += A[rowA][rowA] * x[i_sel] * x[i_sel];
        for (j = i+1; j <= (*nsel)-1; j++) {
            z += 2 * A[rowA][ini + j] * x[i_sel] * x[sel[j]];
        }
    }
    return(z);
}



/*
 * Multiply transpose of matrix A[rowini..rowfi][colini..colfi] by
 * vector x[rowini..rowfi].
 * Store result in vector z.
 */
void Atx(double **A,
         const double *x,
         double *z,
         int rowini,
         int rowfi,
         int colini,
         int colfi)
{
    int i;
    int j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = colini; i <= colfi; i++) {
        z[i] = 0.0;
        for (j = rowini; j <= rowfi; j++) {
            z[i] += A[j][i] * x[j];
        }
    }
}


/*
 * Multiply matrix A[rowiniA..rowfiA][coliniA..colfiA] by
 * matrix B[rowiniB..rowfiB][coliniB..colfiB].
 * Store result in matrix C.
 */
void AB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C) {
  int i,j,k;
    if ((colfiA-coliniA) != (rowfiB-rowiniB)) { errorC("AB", "dimensions don't match", 1); }
    for (i = rowiniA; i <= rowfiA; i++) {
        for (j = coliniB; j <= colfiB; j++) {
            C[i][j] = 0.0;
            for (k = coliniA; k <= colfiA; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


/*
 * Multiply transpose of matrix A[rowiniA..rowfiA][coliniA..colfiA] by
 * matrix B[rowiniB..rowfiB][coliniB..colfiB].
 * Store result in matrix C.
 */
void AtB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C) {
  int i,j,k;

    if ((rowfiA-rowiniA) != (rowfiB-rowiniB)) {
        errorC("AtB", "dimensions don't match", 1);
        /*NOTREACHED*/
    }
    for (i = coliniA; i <= colfiA; i++) {
        for (j = coliniB; j <= colfiB; j++) {
            C[i][j] = 0.0;
            for (k = rowiniA; k <= rowfiA; k++) {
                C[i][j] += A[k][i] * B[k][j];
            }
        }
    }
}


/*
 * Multiply matrix A[rowiniA..rowfiA][coliniA..colfiA] by
 * transpose of matrix B[rowiniB..rowfiB][coliniB..colfiB].
 * Store result in matrix C.
 */
void ABt(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C) {
  int i,j,k;

    if ((colfiA-coliniA) != (colfiB-coliniB)) { errorC("AtB", "dimensions don't match", 1); }
    for (i = rowiniA; i <= rowfiA; i++) {
        for (j = rowiniB; j <= rowfiB; j++) {
            C[i][j] = 0.0;
            for (k = coliniA; k <= colfiA; k++) {
                C[i][j] += A[i][k] * B[j][k];
            }
        }
    }
}


//C= t(A) %*% B, where C[1..ncolA][1..nrowB]
void AvectBvec(double *A, int nrowA, int ncolA, double *B, int nrowB, int ncolB, double **C) {
    int i;
    int j;
    int k;

    //assert(A != NULL);
    //assert(B != NULL);
    //assert(C != NULL);

    if (nrowA != nrowB) {
        errorC("AvectBvec", "dimensions don't match", 1);
        /*NOTREACHED*/
    }
    for (i = 1; i <= ncolA; i++) {
      int offsetA= (i-1)*nrowA;
        for (j = 1; j <= ncolB; j++) {
            C[i][j] = 0.0;
	    int offsetB= (j-1)*nrowB;
            for (k = 1; k <= nrowA; k++) {
	      C[i][j] += A[k-1 + offsetA] * B[k-1 + offsetB];
            }
        }
    }
}


void a_plus_b(const double *a,
              const double *b,
              double *c,
              int ini,
              int fi)
{
    int i;

    //assert(a != NULL);
    //assert(b != NULL);
    //assert(c != NULL);

    for (i = ini; i <= fi; i++) {
        c[i] = a[i] + b[i];
    }
}


void a_prod_b(const double *a,
              const double *b,
              double *c,
              int ini,
              int fi)
{
    int i;

    //assert(a != NULL);
    //assert(b != NULL);
    //assert(c != NULL);

    for (i = ini; i <= fi; i++) {
        c[i] = a[i] * b[i];
    }
}


void a_prod_b_sel(const double *a,
                  const double *b,
                  double *c,
                  const int *lengtha,        /* :TBD: unused? */
                  const int *nsel,
                  const int *sel)
{
    int i;

    //assert(a != NULL);
    //assert(b != NULL);
    //assert(c != NULL);

    for (i = 0; i <= (*nsel-1); i++) {
        int idx;

        idx = sel[i];
        c[idx] = a[idx] * b[idx];
    }
}


void a_zero(double *a,
            int p)
{
    int i;

    //assert(a != NULL);

    for (i = 0; i < p; i++) {
        a[i] = 0.0;
    }
}


void R_zero(double **A,
            int p,
            int q)
{
    int i;
    int j;

    //assert(A != NULL);

    for (i = 0; i < p; i++) {
        for (j = 0; j < q; j++) {
            A[i][j] = 0.0;
        }
    }
}


/* Diagonal matrix */
void ddiag(double **A,
           int ini,
           int fi)
{
    int i;
    int j;

    //assert(A != NULL);

    for (i = ini; i <= fi; i++) {
        for (j = ini; j <= fi; j++) {
            A[i][j] = (i == j) ? 1 : 0;
        }
    }
}


/* Returns the absolute value of parameter x */
/* NOTE: redefines C89 standard library function abs() for no good reason. */
int iabs(int x)
{
    return (x >= 0) ? x : -x;
}


int imax_xy(int x,
            int y)
{
    return (x > y) ? x : y;
}


int imin_xy(int x,
            int y)
{
    return (x < y) ? x : y;
}


double max_xy(double x,
              double y)
{
    return (x > y) ? x : y;
}


double min_xy(double x,
              double y)
{
    return (x < y) ? x : y;
}


/*
 * Minimum of vector x[ini..fi] is returned in xmin.
 * Position at which min occurs is returned in minpos
 */
void minvec(const double *x,
            int ini,
            int fi,
            double *xmin,
            int *minpos)
{
    int i;

    //assert(x != NULL);
    //assert(xmin != NULL);
    //assert(minpos != NULL);

    *xmin = x[ini];
    *minpos = ini;

    for (i = ini+1; i <= fi; i++) {
        if (x[i] < (*xmin)) {
            *xmin = x[i];
            *minpos= i;
        }
    }
}


/*
 * Maximum of vector x[ini..fi] is returned in xmax.
 * Position at which max occurs is returned in maxpos
 */
void maxvec(const double *x,
            int ini,
            int fi,
            double *xmax,
            int *maxpos)
{
    int i;

    *xmax = x[ini];
    *maxpos = ini;
    for (i = ini+1; i <= fi; i++) {
        if (x[i] > (*xmax)) {
            *xmax = x[i];
            *maxpos = i;
        }
    }
}


/* Cumulative sum of elements in x */
void cumsum(double *x, double *cumsumx, int *n) {
  int i;
  cumsumx[0]= x[0];
  for (i=1; i< *n; i++) cumsumx[i]= cumsumx[i-1] + x[i];
} 


/* Make matrix a positive definite by replacing it by a - (lmin + offset) I, where lmin is the smallest eigenvalue of a and I the identity matrix  */
void make_posdef(double **a, int n, double offset) {
  int i;
  double lmin=0, *vals;
  vals= dvector(1,n);
  eigenvals(a,n,vals);
  for (i=1; i<= n; i++) if (vals[i]<lmin) lmin= vals[i];
  lmin = -lmin + offset;
  for (i=1; i<= n; i++) a[i][i] += lmin;
  free_dvector(vals,1,n);
}



void choldc(double **a, int n, double **aout, bool *posdef) {
/* Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
decomposition, A = L * L' . On input, only the upper triangle of a need be given;
 The Cholesky factor L is returned in the lower triangle of aout (upper-diag elem are set to 0) */
  int i,j,k;
  double sum, *p, max_a;

  *posdef= true;
  for (i=1;i<=n;i++) { for (j=i;j<=n;j++) { aout[i][j]= a[i][j]; } }  //copy a into aout
  p= dvector(1,n);
  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      for (sum=aout[i][j],k=i-1;k>=1;k--) sum -= aout[i][k]*aout[j][k];
      if (i == j) {
	if (sum <= 0.0) *posdef= false;
	aout[i][i]=sqrt(sum);
      } else {
	max_a=max_xy(fabs(aout[i][i]), 1e-10);
	aout[j][i]=sum/max_a;
      }
    }
  }
  free_dvector(p,1,n);
  for (i=1;i<=n;i++) { for (j=i+1;j<=n;j++) { aout[i][j]= 0; } }  //set upper-diagonal elem to 0
}


/*Given a positive-definite symmetric matrix a[1..n][1..n], this routine computes the inverse
 of its Cholesky matrix. That is, if A=L * L' it returns the inverse of L
 (note that inv(A)= inv(L)' * inv(L)) */
void choldc_inv(double **a, int n, double **aout, bool *posdef) {
  choldc(a, n, aout, posdef);
  if (*posdef) {
    choldc_inv_internal(aout, n);
  }
}


/*Given the Cholesky decomposition of a matrix S, which we denote cholS, returns the inverse of cholS */
void cholS_inv(double **cholS, int n, double **cholSinv) {
  int i,j;
  for (i=1;i<=n;i++) for (j=1;j<=i;j++) cholSinv[i][j]= cholS[i][j];
  choldc_inv_internal(cholSinv,n);
}


/*Computes inverse of Cholesky matrix cholS and stores the result in cholS*/
void choldc_inv_internal(double **cholS, int n) {
  int i,j,k;
  double sum, max_a;
  for (i=1;i<=n;i++) {
    max_a=max_xy(cholS[i][i], 1e-10);
    cholS[i][i]=1.0/max_a;
    for (j=i+1;j<=n;j++) {
      sum=0.0;
      for (k=i;k<j;k++) sum -= cholS[j][k]*cholS[k][i];
      max_a=max_xy(cholS[j][j], 1e-10);
      cholS[j][i]=sum/max_a;
    }
  }
}


/*
 * Find determinant of the matrix having chols as its Cholesky decomposition.
 *
 * Example:
 *   choldc(S, n, cholS, posdef);
 *   det = choldc_det(cholS, n);
 *
 * Another example:
 *   choldc_inv(S, n, cholSinv, posdef);
 *   det = 1.0 / choldc_det(cholSinv, n);
 */
double choldc_det(double **chols, int n) {
    int i;
    double value, det = 1.0;

    //assert(chols != NULL);
    for (i = 1; i <= n; i++) { value = chols[i][i]; det *= value * value; }
    return(det);
}

double logcholdc_det(double **chols, int n) {
    int i;
    double logdet = 0;

    //assert(chols != NULL);
    for (i = 1; i <= n; i++) { logdet += log(chols[i][i]); }
    return(2.0*logdet);
}


/*
  Inverse of a symmetric, positive definite matrix a[1..n][1..n] using Cholesky decomposition. 

  Input: either a, its Cholesky decomposition chola, or the inverse of its Cholesky decomposition cholainv. If chola is provided then a is ignored. If cholainv is provided, then chola is ignored.

  Output: aout contains the inverse of a, posdef returns if matrix was indeed positive definite

 */
void inv_posdef(double **a, int n, double **aout, bool *posdef, double **chola, double **cholainv) {
  int i, j;
  double **b;

  if (cholainv == NULL) {
    b = dmatrix(1, n, 1, n);
    if (chola == NULL) {
      choldc_inv(a, n, b, posdef); //inverse of chol(a)
    } else {
      cholS_inv(chola, n, b);  //inverse of chola
    }
  } else {
    b= cholainv;
  }
  
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      int k;
      double sum;

      sum = 0.0;
      for (k = 1; k <= n; k++) { sum += b[k][i] * b[k][j]; }
      aout[i][j] = sum;
    }
  }
  
  if (cholainv == NULL) free_dmatrix(b, 1, n, 1, n);

  for (i = 2; i <= n; i++) {
    for (j = 1; j < i; j++) {
      aout[i][j] = aout[j][i];
    }
  }
}


/*
 * Inverse of a symmetric, positive definite matrix a[1..n][1..n] using
 * Cholesky decomposition. Result is returned in aout.
 * Does the same as inv_posdef, except that here only upper triangular
 * elements are returned.
 */
void inv_posdef_upper(double **a,
                      int n,
                      double **aout,
                      bool *posdef)
{
    int i;
    int j;
    double **b;

    //assert(a != NULL);
    //assert(aout != NULL);

    b = dmatrix(1, n, 1, n);
    choldc_inv(a, n, b, posdef);
    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) {
                sum += b[k][i] * b[k][j];
            }
            aout[i][j] = sum;
        }
    }
    free_dmatrix(b, 1, n, 1, n);
}


/*
 * Inverse and determinant of a positive definite matrix a[1..n][1..n] using
 * Cholesky decomposition. Inverse is returned in aout, determinant in det_a.
 */
void invdet_posdef(double **a, int n, double **aout, double *det_a) {
    bool posdef;
    int i;
    int j;
    double **b;

    b = dmatrix(1, n, 1, n);
    choldc_inv(a, n, b, &posdef);
    *det_a = 1.0;
    for (i = 1; i <= n; i++) {
        double value;

        value = b[i][i];
        (*det_a) *= 1 / (value * value);
    }

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) {
                sum += b[k][i] * b[k][j];
            }
            aout[i][j] = sum;
        }
    }
    free_dmatrix(b, 1, n, 1, n);

    for (i = 2; i <= n; i++) {
        for (j = 1; j < i; j++) {
            aout[i][j] = aout[j][i];
        }
    }
}


/*
 * Inverse and determinant of a positive definite matrix a[1..n][1..n] using
 * Cholesky decomposition. Inverse is returned in aout, determinant in det_a.
 */
/*
void invdet_posdef(arma::mat *a, arma::mat *aout, double *det_a) {
    bool posdef;
    int i, j, n= a->n_cols;

    arma::mat b(n,n);
    choldc_inv(a, b, &posdef);
    *det_a = 1.0;
    for (i = 0; i < n; i++) {
        double value;
        value = b->at(i,i);
        (*det_a) *= 1 / (value * value);
    }

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += b->at(k,i) * b->at(k,j);
            }
            aout->at(i,j) = sum;
        }
    }

    for (i = 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            aout->at(i,j) = aout->at(j,i);
        }
    }
}
*/


/*
 * Inverse of a positive definite matrix with inverse of Cholesky decomposition
 * stored in invchol.
 * Result is returned in aout.
 *
 * Example:
 *   choldc_inv(a,n,invchol,posdef);
 *   inv_posdef_chol(invchol,n,ainv);
 */
void inv_posdef_chol(double **invchol,
                     int n,
                     double **aout)
{
    int i, j;

    //assert(invchol != NULL);
    //assert(aout != NULL);

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) {
                sum += invchol[k][i] * invchol[k][j];
            }
            aout[i][j] = sum;
        }
    }

    for (i = 2; i <= n; i++) {
        for (j = 1; j < i; j++) {
            aout[i][j] = aout[j][i];
        }
    }
}

/* Cholesky decomposition

Given a positive-definite symmetric matrix A, this routine constructs its Cholesky
decomposition, A = L * L' . On input, only the upper triangle of a need be given;
The Cholesky factor L is returned in the lower triangle of aout (upper-diag elem are set to 0) 

*/
void choldc(arma::mat *A, arma::mat *cholA, bool *posdef) {

  int i, j, k, n= A->n_cols;
  double sum, max_a;

  *posdef= true;
  for (i=0; i<n; i++) { for (j=i; j<n; j++) { cholA->at(i,j)= A->at(i,j); } }  //copy A into cholA
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) {
      for (sum=cholA->at(i,j),k=i-1; k>=0; k--) sum -= cholA->at(i,k)*cholA->at(j,k);
      if (i == j) {
	      if (sum <= 0.0) *posdef= false;
	      cholA->at(i,i)=sqrt(sum);
      } else {
	      max_a=max_xy(fabs(cholA->at(i,i)), 1e-10);
	      cholA->at(j,i)=sum/max_a;
      }
    }
  }
  for (i=0; i<n; i++) { for (j=i+1; j<n; j++) { cholA->at(i,j)= 0; } }  //set upper-diagonal elem to 0

}

/*Computes inverse of Cholesky matrix cholA and stores the result in cholA */
void choldc_inv_internal(arma::mat *cholA) {
  int i, j, k, n= cholA->n_rows;
  double sum, max_a;
  for (i=0; i<n; i++) {
    max_a=max_xy(cholA->at(i,i), 1e-10);
    cholA->at(i,i)=1.0/max_a;
    for (j=i+1; j<n; j++) {
      sum=0.0;
      for (k=i; k<j; k++) sum -= cholA->at(j,k) * cholA->at(k,i);
      max_a= max_xy(cholA->at(j,j), 1e-10);
      cholA->at(j,i)= sum/max_a;
    }
  }
}

/* Inverse of submatrix cholA[0:n,0:n], where cholA is a Cholesky matrix

  INPUT
  - cholA
  - n: number of rows and columns of cholA used to obtain the submatrix

  OUTPUT
  - cholAinv: inverse of cholA[0:n,0:n]

*/
void choldc_inv_internal(arma::mat *cholA, int n, arma::mat *cholAinv) {

  int i, j, k;
  double sum, max_a;
  for (i=0; i<=n; i++) {
    max_a= max_xy(cholA->at(i,i), 1e-10);
    cholAinv->at(i,i)= 1.0/max_a;
    for (j=i+1; j<=n; j++) {
      sum= 0.0;
      for (k=i; k<j; k++) sum -= cholA->at(j,k) * cholAinv->at(k,i);
      max_a= max_xy(cholA->at(j,j), 1e-10);
      cholAinv->at(j,i)= sum / max_a;
    }
  }

}           



/* Inverse of Cholesky decomposition of A 

That is, if A=L * L' it returns the inverse of L
 (note that inv(A)= inv(L)' * inv(L))

INPUT
  - A: input matrix
OUTPUT
  - cholAinv: inverse of Cholesky decomposition of A
  - posdef: true if A is positive definite, false otherwise

*/
void choldcinv(arma::mat *cholAinv, bool *posdef, arma::mat *A) {
  choldc(A, cholAinv, posdef);
  if (*posdef) {
    choldc_inv_internal(cholAinv);
  }
}

/* Same, but also returns cholA */
void choldcinv(arma::mat *cholAinv, arma::mat *cholA, bool *posdef, arma::mat *A) {
  choldc(A, cholA, posdef);
  if (*posdef) {
    choldc_inv_internal(cholA, cholA->n_cols - 1, cholAinv);           //Inverse of cholA[0:n,0:n], where n= cholA->n_cols - 1
  }
}

/* Inverse of a positive definite matrix A, inverse of its Cholesky decomposition and log(det(Ainv)) */
void choldcinv_det(arma::mat *Ainv, arma::mat *cholAinv, double *logdet_Ainv, arma::mat *A) {
  bool posdef;
  int i, npar= A->n_cols;
  choldcinv(cholAinv, &posdef, A);
  if (Ainv != nullptr) tchol_times_chol(cholAinv, Ainv); // (*Ainv)= cholAinv->t() * (*cholAinv);
  for (i=0, (*logdet_Ainv=0); i<npar; i++) (*logdet_Ainv) += log(cholAinv->at(i,i));
  (*logdet_Ainv) *= 2;

}

/* Inverse, Cholesky decomposition and determinant of A. Also, Cholesky decomposition of A */
void choldcinv_det(arma::mat *Ainv, arma::mat *cholAinv, double *logdet_Ainv, arma::mat *cholA, arma::mat *A) {
  bool posdef;
  int i, npar= A->n_cols;

  choldcinv(cholAinv, cholA, &posdef, A);
  if (Ainv != nullptr) {
    tchol_times_chol(cholAinv, Ainv); // (*Ainv)= cholAinv->t() * (*cholAinv);
  }
  for (i=0, (*logdet_Ainv=0); i<npar; i++) (*logdet_Ainv) += log(cholAinv->at(i,i));
  (*logdet_Ainv) *= 2;

}

//Obtain ans= cholA * b, using that cholA is lower-triangular
void chol_times_vec(arma::mat *cholA, arma::vec *b, arma::vec *ans) {
  int i, j, n= cholA->n_rows;
  double sum;
  for (i=0; i<n; i++) {
    for (j=0, sum=0; j<=i; j++) sum += cholA->at(i,j) * b->at(j);
    ans->at(i)= sum;
  }
}

//Obtain ans= cholA * B, using that cholA is lower-triangular
void chol_times_mat(arma::mat *cholA, arma::mat *B, arma::mat *ans) {
  int i, j, k, nrowA= cholA->n_rows, ncolB= B->n_cols;
  double sum;
  for (i=0; i<nrowA; i++) {
    for (j=0; j<ncolB; j++) {
      for (k=0, sum=0; k<=i; k++) sum += cholA->at(i,k) * B->at(k,j);
      ans->at(i,j)= sum;
    }
  }
}

//Obtain ans= cholA^T cholA, using that cholA is lower-triangular
void tchol_times_chol(arma::mat *cholA, arma::mat *ans) {
  int i, j, k, n= cholA->n_rows;
  double sum;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      for (k= max_xy(i,j), sum=0; k<n; k++) sum += cholA->at(k,i) * cholA->at(k,j);
      ans->at(i,j)= sum;
    }
  }
}



/* Cholesky decomposition of A= L L^T + x x^T. The output is returned in the lower-triangular part of L */
void choldc_rank1_update(arma::mat *L, arma::vec *x) {

  int k, l, p = L->n_cols;
  double c, r, s;

  for (k=0; k<p; k++) {
    r = sqrt(pow(L->at(k,k),2) + pow(x->at(k),2));
    c = r / L->at(k,k);
    s = x->at(k) / L->at(k,k);
    L->at(k,k) = r;
    if (k < p) {
      for (l=k+1; l<p; l++) {
        L->at(l,k) = (L->at(l,k) + s * x->at(l)) / c;
        x->at(l) = c * x->at(l) - s * L->at(l,k);
      }
    }
  }

}   

/* Same as previous choldc_rank1_update, but result is returned in cholA[cholA_rowini:cholA_rowfi, cholA_rowini:cholA_rowfi]

  The Cholesky decomposition being updated is stored in L[L_rowini:L_rowfi, L_rowini:L_rowfi]
  and the rank 1 update is given by x[x_rowini:x_rowfi]

*/
void choldc_rank1_update(arma::mat *cholA, int cholA_rowini, int cholA_rowfi, arma::mat *L, int L_rowini, int L_rowfi, double *x, int x_rowini, int x_rowfi) {

  int k, l, p = cholA_rowfi - cholA_rowini, *L_index, *x_index;
  double c, d, r, s;

  if (p != L_rowfi - L_rowini) Rf_error("In choldc_rank1_update: dimensions of cholA and L do not match");
  if (p != x_rowfi - x_rowini) Rf_error("In choldc_rank1_update: dimensions of cholA and x do not match");

  //Map indexes of row/column in cholA to those in L and x
  L_index= ivector(cholA_rowini, cholA_rowfi);
  x_index= ivector(cholA_rowini, cholA_rowfi);
  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    L_index[k]= k - cholA_rowini + L_rowini;
    x_index[k]= k - cholA_rowini + x_rowini;
  }

  //Copy L[L_rowini,L_rowfi] into cholA[cholA_rowini:cholA_rowfi, cholA_rowini:cholA_rowfi]
  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    for (l= cholA_rowini; l <= k; l++) {
      cholA->at(k,l)= L->at(L_index[k], L_index[l]);
    }
  }

  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    d= cholA->at(k,k);
    r = sqrt(pow(d,2) + pow(x[x_index[k]],2));
    cholA->at(k,k) = r;
    if (k < cholA_rowfi) {
      c = r / d; 
      s = x[x_index[k]] / d;
      for (l= k+1; l <= cholA_rowfi; l++) {
        cholA->at(l,k) = (cholA->at(l,k) + s * x[x_index[l]]) / c;
        x[x_index[l]] = c * x[x_index[l]] - s * cholA->at(l,k);
      }
    }
  }

  free_ivector(L_index, cholA_rowini, cholA_rowfi);
  free_ivector(x_index, cholA_rowini, cholA_rowfi);

}


/* Cholesky decomposition of A= L L^T - x x^T. The output is returned in the lower-triangular part of L */
void choldc_rank1_downdate(arma::mat *L, arma::vec *x) {

  int k, l, p = L->n_cols;
  double c, r, s;

  for (k=0; k<p; k++) {
    r = sqrt(pow(L->at(k,k),2) - pow(x->at(k),2));
    c = r / L->at(k,k);
    s = x->at(k) / L->at(k,k);
    L->at(k,k) = r;
    if (k < p) {
      for (l=k+1; l<p; l++) {
        L->at(l,k) = (L->at(l,k) - s * x->at(l)) / c;
        x->at(l) = c * x->at(l) - s * L->at(l,k);
      }
    }
  }
  
}


/* Same as previous choldc_rank1_downdate, but result is returned in cholA[cholA_rowini:cholA_rowfi, cholA_rowini:cholA_rowfi]

  The Cholesky decomposition being updated is stored in L[L_rowini:L_rowfi, L_rowini:L_rowfi]
  and the rank 1 update is given by x[x_rowini:x_rowfi]

*/
void choldc_rank1_downdate(arma::mat *cholA, int cholA_rowini, int cholA_rowfi, arma::mat *L, int L_rowini, int L_rowfi, double *x, int x_rowini, int x_rowfi) {

  int k, l, p = cholA_rowfi - cholA_rowini, *L_index, *x_index;
  double c, d, r, s;

  if (p != L_rowfi - L_rowini) Rf_error("In choldc_rank1_update: dimensions of cholA and L do not match");
  if (p != x_rowfi - x_rowini) Rf_error("In choldc_rank1_update: dimensions of cholA and x do not match");

  //Map indexes of row/column in cholA to those in L and x
  L_index= ivector(cholA_rowini, cholA_rowfi);
  x_index= ivector(cholA_rowini, cholA_rowfi);
  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    L_index[k]= k - cholA_rowini + L_rowini;
    x_index[k]= k - cholA_rowini + x_rowini;
  }

  //Copy L[L_rowini,L_rowfi] into cholA[cholA_rowini:cholA_rowfi, cholA_rowini:cholA_rowfi]
  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    for (l= cholA_rowini; l <= k; l++) {
      cholA->at(k,l)= L->at(L_index[k], L_index[l]);
    }
  }

  for (k= cholA_rowini; k <= cholA_rowfi; k++) {
    d= cholA->at(k,k);
    r = sqrt(pow(d,2) - pow(x[x_index[k]],2));
    cholA->at(k,k) = r;
    if (k < cholA_rowfi) {
      c = r / d; 
      s = x[x_index[k]] / d;
      for (l= k+1; l <= cholA_rowfi; l++) {
        cholA->at(l,k) = (cholA->at(l,k) - s * x[x_index[l]]) / c;
        x[x_index[l]] = c * x[x_index[l]] - s * cholA->at(l,k);
      }
    }
  }

  free_ivector(L_index, cholA_rowini, cholA_rowfi);
  free_ivector(x_index, cholA_rowini, cholA_rowfi);

}


/* Cholesky decomposition after adding 1 row/column
  This function returns the Cholesky decomposition of B, where
        (A11 A12 A13)
  B =   (A21 A22 A23)
        (A31 A32 A33)
  given cholA, which is the Cholesky decomposition of
  A =   (A11 A13)
        (A31 A33)
 where x= (A12, A22, A23)^T is the d^th column of B

 In the comments below, the submatrices of cholB are denoted as

            (S11,   0,   0)
  cholB =   (S12, S22,   0)
            (S13, S23, S33)
*/
void choldc_addrow(arma::mat *cholA, double *x, int d, arma::mat *cholB) {
  int i, j, p= cholA->n_cols;
  double mysum;

  // Obtain S11, S12, S13, S22
  mysum = 0;
  if (d > 0) { //if d is the first column in B, then (S11, S12, S13) are not needed

    arma::mat L11inv(d,d);
    choldc_inv_internal(cholA, d-1, &L11inv);  //Inverse of cholA[0:d-1,0:d-1]
    for (i=0; i <= d-1; i++) {
      cholB->at(d,i)= 0;
      for (j=0; j<=i; j++) {
        cholB->at(i,j)= cholA->at(i,j); //S11
        cholB->at(d,i) += L11inv.at(i,j) * x[j]; //S12
      }
      for (j=d+1; j <= p; j++) cholB->at(j,i)= cholA->at(j-1,i); //S13
      mysum += pow(cholB->at(d,i), 2);
    }

  }

  cholB->at(d,d)= sqrt(x[d] - mysum); //S22

  if (d+1 <= p) { //if d is the last column in B, then (S13, S23, S33) are not needed

    // Obtain S23
    for (i=d+1; i <= p; i++) {
      mysum = 0;
      for (j=0; j <= d-1; j++) mysum += cholB->at(i,j) * cholB->at(d,j);
      cholB->at(i,d)= (x[i] - mysum) / cholB->at(d,d);
    }

    // Obtain S33=cholB[d+1:p,d+1:p], given by Cholesky decomposition of L^T L - x x^T, where L=cholA[d:p-1,d:p-1],  x= cholB[d+1:p,d]
    double *z, *cold= cholB->colptr(d);
    z= dvector(d+1, p);
    for (j= d+1; j <= p; j++) z[j]= cold[j];
    choldc_rank1_downdate(cholB, d+1, p, cholA, d, p-1, z, d+1, p);
    free_dvector(z, d+1, p);

  }

}


/* Cholesky decomposition of A=B[-d,-d], given cholB = Cholesky decomposition of B

  INPUT
  - cholB: Cholesky decomposition of B
  - d: row/column of B that is dropped to obtain A

  OUTPUT
  - cholA: Cholesky decomposition of A

*/
void choldc_droprow(arma::mat *cholB, int d, arma::mat *cholA) {
  int i, j, p= cholB->n_cols;

  for (i=0; i<d; i++) {
    for (j=0; j<=i; j++) cholA->at(i,j)= cholB->at(i,j);       //L11
    for (j=d; j<p-1; j++) cholA->at(j,i)= cholB->at(j+1,i); //L13 
  }

  //Obtain cholA[d-1:p-2, d-1:p-2], given by Cholesky decomposition of L L^T + x x^T, where L=cholB[d:(p-1),d:(p-1)], x= cholB[d+1:p-1,d]
  double *z, *S23= cholB->colptr(d);
  if (d <= p-2) {
    z= dvector(d+1, p-1);
    for (j= d+1; j <= p-1; j++) z[j]= S23[j];
    choldc_rank1_update(cholA, d, p-2, cholB, d+1, p-1, z, d+1, p-1);
    free_dvector(z, d+1, p-1);
  }

}


/* Inverse, Cholesky decomp & determinant of A=B[-d,-d] 

  INPUT
  - cholB: Cholesky decomposition of B (lower-triangular)
  - d: index of row/column being dropped from B

  OUTPUT
  - Ainv: inverse of A (if nullptr, this is not returned)
  - cholAinv: Cholesky decomposition of Ainv (if nullptr, this is not returned)
  - logdet_Ainv: logarithm of the determinant of Ainv
  - cholA: Cholesky decomposition of A

*/
void choldcinv_det_droprow(arma::mat *Ainv, arma::mat *cholAinv, double *logdet_Ainv, arma::mat *cholA, arma::mat *cholB, int d) {
  int i, ncolA= cholA->n_rows;

  //Obtain cholA
  choldc_droprow(cholB, d, cholA); 

  //Obtain logdet_Ainv
  if (logdet_Ainv != nullptr) { 
    for (i=0, (*logdet_Ainv=0); i < ncolA; i++) (*logdet_Ainv) -= log(cholA->at(i,i));
    (*logdet_Ainv) *= 2;
  }

  //Obtain cholAinv
  if (cholAinv != nullptr) {
    choldc_inv_internal(cholA, ncolA - 1, cholAinv);

    //Obtain Ainv
    if (Ainv != nullptr) (*Ainv)= cholAinv->t() * (*cholAinv);

  }

} 

/* Inverse, Cholesky decomp & det of B, given Cholesky decomp of A= B[-d,-d] 

  INPUT
  - cholA: Cholesky decomposition of A
  - x: column of B that is missing in A. It must have length equal to the number of rows/columns in B
  - d: column in B that is dropped to obtain A

  OUTPUT
  - Binv: inverse of B. If nullptr, this is not computed.
  - cholBinv: Cholesky decomposition of Binv. If nullptr, this is not computed.
  - logdet_Binv: log-determinant of Binv
  - cholB: Cholesky decomposition of B

*/
void choldcinv_det_addrow(arma::mat *Binv, arma::mat *cholBinv, double *logdet_Binv, arma::mat *cholB, arma::mat *cholA, double *x, int d) {
  int i, ncolB= cholB->n_cols;

  //Obtain cholB
  choldc_addrow(cholA, x, d, cholB);

  //Obtain logdet_Binv
  for (i=0, (*logdet_Binv=0); i < ncolB; i++) (*logdet_Binv) -= log(cholB->at(i,i));
  (*logdet_Binv) *= 2;

  //Obtain cholBinv
  if (cholBinv != nullptr) {
    choldc_inv_internal(cholB, ncolB - 1, cholBinv);

    //Obtain Binv
    if (Binv != nullptr) (*Binv)= cholBinv->t() * (*cholBinv);

  }

}



/* Update the inverse of a symmetric matrix A when a column/row newcol is changed

INPUT
  - Ainv: inverse of A before A[,newcol] was changed
  - A_newcol: newvalues for A[,newcol]
  - newcol: index of the column of A being changed

OUTPUT
  - Ainv: inverse of A after A[,newcol] is changed

*/
void update_inverse(arma::mat *Ainv, arma::sp_mat *A_newcol, int *newcol) {

  int i, j, d= Ainv->n_cols - 1;
  double gam = A_newcol->at(*newcol), tmp;

  arma::mat A11inv(d,d), Ainv_newcol;
  Ainv_newcol= Ainv->col(*newcol);
  double Ainv_newcol_diag= Ainv_newcol.at(*newcol,0);
  Ainv_newcol.shed_row(*newcol);

  //Obtain A11inv= Ainv[-newcol, -newcol] - Ainv_newcol * Ainv_newcol.t() / Ainv[newcol, newcol]
  for (i=0; i < *newcol; i++) {
    tmp= Ainv_newcol.at(i) / Ainv_newcol_diag;
    for (j=0; j < *newcol; j++) A11inv.at(i,j)= A11inv.at(j,i)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
    for (j= *newcol + 1; j <= i; j++) A11inv.at(i,j-1)= A11inv.at(j-1,i)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
  }
  for (i= *newcol + 1; i <= d; i++) {
    tmp= Ainv_newcol.at(i-1) / Ainv_newcol_diag;
    for (j=0; j < *newcol; j++) A11inv.at(i-1,j)= A11inv.at(j,i-1)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
    for (j= *newcol + 1; j <= i; j++) A11inv.at(i-1,j-1)= A11inv.at(j-1,i-1)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j-1);
  }

  //omega_12= A_newcol[-newcol,1,drop=FALSE]
  arma::sp_mat omega_12= *A_newcol;
  omega_12.shed_row(*newcol);

  arma::mat v= A11inv * omega_12;

  //gam = A_newcol[newcol,1] - as.double(t(v) %*% omega_12)
  gam -= arma::as_scalar(v.t() * omega_12);

  //RETURN Ainv

  //1. Diagonal entry
  Ainv->at(*newcol, *newcol)= 1 / gam;

  //2. Column newcol
  //Ainv[-newcol,newcol] = - v / gam
  //Ainv[newcol,-newcol] = - v / gam
  arma::vec v_over_gam(d);
  for (i=0; i < *newcol; i++) {
    v_over_gam.at(i)= v.at(i,0) / gam;
    Ainv->at(*newcol, i)= Ainv->at(i, *newcol)= -v_over_gam.at(i);
  }
  for (i= *newcol + 1; i <= d; i++) {
    v_over_gam.at(i-1)= v.at(i-1,0) / gam;
    Ainv->at(*newcol, i)= Ainv->at(i, *newcol)= -v_over_gam.at(i-1);
  }

  //3. Entries outside row/column newcol
  //Ainv[ -newcol, -newcol ] = A11inv + v %*% t( v ) / gam
  for (i=0; i < *newcol; i++) {
    for (j=0; j < *newcol; j++) {
      Ainv->at(i,j)= Ainv->at(j,i)= A11inv.at(i,j) + v.at(i) * v_over_gam.at(j);
    }
  } 
  for (i= *newcol + 1; i <= d; i++) {
    for (j=0; j < *newcol; j++) {
      Ainv->at(i,j)= Ainv->at(j,i)= A11inv.at(i-1,j) + v.at(i-1) * v_over_gam.at(j);
    }
    for (j= *newcol + 1; j <=d; j++) {
      Ainv->at(i,j)= Ainv->at(j,i)= A11inv.at(i-1,j-1) + v.at(i-1) * v_over_gam.at(j-1);
    }

  } 

}


//Same as update_inverse, but result is stored in Ainv[indexes[i],indexes[j]] instead of A[i,j]
void update_inverse_permutedindex(arma::mat *Ainv, int *indexes, arma::sp_mat *A_newcol, int *newcol) {

  int i, j, d= Ainv->n_cols - 1;
  double gam = A_newcol->at(*newcol), tmp;

  arma::mat A11inv(d,d), Ainv_newcol;
  Ainv_newcol= Ainv->col(*newcol);
  double Ainv_newcol_diag= Ainv_newcol.at(*newcol,0);
  Ainv_newcol.shed_row(*newcol);

  //Obtain A11inv= Ainv[-newcol, -newcol] - Ainv_newcol * Ainv_newcol.t() / Ainv[newcol, newcol]
  for (i=0; i < *newcol; i++) {
    tmp= Ainv_newcol.at(i) / Ainv_newcol_diag;
    for (j=0; j < *newcol; j++) A11inv.at(i,j)= A11inv.at(j,i)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
    for (j= *newcol + 1; j <= i; j++) A11inv.at(i,j-1)= A11inv.at(j-1,i)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
  }
  for (i= *newcol + 1; i <= d; i++) {
    tmp= Ainv_newcol.at(i-1) / Ainv_newcol_diag;
    for (j=0; j < *newcol; j++) A11inv.at(i-1,j)= A11inv.at(j,i-1)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j);
    for (j= *newcol + 1; j <= i; j++) A11inv.at(i-1,j-1)= A11inv.at(j-1,i-1)= Ainv->at(i,j) - tmp * Ainv_newcol.at(j-1);
  }

  //omega_12= A_newcol[-newcol,1,drop=FALSE]
  arma::sp_mat omega_12= *A_newcol;
  omega_12.shed_row(*newcol);

  arma::mat v= A11inv * omega_12;

  //gam = A_newcol[newcol,1] - as.double(t(v) %*% omega_12)
  gam -= arma::as_scalar(v.t() * omega_12);

  //RETURN Ainv permuted

  //1. Diagonal entry
  Ainv->at(indexes[*newcol], indexes[*newcol])= 1 / gam;

  //2. Column newcol
  arma::vec v_over_gam(d);
  for (i=0; i < *newcol; i++) {
    v_over_gam.at(i)= v.at(i,0) / gam;
    Ainv->at(indexes[*newcol], indexes[i])= Ainv->at(indexes[i], indexes[*newcol])= -v_over_gam.at(i);
  }
  for (i= *newcol + 1; i <= d; i++) {
    v_over_gam.at(i-1)= v.at(i-1,0) / gam;
    Ainv->at(indexes[*newcol], indexes[i])= Ainv->at(indexes[i], indexes[*newcol])= -v_over_gam.at(i-1);
  }

  //3. Entries outside row/column newcol
  for (i=0; i < *newcol; i++) {
    for (j=0; j < *newcol; j++) {
      Ainv->at(indexes[i],indexes[j])= Ainv->at(indexes[j],indexes[i])= A11inv.at(i,j) + v.at(i) * v_over_gam.at(j);
    }
  } 
  for (i= *newcol + 1; i <= d; i++) {
    for (j=0; j < *newcol; j++) {
      Ainv->at(indexes[i],indexes[j])= Ainv->at(indexes[j],indexes[i])= A11inv.at(i-1,j) + v.at(i-1) * v_over_gam.at(j);
    }
    for (j= *newcol + 1; j <=d; j++) {
      Ainv->at(indexes[i],indexes[j])= Ainv->at(indexes[j],indexes[i])= A11inv.at(i-1,j-1) + v.at(i-1) * v_over_gam.at(j-1);
    }

  } 

}




/* Update the inverse of symmetric matrix A after replacing row/col colid by difvals, given current inverse of A

  The function uses fast rank 1 computations based on the Woodbury formula

  Updated A adds difvals to A[colid,] then adds v2 to A[,colid], where v2= difvals with a zero in position colid

  INPUT
  - Ainv: inverse of A before the update
  - difvals: new values of A[,colid] minus previous values of A[,colid] (and A[colid,], by symmetry)
  - colid: row being updated

  OUPUT: Ainv stores the inverse of A after the update

*/
void symmat_inv_colupdate(arma::mat *Ainv, arma::sp_mat *difvals, int colid) {
  updateinv_rank1(Ainv, difvals, colid); //inverse of A + difvals e^T
  double tmp= difvals->at(colid,0);
  difvals->at(colid,0)= 0;
  updateinv_rank1(Ainv, colid, difvals); //inverse of A + difvals e^T + e v2^T
  difvals->at(colid,0)= tmp;
}



/*Compute inverse of A + u v^T, given Ainv=A^{-1}

 INPUT
 - Ainv: A^{-1}
 - u: column vector giving the rank 1 update
 - v: column vector giving the rank 1 update

 OUTPUT: Ainv stores the inverse of A + u v^T

 Details: (A + u v^T)^{-1}= A^{-1} - A^{-1} u v^T A^{-1} / (1 + v^T A^{-1} u)
*/
void updateinv_rank1(arma::mat *Ainv, arma::sp_mat *u, arma::sp_mat *v) {
  arma::mat Ainvu= (*Ainv) * (*u);
  arma::mat vTAinvu= v->t() * Ainvu;
  double den= 1.0 + vTAinvu.at(0,0);
 
  arma::sp_mat vden= (*v) / den;
  (*Ainv)= (*Ainv) - Ainvu * vden.t() * (*Ainv);
}


/* Particular case of updateinv_rank1 where only 1 row in A is updated by adding v^T

 This is equivalent to setting u^T=(0,..,0,1,0,..,0) where 1 is in entry rowid
*/
void updateinv_rank1(arma::mat *Ainv, int rowid, arma::sp_mat *v) {
  arma::mat Ainvu= Ainv->col(rowid);
  arma::mat vTAinvu= v->t() * Ainvu;
  double den= 1.0 + vTAinvu.at(0,0);

  arma::sp_mat vden= (*v) / den;
  (*Ainv)= (*Ainv) - Ainvu * vden.t() * (*Ainv);
}

/* Particular case of updateinv_rank1 where only column A[,colid] is updated by adding u

 This is equivalent to setting v^T=(0,..,0,1,0,..,0) where 1 is in entry colid
*/
void updateinv_rank1(arma::mat *Ainv, arma::sp_mat *u, int colid) {
  arma::mat vTAinv= Ainv->row(colid);
  arma::mat vTAinvu= vTAinv * (*u);
  double den= 1.0 + vTAinvu.at(0,0);

  arma::sp_mat uden= (*u)/den;
  (*Ainv)= (*Ainv) - (*Ainv) * uden * vTAinv; 
}



/*
 * LU decomposition, Inverse and determinant of a non-singular matrix.
 * Given a matrix a[1..n][1..n], replace it by the LU decomposition of a
 * row-wise permutation of itself. a and n are input. a is output, arranged
 * as in equation (2.3.14) above; indx[1..n] is an output vector that records
 * the row permutation e ected by the partial pivoting; d is output as 1
 * depending on whether the number of row interchanges was even or odd,
 * respectively.
 *
 * Used in combination with lu_solve() to solve linear equations or invert
 * a matrix.
 */
void ludc(double **a,
          int n,
          int *indx,
          double *d)
{
    const double TINY = 1.0e-20;
    int i, j;
    int imax = 1;
    double big;
    double *vv;                       /* Stores implicit scaling of each row */

    //assert(a != NULL);
    //assert(indx != NULL);
    //assert(d != NULL);

    vv = dvector(1, n);
    *d = 1.0;                          /* No row interchanges yet */

    /* Loop over rows to get the implicit scaling information */
    for (i = 1; i <= n; i++) {
        big = 0.0;
        for (j = 1; j <= n; j++) {
            double temp;

            temp = fabs(a[i][j]);
            if (temp > big) {
                big = temp;
            }
        }
        if (big == 0.0) {
            /* No nonzero largest element */
            nrerror("ludc", "", "singular matrix detected");
            /*NOTREACHED*/
        }
        vv[i] = 1.0 / big;             /* Save the scaling */
    }

    /* This is the loop over columns of Crout's method */
    for (j = 1; j <= n; j++) {
        int k;
        double sum;
        double dum;

        /* Equation (2.3.12) except for i = j */
        for (i = 1; i < j; i++) {
            sum = a[i][j];
            for (k = 1; k < i; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
        }

        /* Initialize for the search for largest pivot element */
        big = 0.0;

        /* i = j of equation (2.3.12) and i = j+1 : ::N of equation (2.3.13) */
        for (i = j; i <= n; i++) {
            sum = a[i][j];
            for (k = 1; k < j; k++) {
                sum -= a[i][k] * a[k][j];
            }
            a[i][j] = sum;
            /* Is figure of merit for the pivot better than the best so far? */
            dum = vv[i] * fabs(sum);
            if (dum >= big) {
                big = dum;
                imax = i;
            }
        }

        /* Do we need to interchange rows? */
        if (j != imax) {
            for (k = 1; k <= n; k++) { /* Yes, do so... */
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
            }
            *d = -(*d);                /* Change the parity of d */
            vv[imax] = vv[j];          /* Interchange the scale factor */
        }
        indx[j] = imax;
        /*
         * If the pivot element is zero, the matrix is singular (at least to
         * the precision of the algorithm). For some applications on singular
         * matrices, it is desirable to substitute TINY for zero.
         */
        if (a[j][j] == 0.0) {
            a[j][j] = TINY;
        }

        if (j != n) {
            /* Now, finally, divide by the pivot element */
            dum = 1.0 / a[j][j];
            for (i = j+1; i <= n; i++) {
                a[i][j] *= dum;
            }
        }
        /* Go back for the next column in the reduction... */
    }
    free_dvector(vv, 1, n);
}


/*
 * Solves the set of n linear equations A X = B. Here a[1..n][1..n] is input,
 * not as the matrix A but rather as its LU decomposition, determined by the
 * routine ludc(). indx[1..n] is input as the permutation vector returned by
 * ludcmp(). b[1..n] is input as the right-hand side vector B, and returns
 * with the solution vector X. a, n, and indx are not modified by this routine
 * and can be left in place for successive calls with different right-hand
 * sides b. This routine takes into account the possibility that b will begin
 * with many zero elements, so it is efficient for use in matrix inversion.
 *
 * Usage:
 *   ludc(a,n,indx,&d);
 *   lu_solve(a,n,indx,b);  //answer in b, original matrix A has been destroyed
 *   lu_solve(a,n,indx,c);  //now we solve AX=c,
 *                          //since A hasn't changed we only call lu_solve
 */
void lu_solve(double **a,
              int n,
              const int *indx,
              double b[])
{
    int i, j, ii=0;
    double sum;

    //assert(a != NULL);
    //assert(indx != NULL);

    for (i = 1; i <= n; i++) {
        int ip;

        /*
         * When ii is set to a positive value, it will become the
         * index of the first non-vanishing element of b. We now
         * do the forward substitution, equation (2.3.6).
         * The only new wrinkle is to unscramble the permutation as we go.
         */
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii) {
            for (j = ii; j <= i-1; j++) {
                sum -= a[i][j] * b[j];
            }
        }
        else if (sum != 0.0) {
            /* Nonzero element encountered, from now do sums in above loop */
            ii = i;
        }
        b[i] = sum;
    }

    /* Now we do the backsubstitution, equation (2.3.7) */
    for (i = n; i >= 1; i--) {
        sum = b[i];
        for (j = i+1; j <= n; j++) {
            sum -= a[i][j] * b[j];
        }
        /* Store a component of the solution vector X */
        b[i] = sum / a[i][i];
    }
}


/*
 * Inverse of a non-singular matrix a.
 * Result is stored in aout, original matrix a is destroyed.
 */
void lu_inverse(double **a,
                int n,
                double **aout)
{
    int i, j;
    double *col;
    int *indx;

    //assert(a != NULL);
    //assert(aout != NULL);

    col  = dvector(1, n);
    indx = ivector(1, n);

    /* Decompose the matrix just once */
    {
        double d;

        ludc(a, n, indx, &d);
    }

    /* Find inverse by columns */
    for (j = 1; j <= n; j++) {
        for (i = 1; i <= n; i++) {
            col[i] = 0.0;
        }
        col[j] = 1.0;
        lu_solve(a, n, indx, col);
        for (i = 1; i <= n; i++) {
            aout[i][j] = col[i];
        }
    }
    free_dvector(col,  1, n);
    free_ivector(indx, 1, n);
}


/*
 * Determinant of a matrix a.
 * Original matrix a is destroyed and its LU decomposition is returned.
 */
double lu_det(double **a,
              int n)
{
    double d;
    int j;
    int *indx;

    //assert(a != NULL);

    indx = ivector(1, n);
    ludc(a, n, indx, &d); /* returns d as +/-1 */
    for (j = 1; j <= n; j++) {
        d *= a[j][j];
    }
    free_ivector(indx, 1, n);
    return(d);
}


void eigenvals(double **a, int n, double *vals) {
  //Computes eigenvalues of a real, symmetric matrix a[1..n][1..n]
  int i,j;
  double *e, **z, **b;

  b= dmatrix(1,n,1,n);
  for (i=1; i<=n; i++) for (j=1; j<=n; j++) b[i][j]= a[i][j];

  e= dvector(1,n);
  z= dmatrix(1,n,1,n);

  tred2(b, n, vals, e, false);
  tqli(vals, e, n, z, false);

  free_dmatrix(b,1,n,1,n);
  free_dmatrix(z,1,n,1,n);
  free_dvector(e,1,n);
}




void tred2(double **a, int n, double d[], double e[], bool getVecs) {
  /*  Householder reduction of a real, symmetric matrix a[1..n][1..n]. On output, a is replaced by the orthogonal matrix Q effecting the transformation. d[1..n] returns the diagonal ele- ments of the tridiagonal matrix, and e[1..n] the off-diagonal elements, with e[1]=0.
      Set getVecs=true to obtain eigenvectors, getVecs=false to get only eigenvalues */
  int l,k,j,i;
  double scale,hh,h,g,f;
  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++) scale += fabs(a[i][k]);
      if (scale == 0.0) {  //Skip transformation
	e[i]=a[i][l];
      } else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale; //Use scaled a's for transformation
	  h += a[i][k]*a[i][k]; //Form sigma in h.
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g; //Now h is equation (11.2.4)
	a[i][l]=f-g; //Store u/H in the ith row of a.
	f=0.0;
	for (j=1;j<=l;j++) {
	  a[j][i]=a[i][j]/h;   // This line could be omitted if getVecs=false
	  g=0.0;
	  for (k=1;k<=j;k++) g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++) g += a[k][j]*a[i][k];
	  e[j]=g/h;  //Form element of p in temporarily unused element of e
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++) a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;   // This line could be omitted if getEigen=false
  e[1]=0.0;
  if (getVecs) {
    for (i=1;i<=n;i++) { //Begin accumulation of transformation matrices.
      l=i-1;
      //This block skipped when i=1.
      //Use u and u/H stored in a to form QQQQQQQ.
      if (d[i]) {
        for (j=1;j<=l;j++) {
	  g=0.0;
	  for (k=1;k<=l;k++) g += a[i][k]*a[k][j];
	  for (k=1;k<=l;k++) a[k][j] -= g*a[k][i];
        }
      }
      d[i]=a[i][i];
      a[i][i]=1.0;
      for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
    }
  } else {
    for (i=1;i<=n;i++) d[i]=a[i][i];
  }
}


void tqli(double d[], double e[], int n, double **z, bool getVecs) {
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	iter++;
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SETSIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  /* calc eigenvectors iff getVecs=true */
	  if (getVecs) {
	    for (k=1;k<=n;k++) {
	      f=z[k][i+1];
	      z[k][i+1]=s*z[k][i]+c*f;
	      z[k][i]=c*z[k][i]-s*f;
	     }
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while ((m != l) && (iter<=30));
  }
}

double pythag(double a, double b) {
  //Note: DSQR and dsqrarg definition used to be in cstat.h, moved here to avoid compiler warning of dsqrarg being defined but not used
  #if !defined(DSQR)
  static double dsqrarg;
  #define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
  #endif

  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}



/* Comparison function used by qsort() for doubles */
int dcompare(const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}

int dcompare_decreasing(const void *a, const void *b) {
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}


/* Sorts double vector increasingly */
void dvecsort(double *v, int size) {
    qsort(v, size, sizeof(double), dcompare);
}

/* Sorts double vector decreasingly */
void dvecsort_decreasing(double *v, int size) {
    qsort(v, size, sizeof(double), dcompare_decreasing);
}


/*
 * Sorts vector of doubles x by rearranging values in index with quicksort
 * algorithm, e.g. x[index[ilo]], x[index[ilo+1]]... x[index[ihi]] is ordered
 *
 * Input:
 *   x    : vector of doubles to be ordered from position ilo to ihi
 *   index: vector of integers indexing the values of x
 *   ilo  : first element of x to order
 *   ihi  : last element of x to order
 *   incr : for incr==1 x is returned in increasing order;
 *          incr==-1 in decreasing order
 * Output:
 *   index: rearranged so that x[index[lo]], x[index[lo+1]]...x[index[ihi]]
 *          is ordered
 */
void dindexsort(double *x, int *index, int ilo, int ihi, int incr) {
    int pivot;              /* pivot value for partitioning array      */
    int uhi, ulo;           /* indices at ends of unpartitioned region */
    int tempEntry;          /* temporary entry used for swapping       */
    bool sortup, sortlo;    /* indicate if sub-vectors are sorted so
                             * no further subdivision is needed        */

    //assert(x != NULL);
    //assert(index != NULL);

    if (ilo >= ihi) {
        return;
    }

    sortup = sortlo = true;

    /* Select a pivot value */
    pivot = (ilo + ihi) / 2;

    /* Initialize ends of unpartitioned region */
    ulo = ilo;
    uhi = ihi;

    /* While the unpartitioned region is not empty, try to reduce its size */
    while (ulo < uhi) {
        if ((x[index[uhi]] * incr) > (x[index[pivot]] * incr)) {
            /* Check if upper subvector is ordered */
            if ((uhi < ihi) &&
                ((x[index[uhi]] * incr) > (x[index[uhi+1]] * incr))) {
                sortup = false;
            }

            /* Reduce the size of the unpartitioned region */
            uhi--;
            if ((uhi == pivot) && (ulo < pivot)) {
                tempEntry      = index[pivot];
                index[pivot]   = index[pivot-1];
                index[pivot-1] = tempEntry;
                pivot--;
            }
        }
        else {
            /* Swap entries at indices ulo and uhi */
            tempEntry  = index[ulo];
            index[ulo] = index[uhi];
            index[uhi] = tempEntry;

            if (pivot == ulo) {
                pivot = uhi;
            }

            /* Check if lower subvector is ordered */
            if ((ulo > ilo) &&
                ((x[index[ulo]] * incr) < (x[index[ulo-1]] * incr))) {
                sortlo = false;
            }

            /* Reduce the size of the unpartitioned region */
            ulo++;
            if ((ulo == pivot) && (uhi > (pivot+1))) {
                tempEntry      = index[pivot];
                index[pivot]   = index[pivot+1];
                index[pivot+1] = tempEntry;
                pivot++;
            }
        }
    }

    /*
     * Entries from ilo to pivot-1 are < or > pivot and
     * from pivot+1 to ihi are > or < pivot.
     * The two regions can be sorted recursively.
     */
    if ((sortlo == false) && (ilo < (pivot-1))) {
        dindexsort(x, index, ilo, pivot-1, incr);
    }
    if ((sortup == false) && (ihi > (pivot+1))) {
        dindexsort(x, index, pivot+1, ihi, incr);
    }
}


void iindexsort(int *x, int *index, int ilo, int ihi, int incr) {
/* Sorts vector of integers x by rearranging values in index with quicksort algorithm e.g. x[index[ilo]], x[index[ilo+1]]... x[index[ihi]] is ordered */
/* Input
   - x: vector of doubles that we want to order from position ilo to ihi
   - index: vector of integers indexing the values of x
   - ilo: first element of x we want to order
   - ihi: last element of x we want to order
   - incr: for incr==1 x is returned in increasing order; incr==-1 in decreasing order
   Output: vector index rearranged so that x[index[lo]], x[index[lo+1]]... x[index[ihi]] is ordered
*/

int pivot;              // pivot value for partitioning array
int ulo, uhi;           // indices at ends of unpartitioned region
int tempEntry;          // temporary entry used for swapping
int sortlo, sortup;      // indicate if sub-vectors are sorted so no further subdivision is needed

if (ilo >= ihi) { return; }

sortlo= sortup= 1;
pivot = (ilo + ihi)/2;                                                // Select a pivot value
ulo = ilo; uhi = ihi;                                                 // Initialize ends of unpartitioned region
// While the unpartitioned region is not empty, try to reduce its size.
while (ulo < uhi) {
  if ((x[index[uhi]]*incr) > (x[index[pivot]]*incr)) {                // Here, we can reduce the size of the unpartitioned region and try again.
    if ((uhi<ihi) && ((x[index[uhi]]*incr)>(x[index[uhi+1]]*incr))) sortup= 0;      // Check if upper subvector is ordered
    uhi--;
    if ((uhi==pivot) && (ulo<pivot))  { tempEntry= index[pivot]; index[pivot]= index[pivot-1]; index[pivot-1]= tempEntry; pivot--; }
  } else {                                                            // Here, x[index[uhi]] <= x[index[pivot]], so swap entries at indices ulo and uhi.
    tempEntry = index[ulo]; index[ulo] = index[uhi]; index[uhi] = tempEntry;
    if (pivot==ulo) pivot= uhi;
    if ((ulo>ilo) && ((x[index[ulo]]*incr)<(x[index[ulo-1]]*incr))) sortlo= 0;
    ulo++;                                                            // Reduce the size of the unpartitioned region
    if ((ulo==pivot) && (uhi>(pivot+1)))  { tempEntry= index[pivot]; index[pivot]= index[pivot+1]; index[pivot+1]= tempEntry; pivot++; }
  }
}

// Entries from ilo to pivot - 1 are < or > pivot and from pivot+1 to ihi are > or < pivot. The two regions can be sorted recursively.
if ((sortlo==0) && (ilo<(pivot-1))) iindexsort(x, index, ilo, pivot - 1, incr);
if ((sortup==0) && (ihi>(pivot+1))) iindexsort(x, index, pivot + 1, ihi, incr);

}


std::vector<int> sorted_indexes(const std::vector<double>& input, bool decreasing=false) {
    // Create a vector of indexes
    std::vector<int> indexes(input.size());
    std::iota(indexes.begin(), indexes.end(), 0); // Fill with 0, 1, 2, ..., input.size()-1

    // Sort the indexes based on the corresponding values in the input vector
    if (decreasing) {
      std::sort(indexes.begin(), indexes.end(), [&input](int i, int j) { return input[i] > input[j]; });
    } else {
      std::sort(indexes.begin(), indexes.end(), [&input](int i, int j) { return input[i] < input[j]; });
    }

    return indexes;
}

std::vector<int> sorted_indexes(const std::vector<int>& input, bool decreasing=false) {
    // Create a vector of indexes
    std::vector<int> indexes(input.size());
    std::iota(indexes.begin(), indexes.end(), 0); // Fill with 0, 1, 2, ..., input.size()-1

    // Sort the indexes based on the corresponding values in the input vector
    if (decreasing) {
      std::sort(indexes.begin(), indexes.end(), [&input](int i, int j) { return input[i] > input[j]; });
    } else {
      std::sort(indexes.begin(), indexes.end(), [&input](int i, int j) { return input[i] < input[j]; });
    }

    return indexes;
}



/**************************************************************/
/* Random sampling                                            */
/**************************************************************/

/*
 * Sample from integer vector.
 *
 * Sample of size n without replacement from vector x of length popsize
 * Result is returned in the first n elements of x
 */
void samplei_wr(int *x, int popsize, int n) {
  int i, r, temp;

  for (i = 0; i < n; i++) {

    r = i + floor((popsize - i) * runifC());
    //r = i + (int)((popsize - i - 1) * runifC()); //wrong

    temp = x[i];
    x[i] = x[r];
    x[r] = temp;

  }
}


/*
 * Sample from double vector.
 *
 * Sample of size n without replacement from vector x of length popsize
 * Result is returned in the first n elements of x
 */
void sampled_wr(double *x, int popsize, int n) {
  int i, r;
  double temp;

  for (i = 0; i < n; i++) {

    r = i + floor((popsize - i) * runifC());

    temp = x[i];
    x[i] = x[r];
    x[r] = temp;

  }
}


/*Random draw from a birth-death proposal. With probability pbirth take a random zero entry in model and set it to 1. With probability 1-pbirth take a random non-zero entry and set it to zero

 Input

 - model: input model, structured as a single-column sparse matrix

 - pbirth: birth probability

 Output

 - index: index of the entry that was born / killed. Returns -1 if a birth move was chosen but all entries were already 1, or a death move was chosen but all entries were already 0

 - birth: indicates whether the entry given by index was born / killed, e.g. index= 0, birth=false means that entry 0 in model was killed


*/
void rbirthdeath(int *index, bool *birth, arma::SpMat<short> *model, double pbirth) {

  bool found= false;
  int i, dnonzero= model->n_nonzero, d= model->n_rows;
  arma::SpMat<short>::iterator it;

  (*index)= -1; //return value if proposed update cannot be done
  (*birth)= (runifC() < pbirth);

  if (*birth) {  //propose birth

    if (dnonzero < d) {  //if dnonzero==d, all entries are already non-zero

      int idx_birth= runifdisc(1, d - dnonzero); //choose one of the d - dnonzero entries and set it to 1
      int dnonzero_slice, previous_row= -1;

      for (it= model->begin(); (it != model->end()) && (!found); ++it) {
        dnonzero_slice= it.row() - previous_row - 1;

        if (dnonzero_slice >= idx_birth) {
          (*index)= previous_row + idx_birth;
          found= true;
        } else {
          idx_birth -= dnonzero_slice;
          previous_row= it.row();
        }

      }

      if (!found) (*index)= previous_row + idx_birth;

    }

  } else {  //propose death

    if (dnonzero>0) {  //if all entries are already zero, none can be killed

      int idx_death= runifdisc(0, dnonzero - 1);

      for (it= model->begin(), i=0; (it != model->end()) && (!found); ++it, i++) {
        if (i == idx_death) {
          (*index)= (int) it.row();
          found= true;
        }
      }

    }

  }

}


//Probability mass function of a birth-death proposal
double dbirthdeath(arma::SpMat<short> *modelnew, arma::SpMat<short> *model, double pbirth, bool logscale=true) {

  int d= model->n_nonzero, dnew= modelnew->n_nonzero, nrows= model->n_rows;
  double ans;

  if (dnew == d+1) { //birth was proposed

    ans= pbirth / (double) (model->n_rows - d);

  } else if (dnew == d-1) { //death was proposed

    ans= (1.0 - pbirth) / (double) d;

  } else if ((dnew == d) && (d == nrows)) { //if birth was proposed, but all entries were non-zero

    ans= pbirth; 

  } else if ((dnew == d) && (d == 0)) { //if death was proposed, but all entries were zero

    ans= 1.0 - pbirth;

  } else {

    ans= 0;

  }

  if (logscale) ans= log(ans);

  return ans;

}



/*Random draw from a birth-death-swap proposal

  With probability pbirth take a random zero entry in model and set it to 1
  With probability pdeath take a random non-zero entry and set it to zero
  With probability 1-pbirth-pdeath take a random zero entry in model and set it to 1, and a random non-zero entry and set it to 0

 Input

 - model: input model, structured as a single-column sparse matrix
 - pbirth: birth probability
 - pdeath: death probability

 Output

 - index_birth: index of the entry that was born / killed. Returns -1 if a birth move was chosen but all entries were already 1
 - index_death: index of the dead entry that was borh. Returns -1 if a death move was chosen but all entries were already 0
 - movetype: 1 if a birth move was chosen, 2 if a death move was chosen, 3 if a swap move was chosen

*/
void rbirthdeathswap(int *index_birth, int *index_death, int *movetype, arma::SpMat<short> *model, double pbirth, double pdeath) {

  bool found;
  int i, dnonzero= model->n_nonzero, d= model->n_rows;
  double u;
  arma::SpMat<short>::iterator it;

  (*index_birth)= (*index_death)= -1; //return value if proposed update cannot be done
  u = runifC();
  if (u < pbirth) {

    (*movetype)= 1;
    if (dnonzero < d) {  //if dnonzero==d, all entries are already non-zero

      int idx_birth= runifdisc(1, d - dnonzero); //choose one of the d - dnonzero entries and set it to 1
      int dnonzero_slice, previous_row= -1;

      found= false;
      for (it= model->begin(); (it != model->end()) && (!found); ++it) {
        dnonzero_slice= it.row() - previous_row - 1;

        if (dnonzero_slice >= idx_birth) {
          (*index_birth)= previous_row + idx_birth;
          found= true;
        } else {
          idx_birth -= dnonzero_slice;
          previous_row= it.row();
        }

      }

      if (!found) (*index_birth)= previous_row + idx_birth;

    }

  } else if (u < pbirth + pdeath) {

    (*movetype)= 2;
    if (dnonzero>0) {  //if all entries are already zero, none can be killed

      int idx_death= runifdisc(0, dnonzero - 1);

      found= false;
      for (it= model->begin(), i=0; (it != model->end()) && (!found); ++it, i++) {
        if (i == idx_death) {
          (*index_death)= (int) it.row();
          found= true;
        }
      }

    }

  } else {

    (*movetype)= 3;
    if ((dnonzero < d) && (dnonzero >0)) { //if all entries are one or all are zero, a swap move is not possible

      //Choose entry to be born
      int idx_birth= runifdisc(1, d - dnonzero); //choose one of the d - dnonzero entries and set it to 1
      int dnonzero_slice, previous_row= -1;

      found= false;
      for (it= model->begin(); (it != model->end()) && (!found); ++it) {
        dnonzero_slice= it.row() - previous_row - 1;

        if (dnonzero_slice >= idx_birth) {
          (*index_birth)= previous_row + idx_birth;
          found= true;
        } else {
          idx_birth -= dnonzero_slice;
          previous_row= it.row();
        }

      }

      if (!found) (*index_birth)= previous_row + idx_birth;

      //Choose entry to kill
      int idx_death= runifdisc(0, dnonzero - 1);

      found= false;
      for (it= model->begin(), i=0; (it != model->end()) && (!found); ++it, i++) {
        if (i == idx_death) {
          (*index_death)= (int) it.row();
          found= true;
        }
      }

    }
  }

}


//Probability mass function of a birth-death-swap proposal
double dbirthdeathswap(arma::SpMat<short> *modelnew, arma::SpMat<short> *model, double pbirth, double pdeath, double pswap, bool logscale=true) {

  int d= model->n_nonzero, dnew= modelnew->n_nonzero, nrows= model->n_rows;
  double ans;

  if (dnew == d+1) { //birth was proposed

    ans= pbirth / (double) (model->n_rows - d);

  } else if (dnew == d-1) { //death was proposed

    ans= pdeath / (double) d;

  } else if ((dnew == d) && (d == nrows)) { //either birth or swap were proposed, but all entries were non-zero

    ans= pbirth + pswap;

  } else if ((dnew == d) && (d == 0)) { //either death or swap were proposed, but all entries were zero

    ans= pdeath + pswap;

  } else if (dnew == d) { //if swap move was proposed

    ans= pswap / (d * (model->n_rows - d));

  } else {

    ans= 0;

  }

  if (logscale) ans= log(ans);

  return ans;

}


/************************************************************************
                       RANDOM VARIATE GENERATION
************************************************************************/

/* call setall(is1,is2) */
void setseed(long is1,
             long is2)
{
    cstat_set = 1;
    setall(is1, is2);
}


double runifC(void)
{
    double x;

    if (cstat_set == 0) {
        setall(is1, is2);
        cstat_set = 1;
    }

    /* assign to double x for conversion */
    x = genunf(0.0, 1.0);
    return(x);
}


/* Density of a Unif(a,b) */
double dunifC(double x,
              double a,
              double b)
{
    return ((x > a) && (x < b)) ? (1.0 / (b - a)) : 0.0;
}


/* Returns integer value between min and max (both included) */
int runifdisc(int min,
              int max)
{
    return(min + (int)(runifC()*(max+1-min)));
}


/*
   Random deviates from a discrete distribution with values 0,1...nvals-1
   and probabilities probs[0],probs[1]...probs[nvals-1]
   Returns 0 with probability probs[0], 1 with probability probs[1] etc.

   Note: for efficiency, best to give probs in decreasing order (probs[0] >= probs[1] >= ...)
 */
int rdisc(const double *probs, int nvals) {
    int i;
    double u, pcum;

    u = runifC();
    pcum = probs[0];
    for (i = 1; (pcum < u) && (i < nvals); i++) {
        pcum += probs[i];
    }
    return(i-1);
}

/* Same as rdisc, but given cumulative probabilities (cdf) instead of cdf */
int rdisc_pcum(const double *cdf, int nvals) {
    int i;
    double u, pcum;

    u = runifC();
    pcum = cdf[0];
    for (i = 1; (pcum < u) && (i < nvals); i++) {
        pcum = cdf[i];
    }
    return(i-1);
}

double rbetaC(double alpha,
              double beta)
{
    double x;
    double y;

    x = gamdev(alpha);               /* X ~ gamma(alpha) */
    y = gamdev(beta);                /* Y ~ gamma(beta)  */

    return (x / (x+y));     /* X/(X+y) ~ ebta(alpha, beta) */
}


/* CDF of a Beta distribution */
double pbetaC(double x,
              double a,
              double b)
{
    double bt;
    double c;

    if (x < 0.0 || x > 1.0) {
        nrerror("pbetaC", "", "invalid probability");
        /*NOTREACHED*/
    }

    if (x == 0.0 || x == 1.0) {
        bt = 0.0;
    }
    else {
        /* Factors in front of the continued fraction */
        c = a + b;
        bt = exp(gamln(&c) - gamln(&a) - gamln(&b) + a*log(x) + b*log(1.0-x));
    }

    if (x < (a+1.0) / (a+b+2.0)) {
        /* Use continued fraction directly */
        return bt * betacf(a, b, x) / a;
    }
    else {
        /* Use continued fraction after symmetry transformation */
        return(1.0 - bt * betacf(b, a, 1.0-x) / b);
    }
}


/*
 * Draws from Dirichlet with parameter alpha.
 * The value is saved in w, and p is the dimensionality of the parameter
 */
void rdirichlet(double *w,
                const double *alpha,
                const int *p)
{
    double s = 0.0;
    double W = 1.0;
    double a;
    double b;
    int j;

    //assert(w != NULL);
    //assert(alpha != NULL);
    //assert(p != NULL);

    for (j = 0; j < *p; j++) {
        s += alpha[j];
    }
    b = s;
    for (j = 0; j < *p-1; j++) {
        a = alpha[j];
        b -= alpha[j];
        w[j] = rbetaC(a, b) * W;
        W -= w[j];
    }
    w[*p-1] = W;
    if (W < 0) {
        REprintf("rdirichlet: negative W generated\n");
        /* :TBD: - Should this be considered fatal? */
    }
}


/*
 * Evaluates Dirichlet density at w
 *     alpha: params
 *         p: dimensionality of w and alpha
 */
double ddirichlet(const double *w,
                  double *alpha,    /* maybe const */
                  const int *p)
{
    int i;
    double ans = 0.0;
    double sumalpha = 0.0;

    //assert(w != NULL);
    //assert(alpha != NULL);
    //assert(p != NULL);

    for (i = 0; i < *p; i++) {
        ans += (alpha[i] - 1) * log(w[i]) - gamln(alpha+i);
        sumalpha += alpha[i];
    }
    ans += gamln(&sumalpha);
    return(exp(ans));
}


double gamdev(double alpha)
{
    double value;

    value = gengam(1.0, alpha);
    return value;
}


/* *************************************************
   normal cdf and inv cdf
 ************************************************* */

/* Returns cdf of normal N(0,1) at x */
double pnormC(double y) {
  double cdf, surv;

    if (y < -20.0) {
      cdf= 2.753624e-89;
    }
    else if (y > 20.0) {
      cdf= 1 - 2.753624e-89;
    }
    else {
      cumnor(&y,&cdf,&surv);
    }

    return cdf;
}

/* Returns cdf of normal N(m,s^2) at x */
double pnormC(double y, double m, double s) {
  double cdf, p, mean, sd, bound, x, z;
    /* primitive type conversion */
    x = y;
    mean = m;
    sd = s;
    z = (x - mean) / sd;

    if (z < -20.0) {
      p= 2.753624e-89;
      //p = 2.86e-7F;
    }
    else if (z > 20.0) {
      p= 1 - 2.753624e-89;
      //p = 0.9999997F;
    }
    else {
        double q;
        int status;
        int which = 1;
        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    }

    cdf = p; /* another primitive type conversion */
    return cdf;
}


/* Approximate Normal cdf. 

  Combine Abrawomitz-Stegun's approx 26.2.16 for abs(y) < 3.447088 with 3rd order asymptotic expansion 26.2.12 for abs(y)>= 3.447088

  Max absolute error < 1.153e-05; Max relative error < 0.00623

*/
double apnorm(double y, bool logscale) {
  double ans;

  if (y <= -3.4470887) {
    
    double y2= y*y;
    ans= dnormC(y,1) - log(-y) + log(1.0 - 1.0/y2 + 3.0/(y2*y2));
    if (!logscale) ans= exp(ans);
      
  } else if ((y > -3.4470887) & (y <= 0)) {

    double  a1= 0.4361836, a2= -0.1201676, a3= 0.9372980, p= 0.33267, t= 1.0 / (1.0 + p*(-y)), t2= t*t;
    ans= dnormC(y,1) + log(a1*t + a2*t2 + a3*t2*t);
    if (!logscale) ans= exp(ans);
  
 } else if ((y > 0) & (y <= 3.4470887)) {

    double  a1= 0.4361836, a2= -0.1201676, a3= 0.9372980, p= 0.33267, t= 1.0 / (1.0 + p*y), t2=t*t;
    ans= dnormC(y,1) + log(a1*t + a2*t2 + a3*t2*t);
    if (logscale) { ans= log(1.0 - exp(ans)); } else { ans= 1.0 - exp(ans); }

  } else {
    
    double y2= y*y;
    ans= dnormC(y,1) - log(y) + log(1.0 - 1.0/y2 + 3.0/(y2*y2));
    if (logscale) { ans= log(1.0-exp(ans)); } else { ans= 1.0-exp(ans); }
      
  }
  
  return ans;
}


/* Approximate Normal cdf. Higher accuracy & computational cost than apnorm

  Combine Abrawomitz-Stegun's approx 26.2.17 for abs(z) < 4.056531 with 5th order asymptotic expansion 26.2.12 for For abs(z) > 4.056531

  Reference: Abramowitz & Stegun, p932, http://people.math.sfu.ca/~cbm/aands/frameindex.htm

  Max absolute error <=7.452e-08; Max relative error < 0.00051
*/
double apnorm2(double y, bool logscale) {
  double ans;

  if (y <= -4.056531) {
    
    double y2= y*y, y4= y2*y2, y6= y4*y2;
    ans= dnormC(y,1) - log(-y) + log(1.0 - 1.0/y2 + 3.0/y4 - 15.0/y6 + 105.0/(y6*y2));
    if (!logscale) ans= exp(ans);
      
  } else if ((y > -4.056531) & (y <= 0)) {

    double a1= 0.319381530, a2= -0.356563782, a3= 1.781477937, a4= -1.821255978, a5= 1.330274429, p= 0.2316419, t= 1.0 / (1.0 + p*(-y)), t2= t*t, t3= t2*t, t4=t2*t2;
    ans= dnormC(y,1) + log(a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t4*t);
    if (!logscale) ans= exp(ans);
  
 } else if ((y > 0) & (y <= 4.056531)) {

    double a1= 0.319381530, a2= -0.356563782, a3= 1.781477937, a4= -1.821255978, a5= 1.330274429, p= 0.2316419, t= 1.0 / (1.0 + p*y), t2= t*t, t3= t2*t, t4=t2*t2;
    ans= dnormC(y,1) + log(a1*t + a2*t2 + a3*t3 + a4*t4 + a5*t4*t);
    if (logscale) { ans= log(1.0 - exp(ans)); } else { ans= 1.0 - exp(ans); }

  } else {

    double y2= y*y, y4= y2*y2, y6= y4*y2;
    ans= dnormC(y,1) - log(y) + log(1.0 - 1.0/y2 + 3.0/y4 - 15.0/y6 + 105.0/(y6*y2));
    if (logscale) { ans= log(1.0-exp(ans)); } else { ans= 1.0-exp(ans); }
      
  }
  
  return ans;
}



/*
 * Density of univariate Normal(0,1) evaluated at y.
 * log==1 returns in log-scale.
 */
double dnormC(double y, int logscale) {
    //assert((logscale == 0) || (logscale == 1));

    if (logscale == 1) {
      return(-log(SQ_M_PI_2) - 0.5 * y * y);
    }
    else {
      return(exp(-0.5 * y * y) / SQ_M_PI_2);
    }
}


/*
 * Density of univariate Normal(m,s^2) evaluated at y.
 * log==1 returns in log-scale.
 */
double dnormC(double y,
              double m,
              double s,
              int logscale)
{
    //assert((logscale == 0) || (logscale == 1));

    if (logscale == 1) {
        return(-log(SQ_M_PI_2) - log(s) - 0.5 * (y - m) * (y - m) / (s * s));
    }
    else {
        return(exp(-0.5 * (y - m) * (y - m) / (s * s)) / (SQ_M_PI_2 * s));
    }
}


/* Joint density of y[0]...y[n-1] under Normal(m,s^2) */
double dnormC_jvec(const double *y,
                   int n,
                   double m,
                   double s,
                   int logscale)
{
    int i;
    double ans = 0.0;

    //assert(y != NULL);
    //assert((logscale == 0) || (logscale == 1));

    for (i = 0; i < n; i++) {
        ans += dnormC(y[i], m, s, 1);
    }

    return (logscale == 1) ? ans : exp(ans);
}


/*
 * Density of multivariate Normal evaluated at y[1]...y[p].
 * Input
 * - mu: mean
 * - cholsinv: Cholesky decomposition of the inverse covariance matrix.
 * - det: if logdet==false, det is the determinant of the inverse covariance matrix. if logdet==true, det is the log of the determinant
 * - If transpose==false then Sigma^{-1}= cholsinv * cholsinv', where cholsinv is lower-triangular
 *   If transpose== true then Sigma^{-1}= cholsinv' * cholsinv, where cholsinv is lower-triangular
 * - logscale: if logscale==1 the log of the density is returned
 *
 * Example:
 *   choldc_inv(s,p,cholsinv,posdef);
 *   det= choldc_det(cholsinv,p);
 *   dmvnormC(y,p,mu,cholsinv,det,0,true);
 *
 * Example:
 *   choldc(sinv,p,cholsinv,posdef);
 *   det= choldc_det(cholsinv,p);
 *   dmvnormC(y,p,mu,cholsinv,det,0,false);
 */
double dmvnormC(const double *y, int p, const double *mu, double **cholsinv, double det, bool transpose, int logscale, bool logdet=false) {
    int i;
    double *z, ans;

    z  = dvector(1, p);
    for (i = 1; i <= p; i++) { z[i] = y[i] - mu[i]; }
    ans= dmvnorm0(z, p, cholsinv, det, transpose, true, logdet);
    free_dvector(z, 1, p);

    return (logscale == 1) ? ans : exp(ans);
}

//same as dmvnorm for particular case mean=0
double dmvnorm0(const double *y, int p, double **cholsinv, double det, bool transpose, int logscale, bool logdet=false) {
  int i;
  double *z2, res=0, ans;
    z2 = dvector(1, p);
    if (transpose) {
      Ax(cholsinv, y, z2, 1, p, 1, p);     /* Find y' cholsinv' cholsinv y */
    } else {
      Atx(cholsinv, y, z2, 1, p, 1, p);    /* Find y' cholsinv cholsinv' y */
    }
    for (i = 1; i <= p; i++) res += z2[i] * z2[i];
    free_dvector(z2, 1, p);

    if (logdet) {
      ans = -p * log(SQ_M_PI_2) + 0.5 * det - 0.5 * res;
    } else {
      ans = -p * log(SQ_M_PI_2) + 0.5 * log(det) - 0.5 * res;
    }
    return (logscale == 1) ? ans : exp(ans);
}

//same as dmvnorm for particular case mean=0, and Cholesky decomp of Sinv given as a vector
double dmvnorm0(const double *y, int p, double *cholsinv, double det, int logscale, bool logdet=false) {
  /* cholsinv stores Cholesky decomp of Sinv as a vector, following column order (1st column, 2nd column, etc).

      Element (l,l) of chol(XtX) is stored into cholXtX[ll], where ll= (l-1)*n - (l-1)*(l-2)/2;

      For l>m, element (l,m) of chol(XtX) is stored into  cholXtX[mm + l - m]. Remaining elements are all zero
  */

  int i, j, ii;
  double *z2, res=0, ans;
  z2 = dvector(1, p);

  /* Find y' cholsinv cholsinv' y */
  for (i=1; i<= p; i++) {
    ii= (i-1)*p - (i-1)*(i-2)/2;
    for (j=i, z2[i]=0; j<=p; j++) { z2[i] += cholsinv[ii+j-i] * y[j]; } //z2[i] += cholsinv[j][i] * y[j]; (only for j<=i)
  }

  for (i = 1; i <= p; i++) res += z2[i] * z2[i];
  free_dvector(z2, 1, p);

  if (logdet) {
    ans = -p * log(SQ_M_PI_2) + 0.5 * det - 0.5 * res;
  } else {
    ans = -p * log(SQ_M_PI_2) + 0.5 * log(det) - 0.5 * res;
  }
  return (logscale == 1) ? ans : exp(ans);
}


//same as dmvnorm, but y is a 0-indexed n x p matrix and the density of each row of y is returned in ans[0],...,ans[n-1]
void dmvnormmatC(double *ans, const double *y, int n, int p, const double *mu, double **cholsinv, double det, bool transpose, int logscale) {
  int i, j;
  double *z, *z2, res = 0.0, ct;

  z  = dvector(1, p);
  z2 = dvector(1, p);

  ct= - ((double) p) * log(SQ_M_PI_2) + 0.5 * log(det);
  for (i=0; i<n; i++) {

    for (j = 1; j <= p; j++) { z[j] = y[i + (j-1)*n] - mu[j]; }
    if (transpose) {
      Ax(cholsinv, z, z2, 1, p, 1, p);     /* Find (y-mu)' * cholsinv' * cholsinv * (y-mu) */
    } else {
      Atx(cholsinv, z, z2, 1, p, 1, p);     /* Find (y-mu)' * cholsinv * cholsinv' * (y-mu) */
    }
    res= 0;
    for (j = 1; j <= p; j++) { res += z2[j] * z2[j]; }

    ans[i] = ct - 0.5 * res;
    if (!logscale) { ans[i]= exp(ans[i]); }
  }

  free_dvector(z, 1, p);
  free_dvector(z2, 1, p);
}

//same as dmvnormmat, but the transpose of y is provided instead
// - ty[0],...,ty[p-1] are the values for 1st row of y
// - ty[p],...,ty[2*p-1] are the values for the 2nd row of y
// - etc.
// Output ans[0],...,ans[n-1] contains the density of the rows in y
void dmvnormmat_transC(double *ans, const double *ty, int n, int p, const double *mu, double **cholsinv, double det, bool transpose, int logscale) {
  int i, j;
  double *z, *z2, res = 0.0;

  z  = dvector(1, p);
  z2 = dvector(1, p);

  for (i=0; i<n; i++) {

    for (j = 1; j <= p; j++) { z[j] = ty[j-1 + i*p] - mu[j]; }
    if (transpose) {
      Ax(cholsinv, z, z2, 1, p, 1, p);     /* Find (y-mu)' * cholsinv' * cholsinv * (y-mu) */
    } else {
      Atx(cholsinv, z, z2, 1, p, 1, p);     /* Find (y-mu)' * cholsinv * cholsinv' * (y-mu) */
    }
    for (j = 1; j <= p; j++) {
        res += z2[j] * z2[j];
    }

    ans[i] = -p * log(SQ_M_PI_2) + 0.5 * log(det) - 0.5 * res;
    if (logscale) { ans[i]= exp(ans[i]); }
  }

  free_dvector(z, 1, p);
  free_dvector(z2, 1, p);
}




/* Returns inv cdf of normal N(m,s^2) at p */
double qnormC(double cdf, double m, double s) {
    double y;

    if ((cdf < 0.0) | (cdf > 1.0)) {
        errorC("qnormC", "tried inverse cdf with p<0 or p>1", 1);
        /*NOTREACHED*/
    }

    /* par check */
    if (cdf <= 2.753624e-89) {
	y = -20.0*s + m;
    }
    else if (cdf >= 0.99999999999999989) {
        y = 8.209536 * s + m;
    }
    else {
        /* primitive type conversion */
        double p = cdf;
        double q = 1.0 - p;
        double mean = m;
        double sd = s;
        double bound;
        double x;
        int which = 2;
        int  status;

        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);

        y = x; /* another primitive type conversion */
    }

    return y;
}


/* Returns draw from binomial(n,p) */
int rbinomial(int n,
              double p)
{
    int i;
    int x = 0;

    for (i = 0; i < n; i++) {
        x += (runifC() < p) ? 1 : 0;
    }
    return x;
}


double dbinomial(int x,
                 int n,
                 double p,
                 int logscale)
{
    double ans;

    //assert((logscale == 0) || (logscale == 1));

    ans= lnchoose((double) n, (double) x) + (x+.0)*log(p) + (n-x+.0)*log(1-p);
    return (logscale == 1) ? ans : exp(ans);
}

double dnegbinomial(int x, double r, double p, int logscale) {
  double ans;
 //r: number of failures; p: success prob
  ans= lnchoose(x + r -1.0,(double) x) + (x+.0)*log(p) + r*log(1-p);
  if (logscale==1) return(ans); else return(exp(ans));
}



/*
 * Returns draw from multinomial with cell prob pr.
 *------------------------------------------
 * Value:
 *   x:   vector of indices indicating draws
 *         x in [0..ncells-1]
 *------------------------------------------
 * Input:
 *   ndraws: number of draws
 *   ncells: number of cells
 *   pr:     ncells vector of cell probs (not necessarily standardized)
 * Output:
 *   x:      ndraws vector of indices
 */
void rmultinomial(int ndraws, int ncells, const double *pr, int *x) {
  int i, j;
  double *cum_p;

    cum_p = dvector(0, ncells);

    cum_p[0] = pr[0];
    for (i = 1; i < ncells; i++) {
        cum_p[i] = cum_p[i-1] + pr[i];
    }
    for (j = 0; j < ndraws; j++) {
        double uj;

        uj = runifC() * cum_p[ncells-1];
        for (i = 0; ((uj > cum_p[i]) & (i < ncells)); i++);
        x[j] = i;
    }

    free_dvector(cum_p, 0, ncells);
}


//Single draw from a multinomial with probabilities pr (note: unlike the previous rmultinomial form, pr must add to 1)
int rmultinomial(const std::vector<double> *pr) {
  int i=0;
  double u= runifC(), pcum=pr->at(0);
  while (u > pcum) {
    i++;
    pcum += pr->at(i);
  }
  return i;
}

/*
 * Beta-binomial(alpha,beta) prior probability for a model including
 * k out of p variables.
 */
double bbPrior(int k,
               int p,
               double alpha,
               double beta,
               int logscale)
{
    double ans;

    //assert((logscale == 0) || (logscale == 1));

    ans = lnbeta(alpha + k, beta + p - k) - lnbeta(alpha, beta);
    return (logscale == 1) ? ans : exp(ans);
}


/*
 * Complexity prior probability for a model including
 * k out of p variables.
 */
double complexPrior(int k,
               int p,
               double priorc,
               int logscale)
{
  double priornorm, ans;

    //assert((logscale == 0) || (logscale == 1));
  priornorm= log(1.0 - 1.0/pow((double) p, priorc * ((double) p +1.0))) - log(1.0 - 1.0/pow((double) p, priorc));
  ans= lnbeta(1.0 + (double) k, 1.0 + (double) (p - k)) - (priorc* (double) k) * log((double) p) - priornorm;

  return (logscale == 1) ? ans : exp(ans);
}


/* Draw from univariate Normal(mu,s^2) */
double rnormC(double mu,
              double s)
{
    static bool iset = false;
    static double gset;
    double normdev;

    /* Is a deviate available from a previous invocation? */
    if (iset == false) {
        double fac;
        double rsq;
        double v1;
        double v2;

        do {
            /*
             * Pick two uniform numbers in the square extending
             * from -1 to +1 in each direction
             */
            v1 = 2.0 * runifC() - 1.0;
            v2 = 2.0 * runifC() - 1.0;
            /* See if they are in the unit circle */
            rsq = (v1 * v1) + (v2 * v2);
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        /* Make Box-Muller transformation to get two normal deviates */
        gset = v1 * fac;        /* Save this one for next invocation */
        iset = true;
        normdev = v2 * fac;
    }
    else {
        iset = false;
        normdev = gset;
    }
    return normdev * s + mu;
}


/*
 * Draw from a univariate truncated Normal giving truncation points.
 *    ltrunc - left truncation point
 *    rtrunc - right truncation point
 *    m      - mean
 *    s      - SD
 */
double rnorm_trunc(double ltrunc,
                   double rtrunc,
                   double m,
                   double s)
{
    double lprob;
    double rprob;

    lprob = pnormC(ltrunc, m, s);
    rprob = pnormC(rtrunc, m, s);
    return rnorm_trunc_prob(lprob, rprob, m, s);
}


/*
 * Draw from a univariate truncated Normal giving truncation probabilities.
 *    lprob - prob to the left of the truncation points
 *    rprob - prob to the right of the truncation points
 *    m     - mean
 *    s     - SD
 *
 * For example,
 *   lprob=.05, rprob=.99 means truncate the lower 5% and the upper 1%
 */
double rnorm_trunc_prob(double lprob,
                        double rprob,
                        double m,
                        double s)
{
    double u;

    if (lprob >= rprob) {
        nrerror("rnorm_trunc_prob",
                "",
                "left truncation probability larger than right truncation probability");
        /*NOTREACHED*/
    }
    /* Generate uniform between lprob, rprob */
    u = lprob + runifC() * (rprob-lprob);
    return qnormC(u, m, s);
}

//Obtain n draws from normal given multiple truncation points, i.e. N(th; m,s) prod_i I(ltrunc[i] < th < rtrunc[i])
//Important note: intervals are assumed to be disjoint and ordered.
// Input
// - n: length of output y
// - ltrunc, rtrunc: vectors [0..ntrunc-1] with start / end of truncation intervals.
// - ntrunc: number of truncatio intervals, i.e. length of vectors ltrunc & rtrunc
// - m: mean of underlying normal
// - s: sd
// Output
// - y: random draws
// - pdfy: joint log-density of y[0..n-1]
void rnorm_truncMult(double *y, double *pdfy, int *n, double *ltrunc, double *rtrunc, int ntrunc, double *m, double *s) {
  int i, j;
  double u, **p, *cump;
  //Find quantiles for truncation points
  p = dmatrix(0,ntrunc -1,0,1);
  cump = dvector(0,ntrunc);
  cump[0]= 0;
  for (i=0; i<ntrunc; i++) {
    p[i][0]= pnormC(ltrunc[i], *m, *s);
    p[i][1]= pnormC(rtrunc[i], *m, *s);
    cump[i+1]= max_xy(cump[i] + 1.0e-30, cump[i] + p[i][1] - p[i][0]);
  }
  //Generate random draws
  (*pdfy)= 0;
  for (i=0; i< *n; i++) {
    u= runifC() * cump[ntrunc];
    j= 0;
    while ((u > cump[j+1]) && (j< ntrunc -1)) j++;
    y[i]= qnormC(p[j][0] + u - cump[j], *m, *s);
    (*pdfy) += dnormC(y[i], *m, *s, 1) - log(cump[ntrunc]);
  }
  free_dmatrix(p,0,ntrunc -1,0,1);
  free_dvector(cump,0,ntrunc);
}

// [[Rcpp::export]]
SEXP rnorm_truncMultCI(SEXP n, SEXP ltrunc, SEXP rtrunc, SEXP m, SEXP s) {
  double pdfans;
  SEXP ans;
  PROTECT(ans= Rf_allocVector(REALSXP,INTEGER(n)[0]));
  rnorm_truncMult(REAL(ans),&pdfans,INTEGER(n),REAL(ltrunc),REAL(rtrunc),LENGTH(ltrunc),REAL(m),REAL(s));
  UNPROTECT(1);
  return ans;
}



//Rectangular Truncated Multivariate Normal draws (within==1 for truncation inside the hyper-rectanble, within==0 for outside)
// Input
// - n: number of draws
// - p: number of variables
// - mu: mean [1..p]
// - Sigma: covariance matrix [1..p][1..p]
// - lower: lower truncation points [1..p]
// - upper: upper truncation points [1..p]
// - within: if within==1 random draws are restricted so that lower[i] <= x[i] <= upper[i]. Else, x[i] <= lower[i] or x[i] >= upper[i]
// - method: if method==1 Gibbs sampling is used. Else independent proposal Metropolis-Hastings is used.
void rtmvnorm(double *ans, int n, int p, double *mu, double **Sigma, double *lower, double *upper, int within, int method) {
  bool posdef;
  int i, j, k;
  double **D, **K, *alpha, *ansortho, paccept;
  D= dmatrix(1,p,1,p); K= dmatrix(1,p,1,p);
  alpha= dvector(1,p); ansortho= dvector(0,n*p -1);
  choldc(Sigma, p, D, &posdef);
  choldc_inv(Sigma, p, K, &posdef);
  //Draws from orthogonal transformation
  Ax(K, mu, alpha, 1, p, 1, p); //alpha= K %*% mu
  if (method==1) {
    if (within==1) {
      rtmvnormWithin(ansortho, n, p, alpha, D, lower, upper);
    } else {
      rtmvnormOutside(ansortho, n, p, alpha, D, lower, upper);
    }
  } else {
    double det= choldc_det(K,p);
    rtmvnormMH(ansortho, &paccept, n, p, alpha, D, K, det, lower, upper, within);
  }
  //Map draws back to original variables, i.e. ans= ansortho %*% t(D)
  for (i=0; i<n; i++) {
    for (j=0; j<p; j++) {
      ans[i+j*n]= 0;
      for (k=0; k<p; k++) ans[i+j*n] += ansortho[i+k*n] * D[j+1][k+1];
    }
  }
  free_dmatrix(D,1,p,1,p); free_dmatrix(K,1,p,1,p);
  free_dvector(alpha,1,p); free_dvector(ansortho,0,n*p -1);
}

// [[Rcpp::export]]
SEXP rtmvnormCI(SEXP n, SEXP mu, SEXP Sigma, SEXP lower, SEXP upper, SEXP within, SEXP method) {
  int i, j, p= LENGTH(mu);
  double **S;
  SEXP ans;

  S= dmatrix(1,p,1,p);
  for (i=1; i<=p; i++) S[i][i]= REAL(Sigma)[(i-1)*p+i-1];
  for (i=1; i<=p; i++) {
    int ip= (i-1)*p;
    for (j=1; j<i; j++) S[i][j]= S[j][i]= REAL(Sigma)[ip+j-1];
  }

  PROTECT(ans= Rf_allocVector(REALSXP,INTEGER(n)[0]*p));
  rtmvnorm(REAL(ans), INTEGER(n)[0], p, REAL(mu)-1, S, REAL(lower)-1, REAL(upper)-1, INTEGER(within)[0], INTEGER(method)[0]);
  free_dmatrix(S,1,p,1,p);
  UNPROTECT(1);
  return ans;
}


//Metropolis-Hastings to obtain n multiv draws from z ~ N(alpha, 1) with D[,j] * z[j] restricted either within or outside interval (lower[j] - D[,-j] %*% z[-j], upper[j] - D[,-j] %*% z[-j])
//Input
// - n: number of draws
// - p: dimension of multivariate normal
// - alpha: mean vector 1..p
// - D: matrix [1..p][1..p] with Cholesky decomp of original covariance matrix of x= D %*% z, i.e. Cov(x)= D'D
// - K: matrix [1..p][1..p] with Cholesky decomp of inverse covariance matrix
// - det: determinant of D'D
// - lower: lower truncation points [1..p]
// - upper: upper truncation points [1..p]
// - within: if within==1 draws are restricted to be within interval. Else they are restricted to be outside interval.
//Output
// - ans: n \times p matrix with Gibbs draws, in column order
// - paccept: proportion of accepted proposals
void rtmvnormMH(double *ansortho, double *paccept, int n, int p, double *alpha, double **D, double **K, double det, double *lower, double *upper, int within) {
  int i, j, naccept;
  double *x0, propold, propnew, lold, lnew;
  x0= dvector(1,p);
  rtmvnormProp(x0, &propold, p, alpha, D, lower, upper, within);
  lold= 0;
  for (j=1; j<=p; j++) {
    lold -= 0.5 * (x0[j]-alpha[j]) * (x0[j]-alpha[j]);
    ansortho[n*(j-1)]= x0[j];
  }
  naccept= 1;
  for (i=1; i<n; i++) {
    rtmvnormProp(x0, &propnew, p, alpha, D, lower, upper, within);
    lnew= 0;
    for (j=1; j<=p; j++) lnew -= 0.5 * (x0[j]-alpha[j]) * (x0[j]-alpha[j]);
    if (runifC() <= exp(lnew - lold + propold - propnew)) {
      for (j=1; j<=p; j++) ansortho[i + n*(j-1)]= x0[j];
      lold= lnew; propold= propnew;
      naccept++;
    } else {
      for (j=1; j<=p; j++) ansortho[i + n*(j-1)]= ansortho[i-1 + n*(j-1)];
    }
  }
  (*paccept)= (naccept+.0)/(n+.0);
  free_dvector(x0,1,p);
}


//Gibbs sampling to obtain n multivariate draws from z ~ N(alpha, 1) with D[,j] * z[j] restricted inside interval (lower[j] - D[,-j] %*% z[-j] , upper[j] - D[,-j] %*% z[-j])
//Input
// - n: number of draws
// - p: dimension of multivariate normal
// - alpha: mean vector 1..p
// - D: lower-triangular matrix [1..p][1..p]
// - lower: lower truncation points [1..p]
// - upper: upper truncation points [1..p]
//Output
// - ans: n \times p matrix with Gibbs draws, in column order
void rtmvnormWithin(double *ans, int n, int p, double *alpha, double **D, double *lower, double *upper) {
  int i, j, k;
  double l, u, *Dj, *x0, lprop;
  //Initialize
  Dj= dvector(1,p);
  x0= dvector(1,p);
  rtmvnormProp(x0, &lprop, p, alpha, D, lower, upper, 1);
  for (j=1; j<=p; j++) ans[n*(j-1)]= x0[j];
  //Gibbs
  for (j=1; j<=p; j++) {
    Dj[j]= 0;
    for (k=1; k<=p; k++) Dj[j] += D[j][k] * ans[(k-1)*n];
  }
  i= 1;
  while (i<n) {
    for (j=1; j<=p; j++) {
      //Find truncation points for variable j
      for (k=1; k<=p; k++) Dj[k]= Dj[k] - D[k][j]*ans[i-1 +(j-1)*n];
      l= -1.0e20; u= 1.0e20;
      k= 1;
      while (k <= p) {
        if (D[k][j]>0) {
          l= max_xy(l,(lower[k]-Dj[k])/D[k][j]);
          u= min_xy(u,(upper[k]-Dj[k])/D[k][j]);
        } else if (D[k][j]<0) {
          u= min_xy(u,(lower[k]-Dj[k])/D[k][j]);
          l= max_xy(l,(upper[k]-Dj[k])/D[k][j]);
        }
        k++;
      }
      ans[i +(j-1)*n]= rnorm_trunc(l, u, alpha[j], 1.0);
      for (k=1; k<=p; k++) Dj[k] += ans[i+(j-1)*n] * D[k][j];
    }
    i++;
  }
  free_dvector(x0,1,p);
  free_dvector(Dj,1,p);
}


//Obtain n multivariate draws from z ~ N(alpha, 1) with D[,j] * z[j] restricted outside of interval (lower[j] - D[,-j] %*% z[-j] , upper[j] - D[,-j] %*% z[-j])
//(uses Gibbs sampling)
//Input
// - n: number of draws
// - p: dimension of multivariate normal
// - alpha: mean vector 1..p
// - D: lower-triangular matrix [1..p][1..p]
// - lower: lower truncation points [1..p]
// - upper: upper truncation points [1..p]
//Output
// - ans: n \times p matrix with Gibbs draws, in column order
void rtmvnormOutside(double *ans, int n, int p, double *alpha, double **D, double *lower, double *upper) {
  int i, j, k;
  double *Dj, lprop, *x0;
  //Initialize
  Dj= dvector(1,p);
  x0= dvector(1,p);
  rtmvnormProp(x0, &lprop, p, alpha, D, lower, upper, 0);
  for (j=1; j<=p; j++) ans[n*(j-1)]= x0[j];
  for (j=1; j<=p; j++) {
    Dj[j]= 0;
    for (k=1; k<=p; k++) Dj[j] += D[j][k] * x0[k];
  }
  //Gibbs
  i= 1;
  while (i<n) {
    rtmvnormOutside_Gibbs(x0, Dj, alpha, D, p, lower, upper);
    for (j=1; j<=p; j++) ans[i+n*(j-1)]= x0[j];
    i++;
  }
  free_dvector(Dj,1,p);
  free_dvector(x0,1,p);
}


//Perform 1 Gibbs iteration to sample from z ~ N(alpha, 1) with D[,j] * z[j] restricted outside of interval (lower[j] - D[,-j] %*% z[-j] , upper[j] - D[,-j] %*% z[-j])
//Input:
// - alpha: mean vector [1..p]
// - D: Cholesky decomposition of covariance matrix
// - p: dimensionality of Normal (number of variables)
//Input-Output
// - z: on input contains current draw, on output contains updated draw
// - Dj: on input current value of D %*% z, on output updated according to draw in z
void rtmvnormOutside_Gibbs(double *z, double *Dj, double *alpha, double **D, int p, double *lower, double *upper) {
  int j, k, nrestrict, one=1;
  double *l, *u, oned=1;
  l= dvector(1,p); u= dvector(1,p);
  for (j=1; j<=p; j++) {
    //Find truncation points for variable j
    for (k=1; k<=p; k++) Dj[k]= Dj[k] - D[k][j]*z[j];
    //for (k=1; k<=p; k++) Dj[k]= Dj[k] - D[k][j]*ans[i-1 +(j-1)*n];
    k=1; nrestrict= 0;
    while (k <= p) {
      if (D[k][j]>0) {
	nrestrict++;
	l[nrestrict]= (lower[k]-Dj[k])/D[k][j];
	u[nrestrict]= (upper[k]-Dj[k])/D[k][j];
      } else if (D[k][j]<0) {
	nrestrict++;
	u[nrestrict]= (lower[k]-Dj[k])/D[k][j];
	l[nrestrict]= (upper[k]-Dj[k])/D[k][j];
      }
      k++;
    }
    //Merge excluded regions & define inclusion regions
    if (nrestrict>0) {
      int *o, jj;
      double *lmod, *umod, lprop;
      o= ivector(1,nrestrict);
      lmod= dvector(1,nrestrict+1);
      umod= dvector(1,nrestrict+1);
      for (jj=1; jj<=nrestrict; jj++) o[jj]= jj;
      dindexsort(l, o, 1, nrestrict, 1); //sort indexes according to values in l
      jj=k=2;
      lmod[1]= l[o[1]]; umod[2]= u[o[1]];
      while (jj <= nrestrict) {
        if (u[o[jj]] > umod[k]) {
	  if (l[o[jj]] <= umod[k]) {
	    umod[k]= u[o[jj]];
	  } else {
	    lmod[k]= l[o[jj]];
	    k++;
	    umod[k]= u[o[jj]];
	  }
	}
	jj++;
      }
      umod[1]= -1.0e20; lmod[k]= 1.0e20;
      //Draw random variate
      rnorm_truncMult(z+j, &lprop, &one, umod+1, lmod+1, k, alpha+j, &oned); //note: umod is lower bound for inclusion intervals (lmod is upper bound)
      for (k=1; k<=p; k++) Dj[k] += z[j] * D[k][j];
      //rnorm_truncMult(ans+i+(j-1)*n, &lprop, &one, umod+1, lmod+1, k, alpha+j, &oned); //note: umod is lower bound for inclusion intervals (lmod is upper bound)
      //for (k=1; k<=p; k++) Dj[k] += ans[i+(j-1)*n] * D[k][j];
      free_ivector(o,1,nrestrict);
      free_dvector(lmod,1,nrestrict+1);
      free_dvector(umod,1,nrestrict+1);
    } else {
      z[j]= rnormC(alpha[j], 1);
    }
  }
  free_dvector(l,1,p); free_dvector(u,1,p);
}

//Draw from proposal approximating truncated multivariate normal (rectangular truncation)
// x ~ N(mu,Sigma) restricted to either I(lower[i] <= x[i] <= upper[i]) or I(x[i]<=lower[i]) + I(x[i]>=upper[i])
// Let D be the Cholesky decomp of Sigma, and alpha= D^{-1} %*% mu. rtmvnormProp samples from an approximation to z= D^{-1} %*% x ~ N(alpha, I) * I(lower[i] <= D %*% z <=upper[i])
// Input
// - p: number of variables, i.e. dim(x)
// - alpha: mean vector [1..p]
// - D: Cholesky decomposition of Sigma
// - lower: vector with lower truncation points [1..p]
// - upper: vector with upper truncation points [1..p]
// - within: set within==1 for lower[i] <= x[i] <= upper[i], and within==0 for x[i]<=lower[i] | x[i]>=upper[i]
// Output
// - z: vector [1..p] with draw from proposal distribution (guaranteed to have the same support as target distrib, i.e. to satisfy the constraints)
// - propdens: log-proposal density
void rtmvnormProp(double *z, double *propdens, int p, double *alpha, double **D, double *lower, double *upper, int within) {
  int j, k, one=1;
  double Dj, logp, oned=1.0;
  (*propdens)= 0;
  if (within==1) {
    double l, u;
    l= lower[1]/D[1][1]; u=  upper[1]/D[1][1];
    rnorm_truncMult(z+1, &logp, &one, &l, &u, 1, alpha+1, &oned);
    (*propdens) += logp;
    for (j=2; j<=p; j++) {
      Dj= 0;
      for (k=1; k<=(j-1); k++) Dj += D[j][k] * z[k];
      l= (lower[j]-Dj)/D[j][j]; u= (upper[j]-Dj)/D[j][j];
      rnorm_truncMult(z+j, &logp, &one, &l, &u, 1, alpha+j, &oned);
      (*propdens) += logp;
    }
  } else {
    double *l, *u;
    l= dvector(1,2); u= dvector(1,2);
    l[1]= -1.0e20; l[2]= upper[1]/D[1][1];
    u[1]= lower[1]/D[1][1]; u[2]= 1.0e20;
    rnorm_truncMult(z+1, &logp, &one, l+1, u+1, 2, alpha+1, &oned);
    (*propdens) += logp;
    for (j=2; j<=p; j++) {
      Dj= 0;
      for (k=1; k<=(j-1); k++) Dj += D[j][k] * z[k];
      l[2]= (upper[j]-Dj)/D[j][j]; u[1]= (lower[j]-Dj)/D[j][j];
      rnorm_truncMult(z+j, &logp, &one, l+1, u+1, 2, alpha+j, &oned);
      (*propdens) += logp;
    }
    free_dvector(l,1,2); free_dvector(u,1,2);
  }
}


//Multiv trunc Normal given constraint on product
void rtmvnormProd(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, double upper, int is_low_trunc, int is_up_trunc, int burnin) {
  if (is_low_trunc==1 && is_up_trunc==0) {
    rtmvnormProd_low(ans, n, p, mu, Sinv, k, lower, burnin);
  } else if (is_low_trunc==0 && is_up_trunc==1) {
    rtmvnormProd_up(ans, n, p, mu, Sinv, k, upper, burnin);
  } else if (is_low_trunc==1 && is_up_trunc==1) {
    rtmvnormProd_lowup(ans, n, p, mu, Sinv, k, lower, upper, burnin);
  } else {
    bool posdef;
    double *y, **S, **cholS;
    y= dvector(1,p); S= dmatrix(1,p,1,p); cholS= dmatrix(1,p,1,p);
    inv_posdef(Sinv,p,S,&posdef);
    choldc(S,p,cholS,&posdef);
    rmvnormC(y, p, mu, cholS);
    free_dvector(y,1,p); free_dmatrix(S,1,p,1,p); free_dmatrix(cholS,1,p,1,p);
  }
}

// [[Rcpp::export]]
SEXP rtmvnormProdCI(SEXP n, SEXP mu, SEXP Sigma, SEXP k, SEXP lower, SEXP upper, SEXP is_low_trunc, SEXP is_up_trunc, SEXP burnin) {
  bool posdef;
  int i,j, p=LENGTH(mu), nn=INTEGER(n)[0];
  double **S, **Sinv;

  S= dmatrix(1,p,1,p); Sinv= dmatrix(1,p,1,p);
  for (i=1; i<=p; i++) S[i][i]= REAL(Sigma)[(i-1)*p+i-1];
  for (i=1; i<=p; i++) {
    int ip= (i-1)*p;
    for (j=1; j<i; j++) S[i][j]= S[j][i]= REAL(Sigma)[ip+j-1];
  }
  inv_posdef(S, p, Sinv,&posdef);
  free_dmatrix(S,1,p,1,p);

  SEXP ans;
  PROTECT(ans= Rf_allocVector(REALSXP, nn*p));
  rtmvnormProd(REAL(ans), nn, p, REAL(mu)-1, Sinv, INTEGER(k)[0], REAL(lower)[0], REAL(upper)[0], INTEGER(is_low_trunc)[0], INTEGER(is_up_trunc)[0], INTEGER(burnin)[0]);

  free_dmatrix(Sinv,1,p,1,p);
  UNPROTECT(1);
  return(ans);
}


//Product-Truncated Multivariate Normal draws (via Gibbs), i.e. lower < prod abs(x)^k < upper
// Input
// - n: number of draws
// - p: number of variables
// - mu: vector [1..p] with mean
// - Sinv: precision matrix [1..p][1..p]
// - k: power for the truncation
// - lower: lower truncation point
// - upper: upper truncation point
// - burnin: number of burnin iterations
void rtmvnormProd_lowup(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, double upper, int burnin) {
  int i, j, kk, one=1;
  double *condvar, *condsd, *xcur, l, u, m, pdfy, *ltrunc, *rtrunc;
  //Pre-compute stuff
  condvar= dvector(1,p); condsd= dvector(1,p);
  xcur= dvector(1,p);
  ltrunc= dvector(1,2); rtrunc= dvector(1,2);
  for (i=1; i<=p; i++) { condvar[i]= 1/Sinv[i][i]; condsd[i]= sqrt(condvar[i]); }
  //Initialize
  l= pow(lower, 1.0/(k*p+.0));
  u= pow(upper, 1.0/(k*p+.0));
  for (i=1; i<=p; i++) {
    if (mu[i]>=l && mu[i]<=u) xcur[i]= mu[i]; else if (mu[i]<l) xcur[i]= l + .1*(u-l); else if (mu[i]>u) xcur[i]= u - .1*(u-l);
  }
  //Burn-in
  l= pow(lower,1.0/(k+.0)); u= pow(upper,1.0/(k+.0));
  for (i=1; i<=p; i++) {
    l= l/fabs(xcur[i]);
    u= u/fabs(xcur[i]);
  }
  for (i=0; i<burnin; i++) {
    for (j=1; j<=p; j++) {
      l= l*fabs(xcur[j]); u= u*fabs(xcur[j]);
      m= mu[j];
      for (kk=1; kk<j; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      for (kk=j+1; kk<=p; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      ltrunc[1]= -u; rtrunc[1]= -l;
      ltrunc[2]= l; rtrunc[2]= u;
      rnorm_truncMult(xcur+j, &pdfy, &one, ltrunc+1, rtrunc+1, 2, &m, condsd+j);
      l= l/fabs(xcur[j]); u= u/fabs(xcur[j]);
    }
  }
  //Gibbs
  for (i=0; i<n; i++) {
    for (j=1; j<=p; j++) {
      l= l*fabs(xcur[j]); u= u*fabs(xcur[j]);
      m= mu[j];
      for (k=1; k<j; k++) m -= Sinv[j][k] * (xcur[k]-mu[k]) * condvar[j];
      for (k=j+1; k<=p; k++) m -= Sinv[j][k] * (xcur[k]-mu[k]) * condvar[j];
      ltrunc[1]= -u; rtrunc[1]= -l;
      ltrunc[2]= l; rtrunc[2]= u;
      rnorm_truncMult(xcur+j, &pdfy, &one, ltrunc+1, rtrunc+1, 2, &m, condsd+j);
      l= l/fabs(xcur[j]); u= u/fabs(xcur[j]);
      ans[i+(j-1)*n]= xcur[j];
    }
  }
  free_dvector(condvar,1,p); free_dvector(condsd,1,p);
  free_dvector(xcur,1,p);
  free_dvector(ltrunc,1,2); free_dvector(rtrunc,1,2);
}

//Only lower truncation on product. Input/output parameters as in rtmvnormProd_lowup
void rtmvnormProd_low(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, int burnin) {
  int i, j, kk, one=1;
  double *condvar, *condsd, *xcur, l, m, *ltrunc, *rtrunc, pdfy;
  //Pre-compute stuff & initialize
  condvar= dvector(1,p); condsd= dvector(1,p);
  xcur= dvector(1,p);
  ltrunc= dvector(1,2); rtrunc= dvector(1,2);
  l= pow(lower, 1.0/(k*p+.0));
  for (i=1; i<=p; i++) {
    condvar[i]= 1/Sinv[i][i]; condsd[i]= sqrt(condvar[i]);
    if (mu[i]>l) xcur[i]= mu[i]; else xcur[i]= l + .1*condsd[i];
  }
  //Burn-in
  l= pow(lower,1.0/(k+.0));
  for (i=1; i<=p; i++) { l= l/fabs(xcur[i]); }
  for (i=0; i<burnin; i++) {
    for (j=1; j<=p; j++) {
      l= l*fabs(xcur[j]);
      m= mu[j];
      for (kk=1; kk<j; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      for (kk=j+1; kk<=p; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      ltrunc[1]= -1.0e20; rtrunc[1]= -l;
      ltrunc[2]= l; rtrunc[2]= 1.0e20;
      rnorm_truncMult(xcur+j, &pdfy, &one, ltrunc+1, rtrunc+1, 2, &m, condsd+j);
      l= l/fabs(xcur[j]);
    }
  }
  //Gibbs
  for (i=0; i<n; i++) {
    for (j=1; j<=p; j++) {
      l= l*fabs(xcur[j]);
      m= mu[j];
      for (kk=1; kk<j; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      for (kk=j+1; kk<=p; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      ltrunc[1]= -1.0e20; rtrunc[1]= -l;
      ltrunc[2]= l; rtrunc[2]= 1.0e20;
      rnorm_truncMult(xcur+j, &pdfy, &one, ltrunc+1, rtrunc+1, 2, &m, condsd+j);
      l= l/fabs(xcur[j]);
      ans[i+(j-1)*n]= xcur[j];
    }
  }
  free_dvector(condvar,1,p); free_dvector(condsd,1,p);
  free_dvector(xcur,1,p);
  free_dvector(ltrunc,1,2); free_dvector(rtrunc,1,2);
}

//Only upper truncation on product. Input/output parameters as in rtmvnormProd_lowup
void rtmvnormProd_up(double *ans, int n, int p, double *mu, double **Sinv, int k, double upper, int burnin) {
  int i, j, kk;
  double *condvar, *condsd, *xcur, u, m;
  //Pre-compute stuff
  condvar= dvector(1,p); condsd= dvector(1,p);
  xcur= dvector(1,p);
  u= pow(upper, 1.0/(k*p+.0));
  for (i=1; i<=p; i++) {
    condvar[i]= 1/Sinv[i][i]; condsd[i]= sqrt(condvar[i]);
    if (mu[i]<u) xcur[i]= mu[i]; else if (mu[i]>=u) xcur[i]= u - .1*condsd[i];
  }
  //Burn-in
  u= pow(upper,1.0/(k+.0));
  for (i=1; i<=p; i++) { u= u/fabs(xcur[i]); }
  for (i=0; i<burnin; i++) {
    for (j=1; j<=p; j++) {
      u= u*fabs(xcur[j]);
      m= mu[j];
      for (kk=1; kk<j; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      for (kk=j+1; kk<=p; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      xcur[j]= rnorm_trunc(-u, u, m, condsd[j]);
      u= u/fabs(xcur[j]);
    }
  }
  //Gibbs
  for (i=0; i<n; i++) {
    for (j=1; j<=p; j++) {
      u= u*fabs(xcur[j]);
      m= mu[j];
      for (kk=1; kk<j; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      for (kk=j+1; kk<=p; kk++) m -= Sinv[j][kk] * (xcur[kk]-mu[kk]) * condvar[j];
      xcur[j]= rnorm_trunc(-u, u, m, condsd[j]);
      u= u/fabs(xcur[j]);
      ans[i+(j-1)*n]= xcur[j];
    }
  }
  free_dvector(condvar,1,p); free_dvector(condsd,1,p);
  free_dvector(xcur,1,p);
}




/*
 * Draw from multivar Normal with n dimensions.
 * Result is stored in y[1..n]. mu is the location parameter, chols is the
 * Cholesky decomposition of the covariance matrix. That is, the covariance
 * is s and s=chols*chols'
 * Note: both y and mu should have length n, and s should be an n*n matrix.
 * The routine doesn't check it.
 *
 * Example:
 *   choldc(s,n,chols,posdef); //compute cholesky decomposition
 *   rmvnormC(y,n,mu,chols); //generate random variate
 */
void rmvnormC(double *y,
              int n,
              const double *mu,
              double **chols)
{
    int i;
    double *z;

    //assert(y != NULL);
    //assert(mu != NULL);
    //assert(chols != NULL);

    z = dvector(1, n);
    /* Generate n independent draws from a univariate Normal */
    for (i = 1; i <= n; i++) {
        z[i] = rnormC(0, 1);
    }
    /* Compute mu + chols*z */
    Ax_plus_y(chols, z, mu, y, 1, n);
    free_dvector(z, 1, n);
}


//draw from multivariate Normal for Armadillo input/output arguments. y and mu are taken as a column vector
// - y: output
// - mu: mean
// - chols: Cholesky decomposition of the covariance, i.e. cov = chols.t() * chols
void rmvnormC(arma::mat *y, arma::mat *mu, arma::mat *chols) {
  int i, n= y->n_rows;
  arma::mat z(n,1);

  for (i = 0; i < n; i++) z.at(i,0) = rnormC(0, 1);
  (*y)= chols->t() * z + (*mu);

}

//Mill's ratio (1-pnorm(z))/dnorm(z)
double millsnorm(double z) {
  return (1.0 - pnormC(z)) / dnormC(z,0);
}

//Inverse Mill's ratio dnorm(z)/pnorm(z)
double invmillsnorm(double z) {
  return dnormC(z,0) / pnormC(z);
}


/* Approximate Inverse Mill's ratio dnorm(z)/pnorm(z)

   For z< -1.756506 use optimized Laplace continuous fraction in Expression (5.3), Chu-In Charles Lee, Annals Inst Statis Math 1992, 44, 107-120

   For z>= -1.756506 use apnorm. This means
   - For z in (-1.756506, 3.447088) use approx 26.2.16
   - For z > 3.447088 use 3rd order asymptotic expansion 26.2.12

   Absolute error < 0.000185; Relative error < 0.000102
*/

double ainvmillsnorm(double z) {
  
  double ans;
  if (z< -1.756506) {
    ans= (-z + 1.0/(-z+2.0/(-z+3.0/(-z+4.0/(-z+5.0/(-z+11.5/(-z + 4.890096)))))));
  } else {
    ans= exp(dnormC(z,1) - apnorm(z,true));
  }
  return ans;
  
}


/* Approximate Inverse Mill's ratio dnorm(z)/pnorm(z). Higher accuracy and cost than ainvmillsnorm

   For z< -2.649083 use optimized Laplace continuous fraction in Expression (5.3), Chu-In Charles Lee, Annals Inst Statis Math 1992, 44, 107-120

   For z>= -2.649083 use apnorm2

   Absolute error < 0.000185, Relative error < 0.0000308
*/

double ainvmillsnorm2(double z) {
  double ans;
  if (z< -2.649083) {
    ans= (-z + 1.0/(-z+2.0/(-z+3.0/(-z+4.0/(-z+5.0/(-z+11.5/(-z + 4.890096)))))));
  } else {
    ans= exp(dnormC(z,1) - apnorm2(z,true));
  }
  return(ans);
}




 //Raw moment of N(m,sd) of order "order"
SEXP mnormCI(SEXP order, SEXP m, SEXP sd) {
  SEXP ans;

  ans= PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(ans)[0]= mnorm(REAL(order)[0], REAL(m)[0], REAL(sd)[0]);
  UNPROTECT(1);
  return ans;

}

/*
 * Raw moment of N(m,sd) of order "order"
 * Adapted from that in R package "actuar"
 */
double mnorm(double order, double m, double sd) {
    int n = (int) order;
    double ans;

    if (order == 0.0) {
        ans = 1.0;        /* Trivial case */
    }
    else if ((n % 2) == 1 && m == 0.0) {
        ans = 0.0;        /* Odd moments about 0 are equal to 0 */
    }
    else {
        int i;
        double temp1;
        double temp2;
        double res;

        res = 0.0;
        for (i = 0; i <= n/2; i++) {
            temp1 = i+1;
            temp2 = order-2*i+1;
            res += pow(sd, (double) (2*i)) * pow(m, (double) (n-2*i)) /
	      (pow(2.0, (double) i) * exp(gamln(&temp1)) * exp(gamln(&temp2)));
        }
        temp1 = order+1;
        ans = exp(gamln(&temp1)) * res;
    }
    return(ans);
}



/* Get the next n-tuple in lexicographical order with values 0, 1, ..., base-1.
   Return true if successful, false if unable to increment further. */
int GetNextTuple(int* tuple, int n, int base)
{
  int j = 0;
  while (j < n && tuple[j] == base-1)
    {
      tuple[j] = 0;
      j++;
    }
  if (j < n)
    tuple[j]++;
  return (j < n);
}

int BinomialCoefficient(int power, int nu) {
  int ans;
  /* power can only be 2 or 4 */
  /* 0 <= nu <= power */
  if (power == 2) {
    ans= 1 + nu%2;
  } else if (power == 4) {
      if (nu == 0 || nu == 4)
	ans= 1;
      if (nu == 1 || nu == 3)
	ans= 4;
      if (nu == 2)
	ans= 6;
  } else ans= 0; /* should never happen */
  return ans;
}

/* NB: The equation in the paper is off by one.
   The following code, taken from the paper's definition of kappa,
   acually computes kappa + 1. */
double one_plus_kappa
(
 double dof, /* degrees of freedom */
 int r/* summation index */
 )
{
  double product = 1.0;
  int i;

  if (r == 0)
    return 1.0;

  for (i = 1; i <= r; i++)
    product *= 0.5*dof - i;
  return pow(0.5*dof- 1.0, r) / product;
}


/*
    Expectation of prod (x_i)^(2*power), where x_i ~ multivariate T_dof(mu,sigma). Set dof= -1 for x_i ~ N(mu,sigma)
    mu[1..n]    = mean, a vector of length n
    sigma[1..n][1..n] = covariance, n by n matrix
    dof   = degrees of freedom, -1 signals normal case
    power = 2 or 4, exponent on components
    Note: This function was written by John D. Cook
*/

double mvtexpect(double* mu, double** sigma, int n, int power, double dof)
{
  double product = 1.0;
  double sum = 0.0;
  double covariance_term, mean_term, temp;
  int index_sum;

  int s = power*n;
  int half_s = s/2;
  int j, k, r, half_power;

  int* nu_index; /* indices nu_0 through nu_{n-1} */
  nu_index = ivector(0,n);

  for (r = 0; r <= half_s; r++)
    {

      for (j = 0; j < n; j++)
	nu_index[j] = 0;

      do
	{
	  product = 1.0;

	  index_sum = 0;
	  for (j = 0; j < n; j++)
	    index_sum += nu_index[j];
	  if (index_sum % 2)
	    product *= -1;

	  for (j = 0; j < n; j++)
	    product *= BinomialCoefficient(power, nu_index[j]);

	  if (dof > 0.0)
	    product *= one_plus_kappa(dof, r);
	  /* else normal case and multiplicative term is 1 */

	  covariance_term = 0.0;
	  half_power = power/2;
	  for (j = 0; j < n; j++)
	    {
	      temp = 0.0; /* double_sum accumulates the product sigma*h */
	      for (k = 0; k < n; k++)
		{
		  /* sigma is stored by columns, so sigma_{ij} = sigma[i + j*n] */
		  temp += sigma[j+1][k+1] * (half_power - nu_index[k]);
		}
	      covariance_term += (half_power - nu_index[j])*temp;
	    }
	  product *= pow(0.5*covariance_term, r);

	  mean_term = 0.0;
	  for (j = 0; j < n; j++)
	    mean_term += (half_power - nu_index[j])*mu[j+1];
	  product *= pow(mean_term, s - 2*r);
	  product /= exp(lfact(r) + lfact(s - 2*r));

	  sum += product;

	} while( GetNextTuple(nu_index, n, power+1) );
    }
  free_ivector(nu_index,0,n);
  return sum;
}

//save as mvtexpect but mu is vector [0..n-1] and sigma is 0-indexed vector instead of matrix
double mvtexpect_vec(double* mu, double* sigma, int n, int power, double dof)
{
  double product = 1.0;
  double sum = 0.0;
  double covariance_term, mean_term, temp;
  int index_sum;

  int s = power*n;
  int half_s = s/2;
  int j, k, r, half_power;

  int* nu_index; /* indices nu_0 through nu_{n-1} */
  nu_index = ivector(0,n);

  for (r = 0; r <= half_s; r++)
    {

      for (j = 0; j < n; j++)
	nu_index[j] = 0;

      do
	{
	  product = 1.0;

	  index_sum = 0;
	  for (j = 0; j < n; j++)
	    index_sum += nu_index[j];
	  if (index_sum % 2)
	    product *= -1;

	  for (j = 0; j < n; j++)
	    product *= BinomialCoefficient(power, nu_index[j]);

	  if (dof > 0.0)
	    product *= one_plus_kappa(dof, r);
	  /* else normal case and multiplicative term is 1 */

	  covariance_term = 0.0;
	  half_power = power/2;
	  for (j = 0; j < n; j++)
	    {
	      temp = 0.0; /* double_sum accumulates the product sigma*h */
	      for (k = 0; k < n; k++)
		{
		  /* sigma is stored by columns, so sigma_{ij} = sigma[i + j*n] */
		  temp += sigma[j + k*n]*(half_power - nu_index[k]);
		}
	      covariance_term += (half_power - nu_index[j])*temp;
	    }
	  product *= pow(0.5*covariance_term, r);

	  mean_term = 0.0;
	  for (j = 0; j < n; j++)
	    mean_term += (half_power - nu_index[j])*mu[j];
	  product *= pow(mean_term, s - 2*r);
	  product /= exp(lfact(r) + lfact(s - 2*r));

	  sum += product;

	} while( GetNextTuple(nu_index, n, power+1) );
    }
  free_ivector(nu_index,0,n);
  return sum;
}



/* Density of t with nu df, location mu and scale s^2 */
double dtC(double y,
           double mu,
           double s,
           int nu)
{
    double normk;
    double t1;
    double t2;

    t2 = 0.5 * nu;
    t1 = t2 + 0.5;
    normk = exp(gamln(&t1) - gamln(&t2)) / (sqrt(nu*M_PI) * s);
    return(normk * pow(1.0+(y-mu)*(y-mu) / (s*s*(nu+0.0)), -t1));
}


/*
 * Density of t mixtures with ncomp components,
 * i.e. sum T_nu(y;mu[i],s[i]^2) * probs[i]
 */
double dtmixC(double y,
              const double *mu,
              const double *s,
              const double *probs,
              int nu,
              int ncomp,
              int logscale)
{
    int i;
    double ans = 0.0;

    //assert(mu != NULL);
    //assert(s != NULL);
    //assert(probs != NULL);
    //assert((logscale == 0) || (logscale == 1));

    for (i = 0; i < ncomp; i++) {
        ans += dtC(y, mu[i], s[i], nu) * probs[i];
    }
    return (logscale == 1) ? log(ans) : ans;
}


/*
 * Density of multivariate T with nu df and n dimensions.
 *    mu: location
 *    cholsinv, det: cholesky decomp and determinant of inverse cov matrix
 *    logscale: set to 1 to return density in log-scale
 *
 * Example:
 *     choldc_inv(s, n, cholsinv,posdef);
 *     det = choldc_det(cholsinv, n);
 *     dmvtC(y, n, mu, cholsinv, det, nu, 0);
 */
double dmvtC(const double *y,
             int n,
             const double *mu,
             double **cholsinv,
             double det,
             int nu,
             int logscale)
{
    double res = 0.0;
    double t1;
    double t2;
    double normk;

    //assert(y != NULL);
    //assert(mu != NULL);
    //assert(cholsinv != NULL);
    //assert((logscale == 0) || (logscale == 1));

    /* Find (y-mu)' * cholsinv' * cholsinv * (y-mu) */
    {
        int i;
        double *z;
        double *z2;

        z  = dvector(1, n);
        z2 = dvector(1, n);
        for (i = 1; i <= n; i++) {
            z[i] = y[i] - mu[i];
        }
        Ax(cholsinv, z, z2, 1, n, 1, n);
        for (i = 1; i <= n; i++) {
            res += z2[i] * z2[i];
        }
        free_dvector(z, 1, n);
        free_dvector(z2, 1, n);
    }

    t2 = 0.5 * nu;
    t1 = t2 + 0.5*(n+0.0);
    normk = gamln(&t1) -
            gamln(&t2) -
            0.5*(n+0.0)*(log(nu+0.0)+log(M_PI))+0.5*log(det);

    if (logscale == 1)
        return(normk -t1*log(1+res/(nu+0.0)));
    else
        return(exp(normk) * pow(1.0+res/(nu+0.0), -t1));
}


/* Draw from a univariate standard t with nu degrees of freedom */
double rtC(int nu)
{
    double x;
    double z;

    z = rnormC(0, 1);
    x = gengam(0.5, nu / 2.0);  /* Draw from chi-square with nu DoF */
    return(z * sqrt(nu / x));
}


/* Draw from mixture of t_nu(mu[i],s[i]^2) */
double rtmixC(const double *mu,
              const double *s,
              const double *probs,
              int nu,
              int ncomp)
{
    int comp;

    //assert(mu != NULL);
    //assert(s != NULL);
    //assert(probs != NULL);

    comp = rdisc(probs, ncomp);
    return(mu[comp] + s[comp]*rtC(nu));
}


/*
 * Draw from a univariate truncated t with nu degrees of freedom giving
 * truncation points.
 *     nu     - degrees of freedom
 *     ltrunc - left truncation point
 *     rtrunc - right truncation point
 */
double rt_trunc(int nu,
                double ltrunc,
                double rtrunc)
{
    double lprob;
    double rprob;

    lprob = ptC(ltrunc, nu);
    rprob = ptC(rtrunc, nu);
    return rt_trunc_prob(nu, lprob, rprob);
}


/*
 * Draw from a univariate truncated t with nu degrees of freedom giving
 * truncation probabilities.
 *    nu           - degrees of freedom
 *    lprob, rprob - prob to the left of the truncation points
 *
 * e.g. lprob=.05, rprob=.99 means truncate the lower 5% and the upper 1%
 */
double rt_trunc_prob(int nu,
                     double lprob,
                     double rprob)
{
    double u;

    if (lprob >= rprob) {
        nrerror("rt_trunc_prob",
                "",
                "left truncation probability larger than right truncation probability");
        /*NOTREACHED*/
    }
    /* Generate uniform between lprob, rprob */
    u = lprob + runifC() * (rprob-lprob);
    return(qtC(u, nu));
}


/*
 * Find quantiles of a t-distribution with nu degrees of freedom.
 *    p: probability
 *    nu: degrees of freedom
 *    lower_tail==1 means p indicates left tail area
 *
 * Algorithm 396: Student's t-quantiles by G.W. Hill CACM 13(10), 619-620, October 1970
 *
 * @author Sundar Dorai-Raj
 * See the GNU General Public License for more details at http://www.gnu.org
 */
double qtC(double p,
           int nu)
{
  int neg;
  double ndf, eps, P, q, prob, a, b, c, d, y, x;

  ndf= nu + 0.0;
  if(p<=0 || p>=1 || ndf<1) return(-1);
  eps=1e-12;
  if(p > 0.5) {
    neg = 0;
    P = 2.0 * (1.0 - p);
  }
  else {
    neg = 1;
    P = 2.0 * p;
  }

  if(fabs(ndf - 2.0) < eps) {   /* df ~= 2 */
    q=sqrt(2.0 / (P * (2.0 - P)) - 2.0);
  }
  else if (ndf < 1 + eps) {   /* df ~= 1 */
    prob = P * M_PI_2;
    q = cos(prob)/sin(prob);
  }
  else {      /*-- usual case;  including, e.g.,  df = 1.1 */
    a = 1.0 / (ndf - 0.5);
    b = 48.0 / (a * a);
    c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
    d = ((94.5 / (b + c) - 3.0) / b + 1.0) * sqrt(a * M_PI_2) * ndf;
    y = pow(d * P, 2.0 / ndf);
    if (y > 0.05 + a) {  /* Asymptotic inverse expansion about normal */
      x = qnormC(0.5*P, 0.0, 1.0);
      y = x * x;
      if (ndf < 5)
	c += 0.3 * (ndf - 4.5) * (x + 0.6);
      c = (((0.05 * d * x - 5.0) * x - 7.0) * x - 2.0) * x + b + c;
      y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
      y = a * y * y;
      if (y > 0.002) y = exp(y) - 1.0;
      else { /* Taylor of    e^y -1 : */
	y = (0.5 * y + 1.0) * y;
      }
    } else {
      y = ((1.0 / (((ndf + 6.0) / (ndf * y) - 0.089 * d - 0.822) * (ndf + 2.0) * 3.0) + 0.5 / (ndf + 4.0))* y - 1.0) * (ndf + 1.0) / (ndf + 2.0) + 1.0 / y;
    }
    q = sqrt(ndf * y);
  }
  if (neg==1) q = -q;
  return(q);
}


/* CDF of a t-Student distribution */
double ptC(double x,
           int nu)
{
    if (x > 0) {
        return(1 - 0.5 * pbetaC((nu+0.0) / (x*x+nu), 0.5*nu, 0.5));
    }
    else if (x < 0) {
        return(0.5 * pbetaC((nu+0.0) / (x*x+nu), 0.5*nu, 0.5));
    }
    else {
        return(0.5);
    }
}


/*
 * Draw from multivar T with n dimensions and nu degrees of freedom
 * Result is stored in y[1..n].
 *     mu is the location parameter
 *     chols is the Cholesky decomposition of the covariance matrix.
 * That is, the covariance is s*nu/(nu-2) and s=chols*chols'
 * and nu are the degrees of freedom
 * Note: both y and mu should have length n, and s should be an n*n matrix.
 * The routine doesn't check it.
 *
 * Example:
 *   choldc(s,n,chols,posdef); //compute cholesky decomposition
 *   rmvtC(y,n,mu,chols,nu); //generate random variate
 */
void rmvtC(double *y,
           int n,
           const double *mu,
           double **chols,
           int nu)
{
    int i;
    double x;
    double *z;

    //assert(y != NULL);
    //assert(mu != NULL);
    //assert(chols != NULL);

    /* Draw from chi-square with nu degrees of freedom */
    x = sqrt(nu / gengam(0.5, nu / 2.0));

    /* Multiple n indep normal draws by the common chi-square */
    z = dvector(1, n);
    for (i = 1; i <= n; i++) {
        z[i] = x * rnormC(0, 1);
    }
    /* Compute mu + chols*z */
    Ax_plus_y(chols, z, mu, y, 1, n);

    free_dvector(z, 1, n);
}

/* Generate from a chi-square with df degrees of freedom */
double rchisqC(int df) {
  int i;
  double ans= 0;
  for (i=1; i<=df; i++) ans += pow(rnormC(0,1),2.0);
  return(ans);
}

/* Generate from p x p dimensional Wishart(df,S)
   Input
   - df: degrees of freedom
   - cholS: Cholesky decomposition of the scale matrix S
   - p: number of rows and columns of S
   - returnChol: if true the Cholesky decomposition of the sample matrix is returned (faster)

   Output
   - ans: if returnChol==false, a draw Sigma ~ Wishart(df,S). If returnChol=true, chol(Sigma) is returned instead (Sigma= chol(Sigma) %*% t(chol(Sigma)))

   Example to sample from Wishart(nu,S)

   choldc(S, p, cholS, posdef);
   rwishartC(ans, df, cholS, p, false);

   Example to sample from Inverse Wishart(nu,S)

   inv_posdef_upper(S, p, Sinv, &posdef);
   choldc(Sinv, p, cholSinv, &posdef);
   rwishartC(Sigmainv, df, cholSinv, p, false);
   inv_posdef_upper(Sigmainv, p, Sigma, &posdef);
 */
void rwishartC(double **ans, int df, double **cholS, int p, bool returnChol) {
  int i,j;
  double **Z, **cholans;
  Z= dmatrix(1,p,1,p);
  if (returnChol) { cholans= ans; } else { cholans= dmatrix(1,p,1,p); }
  for (i=1; i<=p; i++) {
    Z[i][i]= sqrt(rchisqC(df-p+i));
    if (p>1) {
      for (j=1; j<i; j++) Z[i][j]= rnormC(0,1);
      for (j=i+1; j<=p; j++) Z[i][j]= 0;
    }
  }
  AB(cholS,1,p,1,p,Z,1,p,1,p,cholans); //cholans= cholS %*% Z
  if (!returnChol) {
    ABt(cholans,1,p,1,p,cholans,1,p,1,p,ans); //ans= cholans %*% t(cholans)
    free_dmatrix(cholans,1,p,1,p);
  }
  free_dmatrix(Z,1,p,1,p);
}



/*
 * Generate from a Gamma(a,b)
 *   (a is shape; b location; mean= a/b)
 */
double rgammaC(double a,
               double b)
{
    return(gengam(b, a));
}


/*
 * Density of a Gamma(a,b)
 *   (a is shape; b location; mean= a/b)
 */
double dgammaC(double x,
               double a,
               double b)
{
    double ans = 0.0;

    if (x != 0.0) {
        ans = exp(a * log(b) - gamln(&a) + (a - 1) * log(x) - x * b);
    }
    else if (a == 1.0) {
        ans = b;
    }
    return ans;
}


/*
 * Density of an Inverse Gamma(a,b)
 *   (a: shape; b: location; mean of 1/x= a/b)
 */
double dinvgammaC(double x,
                  double a,
                  double b,
                  int logscale)
{
    double ans = 0.0;

    if (x != 0.0) {
        ans = a * log(b) - gamln(&a) - (a + 1) * log(x) - b / x;
        if (logscale==0) ans= exp(ans);
    }
    return ans;
}


/************************************************************************
                       NON-LOCAL PRIOR DENSITIES
************************************************************************/

/*
 * Normal MOM prior.
 * Density is proportional to (y-m)^(2*r) N(y;m,tau*phi), where tau*phi
 * is the variance.
 */
double dmom(double y, double m, double tau, double phi, int r, int logscale) {
    double normct[] = {
         0,
         1.098612,
         2.70805,
         4.65396,
         6.851185,
         9.24908,
        11.81403,
        14.52208,
        17.35529,
        20.29973
    };
    double ans;

    //assert((logscale == 0) || (logscale == 1));

    ans = r * log((y - m) * (y - m) / (tau * phi)) +
          dnormC(y, m, sqrt(tau * phi), 1) - normct[r-1];

    return (logscale == 1) ? ans : exp(ans);
}

//Multivariate MOM prior
double dmomvec(double *y, int n, double m, double tau, double phi, int r, int logscale) {
  int i;
  double ans=0;
  for (i=0; i<n; i++) { ans+= dmom(y[i],m,tau,phi,r,1); }
  return (logscale == 1) ? ans : exp(ans);
}

//Univariate iMOM prior
double dimom(double y, double m, double tau, double phi, int logscale) {
  double y2, ans;
  y2= (y-m)*(y-m);
  ans= .5*(log(tau)+log(phi)) - .5*LOG_M_PI - log(y2) - tau*phi/y2;
  if (logscale==0) ans= exp(ans);
  return(ans);
}

//Multivariate iMOM prior
double dimomvec(double *y, int n, double m, double tau, double phi, int logscale) {
  int i;
  double ans=0;
  for (i=0; i<n; i++) { ans+= dimom(y[i],m,tau,phi,logscale); }
  return (logscale == 1) ? ans : exp(ans);
}

//Univariate eMOM prior
double demom(double y, double tau, double phi, int logscale) {
  double pen, ans;
  pen= -tau * phi / (y*y);
  ans= pen + dnormC(y,0,sqrt(tau*phi),1) + sqrt(2.0);
  if (logscale==0) ans= exp(ans);
  return(ans);
}

//Multivariate eMOM prior
double demomvec(double *y, int n, double tau, double phi, int logscale) {
  int i;
  double ans=0;
  for (i=0; i<n; i++) { ans+= demom(y[i],tau,phi,logscale); }
  return (logscale == 1) ? ans : exp(ans);
}


/************************************************************************
               NON-LOCAL PRIOR DENSITY DERIVATIVES
************************************************************************/

//Gradient of log-pMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
void dmomgrad(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i;
  for (i=1; i<=(*n); i++) {
    ans[i]= 2.0/th[i] - th[i]/(exp(*logphi)*(*tau));
  }
}


//Same but univariate
double dmomgraduniv(double *th, double *logphi, double *tau) {
  return (2.0/(*th) - (*th)/(exp(*logphi)*(*tau)));
}


//Hessian of log-pMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
//Output ans is a vector, as non-diagonal elements are 0
void dmomhess(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i;
  for (i=1; i<=(*n); i++) { ans[i]= -2.0/pow(th[i],2) - 1.0/(exp(*logphi)*(*tau)); }
}

//same but univariate
double dmomhessuniv(double *th, double *logphi, double *tau) {
  return (-2.0/pow(*th,2) - 1.0/(exp(*logphi)*(*tau)));
}


//Gradient of log-pMOM(th;0,phi*tau) + log-IG(phi,alpha/2,lambda/2) wrt (th, logphi) where logphi=log(phi)
//Output: ans is vector [1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
void dmomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, p=(*n)-1;
  double sumth2=0;
  if (p>0) {  //th has some elements
    for (i=1; i<=p; i++) {
      ans[i]= 2.0/th[i] - th[i]/(exp(*logphi)*(*tau));  //same as in dmomgrad
      sumth2 += (th[i]*th[i]);
    }
    ans[*n]= -1.5*p - 0.5*(*alpha) -1.0 + 0.5*(sumth2/(*tau) + (*lambda)) * exp(-(*logphi));
    } else {  //th is empty
    ans[1]= -0.5*(*alpha) -1.0 + 0.5*(*lambda)*exp(-(*logphi));
  }
}

//Hessian of log-pMOM(0,phi*tau) + log-IG(phi,alpha,lambda) wrt (th, logphi) where logphi=log(phi)
//Output: ans is matrix [1..n][1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
//Input:
// - th[1..n-1]: vector with th values
// - logphi: log(phi)
// - tau: prior dispersion parameter
// - alpha, lambda: prior parameters for phi
void dmomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, j, p=(*n)-1;
  double sumth2=0;
  for (i=1; i<=p; i++) {
    for (j=1; j<=i; j++) { ans[i][j]= ans[j][i]=0; }
    ans[i][i]= -2.0/pow(th[i],2.0) - 1.0/(exp(*logphi)*(*tau)); //same as in dmomhess
    sumth2 += (th[i]*th[i]);
    for (j=i+1; j<p; j++) { ans[i][j]= ans[j][i]=0; }
    ans[i][*n]= ans[*n][i]= th[i] / (exp(*logphi)*(*tau));
  }
  ans[*n][*n]= -0.5 * exp(-(*logphi)) * (sumth2/(*tau)+(*lambda));
}


//Gradient of log-piMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
void dimomgrad(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i;
  for (i=1; i<=(*n); i++) { ans[i]= 2.0 * (*tau) * exp(*logphi) / (th[i]*th[i]*th[i]) - 2.0 / th[i]; }
}

//same but univariate
double dimomgraduniv(double *th, double *logphi, double *tau) {
  return (2.0 * (*tau) * exp(*logphi) / pow(*th,3.0) - 2.0 / (*th));
}

//Hessian of log-piMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
//Output ans is a vector, as non-diagonal elements are 0
void dimomhess(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i; double th2;
  for (i=1; i<=(*n); i++) {
    th2= pow(th[i],2.0);
    ans[i]= -6.0 * (*tau) * exp(*logphi) / pow(th2,2.0) + 2.0/th2;
  }
}

double dimomhessuniv(double *th, double *logphi, double *tau) {
  double th2;
  th2= pow(*th,2.0);
  return(-6.0 * (*tau) * exp(*logphi) / pow(th2,2.0) + 2.0/th2);
}


//Gradient of log-piMOM(th;0,phi*tau) + log-IG(th;phi,alpha/2,lambda/2) wrt (th, logphi) where logphi=log(phi)
//Output: ans is vector [1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
void dimomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, p=(*n)-1;
  double th2, suminvth2=0;
  if (p>0) {  //th is non-empty
    for (i=1; i<=p; i++) {
      th2= th[i]*th[i];
      ans[i]= 2.0 * (*tau) * exp(*logphi) / (th2*th[i]) - 2.0 / th[i]; //same as dimomgrad
      suminvth2+= 1.0/th2;
    }
    ans[*n]= 0.5*p - 0.5*(*alpha) -1.0 + 0.5*(*lambda)*exp(-(*logphi)) - exp(*logphi)*(*tau) * suminvth2;
  } else {  //th is empty
    ans[1]= -0.5*(*alpha) -1.0 + 0.5*(*lambda)*exp(-(*logphi));
  }
}

//Hessian of log-piMOM(0,phi*tau) + log-IG(phi,alpha,lambda) wrt (th, logphi) where logphi=log(phi)
//Output: ans is matrix [1..n][1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
void dimomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, j, p=(*n)-1;
  double th2, suminvth2=0;
  for (i=1; i<=p; i++) {
    for (j=1; j<i; j++) { ans[i][j]= ans[j][i]=0; }
    th2= th[i]*th[i];
    suminvth2+= 1.0/th2;
    ans[i][i]= -6.0 * (*tau) * exp(*logphi) / (th2 * th2) + 2.0/th2; //same as in dimomhess
    for (j=i+1; j<=p; j++) { ans[i][j]= ans[j][i]=0; }
    ans[i][*n]= ans[*n][i]= 2.0 * (*tau) * exp(*logphi) / (th2 * th[i]);
  }
  ans[*n][*n]= -0.5*exp(-(*logphi))*(*lambda) - (*tau) * exp(*logphi) * suminvth2;
}

//Gradient of log-peMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
void demomgrad(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i;
  for (i=1; i<=(*n); i++) { ans[i]= 2.0 * (*tau) * exp(*logphi) / pow(th[i],3.0) - th[i]*exp(-(*logphi))/(*tau); }
}

//same but univariate
double demomgraduniv(double *th, double *logphi, double *tau) {
  return(2.0 * (*tau) * exp(*logphi) / pow(*th,3.0) - (*th)*exp(-(*logphi))/(*tau));
}

//Hessian of log-peMOM(th;0,phi*tau) density wrt th. Note: n=dim(th)
//Output ans is a vector, as non-diagonal elements are 0
void demomhess(double *ans, int *n, double *th, double *logphi, double *tau) {
  int i;
  for (i=1; i<=(*n); i++) {
    ans[i]= -6.0 * (*tau) * exp(*logphi)/pow(th[i],4.0) - exp(-(*logphi))/(*tau);
  }
}

//same but univariate
double demomhessuniv(double *th, double *logphi, double *tau) {
  return(-6.0 * (*tau) * exp(*logphi)/pow(*th,4.0) - exp(-(*logphi))/(*tau));
}


//Gradient of log-peMOM(th;0,phi*tau) + log-IG(th;phi,alpha/2,lambda/2) wrt (th, logphi) where logphi=log(phi)
//Output: ans is vector [1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
void demomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, p=(*n)-1;
  double th2, suminvth2=0, sumth2=0;
  if (p>0) {  //th is non-empty
    for (i=1; i<=p; i++) {
      th2= th[i]*th[i];
      sumth2+= th2;
      suminvth2+= 1.0/th2;
      ans[i]= 2.0 * (*tau) * exp(*logphi) /(th2*th[i]) - th[i]*exp(-(*logphi))/(*tau);
    }
    ans[*n]= -0.5*p - 0.5*(*alpha) -1.0 + 0.5*(sumth2/(*tau) + (*lambda)) * exp(-(*logphi)) - exp(*logphi)*(*tau)*suminvth2;
  } else {  //th is empty
    ans[1]= -0.5*(*alpha) -1.0 + 0.5*(*lambda) * exp(-(*logphi));
  }
}

//Hessian of log-peMOM(0,phi*tau) + log-IG(phi,alpha,lambda) wrt (th, logphi) where logphi=log(phi)
//Output: ans is matrix [1..n][1..n] where n=dim(th)+1, i.e. n==1 indicates dim(th)=0
void demomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda) {
  int i, j, p=(*n)-1;
  double th2, sumth2=0, suminvth2=0;
  for (i=1; i<=p; i++) {
    for (j=1; j<i; j++) { ans[i][j]= ans[j][i]=0; }
    th2= th[i]*th[i];
    sumth2+= th2;
    suminvth2+= 1.0/th2;
    ans[i][i]= -6.0 * (*tau) * exp(*logphi)/(th2*th2) - exp(-(*logphi))/(*tau); //same as in demomhess
    for (j=i+1; j<=p; j++) { ans[i][j]= ans[j][i]=0; }
    ans[i][*n]= ans[*n][i]= th[i]/(exp(*logphi)*(*tau)) + 2.0*(*tau)*exp(*logphi)/(th2*th[i]);
  }
  ans[*n][*n]= -0.5 * exp(-(*logphi)) * (sumth2/(*tau)+(*lambda)) - (*tau) * exp(*logphi) * suminvth2;
}


/************************************************************************
               POSTERIOR SAMPLING UNDER NON-LOCAL PRIORS
************************************************************************/


//Sample from posterior for several non-local priors under a linear model
// Input:
// - niter: number of Gibbs iterations
// - burnin: number of burn-in iterations
// - thinning: only 1 out of each thinning iterations are saved in ans
// - y: y[0..n-1] contains observed response
// - x: covariates stored as a vector (in column order)
// - n: number of observations
// - p: number of variables
// - r: MOM power parameter, i.e. penalty is th[i]^(2*r)
// - tau: prior dispersion
// - a_phi, b_phi: prior for residual variance phi ~ IG(a_phi/2, b_phi/2)
// - prior: prior==0 for MOM, prior==1 for iMOM, prior==2 for eMOM

// [[Rcpp::export]]
SEXP rnlpPostCI_lm(SEXP niter, SEXP burnin, SEXP thinning, SEXP y, SEXP x, SEXP p, SEXP r, SEXP tau, SEXP a_phi, SEXP b_phi, SEXP prior) {
  int n= LENGTH(y), nsave;
  SEXP ans;
  nsave= (int) (floor(INTEGER(niter)[0] - INTEGER(burnin)[0] +.0)/(INTEGER(thinning)[0] +.0));
  ans= PROTECT(Rf_allocVector(REALSXP, nsave * (INTEGER(p)[0]+1)));
  rnlpPost_lm(REAL(ans), INTEGER(niter)[0], INTEGER(burnin)[0], INTEGER(thinning)[0], REAL(y), REAL(x), n, INTEGER(p)[0], INTEGER(r)[0], REAL(tau)[0], REAL(a_phi)[0], REAL(b_phi)[0], INTEGER(prior)[0]);
  UNPROTECT(1);
  return ans;
}

void rnlpPost_lm(double *ans, int niter, int burnin, int thinning, double *y, double *x, int n, int p, int r, double tau, double a_phi, double b_phi, int prior) {
  bool posdef;
  int i, j, k, isave, nsave;
  double *m, *mortho, *alpha, **S, **Sinv, **cholSinv, **inv_cholSinv, **K, **D, tauinv= 1.0/tau, *Xty, *thcur, phicur, phinew, sqrtphi, th2sum, th2invsum, apost, bpost, *linpred, ssr;
  //Pre-compute stuff
  nsave= (int) floor((niter - burnin +.0)/(thinning +.0));
  m= dvector(1,p); mortho= dvector(1,p); alpha= dvector(1,p); thcur= dvector(1,p); linpred= dvector(0,n-1);
  S= dmatrix(1,p,1,p); Sinv= dmatrix(1,p,1,p); cholSinv= dmatrix(1,p,1,p); inv_cholSinv= dmatrix(1,p,1,p); K= dmatrix(1,p,1,p);
  D= dmatrix(1,p,1,p);

  AvectBvec(x, n, p, x, n, p, S); //S= t(x) %*% x + 1/tau
  for (i=1; i<=p; i++) S[i][i] += tauinv;
  inv_posdef(S, p, Sinv, &posdef);
  choldc(Sinv, p, cholSinv, &posdef);
  choldc_inv(Sinv, p, inv_cholSinv, &posdef); //inverse of chol(Sinv)

  Xty= dvector(1,p);
  Atvecx(x, y, Xty+1, 0, p-1, 0, n-1); //m= solve(S) %*% t(x) %*% y
  Ax(Sinv, Xty, m, 1, p, 1, p);
  Ax(inv_cholSinv, m, mortho, 1, p, 1, p);
  free_dvector(Xty,1,p);

  if (prior==0) apost= .5*(a_phi+n+3*p); else if (prior==1) apost= .5*(a_phi+n-p); else apost= .5*(a_phi+n+p);
  //Initialize
  th2sum= 0; phicur= sqrtphi= 1.0;
  for (j=1; j<=p; j++) { thcur[j]= m[j]; th2sum += thcur[j]*thcur[j]; }
  //Ax(cholSinv, thcur, Dthcur, 1, p, 1, p);
  isave= 0;
  for (i=1; i<=niter; i++) {
    //for (j=1; j<=p; j++) Dthcur[j] = Dthcur[j] / sqrtphi;
    Avecx(x, thcur+1, linpred, 0, n-1, 0, p-1);
    ssr=0;
    for (j=0; j<n; j++) ssr += pow(y[j]-linpred[j],2.0);
    if (prior==0) {
      bpost= .5*(b_phi + th2sum/tau + ssr);
      phicur= 1.0/rgammaC(apost, bpost);
      sqrtphi= sqrt(phicur);
    } else {
      if (prior==1) bpost= .5*(b_phi + ssr); else bpost= .5*(b_phi + th2sum/tau + ssr);
      phinew= 1.0/rgammaC(apost, bpost);
      th2invsum= 0;
      for (j=1; j<=p; j++) th2invsum += 1/(thcur[j]*thcur[j]);
      if (runifC() < exp((phicur-phinew)*tau*th2invsum)) { phicur= phinew; sqrtphi= sqrt(phicur); }
    }
    for (j=1; j<=p; j++) {
      alpha[j]= mortho[j]/sqrtphi;
      //Dthcur[j]= Dthcur[j] * sqrtphi;
      for (k=1; k<=j; k++) { D[j][k]= cholSinv[j][k] * sqrtphi; K[j][k]= inv_cholSinv[j][k] / sqrtphi; }
    }
    rnlp_Gibbs(thcur, p, alpha, D, K, &tau, &phicur, r, prior);
    if (i>burnin && ((i-burnin) % thinning)==0) {
      for (j=1; j<=p; j++) ans[isave + (j-1)*nsave]= thcur[j];
      ans[isave + p*nsave]= phicur;
      isave++;
    }
  }
  free_dvector(m,1,p);
  free_dvector(mortho,1,p);
  free_dvector(alpha,1,p);
  free_dvector(thcur,1,p);
  free_dvector(linpred,0,n-1);
  free_dmatrix(S,1,p,1,p); free_dmatrix(Sinv,1,p,1,p); free_dmatrix(cholSinv,1,p,1,p); free_dmatrix(inv_cholSinv,1,p,1,p); free_dmatrix(K,1,p,1,p);
  free_dmatrix(D,1,p,1,p);
}

//Single Gibbs update (th,l) ~ N(th;m,S) * prod g(th[i]) >l[i]
//Input
// - p: dimensionality of th (number of variables)
// - m: m[1..p] is mean of Normal kernel
// - cholS: cholS[1..p][1..p] is Cholesky decomp of covariance
// - K: inverse of cholS
// - tau: value of tau (prior dispersion)
// - phi: value of phi (residual variance)
// - r: power parameter is 2*r
// - prior: prior==0 for MOM, prior==1 for iMOM, prior==2 for eMOM
//Input-Output
// - th: at input th[1..p] is current value of th; at output the updated value
void rnlp_Gibbs(double *th, int p, double *m, double **cholS, double **K, double *tau, double *phi, int r, int prior) {
  int i;
  double *lower, *upper, *l, *z, upperb;
  lower= dvector(1,p); upper= dvector(1,p); l= dvector(1,p); z= dvector(1,p);
  //Update l
  if (prior==0) {
    for (i=1; i<=p; i++) {
      upperb= pen_mom(th+i, phi, tau, r);
      l[i]= runifC() * upperb;  //l[i]= runifC() * pow(th[i]*th[i] / (phi*tau), r + .0);
      if (r==1) { upper[i]= sqrt(l[i] * (*tau) * (*phi)); } else { upper[i]= pow(l[i] * (*tau) * (*phi), 1.0/(2.0*r)); }
      lower[i]= -upper[i];
    }
  } else if (prior==1) {
    for (i=1; i<=p; i++) {
      upperb= pen_imom(th+i, phi, tau, 1);
      l[i]= log(runifC()) + upperb;
      upper[i]= invpen_imom_sandwich(l+i, phi, tau);
      lower[i]= -upper[i];
    }
  } else if (prior==2) {
    for (i=1; i<=p; i++) {
      upperb= pen_emom(th+i, phi, tau, 1);
      l[i]= runifC() * exp(upperb);  //l[i]= runifC() * exp(sqrt(2) + tau*phi/th[i]^2)
      upper[i]= sqrt(fabs((*tau) * (*phi)/(log(l[i])-sqrt(2.0))));
      lower[i]= -upper[i];
    }
  }
  //Update th, cholSth
  Ax(K, th, z, 1, p, 1, p); //z= K th;
  rtmvnormOutside_Gibbs(z, th, m, cholS, p, lower, upper);
  Ax(cholS, z, th, 1, p, 1, p); //th= D z
  free_dvector(lower,1,p); free_dvector(upper,1,p); free_dvector(l,1,p); free_dvector(z,1,p);
}



//Sample from posterior for several non-local priors given mean and variance. Density prop to d(theta) N(theta;m,V)
// Input:
// - niter: number of Gibbs iterations
// - burnin: number of burn-in iterations
// - thinning: only 1 out of each thinning iterations are saved in ans
// - m: mean
// - V: covariance
// - p: number of variables
// - r: MOM power parameter, i.e. penalty is th[i]^(2*r)
// - tau: prior dispersion
// - prior: prior==0 for MOM, prior==1 for iMOM, prior==2 for eMOM

// [[Rcpp::export]]
SEXP rnlpCI(SEXP niter, SEXP burnin, SEXP thinning, SEXP m, SEXP V, SEXP p, SEXP r, SEXP tau, SEXP prior) {
  int nsave;
  SEXP ans;
  nsave= (int) (floor(INTEGER(niter)[0] - INTEGER(burnin)[0] +.0)/(INTEGER(thinning)[0] +.0));
  ans= PROTECT(Rf_allocVector(REALSXP, nsave * (INTEGER(p)[0])));
  rnlp(REAL(ans), INTEGER(niter)[0], INTEGER(burnin)[0], INTEGER(thinning)[0], REAL(m), REAL(V), INTEGER(p)[0], INTEGER(r)[0], REAL(tau)[0], INTEGER(prior)[0]);
  UNPROTECT(1);
  return ans;
}

//Draw from density proportional to d(theta) * N(theta; m, V), where d(theta) is the non-local prior penalty
void rnlp(double *ans, int niter, int burnin, int thinning, double *m, double *Vvec, int p, int r, double tau, int prior) {
  bool posdef;
  int i, j;
  double **S, **Sinv, **cholSinv, **inv_cholSinv, **K, **D;
  //Pre-compute stuff
  S= dmatrix(1,p,1,p); Sinv= dmatrix(1,p,1,p); cholSinv= dmatrix(1,p,1,p); inv_cholSinv= dmatrix(1,p,1,p); K= dmatrix(1,p,1,p);
  D= dmatrix(1,p,1,p);

  for (i=1; i<=p; i++) {
    Sinv[i][i]= Vvec[i-1+p*(i-1)];
    for (j=1; j<i; j++) { Sinv[i][j]= Sinv[j][i]= Vvec[i-1+p*(j-1)]; }
  }
  inv_posdef(Sinv, p, S, &posdef);
  choldc(Sinv, p, cholSinv, &posdef);
  choldc_inv(Sinv, p, inv_cholSinv, &posdef); //inverse of chol(Sinv)

  rnlp_Gibbs_multiple(ans, m-1, p, m-1, cholSinv, inv_cholSinv, &tau, r, prior, niter, burnin, thinning);

  free_dmatrix(S,1,p,1,p); free_dmatrix(Sinv,1,p,1,p); free_dmatrix(cholSinv,1,p,1,p); free_dmatrix(inv_cholSinv,1,p,1,p); free_dmatrix(K,1,p,1,p);
  free_dmatrix(D,1,p,1,p);
}



//Multiple Gibbs draws (th,l) ~ N(th;m,S) * prod g(th[i]) >l[i]
//Input
// - thini: thini[1..p] contains initial value for th
// - B: number of draws
// - p: dimensionality of th (number of variables)
// - m: m[1..p] is mean of Normal kernel
// - cholS: cholS[1..p][1..p] is Cholesky decomp of covariance S
// - K: inverse of cholS
// - tau: value of tau (prior dispersion)
// - r: power parameter for MOM is 2*r
// - prior: prior==0 for MOM, prior==1 for iMOM, prior==2 for eMOM
//Output
// - th: matrix with B rows and p columns (in th[0..(B*p-1)] vector format) with drawn values of th
void rnlp_Gibbs_multiple(double *th, double *thini, int p, double *m, double **cholS, double **K, double *tau, int r, int prior, int niter, int burnin, int thinning) {
  int i, isave, nsave, b;
  double *thcur, *lower, *upper, *l, *z, upperb, phi=1, *alpha;
  nsave= (int) floor((niter - burnin +.0)/(thinning +.0));
  thcur= dvector(1,p); alpha= dvector(1,p);
  lower= dvector(1,p); upper= dvector(1,p); l= dvector(1,p); z= dvector(1,p);
  Ax(K, m, alpha, 1, p, 1, p); //alpha= K %*% mu
  for (i=1; i<=p; i++) { thcur[i]= thini[i]; }
  isave= 0;
  for (b=1; b<=niter; b++) {
    //Update l
    if (prior==0) {
      for (i=1; i<=p; i++) {
        upperb= pen_mom(thcur+i, &phi, tau, r);
        l[i]= runifC() * upperb;
        if (r==1) { upper[i]= sqrt(l[i] * (*tau) * phi); } else { upper[i]= pow(l[i] * (*tau) * phi, 1.0/(2.0*r)); }
        lower[i]= -upper[i];
      }
    } else if (prior==1) {
      for (i=1; i<=p; i++) {
        upperb= pen_imom(thcur+i, &phi, tau, 1);
        l[i]= log(runifC()) + upperb;
        upper[i]= invpen_imom_sandwich(l+i, &phi, tau);
        lower[i]= -upper[i];
      }
    } else if (prior==2) {
      for (i=1; i<=p; i++) {
        upperb= pen_emom(thcur+i, &phi, tau, 1);
        l[i]= runifC() * exp(upperb);
        upper[i]= sqrt(fabs((*tau) * phi/(log(l[i])-sqrt(2.0))));
        lower[i]= -upper[i];
      }
    }
    //Update th, cholSth
    Ax(K, thcur, z, 1, p, 1, p); //z= K th;
    rtmvnormOutside_Gibbs(z, thcur, alpha, cholS, p, lower, upper);
    Ax(cholS, z, thcur, 1, p, 1, p); //th= D z
    if (b>burnin && ((b-burnin) % thinning)==0) {
      for (i=1; i<=p; i++) th[isave + (i-1)*nsave]= thcur[i];
      isave++;
    }
  }
  free_dvector(thcur,1,p); free_dvector(alpha,1,p); free_dvector(lower,1,p); free_dvector(upper,1,p); free_dvector(l,1,p); free_dvector(z,1,p);
}


//Evaluates MOM prior penalty, i.e. (th^2 / (phi*tau))^r
double pen_mom(double *th, double *phi, double *tau, int r) {
  double ans;
  ans= pow(th[0]*th[0] / ((*phi) * (*tau)), r + .0);
  return ans;
}

//Evaluates eMOM prior penalty, i.e. exp(-sqrt(2)*tau*phi/th^2)
double pen_emom(double *th, double *phi, double *tau, int logscale) {
  double ans;
  ans= sqrt(2.0) - (*tau) * (*phi) / (th[0]*th[0]);
  if (logscale==0) ans= exp(ans);
  return ans;
}

//Evaluates iMOM prior penalty, i.e. dimom(th) / dnorm(th)
double pen_imom(double *th, double *phi, double *tau, int logscale) {
  double ans;
  ans= dimom(*th, 0, *tau, *phi, 1) - dnormC(*th, 0, sqrt((*tau)*(*phi)), 1);
  if (logscale==0) ans= exp(ans);
  return(ans);
}

//Evaluates inverse of iMOM penalty, i.e. find th s.t. pen_imom(th,phi,tau) = lambda
//Refines initial search using Newton's algorithm
double invpen_imom_newton(double *loglambda, double *phi, double *tau) {
  int i, maxiter=50;
  double b, d, zcur, thcur, fcur, fpcur, err, ftol=1.0e-5, tauphi= (*tau)*(*phi), halftauphi;
  //Initial guess
  halftauphi= .5*tauphi;
  b= .5*(log((*tau)*(*tau)) + 2.0*log(*phi) + log(2.0)) - (*loglambda);
  d= sqrt(b*b + 2.0);
  zcur= (*tau)*(*phi)*(-b+d);
  thcur= sqrt(zcur);
  fcur= pen_imom(&thcur, phi, tau, 1);
  //Newton search
  err= *loglambda - fcur; i= 1;
  while ((i<maxiter) && (fabs(err) > ftol)) {
    fpcur= -1/zcur + tauphi/(zcur*zcur) + halftauphi;
    zcur += err/fpcur;
    thcur= sqrt(zcur);
    fcur= pen_imom(&thcur, phi, tau, 1);
    err= *loglambda - fcur;
    i++;
  }
  return thcur;
}


// Uses an initial guess to bound the solution. Then uses a sandwhich approach based on recursive linear interpolation
double invpen_imom_sandwich(double *loglambda, double *phi, double *tau) {
  int i, maxiter=50;
  double b, d, zcur, thcur, fcur, zlow, thlow, flow, zup, thup, fup, err, ftol=1.0e-5;
  //Initial guess
  b= .5*(log((*tau)*(*tau)) + 2.0*log(*phi) + log(2.0)) - (*loglambda);
  d= sqrt(b*b + 2.0);
  zcur= (*tau)*(*phi)*(-b+d);
  thcur= sqrt(zcur);
  fcur= pen_imom(&thcur, phi, tau, 1);  //log-penalty at zcur
  //Lower & upper bound
  if (fcur >= (*loglambda)) {
    zlow= .8*.8*zcur; thlow= sqrt(zlow); flow= pen_imom(&thlow, phi, tau, 1);
    while (flow >= (*loglambda)) {
      zcur= zlow; thcur= thlow; fcur= flow;
      zlow= .8*.8*zcur; thlow= sqrt(zlow); flow= pen_imom(&thlow, phi, tau, 1);
    }
    zup= zcur; thup= thcur; fup= fcur;
  } else {
    zup= 1.2*1.2*zcur; thup= sqrt(zup); fup= pen_imom(&thup, phi, tau, 1);
    while (fup <= (*loglambda)) {
      zcur= zup; thcur= thup; fcur= fup;
      zup= 1.2*1.2*zcur; thup= sqrt(zup); fup= pen_imom(&thup, phi, tau, 1);
    }
    zlow= zcur; thlow= thcur; flow= fcur;
  }
  //Search by sandwich linear interpolation
  err= fcur - *loglambda; i= 1;
  while ((i<maxiter) && (fabs(err) > ftol)) {
    b= (fup-flow)/(zup-zlow);
    zcur= zlow + ((*loglambda)-flow)/b; //approx is flow + b*(z-zlow)= loglambda
    thcur= sqrt(zcur);
    fcur= pen_imom(&thcur, phi, tau, 1);
    err= fcur - *loglambda;
    if (err > 0) { zup= zcur; fup= fcur; } else { zlow= zcur; flow= fcur; }
    i++;
  }
  return thcur;
}


/************************************************************************
                       MORE RANDOM VARIATE STUFF
************************************************************************/

/*
 * Truncates a double precision number to an integer.
 *     a - number to be truncated
 */
double fifdint(double a)
{
    long temp = (long)(a);
    return (double)(temp);
}


/**********************************************************************

      void cdfnor(int *which,double *p,double *q,double *x,double *mean,
            double *sd,int *status,double *bound)

               Cumulative Distribution Function
               NORmal distribution


                              Function


     Calculates any one parameter of the normal
     distribution given values for the others.


                              Arguments


     WHICH  --> Integer indicating  which of the  next  parameter
     values is to be calculated using values  of the others.
     Legal range: 1..4
               iwhich = 1 : Calculate P and Q from X,MEAN and SD
               iwhich = 2 : Calculate X from P,Q,MEAN and SD
               iwhich = 3 : Calculate MEAN from P,Q,X and SD
               iwhich = 4 : Calculate SD from P,Q,X and MEAN

     P <--> The integral from -infinity to X of the normal density.
            Input range: (0,1].

     Q <--> 1-P.
            Input range: (0, 1].
            P + Q = 1.0.

     X < --> Upper limit of integration of the normal-density.
             Input range: ( -infinity, +infinity)

     MEAN <--> The mean of the normal density.
               Input range: (-infinity, +infinity)

     SD <--> Standard Deviation of the normal density.
             Input range: (0, +infinity).

     STATUS <-- 0 if calculation completed correctly
               -I if input parameter number I is out of range
                1 if answer appears to be lower than lowest
                  search bound
                2 if answer appears to be higher than greatest
                  search bound
                3 if P + Q .ne. 1

     BOUND <-- Undefined if STATUS is 0

               Bound exceeded by parameter number I if STATUS
               is negative.

               Lower search bound if STATUS is 1.

               Upper search bound if STATUS is 2.


                              Method




     A slightly modified version of ANORM from

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     is used to calulate the  cumulative standard normal distribution.

     The rational functions from pages  90-95  of Kennedy and Gentle,
     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
     starting values to Newton's Iterations which compute the inverse
     standard normal.  Therefore no  searches  are necessary for  any
     parameter.

     For X < -15, the asymptotic expansion for the normal is used  as
     the starting value in finding the inverse standard normal.
     This is formula 26.2.12 of Abramowitz and Stegun.


                              Note


      The normal density is proportional to
      exp( - 0.5 * (( X - MEAN)/SD)**2)

**********************************************************************/
void cdfnor(int *which,
            double *p,
            double *q,
            double *x,
            double *mean,
            double *sd,
            int *status,
            double *bound)
{
static int K1 = 1;
static double z,pq;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    *status = 0;
    if(!(*which < 1 || *which > 4)) goto S30;
    if(!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
S10:
    *bound = 4.0e0;
S20:
    *status = -1;
    return;
S30:
    if(*which == 1) goto S70;
/*
     P
*/
    if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if(!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
S40:
    *bound = 1.0e0;
S50:
    *status = -2;
    return;
S70:
S60:
    if(*which == 1) goto S110;
/*
     Q
*/
    if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if(!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
S80:
    *bound = 1.0e0;
S90:
    *status = -3;
    return;
S110:
S100:
    if(*which == 1) goto S150;
/*
     P + Q
*/
    pq = *p+*q;
    if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*spmpar(&K1))) goto S140;
    if(!(pq < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
S120:
    *bound = 1.0e0;
S130:
    *status = 3;
    return;
S150:
S140:
    if(*which == 4) goto S170;
/*
     SD
*/
    if(!(*sd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
S170:
S160:
/*
     Calculate ANSWERS
*/
    if(1 == *which) {
/*
     Computing P
*/
        z = (*x-*mean)/ *sd;
        cumnor(&z,p,q);
    }
    else if(2 == *which) {
/*
     Computing X
*/
        z = dinvnr(p,q);
        *x = *sd*z+*mean;
    }
    else if(3 == *which) {
/*
     Computing the MEAN
*/
        z = dinvnr(p,q);
        *mean = *x-*sd*z;
    }
    else if(4 == *which) {
/*
     Computing SD
*/
        z = dinvnr(p,q);
        *sd = (*x-*mean)/z;
    }
    return;
}


/*
-----------------------------------------------------------------------

     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.

-----------------------------------------------------------------------
     WRITTEN BY
        ALFRED H. MORRIS, JR.
        NAVAL SURFACE WARFARE CENTER
        DAHLGREN VIRGINIA
-----------------------------------------------------------------------
-----------------------------------------------------------------------
     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
-----------------------------------------------------------------------
*/
double spmpar(int *i)
{
static int K1 = 4;
static int K2 = 8;
static int K3 = 9;
static int K4 = 10;
static double spmpar,b,binv,bm1,one,w,z;
static int emax,emin,ibeta,m;
/*
     ..
     .. Executable Statements ..
*/
    if(*i > 1) goto S10;
    b = ipmpar(&K1);
    m = ipmpar(&K2);
    spmpar = pow(b,(double)(1-m));
    return spmpar;
S10:
    if(*i > 2) goto S20;
    b = ipmpar(&K1);
    emin = ipmpar(&K3);
    one = 1.0;
    binv = one/b;
    w = pow(b,(double)(emin+2));
    spmpar = w*binv*binv*binv;
    return spmpar;
S20:
    ibeta = ipmpar(&K1);
    m = ipmpar(&K2);
    emax = ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta-1;
    one = 1.0;
    z = pow(b,(double)(m-1));
    w = ((z-one)*b+bm1)/(b*z);
    z = pow(b,(double)(emax-2));
    spmpar = w*z*b*b;
    return spmpar;
}


/*
**********************************************************************

     void cumnor(double *arg,double *result,double *ccum)


                              Function


     Computes the cumulative  of    the  normal   distribution,   i.e.,
     the integral from -infinity to x of
          (1/sqrt(2*pi)) exp(-u*u/2) du

     X --> Upper limit of integration.
                                        X is DOUBLE PRECISION

     RESULT <-- Cumulative normal distribution.
                                        RESULT is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative normal distribution.
                                        CCUM is DOUBLE PRECISION

     Renaming of function ANORM from:

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     with slight modifications to return ccum and to deal with
     machine constants.

**********************************************************************
  Original Comments:
------------------------------------------------------------------

 This function evaluates the normal distribution function:

                              / x
                     1       |       -t*t/2
          P(x) = ----------- |      e       dt
                 sqrt(2 pi)  |
                             /-oo

   The main computation evaluates near-minimax approximations
   derived from those in "Rational Chebyshev approximations for
   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
   This transportable program uses rational functions that
   theoretically approximate the normal distribution function to
   at least 18 significant decimal digits.  The accuracy achieved
   depends on the arithmetic system, the compiler, the intrinsic
   functions, and proper selection of the machine-dependent
   constants.

*******************************************************************
*******************************************************************

 Explanation of machine-dependent constants.

   MIN   = smallest machine representable number.

   EPS   = argument below which anorm(x) may be represented by
           0.5  and above which  x*x  will not underflow.
           A conservative value is the largest machine number X
           such that   1.0 + X = 1.0   to machine precision.
*******************************************************************
*******************************************************************

 Error returns

  The program returns  ANORM = 0     for  ARG .LE. XLOW.


 Intrinsic functions required are:

     ABS, AINT, EXP


  Author: W. J. Cody
          Mathematics and Computer Science Division
          Argonne National Laboratory
          Argonne, IL 60439

  Latest modification: March 15, 1992

------------------------------------------------------------------
*/
void cumnor(double *arg,
            double *result,
            double *ccum)
{
static double a[5] = {
    2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
    1.8154981253343561249e04,6.5682337918207449113e-2
};
static double b[4] = {
    4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
    4.5507789335026729956e04
};
static double c[9] = {
    3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
    5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
    1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
};
static double d[8] = {
    2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
    6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
    3.8912003286093271411e04,1.9685429676859990727e04
};
static double half = 0.5e0;
static double p[6] = {
    2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
    1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
};
static double one = 1.0e0;
static double q[5] = {
    1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
    3.78239633202758244e-3,7.29751555083966205e-5
};
static double sixten = 1.60e0;
static double sqrpi = 3.9894228040143267794e-1;
static double thrsh = 0.66291e0;
static double root32 = 5.656854248e0;
static double zero = 0.0e0;
static int K1 = 1;
static int K2 = 2;
static int i;
static double del,eps,temp,x,xden,xnum,y,xsq,min;
/*
------------------------------------------------------------------
  Machine dependent constants
------------------------------------------------------------------
*/
    eps = spmpar(&K1)*0.5e0;
    min = spmpar(&K2);
    x = *arg;
    y = fabs(x);
    if(y <= thrsh) {
/*
------------------------------------------------------------------
  Evaluate  anorm  for  |X| <= 0.66291
------------------------------------------------------------------
*/
        xsq = zero;
        if(y > eps) xsq = x*x;
        xnum = a[4]*xsq;
        xden = xsq;
        for(i=0; i<3; i++) {
            xnum = (xnum+a[i])*xsq;
            xden = (xden+b[i])*xsq;
        }
        *result = x*(xnum+a[3])/(xden+b[3]);
        temp = *result;
        *result = half+temp;
        *ccum = half-temp;
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
------------------------------------------------------------------
*/
    else if(y <= root32) {
        xnum = c[8]*y;
        xden = y;
        for(i=0; i<7; i++) {
            xnum = (xnum+c[i])*y;
            xden = (xden+d[i])*y;
        }
        *result = (xnum+c[7])/(xden+d[7]);
        xsq = fifdint(y*sixten)/sixten;
        del = (y-xsq)*(y+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for |X| > sqrt(32)
------------------------------------------------------------------
*/
    else  {
        *result = zero;
        xsq = one/(x*x);
        xnum = p[5]*xsq;
        xden = xsq;
        for(i=0; i<4; i++) {
            xnum = (xnum+p[i])*xsq;
            xden = (xden+q[i])*xsq;
        }
        *result = xsq*(xnum+p[4])/(xden+q[4]);
        *result = (sqrpi-*result)/y;
        xsq = fifdint(x*sixten)/sixten;
        del = (x-xsq)*(x+xsq);
        *result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
        *ccum = one-*result;
        if(x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    if(*result < min) *result = 0.0e0;
/*
------------------------------------------------------------------
  Fix up for negative argument, erf, etc.
------------------------------------------------------------------
----------Last card of ANORM ----------
*/
    if(*ccum < min) *ccum = 0.0e0;
}


/*
**********************************************************************

     double dinvnr(double *p,double *q)
     Double precision NoRmal distribution INVerse


                              Function


     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


                              Arguments


     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION

     Q --> 1-P
                    P is DOUBLE PRECISION


                              Method


     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
     value for the Newton method of finding roots.


                              Note


     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)

**********************************************************************
*/
double dinvnr(double *p,
              double *q)
{
#define maxit 100
#define eps 1.0e-13
#define r2pi 0.3989422804014326e0
#define nhalf -0.5e0
#define dennor(x) (r2pi*exp(nhalf*(x)*(x)))
static double dinvnr,strtx,xcur,cum,ccum,pp,dx;
static int i;
static unsigned long qporq;
/*
     ..
     .. Executable Statements ..
*/
/*
     FIND MINIMUM OF P AND Q
*/
    qporq = *p <= *q;
    if(!qporq) goto S10;
    pp = *p;
    goto S20;
S10:
    pp = *q;
S20:
/*
     INITIALIZATION STEP
*/
    strtx = stvaln(&pp);
    xcur = strtx;
/*
     NEWTON INTERATIONS
*/
    for(i=1; i<=maxit; i++) {
        cumnor(&xcur,&cum,&ccum);
        dx = (cum-pp)/dennor(xcur);
        xcur -= dx;
        if(fabs(dx/xcur) < eps) goto S40;
    }
    dinvnr = strtx;
/*
     IF WE GET HERE, NEWTON HAS FAILED
*/
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
S40:
/*
     IF WE GET HERE, NEWTON HAS SUCCEDED
*/
    dinvnr = xcur;
    if(!qporq) dinvnr = -dinvnr;
    return dinvnr;
#undef maxit
#undef eps
#undef r2pi
#undef nhalf
#undef dennor
}


/*
**********************************************************************

     double stvaln(double *p)
                    STarting VALue for Neton-Raphon
                calculation of Normal distribution Inverse


                              Function


     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


                              Arguments


     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION


                              Method


     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980.

**********************************************************************
*/
double stvaln(double *p)
{
static double xden[5] = {
    0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
    0.38560700634e-2
};
static double xnum[5] = {
    -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
    -0.453642210148e-4
};
static int K1 = 5;
static double stvaln,sign,y,z;
/*
     ..
     .. Executable Statements ..
*/
    if(!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
S10:
    sign = 1.0e0;
    z = 1.0e0-*p;
S20:
    y = sqrt(-(2.0e0*log(z)));
    stvaln = y+devlpl(xnum,&K1,&y)/devlpl(xden,&K1,&y);
    stvaln = sign*stvaln;
    return stvaln;
}


/*
**********************************************************************

     double devlpl(double a[],int *n,double *x)
              Double precision EVALuate a PoLynomial at X


                              Function


     returns
          A(1) + A(2)*X + ... + A(N)*X**(N-1)


                              Arguments


     A --> Array of coefficients of the polynomial.
                                        A is DOUBLE PRECISION(N)

     N --> Length of A, also degree of polynomial - 1.
                                        N is INTEGER

     X --> Point at which the polynomial is to be evaluated.
                                        X is DOUBLE PRECISION

**********************************************************************
*/
double devlpl(double a[],
              int *n,
              double *x)
{
static double devlpl,term;
static int i;
/*
     ..
     .. Executable Statements ..
*/
    term = a[*n-1];
    for(i= *n-1-1; i>=0; i--) term = a[i]+term**x;
    devlpl = term;
    return devlpl;
}


/*
-----------------------------------------------------------------------

     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

  INTEGERS.

     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

     IPMPAR(1) = A, THE BASE.

     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.

     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.

  FLOATING-POINT NUMBERS.

     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

     IPMPAR(4) = B, THE BASE.

  SINGLE-PRECISION

     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

  DOUBLE-PRECISION

     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

-----------------------------------------------------------------------

     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED REMOVE
     THE COMMENT DELIMITORS FROM THE DEFINITIONS DIRECTLY BELOW THE NAME
     OF THE MACHINE

-----------------------------------------------------------------------

     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.

-----------------------------------------------------------------------
     .. Scalar Arguments ..
*/
int ipmpar(int *i)
{
static int imach[11];
static int ipmpar;
/*     MACHINE CONSTANTS FOR AMDAHL MACHINES. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
       PC 7300, AND AT&T 6300. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 33;
   imach[3] = 8589934591;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -256;
   imach[7] = 255;
   imach[8] = 60;
   imach[9] = -256;
   imach[10] = 255;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -50;
   imach[10] = 76;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -32754;
   imach[10] = 32780;
*/
/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 48;
   imach[3] = 281474976710655;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -974;
   imach[7] = 1070;
   imach[8] = 95;
   imach[9] = -926;
   imach[10] = 1070;
*/
/*     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS/VE OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -4096;
   imach[7] = 4095;
   imach[8] = 96;
   imach[9] = -4096;
   imach[10] = 4095;
*/
/*     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3. */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 47;
   imach[6] = -8189;
   imach[7] = 8190;
   imach[8] = 94;
   imach[9] = -8099;
   imach[10] = 8190;
*/
/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200. */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE HARRIS 220. */
/*
   imach[1] = 2;
   imach[2] = 23;
   imach[3] = 8388607;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 38;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
       AND DPS 8/70 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 63;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       3 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 39;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       4 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 55;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 9000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -126;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
       5/7/9 AND THE SEL SYSTEMS 85/86. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE IBM PC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
       MACFORTRAN II. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 54;
   imach[9] = -101;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 62;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
       32-BIT INTEGER ARITHMETIC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
       SERIES (MIPS R3000 PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). */

   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 60;
   imach[9] = -1024;
   imach[10] = 1023;
*/
/*     MACHINE CONSTANTS FOR THE VAX 11/780. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
    ipmpar = imach[*i];
    return ipmpar;
}


/************************************************************************
                        EVEN MORE STUFF
************************************************************************/

/*
**********************************************************************
     double genunf(double low,double high)
               GeNerate Uniform Real between LOW and HIGH
                              Function
     Generates a real uniformly distributed between LOW and HIGH.
                              Arguments
     low --> Low bound (exclusive) on real value to be generated
     high --> High bound (exclusive) on real value to be generated
**********************************************************************
*/
double genunf(double low,
              double high)
{
static double genunf;

    if (!(low > high)) goto S10;
    REprintf("genunf: low > high: low=%16.6E, high=%16.6E\n", low, high);
    _cstatfatal();
    /*NOTREACHED*/
S10:
    genunf = low+(high-low)*ranf();
    return genunf;
}


/*
**********************************************************************
     double gengam(double a,double r)
           GENerates random deviates from GAMma distribution
                              Function
     Generates random deviates from the gamma distribution whose
     density is
          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
                              Arguments
     a --> Location parameter of Gamma distribution
     r --> Shape parameter of Gamma distribution
     CAREFUL: order of parameters is reversed wrt usual notation. mean is r/a; var is r/a^2
                              Method
     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               (Case R >= 1.0)
               Ahrens, J.H. and Dieter, U.
               Generating Gamma Variates by a
               Modified Rejection Technique.
               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
     Algorithm GD
               (Case 0.0 <= R <= 1.0)
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling from Gamma,
               Beta, Poisson and Binomial Distributions.
               Computing, 12 (1974), 223-246/
     Adapted algorithm GS.
**********************************************************************
*/
double gengam(double a,
              double r)
{
static double gengam;

    gengam = sgamma(r);
    gengam /= a;
    return gengam;
}


/*
**********************************************************************


     (STANDARD-)  G A M M A  DISTRIBUTION


**********************************************************************
**********************************************************************

               PARAMETER  A >= 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               GENERATING GAMMA VARIATES BY A
               MODIFIED REJECTION TECHNIQUE.
               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.

     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER
                                 (STRAIGHTFORWARD IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************

               PARAMETER  0.0 < A < 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               COMPUTER METHODS FOR SAMPLING FROM GAMMA,
               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.
               COMPUTING, 12 (1974), 223 - 246.

     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)

**********************************************************************
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
double sgamma(double a)
{
static double q1 = 4.166669E-2;
static double q2 = 2.083148E-2;
static double q3 = 8.01191E-3;
static double q4 = 1.44121E-3;
static double q5 = -7.388E-5;
static double q6 = 2.4511E-4;
static double q7 = 2.424E-4;
static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 4.07753E-2;
static double e5 = 1.0293E-2;
static double aa = 0.0;
static double aaa = 0.0;
static double sqrt32 = 5.656854;
static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;
    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm();
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = ranf();
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s-0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-2*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = sexpo();
    u = ranf();
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0+0.3678794*a;
S130:
    p = b*ranf();
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(sexpo() < sgamma) goto S130;
    return sgamma;
S140:
    sgamma = -log((b-p)/ a);
    if(sexpo() < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;
}


/*
**********************************************************************


     (STANDARD-)  N O R M A L  DISTRIBUTION


**********************************************************************
**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
               SAMPLING FROM THE NORMAL DISTRIBUTION.
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************
     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
*/
double snorm(void)
{
static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
};
static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
};
static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
};
static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
};
static long i;
static double snorm,u,s,ustar,aa,w,y,tt;
    u = ranf();
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long) (u);
    if(i == 32) i = 31;
    if(i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = ranf();
S80:
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = ranf();
    if(ustar > tt) goto S50;
    u = ranf();
    if(ustar >= u) goto S150;
    u = ranf();
    goto S140;
}


/* Transfers sign of argument sign to argument num */
double fsign(double num,
             double sign)
{
    return ((sign > 0.0 && num < 0.0) ||
            (sign < 0.0 && num > 0.0)) ? -num : num;
}


/*
**********************************************************************


     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION


**********************************************************************
**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               COMPUTER METHODS FOR SAMPLING FROM THE
               EXPONENTIAL AND NORMAL DISTRIBUTIONS.
               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.

     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM
     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************
     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
*/
double sexpo(void)
{
static double q[8] = {
    0.6931472,0.9333737,0.9888778,0.9984959,0.9998293,0.9999833,0.9999986,1.0
};
static long i;
static double sexpo,a,u,ustar,umin;
static double *q1 = q;
    a = 0.0;
    u = ranf();
    goto S30;
S20:
    a += *q1;
S30:
    u += u;
    if(u <= 1.0) goto S20;
    u -= 1.0;
    if(u > *q1) goto S60;
    sexpo = a+u;
    return sexpo;
S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
S70:
    ustar = ranf();
    if(ustar < umin) umin = ustar;
    i += 1;
    if(u > *(q+i-1)) goto S70;
    sexpo = a+umin**q1;
    return sexpo;
}


/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
long mltmod(long a,
            long s,
            long m)
{
#define h 32768L
static long mltmod,a0,a1,k,p,q,qh,rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if(!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
    REprintf("mltmod: requires (0 < a < m); (0 < s < m): ");
    REprintf("a = %12ld, s = %12ld, m = %12ld\n", a, s, m);
    _cstatfatal();
    /*NOTREACHED*/
S10:
    if(!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
S20:
    a1 = a/h;
    a0 = a-h*a1;
    qh = m/h;
    rh = m-h*qh;
    if(!(a1 >= h)) goto S50;
    a1 -= h;
    k = s/qh;
    p = h*(s-k*qh)-k*rh;
S30:
    if(!(p < 0)) goto S40;
    p += m;
    goto S30;
S40:
    goto S60;
S50:
    p = 0;
S60:
/*
     P = (A2*S*H)MOD M
*/
    if(!(a1 != 0)) goto S90;
    q = m/a1;
    k = s/q;
    p -= (k*(m-a1*q));
    if(p > 0) p -= m;
    p += (a1*(s-k*q));
S70:
    if(!(p < 0)) goto S80;
    p += m;
    goto S70;
S90:
S80:
    k = p/qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = h*(p-k*qh)-k*rh;
S100:
    if(!(p < 0)) goto S110;
    p += m;
    goto S100;
S120:
S110:
    if(!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m/a0;
    k = s/q;
    p -= (k*(m-a0*q));
    if(p > 0) p -= m;
    p += (a0*(s-k*q));
S130:
    if(!(p < 0)) goto S140;
    p += m;
    goto S130;
S150:
S140:
    mltmod = p;
    return mltmod;
#undef h
}


/*
**********************************************************************
     double ranf(void)
                RANDom number generator as a Function
     Returns a random doubleing point number from a uniform distribution
     over 0 - 1 (endpoints of this interval are not returned) using the
     current generator
     This is a transcription from Pascal to Fortran of routine
     Uniform_01 from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
double ranf(void)
{
static double ranf;
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
    ranf = ignlgi()*4.656613057E-10;
    return ranf;
}


/*
**********************************************************************
     void gscgn(long getset,long *g)
                         Get/Set GeNerator
     Gets or returns in G the number of the current generator
                              Arguments
     getset --> 0 Get
                1 Set
     g <-- Number of the current random number generator (1..32)
**********************************************************************
*/
void gscgn(long getset,
           long *g)
{
#define numg 32L
static long curntg = 1;
    if(getset == 0) *g = curntg;
    else  {
        if(*g < 0 || *g > numg) {
            REprintf("gscgn: generator number out of range\n");
            _cstatfatal();
            /*NOTREACHED*/
        }
        curntg = *g;
    }
#undef numg
}


/*
**********************************************************************
     void gsrgs(long getset,long *qvalue)
               Get/Set Random Generators Set
     Gets or sets whether random generators set (initialized).
     Initially (data statement) state is not set
     If getset is 1 state is set to qvalue
     If getset is 0 state returned in qvalue
**********************************************************************
*/
void gsrgs(long getset,
           long *qvalue)
{
static long qinit = 0;

    if(getset == 0) *qvalue = qinit;
    else qinit = *qvalue;
}


/*
**********************************************************************
     void gssst(long getset,long *qset)
          Get or Set whether Seed is Set
     Initialize to Seed not Set
     If getset is 1 sets state to Seed Set
     If getset is 0 returns T in qset if Seed Set
     Else returns F in qset
**********************************************************************
*/
void gssst(long getset,
           long *qset)
{
static long qstate = 0;
    if(getset != 0) qstate = 1;
    else  *qset = qstate;
}


/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
void setall(long iseed1,
            long iseed2)
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1vw,Xa2vw,Xig1[32],Xig2[32];
static long T1;
static long g,ocgn;
static long qrgnin;
    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1,&T1);
    gscgn(0L,&ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    *Xig1 = iseed1;
    *Xig2 = iseed2;
    initgn(-1L);
    for(g=2; g<=numg; g++) {
        *(Xig1+g-1) = mltmod(Xa1vw,*(Xig1+g-2),Xm1);
        *(Xig2+g-1) = mltmod(Xa2vw,*(Xig2+g-2),Xm2);
        gscgn(1L,&g);
        initgn(-1L);
    }
    gscgn(1L,&ocgn);
#undef numg
}


/*
**********************************************************************
     void initgn(long isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
void initgn(long isdtyp)
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1w,Xa2w,Xig1[32],Xig2[32],Xlg1[32],Xlg2[32],Xcg1[32],Xcg2[32];
static long g;
static long qrgnin;
/*
     Abort unless random number generator initialized
*/
    gsrgs(0L,&qrgnin);
    if (qrgnin) goto S10;
        REprintf("initgn: random number generator not initialized\n");
        _cstatfatal();
        /*NOTREACHED*/
S10:
    gscgn(0L,&g);
    if(-1 != isdtyp) goto S20;
    *(Xlg1+g-1) = *(Xig1+g-1);
    *(Xlg2+g-1) = *(Xig2+g-1);
    goto S50;
S20:
    if(0 != isdtyp) goto S30;
    goto S50;
S30:
/*
     do nothing
*/
    if(1 != isdtyp) goto S40;
    *(Xlg1+g-1) = mltmod(Xa1w,*(Xlg1+g-1),Xm1);
    *(Xlg2+g-1) = mltmod(Xa2w,*(Xlg2+g-1),Xm2);
    goto S50;
S40:
    REprintf("initgn: isdtyp not in range\n");
    _cstatfatal();
    /*NOTREACHED*/
S50:
    *(Xcg1+g-1) = *(Xlg1+g-1);
    *(Xcg2+g-1) = *(Xlg2+g-1);
#undef numg
}


/*
**********************************************************************
     long ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
long ignlgi(void)
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern void gssst(long getset,long *qset);
extern void gscgn(long getset,long *g);
extern long Xm1,Xm2,Xa1,Xa2,Xcg1[32],Xcg2[32];
extern long Xqanti[32];
static long ignlgi,curntg,k,s1,s2,z;
static long qqssd,qrgnin;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L,&qrgnin);
    if(!qrgnin) inrgcm();
    gssst(0,&qqssd);
    if(!qqssd) setall(1234567890L,123456789L);
/*
     Get Current Generator
*/
    gscgn(0L,&curntg);
    s1 = *(Xcg1+curntg-1);
    s2 = *(Xcg2+curntg-1);
    k = s1/53668L;
    s1 = Xa1*(s1-k*53668L)-k*12211;
    if(s1 < 0) s1 += Xm1;
    k = s2/52774L;
    s2 = Xa2*(s2-k*52774L)-k*3791;
    if(s2 < 0) s2 += Xm2;
    *(Xcg1+curntg-1) = s1;
    *(Xcg2+curntg-1) = s2;
    z = s1-s2;
    if(z < 1) z += (Xm1-1);
    if(*(Xqanti+curntg-1)) z = Xm1-z;
    ignlgi = z;
    return ignlgi;
#undef numg
}


/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
void inrgcm(void)
{
#define numg 32L
extern void gsrgs(long getset,long *qvalue);
extern long Xm1,Xm2,Xa1,Xa2,Xa1w,Xa2w,Xa1vw,Xa2vw;
extern long Xqanti[32];
static long T1;
static long i;
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1w = 1033780774L;
    Xa2w = 1494757890L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    for(i=0; i<numg; i++) *(Xqanti+i) = 0;
    T1 = 1;
/*
     Tell the world that common has been initialized
*/
    gsrgs(1L,&T1);
#undef numg
}


/************************************************************************
                      INTEGRATION
************************************************************************/

/*
 * Computes the nth stage of refinement of an extended midpoint rule.
 * func is input as a pointer to the function to be integrated between
 * limits a and b. When called with n=1, the routine return the crudest
 * estimate. Subsequent calls with n=2,3... (in that order) will improve
 * the accuracy of s by adding (2/3)*3^(n-1) additional interior points.
 * s should not be modified between sequential calls.
 */
double midpnt(double (* const func)(double),
              double a,
              double b,
              int n)
{
#define FUNC(x) ((*func)(x))

    static double s;

    //assert(func != NULL);

    if (n == 1) {
        s = (b - a) * FUNC(0.5 * (a + b));
    }
    else {
        int j;
        double x;
        double tnm;
        double del;
        double ddel;
        double sum = 0.0;
        int it = 1;

        for (j = 1; j < n-1; j++) {
            it *= 3;
        }
        tnm = it;
        del = (b - a) / (3.0 * tnm);
        ddel = del + del;
        x = a + 0.5 * del;
        for (j = 1; j <= it; j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s = (s + (b - a) * sum / tnm) / 3.0;
    }
    return s;
#undef FUNC
}


/*
 * This routine is an exact replacement for midpnt(), i.e. returns the nth
 * stage of refinement of the integral of funk from aa to bb, except that the
 * function is evaluated at evenly spaced points in 1/x rather than in x.
 * This allows the upper limit bb to be as large and positive as the computer
 * allows, or the lower limit aa to be as large and negative, but not both.
 * aa and bb must have the same sign, and they cannot be equal to zero.
 */
double midinf(double (* const func)(double),
              double aa,
              double bb,
              int n)
{
#define FUNC(x) ((*func)(1.0 / (x)) / ((x) * (x)))

    static double s;
    double a;
    double b;

    //assert(func != NULL);

    if (SIGN(aa) != SIGN(bb)) {
        nrerror("midinf", "", "aa and bb must have same sign");
        /*NOTREACHED*/
    }
    if ((aa == 0.0) || (bb == 0.0)) {
        nrerror("midinf", "", "aa and/or bb is zero");
        /*NOTREACHED*/
    }

    /* Change the limits of integration */
    b = 1.0 / aa;
    a = 1.0 / bb;

    /* From here on, the routine is identical to midpnt() */
    if (n == 1) {
        s = (b - a) * FUNC(0.5 * (a + b));
    }
    else {
        int j;
        double x;
        double tnm;
        double del;
        double ddel;
        double sum = 0.0;
        int it = 1;

        for (j = 1; j < n-1; j++) {
            it *= 3;
        }
        tnm = it;
        del = (b - a) / (3.0 * tnm);
        ddel = del + del;
        x = a + 0.5 * del;
        for (j = 1; j <= it; j++) {
            sum += FUNC(x);
            x += ddel;
            sum += FUNC(x);
            x += del;
        }
        s = (s + (b - a) * sum / tnm) / 3.0;
    }
    return s;
#undef FUNC
}


/*
 * Romberg integration on an open interval. Returns the integral of the
 * function func from a to b, using any specified integrating function choose
 * and Romberg's method. Normally choose will be an open formula like midpnt
 * or midinf, not evaluating the function at the endpoints. It is assumed that
 * choose triples the number of steps on each call, and that its error series
 * contains only even powers of the number of steps. Integration is performed
 * by Romberg's method of order 2K (K=2 is Simpson's rule)
 *
 * Examples:
 *
 * Integrate f from 0 to 2
 * answer = qromo(f, 0.0, 2.0, midpnt);
 *
 * Integrate f from 0 to 1.0e30 (cutpoint should be chosen in the tail of fx)
 * answer = qromo(f, 0.0, 2.0, midpnt) + qromo(f, 2.0, 1.0e30, midinf)
 */
double qromo(double (* const func)(double),
             double a,
             double b,
             double (* const choose)(double(* const)(double),
                                     double,
                                     double,
                                     int))
{
    const double EPS = 1.0e-6;
    const int JMAX = 14;
    const int JMAXP = (JMAX+1);
    const int K = 5;

    int j;
    double ss;
    double dss;
    double h[JMAXP+1];
    double s[JMAXP];

    //assert(func != NULL);
    //assert(choose != NULL);

    h[1] = 1.0;
    for (j = 1; j <= JMAX; j++) {
        s[j] = (*choose)(func, a, b, j);
        if (j >= K) {
            polint(&h[j-K], &s[j-K], K, 0.0, &ss, &dss);
            if (fabs(dss) <= EPS*fabs(ss)) {
                return ss;
            }
        }
        /* Assumes step tripling and even error series */
        h[j+1] = h[j] / 9.0;
    }
    nrerror("qromo", "integrate a function", "");
    /*NOTREACHED*/
    return(0.0);    /* make compiler happy */
}


/************************************************************************
          INTERPOLATION, EXTRAPOLATION AND SPLINES
************************************************************************/

/*
 * Given arrays xa[1..n] and ya[1..n], and given a value x, this routine
 * returns a value y, and and error estimate dy. If P(x) is the polynomial of
 * degree N-1 such that P(xa[i])=ya[i] i=1..n then the returned value y=P(x).
 */
void polint(double xa[],
            double ya[],
            int n,
            double x,
            double *y,
            double *dy)
{
    int i;
    int m;
    int ns = 1;
    double dif;
    double *c;
    double *d;

    //assert(y != NULL);
    //assert(dy != NULL);

    dif = fabs(x - xa[1]);
    c = dvector(1, n);
    d = dvector(1, n);

    /* Find the index ns of the closest table entry */
    for (i = 1; i <= n; i++) {
        double dift;

        dift = fabs(x - xa[i]);
        if (dift < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    *y = ya[ns--];          /* initial approximation to y */

    /* For each column of the tableau, update the c's and d's */
    for (m = 1; m < n; m++) {
        for (i = 1; i <= (n - m); i++) {
            double den;
            double ho;
            double hp;
            double w;

            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            den = ho - hp;
            if (den == 0.0) {
                nrerror("polint",
                        "",
                 "increment x axis in 0 units (two identical input x values)");
                /*NOTREACHED*/
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }

        /*
         * After each column in the tableau is completed, we decide which
         * correction (c or d) to add to our accumulating value of y.
         * i.e., which path to take through the tableau: forking up or down.
         * We do this in such a way as to take the most "straight line"
         * route to its apex, updating ns accordingly to keep track of where
         * we are. This route keeps the partial approximations centered on
         * the target x. The last dy added is thus the error indication.
         */

        *dy = (2 * ns) < (n - m) ? c[ns + 1] : d[ns--];
        *y += *dy;
    }

    free_dvector(d, 1, n);
    free_dvector(c, 1, n);
}


/*
 * Returns the jth B-spline basis evaluated at single value x
 *    x: value at which to evaluate the B-spline basis
 *    j: basis
 *    degree: degree of the B-spline (0: piecewise constant, 1:linear etc.)
 *    knots: sequence of knots
 */
double bspline_singlex(double x,
                       int j,
                       int degree,
                       const double *knots)
{
    double ans;

    //assert(knots != NULL);

    if (degree == 0) {
        ans = (knots[j] <= x && x < knots[j+1]) ? 1.0 : 0.0;
    } else {
        ans = bspline_singlex(x, j, degree-1, knots) *
              (x - knots[j]) / (knots[j+degree] - knots[j]) +
              bspline_singlex(x, j+1, degree-1, knots) *
              (knots[j+degree+1] - x) / (knots[j+degree+1] - knots[j+1]);
    }
    return(ans);
}


/*
 * B-spline basis eval at vector of values x.
 * Normalized to sum to 1 at any x value.
 *
 * Input:
 *    x: vector of values at which to evaluate the B-spline basis
 *    nx: length of x
 *    degree: degree of the spline (0: piecewise constant, 1: linear etc.)
 *    knots: vector with positions of the knots
 *    nknots: length of knots
 *  Output:
 *    W: matrix[nx][nknots-degree-1] containing the B-spline basis
 */
void bspline(double **W,
             const double *x,
             const int *nx,
             const int *degree,
             const double *knots,
             const int *nknots)
{
    //assert(W != NULL);
    //assert(x != NULL);
    //assert(nx != NULL);
    //assert(degree != NULL);
    //assert(knots != NULL);
    //assert(nknots != NULL);

    if (*nknots < (*degree + 2)) {
        REprintf("bspline: number of knots must be >= degree+2\n");
        /* :TBD: - Should this be fatal? */
    }
    else {
        int i, j, ncols = (*nknots - *degree - 1);

        for (i = 0; i < (*nx); i++) {
            for (j = 0; j < ncols; j++) {
                W[i][j] = bspline_singlex(x[i], j, *degree, knots);
            }
        }
    }
}



/*
 * Same as bspline() but uses a vector as its first argument so that it can be called from R
 */
void bspline_vec(double *W,
                 const double *x,
                 const int *nx,
                 const int *degree,
                 const double *knots,
                 const int *nknots)
{
    int i;
    int j;
    double **Wtemp;
    int ncols;

    //assert(W != NULL);
    //assert(x != NULL);
    //assert(nx != NULL);
    //assert(degree != NULL);
    //assert(knots != NULL);
    //assert(nknots != NULL);

    ncols = *nknots - *degree - 1;
    Wtemp = dmatrix(0, *nx, 0, ncols);
    bspline(Wtemp, x, nx, degree, knots, nknots);
    for (i = 0; i < (*nx); i++) {
        for (j = 0; j < ncols; j++) {
            W[(i * ncols) + j] = Wtemp[i][j];
        }
    }
    free_dmatrix(Wtemp, 0, *nx, 0, *nknots- *degree -1);
}


/*
 * M-spline basis eval at vector of values x.
 * Normalized to integrate to 1 wrt x
 */
void mspline(double **W,
             const double *x,
             const int *nx,
             const int *degree,
             const double *knots,
             const int *nknots)
{
    //assert(W != NULL);
    //assert(x != NULL);
    //assert(nx != NULL);
    //assert(degree != NULL);
    //assert(knots != NULL);
    //assert(nknots != NULL);

    if (*nknots < (*degree + 2)) {
        REprintf("mspline: number of knots must be >= degree+2\n");
        /* :TBD: - Should this be fatal? */
    }
    else {
        int i;
        int j;
        int ncols = (*nknots - *degree - 1);

        for (i = 0; i < (*nx); i++) {
            for (j = 0; j < ncols; j++) {
                W[i][j] = bspline_singlex(x[i], j, *degree, knots) *
                          (*degree+1.0) / (knots[j+ *degree +1] - knots[j]);
            }
        }
    }
}


/*
 * Same as mspline() but uses a vector as its first argument so that it can
 * be called from R
 */
void mspline_vec(double *W,
                 const double *x,
                 const int *nx,
                 const int *degree,
                 const double *knots,
                 const int *nknots)
{
    int i;
    int j;
    double **Wtemp;
    int ncols;

    //assert(W != NULL);
    //assert(x != NULL);
    //assert(nx != NULL);
    //assert(degree != NULL);
    //assert(knots != NULL);
    //assert(nknots != NULL);

    ncols = *nknots - *degree - 1;
    Wtemp = dmatrix(0, *nx, 0, ncols);
    mspline(Wtemp, x, nx, degree, knots, nknots);
    for (i = 0; i < (*nx); i++) {
        for (j = 0; j < ncols; j++) {
            W[(i * ncols) + j] = Wtemp[i][j];
        }
    }
    free_dmatrix(Wtemp, 0, *nx, 0, ncols);
}


/************************************************************************
            FUNCTION OPTIMIZATION
************************************************************************/

/*
 * PURPOSE: UNIVARIATE OPTIMIZATION WITHOUT DERIVATIVE INFORMATION
 *
 * Given a function f and a bracketing triplet of abscissas ax, bx, cx
 * (such that bx is between ax & cx and f(bx) is less than f(ax) and f(cx)),
 * this routine isolates the minimum to a fractional precision of about eps
 * using Brent's method. The abscissa of the minimum is returned as xmin and
 * the value of the function at the minimum is the returned value.
 */
double univmin(double ax,
               double bx,
               double cx,
               double (* const f)(double),
               double eps,
               double *xmin,
               int itmax)
{
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

  const double CGOLD = 0.3819660;  //golden ratio
  const double ZEPS = 1.0e-10;   //protects against trying to achieve fractional accuracy when min is exactly zero
  int iter;
  double a,b,d=1,etemp,fu,fv,fw,fx,p,q,r,eps1,eps2,u,v,w,x,xm;
  double e=0.0;

  //assert(f != NULL);
  //assert(xmin != NULL);

  a=(ax < cx ? ax : cx);     //a,b must be in ascending order but input abscissas need not be
  b=(ax > cx ? ax : cx);

  x=w=v=bx;                  //initializations
  fw=fv=fx=(*f)(x);

  for (iter=1; iter<=itmax; iter++) {      //main program loop
    xm=0.5*(a+b);
    eps2= 2.0*(eps1=eps*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (eps2-0.5*(b-a))) {  //test for done here
      *xmin= x;
      return fx;
    }
    if (fabs(e) > eps1) {                  //construct a trial parabolic fit
      r= (x-w)*(fx-fv);
      q= (x-v)*(fx-fw);
      p= (x-v)*q - (x-w)*r;
      q= 2.0*(q-r);
      if (q>0.0) p= -p;
      q= fabs(q);
      etemp= e;
      e= d;
      //determine acceptability of parabolic fit & take the golden section step into the larger of the two segments
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
	d= CGOLD*(e=(x >= xm ? a-x : b-x));
      } else {
	d=p/q;                            //take the parabolic step
	u= x+d;
	if (u-a < eps2 || b-u < eps2) d=SETSIGN(eps1,xm-x);
      }
    } else {
      d= CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u= (fabs(d) >= eps1 ? x+d : x+SETSIGN(eps1,d));
    fu= (*f)(u);                         //this is the one function evaluation per iteration
    if (fu <= fx) {                      //now decide what to do with our function evaluation
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);                     //housekeeping
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }                                   //done with housekeeping. Back for another iteration
  }

  *xmin= x;                             //only get here if iteration limit is reached
  return fx;
#undef SHFT
}


/*
 * Given a function f and its derivative function df, and given a bracketing
 * triplet of abscissas ax, bx, cx [such that bx is between ax and cx, and
 * f(bx) is less than both f(ax) and f(cx)], this routine isolates the minimum
 * to a fractional precision of about eps using a modification of Brent's
 * method that uses derivatives. The abscissa of the minimum is returned as
 * xmin and the value of the function at the minimum is the returned value.
 */
double dunivmin(double ax,
                double bx,
                double cx,
                double (* const f)(double),
                double (* const df)(double),
                double eps,
                double *xmin,
                int itmax)
{
#define MOV3(a,b,c,d,e,f) (a)=(d);(b)=(e);(c)=(f)

  const double ZEPS = 1.0e-10;   //protects against trying to achieve fractional accuracy when min is exactly zero
  int iter,ok1,ok2; //Will be used as flags for whether proposed steps are acceptable or not.
  double a,b,d=1,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,eps1,eps2,u,u1,u2,v,w,x,xm;

  //assert(f != NULL);
  //assert(df != NULL);
  //assert(xmin != NULL);

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x);
  for (iter=1;iter<=itmax;iter++) {
    xm=0.5*(a+b);
    eps1=eps*fabs(x)+ZEPS;
    eps2=2.0*eps1;
    if (fabs(x-xm) <= (eps2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > eps1) {
      d1=2.0*(b-a); //Initialize these d's to an out-of-bracket value
      d2=d1;
      if (dw != dx) d1=(w-x)*dx/(dx-dw); //Secant method with one point.
      if (dv != dx) d2=(v-x)*dx/(dx-dv); //And the other.
      //Which of these two estimates of d shall we take? We will insist that they be within
      //the bracket, and on the side pointed to by the derivative at x:
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e; //Movement on the step before last.
      e=d;
      if (ok1 || ok2) { //Take only an acceptable d, and if both are acceptable, then take the smallest one.
	if (ok1 && ok2)
	  d=(fabs(d1) < fabs(d2) ? d1 : d2);
	else if (ok1)
	  d=d1;
	else
	  d=d2;
	if (fabs(d) <= fabs(0.5*olde)) {
	  u=x+d;
	  if (u-a < eps2 || b-u < eps2)
	    d=SETSIGN(eps1,xm-x);
	} else { //Bisect, not golden section.
	  d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	  //Decide which segment by the sign of the derivative.
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= eps1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SETSIGN(eps1,d);
      fu=(*f)(u);
      if (fu > fx) { //If the minimum step in the downhill direction takes us uphill, then we are done.
	*xmin=x;
	return fx;
      }
    }
    du=(*df)(u); //Now all the housekeeping, sigh.
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw);
      MOV3(w,fw,dw, x,fx,dx);
      MOV3(x,fx,dx, u,fu,du);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	MOV3(v,fv,dv, w,fw,dw);
	MOV3(w,fw,dw, u,fu,du);
      } else if (fu < fv || v == x || v == w) {
	MOV3(v,fv,dv, u,fu,du);
      }
    }
  }

  *xmin= x;                             //only get here if iteration limit is reached
  return fx;
#undef MOV3
}


/*
 * Multivariate function minimization.
 *
 * Input/Output
 *   th[1..n]: initial parameter values.
 *   dirini[1..n][1..n]: matrix with initial directions,
 *                       typically the n unit vectors
 *   n: length of th
 *   eps: relative tolerance to achieve convergence
 *   f: function to minimize (must take a vector of doubles as input)
 * Output
 *   th[1..n]: optimal parameter values
 *   dirini[1..n][1..n]: directions used in the last iteration
 *   iter: number of iterations used
 *   fret: value of the function at the optimum
 */
void minimize(double th[],
              double **dirini,
              int n,
              double eps,
              int *iter,
              double *fret,
              double (* const f)(double []),
              int itmax)
{
    int i;
    int j;
    int ibig;
    double del;
    double fth;
    double fthtt;
    double t;
    double *tht;
    double *thtt;
    double *dirinit;
    bool converged = false;

    //assert(dirini != NULL);
    //assert(iter != NULL);
    //assert(fret != NULL);
    //assert(f != NULL);

    tht     = dvector(1, n);
    thtt    = dvector(1, n);
    dirinit = dvector(1, n);

    /* Initial parameter and function values */
    *fret = (*f)(th);
    for (j = 1; j <= n; j++) {
        tht[j] = th[j];
    }

    for (*iter = 1; (*iter < itmax) && (converged == false); ++(*iter)) {
        fth = (*fret);
        ibig = 0;
        del = 0.0;

        /* Minimize along all directions */
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                dirinit[j] = dirini[j][i];
            }
            fthtt = (*fret);
            dirmin(th, dirinit, n, fret, f, itmax, eps);
            if (fabs(fthtt - (*fret)) > del) {
                del = fabs(fthtt - (*fret));
                ibig = i;
            }
        }

        /*
         * extrapolated point, average direction and function value
         * at extrapolated point
         */
        for (j = 1; j <= n; j++) {
            thtt[j] = 2.0 * th[j] - tht[j];
            dirinit[j] = th[j] - tht[j];
            tht[j] = th[j];
        }
        fthtt = (*f)(thtt);
        if (fthtt < fth) {
            t = 2.0 * (fth - 2.0 * (*fret) + fthtt) *
                sqrt(fth - (*fret) - del) - del * sqrt(fth - fthtt);
            if (t < 0.0) {
                dirmin(th, dirinit, n, fret, f, itmax, eps);
                for (j = 1; j <= n; j++) {
                    dirini[j][ibig] = dirini[j][n];
                    dirini[j][n] = dirinit[j];
                }
            }
        }

        if (2.0 * fabs(fth - (*fret)) <= eps * (fabs(fth) + fabs(*fret))) {
            converged = true;
        }
    }

    free_dvector(dirinit, 1, n);
    free_dvector(thtt, 1, n);
    free_dvector(tht, 1, n);
}


/*
 * Global variables communicate with f1dim.
 */
int ncom;
double *pcom;
double *xicom;
double (*nrfunc)(double []);


/* Must accompany dirmin. */
double f1dim(double x)
{
    int j;
    double f;
    double *xt;

    xt = dvector(1, ncom);
    for (j = 1; j <= ncom; j++) {
        xt[j] = pcom[j] + x * xicom[j];
    }
    f = (*nrfunc)(xt);
    free_dvector(xt, 1, ncom);

    return f;
}


/*
 * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n],
 * moves and resets p to where the function func(p) takes on a minimum along
 * the direction xi from p, and replaces xi by the actual vector displacement
 * that p was moved. Also returns as fret the value of func at the returned
 * location p. This is actually all accomplished by calling the routines
 * mnbrak() and univmin().
 */
void dirmin(double p[],
            double xi[],
            int n,
            double *fret,
            double (* const func)(double []),
            int itmax,
            double dirminEPS)
{
    int j;
    double xmin;
    double fx;
    double fb;
    double fa;
    double bx;
    double ax = 0.0;
    double xx = 1.0;

    //assert(fret != NULL);
    //assert(func != NULL);

    /* Set the global variables */
    ncom = n;
    pcom = dvector(1, n);
    xicom = dvector(1, n);
    nrfunc = func;
    for (j = 1; j <= n; j++) {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }

    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
    *fret = univmin(ax, xx, bx, f1dim, dirminEPS, &xmin, itmax);

    /* Construct the vector results to return */
    for (j = 1; j <= n; j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_dvector(xicom, 1, n);
    free_dvector(pcom, 1, n);
}


static double maxarg1;
static double maxarg2;

#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

/*
 * Given a function func, and given distinct initial points ax and bx, search
 * in the downhill direction (defined by the function as evaluated at the
 * initial points) and returns new points ax, bx, cx that bracket a minimum
 * of the function. Also returned are the function values at the three points,
 * fa, fb, and fc.
 */
void mnbrak(double *ax,
            double *bx,
            double *cx,
            double *fa,
            double *fb,
            double *fc,
            double (* const func)(double))
{
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

    const double GOLD = 1.618034;  /* default ratio by which successive intervals are magnified */
    const double GLIMIT = 100.0;   /* maximum magnification allowed for a parabolic-fit step */
    const double TINY = 1.0e-25;

    //assert(ax != NULL);
    //assert(bx != NULL);
    //assert(cx != NULL);
    //assert(fa != NULL);
    //assert(fb != NULL);
    //assert(fc != NULL);
    //assert(func != NULL);

    *fa = (*func)(*ax);
    *fb = (*func)(*bx);
    if (*fb > *fa) {
        double dum;

        /*
         * Switch roles of a and b so that we can go downhill
         * in the direction from a to b.
         */
        SHFT(dum, *ax, *bx, dum);
        SHFT(dum, *fb, *fa, dum);
    }
    *cx = (*bx) + GOLD * (*bx - *ax);  /* First guess for c */
    *fc = (*func)(*cx);

    /* Keep returning here until we bracket... */
    while (*fb > *fc) {
        double r;
        double q;
        double u;
        double ulim;
        double fu;

        /*
         * Compute u by parabolic extrapolation from a, b, c.
         * TINY prevents any division by zero.
         */
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
                (2.0 * SETSIGN(FMAX(fabs(q-r), TINY), q-r));
        ulim = (*bx) + GLIMIT * (*cx - *bx); /* Go no farther than this. */

        /* Test various possibilities */
        if ((*bx - u) * (u - *cx) > 0.0) {
            /* Parabolic u is between b and c: try it */
            fu = (*func)(u);
            if (fu < *fc) {
                /* Got a minimum between b and c */
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb) {
                /* Got a minimum between between a and u */
                *cx = u;
                *fc = fu;
                return;
            }
            /* Parabolic fit was no use, use default magnification */
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        }
        else if ((*cx - u) * (u - ulim) > 0.0) {
            /* Parabolic fit is between c and its allowed limit */
            fu = (*func)(u);
            if (fu < *fc) {
                SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx));
                SHFT(*fb, *fc, fu, (*func)(u));
            }
        }
        else if ((u - ulim) * (ulim - *cx) >= 0.0) {
            /* Limit parabolic u to maximum allowed value */
            u = ulim;
            fu = (*func)(u);
        }
        else {
            /* Reject parabolic u, use default magnification */
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        }

        /* Eliminate oldest point and continue */
        SHFT(*ax, *bx, *cx, u);
        SHFT(*fa, *fb, *fc, fu);
    }
#undef SHFT
}

