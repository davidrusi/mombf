/***********************************************************
 Basic statistical, input-output and matrix manipulation

 Authors: Peter Mueller, Stephen Morris, David Rossell
***********************************************************/

#ifndef CSTAT_H
#define CSTAT_H 1

#include <math.h>
#include <stdio.h>
#include <stdbool.h>


#if !defined(M_PI)
#define M_PI (3.1415926535897932385)
#endif

#if !defined(M_PI_2)
#define M_PI_2 (1.570796326794896619231321691640)
#endif

#if !defined(LOG_M_2PI)
#define LOG_M_2PI (1.8378770664093453)
#endif

#if !defined(SQ_M_PI_2)
#define SQ_M_PI_2 (2.5066282746310002416123552)
#endif

#if !defined(LOG_M_PI)
#define LOG_M_PI (1.1447298858494)
#endif



/**************************************************************/
/* Functions interfacing with R                               */
/**************************************************************/

extern "C" {

  //Non-local prior sampling
  SEXP rnlpPostCI_lm(SEXP niter, SEXP burnin, SEXP thinning, SEXP y, SEXP x, SEXP p, SEXP r, SEXP tau, SEXP a_phi, SEXP b_phi, SEXP prior);
  SEXP rnlpCI(SEXP niter, SEXP burnin, SEXP thinning, SEXP m, SEXP V, SEXP p, SEXP r, SEXP tau, SEXP prior);

  //Truncated multivariate Normal sampling
  SEXP rnorm_truncMultCI(SEXP n, SEXP ltrunc, SEXP rtrunc, SEXP m, SEXP s);  //R interface for rnorm_truncMult
  SEXP rtmvnormCI(SEXP n, SEXP mu, SEXP Sigma, SEXP lower, SEXP upper, SEXP within, SEXP method); //R interface for rtmvnorm
  SEXP rtmvnormProdCI(SEXP n, SEXP mu, SEXP Sigma, SEXP k, SEXP lower, SEXP upper, SEXP is_low_trunc, SEXP is_up_trunc, SEXP burnin); //R interface for rtmvnormProd

}

/**************************************************************/
/* Functions to compute means & variances                     */
/**************************************************************/

double meani(const int *x, int lim);
double wmeani(const int *x, int lim, const double *w);
double meanx(const double *x, int lim); //mean of x[0..lim]
double wmeanx(const double *x, int lim, const double *w);
double vari(const int *x, int lim, bool unbiased);
double wvari(const int *x, int lim, const double *w);
double varx(const double *x, int lim, bool unbiased);
double wvarx(const double *x, int lim, const double *w);
double cv(const double *x, int ini, int fi); //coefficient of variation i.e. SD/mean
double cvinv(const double *x, int ini, int fi); //coefficient of variation of 1/x

void colMeans(double *m, const double *x, int nrow, int ncol);
void colVar(double *m, const double *x, int nrow, int ncol);
void colCV(double *cv, const double *x, int nrow, int ncol);
void colCVinv(double *cv, const double *x, int nrow, int ncol); //CV of 1/x


/************************************************************************
                         BASIC BAYESIAN MODELS
************************************************************************/

void nn_bayes(double *mpo, double **Spo, double **Spo_inv, int p, double r1, double *mpr, double **Spr_inv, double r2, double *y, double **Slik_inv);  //Posterior of multiv normal mean with normal prior
void nn_bayes_rand(double *theta, int p, double r1, double **Spr_inv, double *mpr, double r2, double **Slik_inv, double *y); //Single draw from posterior of multiv normal mean with normal prior
double nn_integral(const double *x, const double *rx, double **Vxinv, const double *detVx, const double *mpr, const double *rpr, double **Vprinv, const double *detVpr, const int *p, const int *logscale); //Normal-Normal integral (useful to compute Bayes factors etc.)

void lm(double *b, double **XtX, double **invXtX, double *Xty, double *s, double *ypred, const double *y, double **X, const int *n, const int *p, const int *useXtX); //classical multiple linear regression
void lmbayes(double *bpost, double *spost, double *b, double **Vb, double *a_s, double *b_s, double **XtX, double **invXtX, double *Xty, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior, double *nu0, double *s0); //Bayesian multiple linear regression
void lmbayes_knownvar(double *bpost, double *b, double **Vb, double **XtX, double **invXtX, double *Xty, double *sigma, int *B, double *y, double **X, int *n, int *p, int *useXtX, double *mpr, double **Spr_inv, double *tauprior); //same as lmbayes with known variance sigma^2


/**************************************************************/
/* Input/output functions (interface)                         */
/**************************************************************/

FILE *openIn(const char *filename);
FILE *openOut(const char *filename);
//void scanFloat(char *, float *);
//void scanDouble(char *, double *);
//void scanInt(char *, int *);
//void fscanDouble(FILE *, char *, double *);
//void fscanInt(FILE *, char *, int *);
//void scanLong(char *, long *);
// 
//void scanFloatArray(char *, float *, int);
//void scanArray(char *, float *, int);
//void scanDoubleArray(char *, double *, int);
//void scanString(char *txt, char *s, int n);
//void fscanString(FILE *, char *txt, char *s, int n);
//void fscanDoubleArray(FILE *, double *, int);
//void scanDoubleMatrix(char *, double **, int, int);
//void fscanDoubleMatrix(FILE *ifile, double **x, int r, int c);
//void scanIntArray(char *, int *, int);
//void fscanIntArray(FILE *ifile, int *x, int n);

void writeInt(int);
void writeLong(long i);
void writeFloat(float);
void writeDouble(double);

void writeIntArray(int *, int, int);
void fwriteIntArray(FILE *, int *, int, int);
void fwriteIntMatrix(FILE *f, int **x, int rows, int cols);
void writeIntMatrix(int **x, int rows, int cols);
void writeDoubleArray(double *, int, int);
void writeDoubleMatrix2(double **, int , int);
void fwriteDoubleArray(FILE *, double *, int, int);
void fwriteDoubleMatrix2(FILE *, double **, int , int);
void writeDoubleMatrix(double **, int, int);
void writeFloatArray(float *, int, int);
void writeArray(float *, int, int); 


/**************************************************************/
/* Debug messages etc. (mess)                                 */
/**************************************************************/

void errorC(const char *module, const char *msg, int nr);
void err_msg(const char *fct, const char *txt, int n1, int n2, int n3);
void fserror(const char *proc, const char *act, const char *what);
void nrerror(const char *proc, const char *act, const char *what);

/**************************************************************/
/* Memory allocation                                          */
/**************************************************************/



char *charvector(int ,int);
float   *vector(int, int);
double  *dvector(int, int);
double  **dmatrix(int, int, int, int);
//double  ***darray_3(int, int);
//double ***darray3(int n, int p, int q);
double ***darray3(int n1, int n2, int n3);   //allocate 3-way double array [0..n1-1][0..n2-1][0..n3-1]
int     *ivector(int, int);
int     **imatrix(int, int, int, int);
//int ***iarray_3(int lo, int hi);
//int ***iarray3(int p1, int p2, int p3);
int ***iarray3(int n1, int n2, int n3);   //allocate 3-way int array [0..n1-1][0..n2-1][0..n3-1]

void free_charvector(char *v,int nl,int nh);
void free_vector(float *, int, int);
void free_dvector(double *, int, int);
void free_ivector(int *, int, int);
void free_dmatrix(double **, int, int, int, int);
void free_imatrix(int **, int, int, int, int);
void free_darray3(double ***a, int n1, int n2, int n3);
void free_iarray3(int ***a, int n1, int n2, int n3);

/**************************************************************/
/* Mathematical functions                                     */
/**************************************************************/

double gamln(double*);  //log-Gamma function
double gamln1(double*);  //auxiliary function called by gamln
double lfact(int);  //log-factorial
double ldoublefact(double x);  //log-double factorial(x)
double digamma(double x);                          //from S Poetry (by Patrick J. Burns)
double trigamma(double x);
double polygamma(double x, long n, double low, double high, long terms, double nfact); //from S Poetry
double lnbeta(double a, double b); //log of Beta function
double betacf(double a, double b, double x); //continued fraction for incomplete Beta function
double lnchoose(double n, double k);
double choose(double n, double k);

double logit(double x);
double ilogit(double x);

double dsign(double x);  //returns 1.0 if x>=0, -1.0 if x<0
double isign(int x); //returns 1.0 if x>0, 0 if x==0, -1.0 if x<0

/**************************************************************/
/* Vector algebra (vector)                                    */
/**************************************************************/

void grid(double x0, double xn, int n, double *x);
void rA(double r, double **A, double **B, int rowini, int rowfi, int colini, int colfi);  //matrix*scalar
void A_plus_B(double **A, double **B, double **C, int rowini, int rowfi, int colini, int colfi); //matrix + matrix
void rA_plus_sB(double r, double **A, double s, double **B, double **C, int rowini, int rowfi, int colini, int colfi); //matrix*scalar + matrix*scalar
void rAx_plus_sBy(double r, double **A, const double *x, double s, double **B, const double *y, double *z, int rowini, int rowfi, int colini, int colfi); //scalar*matrix*vector + scalar*matrix*vector
void Ax_plus_y(double **A, const double *x, const double *y, double *z, int ini, int fi); //matrix*vector+vector
void xA(const double *x, double **A, double *z, int ini, int fi);  //Multiply vector * matrix
void Ax(double **A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi);  //matrix * vector
void Avecx(const double *A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi); //same but A is in vector format
void Atvecx(const double *A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi); //same for A' (row/col indexes refer to A')
double xtAy(const double *x, double **A, const double *y, int ini, int fi); //t(vector)*matrix*vector

double quadratic_xtAx(const double *x, double **A, int ini, int fi); //t(vector)*matrix*vector for quadratic forms (A symmetric)
double quadratic_xseltAselxsel(const double *x, const double *A, const int *ncol, const int *nsel, const int *sel); // same but A is formatted as vector & only a subset of x, A is to be used
double quadratic_xtAselx(const double *x, const double *A, const int *ncolA, const int *nsel, const int *sel); //same but subset is only for A
double quadratic_xseltAxsel(const double *x, double **A, int ini, const int *nsel, const int *sel); //same but subset is only for x

void Atx(double **A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi); //t(matrix)*vector
void AtB(double **A, int rowiniA, int rowfiA, int coliniA, int colfiA, double **B, int rowiniB, int rowfiB, int coliniB, int colfiB, double **C); //t(matrix)*matrix, stored in C
void AvectBvec(double *A, int nrowA, int ncolA, double *B, int nrowB, int ncolB, double **C); //same but input as 0-indexed vector
void a_plus_b(const double *a, const double *b, double *c, int ini, int fi); //Vector sum i.e. c[i]=a[i]+b[i]
void a_prod_b(const double *a, const double *b, double *c, int ini, int fi); //Vector prod i.e. c[i]=a[i]*b[i]
void a_prod_b_sel(const double *a, const double *b, double *c, const int *lengtha, const int *nsel, const int *sel); //same but only using indexes in sel
void a_zero(double *, int); //Set vector to zero
void R_zero(double **, int, int); //Set matrix to zero
void ddiag(double **A, int ini, int fi); //Diagonal matrix

int iabs(int x);   //absolute value of an integer
int imax_xy(int x, int y);
int imin_xy(int x, int y);
double max_xy(double x, double y);
double min_xy(double x, double y);
void minvec(const double *x, int ini, int fi, double *xmin, int *minpos); //min of a vector and position at which min occurs
void maxvec(const double *x, int ini, int fi, double *xmax, int *maxpos); //max of a vector and position at which max occurs

void choldc(double **a, int n, double **aout);   //Cholesky decomposition
void choldc_inv(double **a, int n, double **aout); //Inverse of Cholesky decomposition
double choldc_det(double **chols, int n); //Determinant of a symmetric def+ using its Cholesky decomp
void inv_posdef(double **a, int n, double **aout); //Inverse of a symmetric, positive definite matrix
void inv_posdef_upper(double **a, int n, double **aout); //Same but only returns upper triangular elements
void invdet_posdef(double **a, int n, double **aout, double *det_a); //Inverse and determinant of positive def matrix
void inv_posdef_chol(double **invchol, int n, double **aout); //Inverse given cholesky decomposition

void ludc(double **a, int n, int *indx, double *d); //LU decomposition (renamed routine ludcmp from NR)
void lu_solve(double **a, int n, const int *indx, double b[]); //Solve A*x=b (renamed routine lubksb from NR)
void lu_inverse(double **a, int n, double **aout); //Inverse of A[1..n][1..n]
double lu_det(double **a, int n); //Determinant of A[1..n][1..n]

int dcompare(const void *a, const void *b);               
void dvecsort(double *v, int size);                           //sort a vector using qsort from stdlib
void dindexsort(double *x, int *index, int ilo, int ihi, int incr); //sort a vector of indexes using self-written quicksort routine
void iindexsort(int *x, int *index, int ilo, int ihi, int incr); //like dindexsort but for integers


/**************************************************************/
/* Random sampling                                            */
/**************************************************************/

void samplei_wr(int *x, int popsize, int n); //sample wo replacement from a vector of integers
void sampled_wr(double *x, int popsize, int n); //same for vector of doubles

/**************************************************************/
/* Probability distributions                                  */
/**************************************************************/

// Several
void setseed(long, long);
int rdisc(const double *probs, int nvals);
double gamdev(double);
int rbinomial(int , double );
double dbinomial(int x, int n, double p, int logscale);
double dnegbinomial(int x, double r, double p, int logscale);
void rmultinomial(int ndraws, int ncells, const double *pr, int *x);
double bbPrior(int k, int p, double alpha, double beta, int logscale);

// Uniform
double runif();
double dunifC(double x, double a, double b);
int runifdisc(int min, int max);

// Beta-Dirichlet
double rbetaC(double , double );
double pbetaC(double x, double pin, double qin); //quantile from a Beta(pin,qin)
void rdirichlet(double *w, const double *alpha, const int *p);
double ddirichlet(const double *w, double *alpha, const int *p); //Dirichlet density

// Normal
double dnormC(double y, double m, double s, int logscale); //density of Normal(m,s^2)
double dnormC_jvec(const double *y, int n, double m, double s, int logscale); //joint density of y[0]...y[n-1] under Normal(m,s^2), i.e. returns scalar
double dmvnormC(const double *y, int n, const double *mu, double **cholsinv, double det, int logscale); //density of multivariate Normal
double	qnormC(double cdf, double m, double s);  //quantile from Normal(m,s^2)
double	pnormC(double y, double m, double s);  //cdf of Normal(m,s^2)
double rnormC(double mu, double s); //draw from univariate Normal(mu,s^2)
void rmvnormC(double *y, int n, const double *mu, double **chols); //draw from multivariate Normal

// Truncated Normal
double rnorm_trunc(double ltrunc, double rtrunc, double m, double s); //draw trunc Normal given two truncation points
double rnorm_trunc_prob(double lprob, double rprob, double m, double s); //draw trunc Normal given trunc probs
void rnorm_truncMult(double *y, double *pdfy, int *n, double *ltrunc, double *rtrunc, int ntrunc, double *m, double *s); //n draws univ trunc Normal given multiple truncation points

void rtmvnorm(double *ans, int n, int p, double *mu, double **Sigma, double *lower, double *upper, int within, int method); //n draws multiv trunc Normal given rectangular constraints
void rtmvnormMH(double *ansortho, double *paccept, int n, int p, double *alpha, double **D, double **K, double det, double *lower, double *upper, int within); //Indep prop MH
void rtmvnormOutside(double *ans, int n, int p, double *alpha, double **D, double *lower, double *upper); //Gibbs sampling, truncation outside interval
void rtmvnormOutside_Gibbs(double *z, double *Dj, double *alpha, double **D, int p, double *lower, double *upper); //Single Gibbs sampling step
void rtmvnormWithin(double *ans, int n, int p, double *alpha, double **D, double *lower, double *upper);  //Gibbs sampling, truncation within interval
void rtmvnormProp(double *z, double *propdens, int p, double *alpha, double **D, double *lower, double *upper, int within); //proposal for multiv trunc Normal

void rtmvnormProd(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, double upper, int is_low_trunc, int is_up_trunc, int burnin); //multiv trunc Normal given constraint on product
void rtmvnormProd_lowup(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, double upper, int burnin); //lower and upper truncation on product
void rtmvnormProd_low(double *ans, int n, int p, double *mu, double **Sinv, int k, double lower, int burnin); //only lower truncation on product
void rtmvnormProd_up(double *ans, int n, int p, double *mu, double **Sinv, int k, double upper, int burnin); //only upper truncation on product

// Moments
double mnorm(double order, double m, double sd); //raw moment of N(m,sd) of order "order"

double mvtexpect(double* mu, double** sigma, int n, int power, double dof); //mean of prod (x_i)^(2*power) when x_i ~ T_dof(mu,sigma). Set dof=-1 for N(mu,sigma). Written by John Cook
double mvtexpect_vec(double* mu, double* sigma, int n, int power, double dof); //same but input args are 0-indexed vectors

// T Student
double dtC(double y, double mu, double s, int nu); //density of t with nu df
double dtmixC(double y, const double *mu, const double *s, const double *probs, int nu, int ncomp, int logscale); //density of t_nu(mu[i],s[i]^2) mixtures with ncomp components
double dmvtC(const double *y, int n, const double *mu, double **cholsinv, double det, int nu, int logscale); //density of multivariate t
double rtC(int nu); //draw from univariate t with nu degrees of freedom
double rtmixC(const double *mu, const double *s, const double *probs, int nu, int ncomp); //draw from mixture of t_nu(mu[i],s[i]^2)
double rt_trunc(int nu, double ltrunc, double rtrunc); //draw from truncated t given trunc points
double rt_trunc_prob(int nu, double lprob, double rprob);  //draw from truncated t given trunc probs
double qtC(double p, int nu);  //quantile from t-Student with nu degrees of freedom
double ptC(double x, int nu);  //CDF of t-Student with nu degrees of freedom
void rmvtC(double *y, int n, const double *mu, double **chols, int nu); //draw from multivar T with nu degrees of freedom

// Gamma & Inverse gamma
double rgammaC(double a, double b); //a: shape; b: location; mean=a/b
double dgammaC(double x, double a, double b); //a: shape; b: location; mean=a/b
double dinvgammaC(double x, double a, double b, int logscale); //a: shape; b: location; mean of x= b/(a-1)

// Non-local prior densities
double dmom(double y, double m, double tau, double phi, int r, int logscale); //Univariate MOM prior (power is 2*r)
void dmomvec(double *y, int n, double m, double tau, double phi, int r, int logscale); //Multivariate MOM prior
double dimom(double y, double m, double tau, double phi, int logscale); //Univariate iMOM prior
void dimomvec(double *y, int n, double m, double tau, double phi, int logscale); //Multivariate iMOM prior
double demom(double y, double tau, double phi, int logscale); //Univariate eMOM prior
void demomvec(double *y, int n, double tau, double phi, int logscale); //Multivariate eMOM prior

// Non-local prior density derivatives
void dmomgrad(double *ans, int *n, double *th, double *logphi, double *tau); //Gradient of log-pMOM density wrt th
void dmomhess(double *ans, int *n, double *th, double *logphi, double *tau); //Hessian of log-pMOM density wrt th
void dmomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda); //Grad log-pMOM + log-IG wrt (th,logphi)
void dmomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda); //Hess log-pMOM + log-IG wrt (th,logphi)

void dimomgrad(double *ans, int *n, double *th, double *logphi, double *tau); //Gradient of log-piMOM density wrt th
void dimomhess(double *ans, int *n, double *th, double *logphi, double *tau); //Hessian of log-piMOM density wrt th
void dimomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda); //Grad log-piMOM + log-IG wrt (th,logphi)
void dimomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda);//Hess log-piMOM + log-IG wrt (th,logphi)

void demomgrad(double *ans, int *n, double *th, double *logphi, double *tau); //Gradient of log-peMOM density wrt th
void demomhess(double *ans, int *n, double *th, double *logphi, double *tau); //Hessian of log-peMOM density wrt th
void demomiggrad(double *ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda); //Grad log-peMOM + log-IG wrt (th,logphi)
void demomighess(double **ans, int *n, double *th, double *logphi, double *tau, double *alpha, double *lambda);//Hess log-peMOM + log-IG wrt (th,logphi)


// Posterior sampling under non-local priors
void rnlpPost_lm(double *ans, int niter, int burnin, int thinning, double *y, double *x, int n, int p, int r, double tau, double a_phi, double b_phi, int prior); //NLP posterior samples under linear model
void rnlp(double *ans, int niter, int burnin, int thinning, double *m, double *Vvec, int p, int r, double tau, int prior); //NLP samples based on mean & covar
void rnlp_Gibbs(double *th, int p, double *m, double **cholS, double **K, double *tau, double *phi, int r, int prior); //single Gibbs update 
void rnlp_Gibbs_multiple(double *th, double *thini, int p, double *m, double **cholS, double **K, double *tau, int r, int prior, int niter, int burnin, int thinning); //multiple updates
double pen_mom(double *th, double *phi, double *tau, int r); //MOM penalty (th^2 / (phi*tau))^r
double pen_emom(double *th, double *phi, double *tau, int logscale); //eMOM penalty exp(-sqrt(2)*tau*phi/th^2)
double pen_imom(double *th, double *phi, double *tau, int logscale); //iMOM penalty dimom(th,tau*phi) / dnorm(th,0,tau*phi)
double invpen_imom_newton(double *lambda, double *phi, double *tau); //inverse of iMOM penalty using Newton's method
double invpen_imom_sandwich(double *lambda, double *phi, double *tau); //inverse of iMOM penalty using Sandwich method

/* More random variate stuff (dcdflib, from CMU statlib "www.stat.cmu.edu") */
double fifdint(double);
void cdfnor(int*, double*, double*, double*, double*, double*, int*, double*);
double spmpar(int*);
void cumnor(double*, double*, double*);
double dinvnr(double *p, double *q);
double stvaln(double*);
double devlpl(double [], int*, double*);
extern int ipmpar(int*);                      /* code in ipmpar.c */

/*even more stuff (ranlib) */
extern double genunf(double low, double high);
extern double gengam(double a, double r);
extern double sgamma(double a);
extern double snorm(void);
double fsign(double num, double sign);
extern double sexpo(void);
extern long mltmod(long a, long s, long m);
extern double ranf(void);
extern void gscgn(long getset, long *g);
extern void setall(long iseed1, long iseed2);  /* code in com.c */
extern void initgn(long isdtyp);               /* code in com.c */
extern long ignlgi(void);                      /* code in com.c */
extern void inrgcm(void);                      /* code in com.c */


/**************************************************************/
/* Integration                                                */
/**************************************************************/

double midpnt(double (* const func)(double), double a, double b, int n); //nth stage refinement of integral of func from a to b (evenly spaced in x)
double midinf(double (* const func)(double), double aa, double bb, int n); //nth stage refinement of integral of func from aa to bb (evenly spaced in 1/x)
double qromo(double (* const func)(double), double a, double b, double (* const choose)(double(* const)(double), double, double, int)); //Romberg integr on open interval (a,b)


/**************************************************************/
/* Interpolation, extrapolation and splines                   */
/**************************************************************/

void polint(double xa[], double ya[], int n, double x, double *y, double *dy); //interpolates via polynomials
double bspline_singlex(double x, int j, int degree, const double *knots); //jth B-spline basis eval at single value x
void bspline(double **W, const double *x, const int *nx, const int *degree, const double *knots, const int *nknots); //B-spline basis eval at vector of values x
void bspline_vec(double *W, const double *x, const int *nx, const int *degree, const double *knots, const int *nknots); //same as bspline but returns a vector, so that it can be called from R
void mspline(double **W, const double *x, const int *nx, const int *degree, const double *knots, const int *nknots); //M-spline basis eval at vector of values x
void mspline_vec(double *W, const double *x, const int *nx, const int *degree, const double *knots, const int *nknots); //same as mspline but returns a vector, so that it can be called from R


/**************************************************************/
/* Function optimization                                      */
/**************************************************************/

double univmin(double ax, double bx, double cx, double (* const f)(double), double tol, double *xmin, int itmax); //univariate minim
double dunivmin(double ax, double bx, double cx, double (* const f)(double), double (* const df)(double), double tol, double *xmin, int itmax);
void minimize(double th[], double **dirini, int n, double ftol, int *iter, double *fret, double (* const f)(double []), int itmax);//multivar minim
void dirmin(double p[], double xi[], int n, double *fret, double (* const func)(double []), int itmax, double dirminEPS); //minim in 1 direction
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double (* const func)(double)); //find bracketing triplets

#endif /* CSTAT_H */

