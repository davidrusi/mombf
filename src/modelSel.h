#ifndef MODELSEL_H
#define MODELSEL_H 1

#include <R.h>
#include <Rinternals.h>


//Define structures
struct marginalPars {
  int *sel;
  int *nsel;
  int *n;
  int *p;
  double *y;
  double *sumy2;
  double *x;
  double *XtX;
  double *ytX;
  double *m;  //Sinv * Xty   (needed by mom and emom)
  double **S;  //XtX + I/tau  (needed by mom and emom)
  int *method;
  int *B;
  double *alpha;    //prior for residual variance is IG(.5*alpha,.5*lambda)
  double *lambda;
  double *phi;      //residual variance
  double *tau;      //coefficients prior dispersion parameter
  int *r;           //MOM power parameter for prior on coefficients
  double *prDeltap; //For Binomial prior on model space, prDeltap is the prob of success
  double *parprDeltap; //For Beta-Binomial prior on model space, parprDeltap[0],parprDeltap[1] are the prior parameters
  int *logscale;
};

struct modavgPars {
  int *n;
  int *p1;
  int *p2;
  int *isbinary;  //isbinary==1 for probit regression (outcome stored in ybinary). isbinary==0 for linear model (outcome in y)
  int *ybinary;
  double *y;
  double *sumy2;
  double *x1;
  double *x2;
  double *XtX;   // t(x1) %*% x1
  double *ytX;   // t(y) %*% x1
  double *cholS2;
  double *S2inv;
  double *cholS2inv;
  double *colsumx1sq; //column sums for x1^2
  double *alpha;  //prior for resiual variance is IG(.5*alpha,.5*lambda)
  double *lambda;
  int *priorCoef; //1: pMOM prior; 2: peMOM prior
  int *r; //pMOM prior power parameter is 2*r
  double *tau1;
  double *tau2;
  int *priorTau1; //0: known; 1: IG(.5*atau1,.5*btau1)
  double *atau1;
  double *btau1;
  int *priorModel; //0 for uniform, 1 for binomial, 2 for Beta-binomial prior
  double *prModelpar; //For priorModel==1, 1st elem is prob of success. For priorModel==2, 1st and 2nd elem are Beta hyper-parameters
};

typedef double(*pt2margFun)(int *, int *, struct marginalPars *);  //pointer to function to compute marginal densities & prior prob (used for model selection)
typedef double(*pt2modavgPrior)(int *, int *, struct modavgPars *);  //pointer to function to compute prior model prob (used in model averaging routines)

//*************************************************************************************
//Setting prior & marginals
//*************************************************************************************

pt2margFun set_marginalFunction(int *prCoef, int *knownphi);
pt2margFun set_priorFunction(int *prDelta);
pt2modavgPrior set_priorFunction_modavg(int *priorModel);

//*************************************************************************************
//General Algebra
//*************************************************************************************

void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans);  //multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]; Return in ans[1..fi]
void Asymsel_x(double *A, int ncolA, double *x, int nsel, int *sel, double *ans); //similar but A formatted as vector. Use only selected elems in A and all elems in x[1..nsel]

void addct2XtX(double *ct, double *XtX, int *sel, int *nsel, int *p, double **V); //add constant to diagonal elem of XtX


//*************************************************************************************
// Model Averaging Routines
//*************************************************************************************

void set_modavgPars(struct modavgPars *pars, int *n, int *p1, int *p2, int *isbinary, int *ybinary, double *y, double *sumy2, double *x1, double *x2, double *XtX, double *ytX, double *cholS2, double *S2inv, double *cholS2inv, double *colsumx1sq, double *alpha, double *lambda, int *priorCoef, int *r, double *tau1, double *tau2, int *priorTau1, double *atau1, double *btau1, int *priorModel, double *prModelpar);

SEXP pmomLM_I(SEXP postModel, SEXP margpp, SEXP postCoef1, SEXP postCoef2, SEXP postPhi, SEXP postOther, SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP XtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar);
void pmomLM(int *postModel, double *margpp, double *postCoef1, double *postCoef2, double *postPhi, double *postOther, struct modavgPars *pars, int *niter, int *thinning, int *burnin, int *niniModel, int *iniModel, double *iniCoef1, double *iniCoef2, double *iniPhi, double *iniOthers, int *verbose);

void sample_latentProbit(double *y, double *res, double *sumres2, int *ybinary, double *linpred1, double *linpred2, struct modavgPars *pars);
void MHTheta1pmom(int *newdelta, double *newcoef, double *pinclude, int *resupdate, double *res, double *partialres, double *sumres2, double *sumpartialres2, int j, int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars);
void proposalpmom(double *propPars, double *m, double *S, double *phi, int *r, double *tau1, int *n, double *e, double *xj, double *m1, int *nu);
double pmomMargKuniv(double *y, double *x, double *sumy2, double *sumxsq, int *n, double *phi, double *tau, int *r, int *logscale);

void simTheta2(double *theta2, double *res, double *phi, struct modavgPars *pars);
double simPhipmom(int *nsel, int *curModel, double *curCoef1, double *curCoef2, double *ssr, struct modavgPars *pars);
double simTaupmom(int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars);


//*************************************************************************************
//General marginal density calculation routines
//*************************************************************************************

void set_marginalPars(struct marginalPars *pars, int *n,int *p,double *y,double *sumy2,double *x,double *XtX,double *ytX,int *method,int *B,double *alpha,double *lambda,double *phi,double *tau,int *r,double *prDeltap,double *parprDeltap, int *logscale);
void set_f2opt_pars(double *m, double **S, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel);
void set_f2int_pars(double *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale);

//Function by John D. Cook
double mvtexpect(const double* mu, const double* sigma, int n, int power, double dof); //mean of prod (x_i)^(2*power) when x_i ~ T_dof(mu,sigma). Set dof=-1 for N(mu,sigma)


//*************************************************************************************
// Model Selection Routines
//*************************************************************************************

void modelSelectionGibbs(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars);
void modelSelectionGibbs2(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars);
void greedyVarSelC(int *postMode, double *postModeProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars);
void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew);

// Priors on Model Space (always return on log scale)
double unifPrior(int *sel, int *nsel, struct marginalPars *pars);
double unifPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);
double binomPrior(int *sel, int *nsel, struct marginalPars *pars);
double binomPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);
double betabinPrior(int *sel, int *nsel, struct marginalPars *pars);
double betabinPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);

//*************************************************************************************
// Product MOM routines
//*************************************************************************************

double f2opt_mom(double *th);
double fmomNegC_non0(double *th, double *m, double **S, double *phi, double *tau, int *r, int *n, int *nsel);
void fppmomNegC_non0(double **ans, double *th, double **S, double *phi, double *tau, int *r, int *n, int *nsel);
void momIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *n, int *nsel, double *m, double **S, double *detS, double *phi, double *tau, int *r, int *logscale);

double rsumlogsq(double *th, int *r, int *nsel);  //compute r*sum(log(th^2))
double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double MC_mom(double *m,double **Sinv,int *r,int *nsel, int *B);  //MC evaluation of E(prod(th^2r)) for th ~ N(m,Sinv)
double MC_mom_T(double *m,double **Sinv,int *nu,int *r,int *nsel, int *B); //MC evaluation of E(prod(th^2r)) for th ~ T_nu(m,Sinv)

//*************************************************************************************
// Product iMOM routines
//*************************************************************************************

double f2opt_imom(double *th);
double fimomNegC(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
double fimomNegC_non0(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void fppimomNegC_non0(double **ans, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *logscale);

double f2opt_imomU(double *th);
double fimomUNegC_non0(double *th, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel);
void imomUIntegralApproxC(double *ILaplace, double *thopt, int *sel, int *nsel, int *n, int *p, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *logscale);

SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale);
double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *B);

double f2int_imom(double phi);
SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda);
double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// Product MOM routines
//*************************************************************************************

double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// Product eMOM routines
//*************************************************************************************

double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

#endif /* MODELSEL_H */

