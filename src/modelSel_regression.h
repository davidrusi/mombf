#ifndef MODELSEL_H
#define MODELSEL_H 1

//#include <R.h>
//#include <Rinternals.h>
#include <list>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "crossprodmat.h"
#include "covariancemat.h"
#include "glm.h"
#include "survival.h"
#include "modselFunction.h"
#include "Polynomial.h"


/*
 * Function Prototypes
 */


//test function for debugging
SEXP testfunctionCI(SEXP x);

//Auxiliary functions
SEXP bsplineCI(SEXP x, SEXP degree, SEXP knots);

SEXP eprod_I(SEXP m, SEXP S, SEXP n, SEXP power, SEXP dof);

  //Posterior sampling for parameters
SEXP pmomLM_I(SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP SXtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar);

  //Model selection
SEXP modelSelectionEnumCI(SEXP Snmodels, SEXP Smodels, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose);


SEXP modelSelectionGibbsCI(SEXP SpostModeini, SEXP SpostModeiniProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose);


SEXP greedyVarSelCI(SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose);

  //Non-local prior marginal likelihoods
SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Sngroups, SEXP Snvaringroup);
SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP Sngroups, SEXP Snvaringroup);

SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Sngroups, SEXP Snvaringroup);
SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP Sngroups, SEXP Snvaringroup);



/*****************************************************************************************************************
  typedefs
*****************************************************************************************************************/

typedef double(*pt2margFun)(int *, int *, struct marginalPars *);  //pointer to function to compute marginal densities & prior prob (used for model selection)
typedef double(*pt2modavgPrior)(int *, int *, struct modavgPars *);  //pointer to function to compute prior model prob (used in model averaging routines)
typedef std::list<int*> intptrlist; //list where each element is a pointer to an integer
typedef std::vector<int*> intptrvec; //vector where each element is a pointer to an integer





/*****************************************************************************************************************
  Define structures
*****************************************************************************************************************/


// Structure containing parameters needed to calculate the integrated likelihood
struct marginalPars {
  int *family;
  int *priorcode;
  int *sel;
  int *nsel;
  int *n;        //number of observations
  int *nuncens;  //number of uncensored observations
  int *p;
  double *y;
  int *uncens;
  double *sumy2;
  double *sumy;
  double *sumlogyfact; //sum(log(y!)), used in Poisson regression
  double *x;
  double *colsumsx;   //column sums of x
  crossprodmat *XtX;  //t(x) %*% x using all observations
  crossprodmat *XtXuncens; //t(x) %*% x using uncensored observations
  covariancemat *V0inv;  // covariance matrix for coef and groups priors
  double *ytX;             //t(x) %*% y using all observations
  double *ytXuncens;       //t(x) %*% y using uncensored observations
  double *m;  //Sinv * Xty   (needed by mom and emom)
  double **S;  //XtX + I/tau  (needed by mom and emom)
  int *method; //method==0 for Laplace; method==1 for Monte Carlo; method==2 for ALA (method== -1 for automatic choice)
  int *adjoverdisp; //Used only for ALA. 0 for no adjustment; 1 to estimate overdispersion from intercept-only model, 2 to estimate from residuals
  int *hesstype; //for asymmetric Laplace residuals hess=1 means using asymptotic hessian, hess=2 means using diagonal adjustment to asymp hessian
  int *optimMethod; //optimization method to find mode
  int *optim_maxit; //maximum number of iterations
  int *usethinit; //usethinit==1 tells optimization algorithms to store the optimal model parameters at thinit; usethinit==2 to initialize at thinit upon entry and store optimal value at thinit upon exit; usethinit==3 to initialize at thinit but not to change thinit; usethinit==0 to ignore thinit
  double *thinit; //thinit[sel[j]] stores initial values for model parameters to be used by optimization algorithms
  int *B;      //number of Monte Carlo samples
  double *alpha;    //prior for residual variance is IG(.5*alpha,.5*lambda)
  double *lambda;
  int *knownphi; //should dispersion parameter be considered known, e.g. error var in Gaussian regression, or phi=1 in logistic/poisson regression
  double *phi;      //residual variance
  double *tau;      //dispersion parameter in prior for regression coefficients. Also used to store the penalty parameter when using an info criteria rather than a prior
  double *taugroup; //prior dispersion parameter on grouped coefficients, e.g. the block Zellner prior is prod_j N(delta_j; 0, (taugroup/ncol(X_j)) (X_j'X_j)^{-1})
  double *taualpha; //dispersion parameter in prior for asymmetry parameter in two-piece Normal or two-piece Laplace residuals
  double *fixatanhalpha; //fixed value for asymmetry parameter (usedful for quantile regression at fixed quantile levels)
  int *r;           //MOM power parameter for prior on coefficients
  double *prDeltap; //For Binomial prior on model space, prDeltap is the prob of success. For complexity prior, the power parameter in the exponential
  double *parprDeltap; //For Beta-Binomial prior on model space, parprDeltap[0],parprDeltap[1] are the prior parameters
  double *prConstrp; //idem for prior on number of included groups under hierarchical constraints
  double *parprConstrp;
  int *logscale;
  double *offset;
  int *groups;  //group that each variable belongs to
  int *isgroup; //isgroup[j]==1 indicates that variable j is in a group
  int *ngroups; //total number of groups
  int *ngroupsconstr; //number of groups that have a hierarchical constraint
  int *nvaringroup; //number of coefficients in group[0],...,group[ngroups-1]
  int *nconstraints; //number of constraints in group[0],...,group[ngroups-1]
  int *ninvconstraints; //number of inverse constraints (number of groups depending on group[0],...,group[ngroups-1])
};



//Structure containing parameters required to compute Bayesian model averaging
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
  crossprodmat *XtX; // t(x1) %*% x1
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



void testfunction(double *x);




//*************************************************************************************
//Setting prior & marginals
//*************************************************************************************

int mspriorCode(int *prCoef, int *prGroup, struct marginalPars *pars);
pt2margFun set_marginalFunction(struct marginalPars *pars);
pt2margFun set_priorFunction(int *prDelta, int *prConstr, int *family);
pt2modavgPrior set_priorFunction_modavg(int *priorModel);

//*************************************************************************************
//General Algebra
//*************************************************************************************

void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans);  //multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]; Return in ans[1..fi]
void Asel_x(double *A, int ncolA, double *x, int nsel, int *sel, double *ans); //multiply symmetric A[1..ncolA][1..ncolA] (formatted as vector) with x[1..nsel].

void addct2XtX(double *ct, crossprodmat *XtX, int *sel, int *nsel, int *p, double **V); //add constant to diagonal elem of XtX



//*************************************************************************************
// Model Averaging Routines
//*************************************************************************************

void set_modavgPars(struct modavgPars *pars, int *n, int *p1, int *p2, int *isbinary, int *ybinary, double *y, double *sumy2, double *x1, double *x2, crossprodmat *XtX, double *ytX, double *cholS2, double *S2inv, double *cholS2inv, double *colsumx1sq, double *alpha, double *lambda, int *priorCoef, int *r, double *tau1, double *tau2, int *priorTau1, double *atau1, double *btau1, int *priorModel, double *prModelpar);

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

void set_marginalPars(struct marginalPars *pars, int *family, int *n,int *nuncens,int *p,double *y,int *uncens,double *sumy2,double *sumy,double *sumlogyfact,double *x,double *colsumsx,crossprodmat *XtX,double *ytX,int *method,int *adjoverdisp,int *hesstype,int *optimMethod,int *optim_maxit,int *usethinit,double *thinit,int *B,double *alpha,double *lambda,int *knownphi,double *phi,double *tau,double *taugroup,double *taualpha, double *fixatanhalpha, int *r,double *prDeltap,double *parprDeltap, double *prConstrp,double *parprConstrp, int *logscale, double *offset, int *groups, int *isgroup, int *ngroups, int *ngroupsconstr, int *nvaringroup, int *nconstraints, int *ninvconstraints, crossprodmat *XtXuncens, double *ytXuncens);
void set_f2opt_pars(double *m, double **S, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel);
void set_f2int_pars(crossprodmat *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale);



//*************************************************************************************
// Model Selection Routines
//*************************************************************************************

void modelSelectionEnum(int *postMode, double *postModeProb, double *postProb, int *nmodels, int *models, int *prDelta, int *prConstr, int *verbose, struct marginalPars *pars);
void modelSelectionGibbs(int *postSample, double *margpp, int *postMode, double *postModeProb, double *postProb, int *prDelta, int *prConstr, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *includevars, intptrvec *constraints, intptrvec *invconstraints, int *verbose, struct marginalPars *pars);
void greedyVarSelC(int *postMode, double *postModeProb, int *prDelta, int *prConstr, int *niter, int *ndeltaini, int *deltaini, int *includevars, intptrvec *constraints, intptrvec *invconstraints, int *verbose, struct marginalPars *pars);

void update_postMode(int *postMode, int nselnew, int *selnew, int p, int family);
bool checkConstraints(int *addgroups, int *naddgroups, int *dropgroups, int *ndropgroups, intptrvec *constraints, int *nconstraints, intptrvec *invconstraints, int *ninvconstraints, int *groups, int *nvaringroup, int *sel, int *nsel);
void sel2selnew(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups, int *nvaringroup, int *firstingroup);
void findselgroups(double *nvarinselgroups, double *firstingroup, double *nselgroups, double *selgroups, int *sel, int *nsel, int *nvaringroup, int *ngroups);
void nselConstraints(int *ngroups0, int *ngroups1, int *sel, int *nsel, int *group, int *nconstraints, int *nvaringroup);
void countConstraints(int *nconstraints, intptrvec *constraints, int *ninvconstraints, intptrvec *invconstraints, int *ngroupsconstr, int *isgroup, int *ngroups, int *nvaringroup, SEXP Sconstraints, SEXP Sinvconstraints);

double pmompenalty_approx(double *thopt, double **Hinv, double *tau, int thlength, double *nvaringroups, double *firstingroup);
double gmompenalty_approx(bool momsingle, bool momgroup, double *thopt, double **Hinv, double *Sinv, double phi, int thlength, int nsel, int nselgroupsint, double *nvarinselgroups, double *firstingroup, double *cholSini);


// Priors on Model Space (always return on log scale)
double unifPrior(int *sel, int *nsel, struct marginalPars *pars);
double unifPriorTP(int *sel, int *nsel, struct marginalPars *pars);
double unifPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);
double vectBinom(int *sel, int *nsel, int len_prDeltap, int len_prConstrp, struct marginalPars *pars);
double binomPrior(int *sel, int *nsel, struct marginalPars *pars);
double binomPriorTP(int *sel, int *nsel, struct marginalPars *pars);
double binomPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);
double betabinPrior(int *sel, int *nsel, struct marginalPars *pars);
double betabinPriorTP(int *sel, int *nsel, struct marginalPars *pars);
double betabinPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);

double complexityPrior(int *sel, int *nsel, struct marginalPars *pars);
double complexityPriorTP(int *sel, int *nsel, struct marginalPars *pars);
double complexityPrior_modavg(int *sel, int *nsel, struct modavgPars *pars);



//*************************************************************************************
// PRIOR DENSITIES ON PARAMETERS
//*************************************************************************************

void dmomgzell(double *ans, double *th, double *tau, double *nvaringroup, double *ngroups, double *detSinv, double *cholSinv, double *cholSini, bool logscale);
void demomgzell(double *ans, double *th, double *tau, double *nvaringroup, double *ngroups, double *detSinv, double *cholSinv, double *cholSini, bool logscale);
void dgzellgzell(double *ans, double *th, double *nvaringroup, double *ngroups, double *ldetSinv, double *cholSinv, double *cholSini, bool logscale);

void gzell_Sinv(double *Sinv, double *cholSinv, double *ldetSinv, int *ngroups, double *nvaringroups, int *sel, double *cholSini, crossprodmat *XtX, double *tau, double *taugroup, bool orthoapprox);
void gzell_Sinv_byprior(double *Sinv, double *cholSinv, double *ldetSinv, int *ngroups, double *nvaringroups, int *sel, double *cholSini, crossprodmat *XtX, int *n, double *tau, double *taugroup, int *priorcode);
double getelem_Sinv(int groupid, int k, int l, double *Sinv, double *cholSini, int ningroup);
void cholSini_indexes(double *cholSini, int *cholSsize, int ngroups, double *nvaringroups);


//*************************************************************************************
// PRIORS, GRADIENTS AND HESSIANS IN FORMAT REQUIRED BY TYPEDEF pt2fun, pt2funupdate, pt2grad
//*************************************************************************************

// pMOM + group Zellner
void pmomgzell_log (double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void pmomgzell_gradhess (double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pmomgzell_grad (double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pmomgzell_hess (double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);

// peMOM + group Zellner
void pemomgzell_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void pemomgzell_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pemomgzell_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pemomgzell_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);

// group Zellner + group Zellner
void gzellgzell_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void gzellgzell_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void gzellgzell_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void gzellgzell_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);

// pMOM + group Zellner + inverse gamma
void pmomgzellig_log (double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void pmomgzellig_gradhess (double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pmomgzellig_grad (double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pmomgzellig_hess (double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);

// peMOM + group Zellner + inverse gamma
void pemomgzellig_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void pemomgzellig_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pemomgzellig_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void pemomgzellig_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);

// group Zellner + group Zellner + inverse gamma
void gzellgzellig_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);
void gzellgzellig_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void gzellgzellig_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);
void gzellgzellig_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs);





//*************************************************************************************
// LEAST SQUARES
//*************************************************************************************

void leastsquares(double *theta, double *phi, double *ypred, double *y, double *x, crossprodmat *XtX, double *ytX, int *n, int *p, int *sel, int *nsel);


//*************************************************************************************
// MARGINAL LIKELIHOOD UNDER NORMAL ERRORS
//*************************************************************************************

// pMOM on all coef
double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

// piMOM on all coef
double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

// peMOM on all coef
double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

// Zellner on all coef
double zellnerMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double zellnerMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

// Normal on all coef
double normalidMarginalKC(int *sel, int *nsel, struct marginalPars *pars);
double normalidMarginalUC(int *sel, int *nsel, struct marginalPars *pars);

// pMOM on individual coef, block Zellner on groups
double pmomgzellMarg(int *sel, int *nsel, struct marginalPars *pars);

// pMOM on individual coef, group MOM on groups
double pmomgmomMarg(int *sel, int *nsel, struct marginalPars *pars);

// peMOM on individual coef, group eMOM on groups
double pemomgemomMarg(int *sel, int *nsel, struct marginalPars *pars);

// peMOM on individual coef, block Zellner on groups
double pemomgzellMarg(int *sel, int *nsel, struct marginalPars *pars);

// Zellner on individual coef, block Zellner on groups
double zellgzellMarg (int *sel, int *nsel, struct marginalPars *pars);

// Zellner on individual coef, normalid on groups
double zellnormidMarg (int *sel, int *nsel, struct marginalPars *pars);

// normalid on individual coef, group Zellner on groups
double normidgzellMarg (int *sel, int *nsel, struct marginalPars *pars);

// BIC criterion for Gaussian linear model
double bic_lm (int *sel, int *nsel, struct marginalPars *pars);




//*************************************************************************************
// MARGINAL LIKELIHOOD UNDER NORMAL / TWO-PIECE NORMAL / LAPLACE / TWO-PIECE LAPLACE RESIDUALS
//*************************************************************************************

double pmomMargTP(int *sel, int *nsel, struct marginalPars *pars);
double pimomMargTP(int *sel, int *nsel, struct marginalPars *pars);
double pemomMargTP(int *sel, int *nsel, struct marginalPars *pars);


//*************************************************************************************
// TWO-PIECE LAPLACE ROUTINES
//*************************************************************************************

double pmomMargLaplU(int *sel, int *nsel, struct marginalPars *pars);
double pimomMargLaplU(int *sel, int *nsel, struct marginalPars *pars);
double pemomMargLaplU(int *sel, int *nsel, struct marginalPars *pars);
double pmomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars);
double pimomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars);
double pemomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars);
double nlpMargAlapl(int *sel, int *nsel, struct marginalPars *pars, int *prior, int *symmetric);

void postmodeAlaplCDA(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *ftol, double *thtol, double *tau, double *taualpha, double *fixatanhalpha, double *alphaphi, double *lambdaphi, int *prior, int *hesstype, int *symmetric);
void mleAlaplCDA(double *thmode, double *fmode, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, bool useinit, int *symmetric, double *fixatanhalpha);

void fnegAlapl(double *ans, double *ypred, double *th, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, bool logscale, int *symmetric, int fixedalpha);
void fpnegAlaplUniv(int j, double *g, double *H, double *th, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric);
void fppnegAlapl(double **H, double *th, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric, int *hesstype);

void loglAlapl(double *ans, double *ypred, double *th, int *nsel, int *sel, int *n, double *scale, double *alpha, double *y, double *x, int *symmetric);
void loglnegGradHessAlaplUniv(int j, double *g, double *H, double *th, int *nsel, int *sel, int *n, int *p, double *y, double *ypred, double *x, crossprodmat *XtX, int *symmetric);
void loglnegHessAlapl(double **H, double *th, int *nsel, int *sel, int *n, int *p, double *y, double *ypred, double *x, crossprodmat *XtX, int *symmetric, int *hesstype);
void quadapproxALaplace(double *hdiag, double **H, int *nsel, int *sel, int *n, double *y0, double *x, double *th, double *vartheta, double *alpha, double *wy0, int *symmetric, double *w1, double *w2);



//*************************************************************************************
// TWO-PIECE NORMAL ROUTINES
//*************************************************************************************

double pmomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars);
double pimomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars);
double pemomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars);
double nlpMargSkewNorm(int *sel, int *nsel, struct marginalPars *pars, int *prior, int *symmetric);

void postmodeSkewNorm(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *tau, double *taualpha, double *alpha, double *lambda, bool *initmle, int *prior);
void postmodeSkewNormCDA(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *ftol, double *thtol, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric);

void fnegSkewnorm(double *ans, double *ypred, double *th, int *sel, int *nsel, int *n, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, bool logscale, int *symmetric);
void fpnegSkewnorm(double *g, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior);
void fpnegSkewnormUniv(int j, double *g, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmtric);
void fppnegSkewnorm(double **H, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric);
void fppnegSkewnormUniv(int j, double *H, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric);

void loglSkewnorm(double *ans, double *ypred, double *th, int *nsel, int *sel, int *n, double *scale, double *alpha, double *y, double *x, crossprodmat *XtX);
void loglnegGradSkewNorm(double *g, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x);
void loglnegGradSkewNormUniv(int j, double *g, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric);
void loglnegHessSkewNorm(double **H, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric);
void loglnegHessSkewNormUniv(int jj, double *H, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric);

void mleSkewnorm(double *thmode, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, bool useinit);



//*************************************************************************************
// Product MOM routines under Normal residuals
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
// Product iMOM routines under Normal residuals
//*************************************************************************************

double f2opt_imom(double *th);
double fimomNegC(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
double fimomNegC_non0(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void fppimomNegC_non0(double **ans, double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel);
void imomModeK(double *th, PolynomialRootFinder::RootStatus_T *status, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *sel, int *nsel, int *p);
void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *logscale);

double f2opt_imomU(double *th);
double fimomUNegC_non0(double *th, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel);
void fppimomUNegC_non0(double **ans, double *th, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel);
void imomModeU(double *th, PolynomialRootFinder::RootStatus_T *status, double *sumy2, crossprodmat *XtX, double *ytX, double *tau, double *alpha, double *lambda, int *sel, int *nsel, int *n, int *p);

void imomUIntegralApproxC(double *ILaplace, double *thopt, int *sel, int *nsel, int *n, int *p, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *logscale);

double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *B);

double f2int_imom(double phi);



#endif /* MODELSEL_H */

