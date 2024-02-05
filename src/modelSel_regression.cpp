// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//Include other headers
#include "cstat.h"
#include "crossprodmat.h"
#include "covariancemat.h"
#include "glm.h"
#include "survival.h"
#include "modelSel_regression.h"
#include "modselIntegrals.h"
#include "modselFunction.h"
#include "Polynomial.h"

#include <csignal>
#include <map>
#include <string>

//Global variables defined for minimization/integration routines
struct marginalPars f2opt_pars, f2int_pars;



//*************************************************************************************
//EXAMPLE FUNCTIONS CREATED BY RcppArmadillo
//*************************************************************************************

// another simple example: outer product of a vector, returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op, Rcpp::Named("inner")=ip);
}



//*************************************************************************************
//TEST FUNCTION TO HELP IN DEBUGGING
//*************************************************************************************


// [[Rcpp::export]]
SEXP testfunctionCI(SEXP x) {
    SEXP ans;

    testfunction(REAL(x)); //run whatever code

    PROTECT(ans = Rf_allocVector(REALSXP, 1)); //return 0
    *REAL(ans)= 0;
    UNPROTECT(1);
    return ans;
}


void testfunction(double *x) {

  Rprintf("apnorm  %f\n", apnorm(*x,false));
  Rprintf("apnorm2 %f\n", apnorm2(*x,false));

  Rprintf("ainvmillsnorm  %f\n", ainvmillsnorm(*x));
  Rprintf("ainvmillsnorm2 %f\n", ainvmillsnorm2(*x));

}

/*
// FUNCTION: sum_k (sel[k]+1) th[k]^2 + sum_{l>k} th[k] th[l], for sel[k]>=0. The minimum is trivially at 0
//
// GRADIENT WRT th[j]: 2 (sel[j]+1) th[j] + sum_{l \neq j} th[l]
// HESSIAN WRT th[j]:  2 (sel[j]+1)
//
// Optionally we can store into "funargs" the following info: sumth= sum_k th[k]; sumth2= sum_k (sel[k]+1) th[k]^2; sumcrossprod= sum_{l>k} th[k] th[l]
//
// Then fun(th)= sumth2 + sumcrossprod, hence changing th[j] to thjnew gives
//
// fun(thnew)= f(th) + (sel[j]+1) (thjnew - th[j])^2 + (thjnew-th[j]) (sumth - th[j])
//
// Also grad(th)= 2 (sel[j]+1) th[j] + sumth - th[j]


//Evaluate function but not funargs
void foo(double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs) {
  int k, l;
  for (k=0, (*f)=0; k< *nsel; k++) {
    (*f) += (double)(sel[k]+1) * th[k] * th[k];
    for (l=k+1; l< *nsel; l++) { (*f) += th[k] * th[l]; }
  }
}

//Compute gradient and hessian wrt th[j], not using funargs
void foogradhess(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int l;
  (*hess)= 2.0 * (double)(sel[j]+1);
  (*grad)= (*hess) * th[j];
  for (l=0; l< j; l++) { (*grad)+= th[l]; }
  for (l=j+1; l< *nsel; l++) { (*grad)+= th[l]; }
}

//Return univariate optimum for th[j], that is thnew= -0.5/(sel[j]+1) * sum_{l \neq j} th[l]. Not using funargs
void fooupdateUniv(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int l;
  *thnew= 0;
  for (l=0; l< j; l++) (*thnew)-= th[l];
  for (l=j+1; l< *nsel; l++) (*thnew)-= th[l];
  (*thnew) *= 0.5/((double)(sel[j]+1));
}


//Evaluate function and funargs
void fooargs(double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs) {
  int k, l;
  double sumth=0, sumth2=0, sumcrossprod=0;
  for (k=0, (*f)=0; k< *nsel; k++) {
    sumth += th[k];
    sumth2 += (double)(sel[k]+1) * th[k] * th[k];
    for (l=k+1; l< *nsel; l++) { sumcrossprod += th[k] * th[l]; }
  }
  (*f)= sumth2 + sumth;
  *(*funargs)["sumth"]= sumth;
  *(*funargs)["sumth2"]= sumth2;
  *(*funargs)["sumcrossprod"]= sumcrossprod;
}

//Update function and funargs from changing th[j] to thjnew
void fooupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs) {
  double thdif= *thjnew - th[j];
  *(*funargs)["sumth"] += thdif;
  *(*funargs)["sumth2"] += (double)(sel[j]+1) * (pow(*thjnew,2) - pow(th[j],2));
  *(*funargs)["sumcrossprod"] += (*thjnew - th[j]) * (*(*funargs)["sumth"] - *thjnew);
  (*fnew)= *(*funargs)["sumth2"] + *(*funargs)["sumcrossprod"];
}

//Compute gradient and hessian wrt th[j], using funargs
void foogradhessargs(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double*> *funargs) {
  (*hess)= 2.0 * (double)(sel[j]+1);
  (*grad)= (*hess) * th[j] + (*(*funargs)["sumth"]) - th[j];
}



void testfunction() {

  int nsel=2, *sel;
  double *thini, *thopt, fopt;
  struct marginalPars *pars= NULL;
  modselFunction *msfun;
  std::map<string, double *> funargs;

  //Alloc memory for elements in funargs. For vector arguments use dvector
  double sumth= 0, sumth2= 0, sumcrossprod= 0;
  funargs["sumth"]= &sumth; funargs["sumth2"]= &sumth2; funargs["sumcrossprod"]= &sumcrossprod;

  sel= ivector(0,nsel); thini= dvector(0,nsel); thopt= dvector(0,nsel);
  sel[0]= 0; sel[1]= 2;
  thini[0]= 1; thini[1]= 1;
  msfun= new modselFunction(sel, &nsel, pars, NULL);

  //Option 1. CDA
  msfun->updateUniv= &fooupdateUniv;
  msfun->cda(thopt, thini);
  Rprintf("cda.               thopt= %f %f\n", thopt[0], thopt[1]);

  //Option 2. CDA providing foo
  msfun->fun= &foo;
  msfun->updateUniv= &fooupdateUniv;
  msfun->cda(thopt, &fopt, thini);
  Rprintf("cda.               thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 3. CDA providing foo and funargs
  msfun->fun= &fooargs;
  msfun->funupdate= &fooupdate;
  msfun->updateUniv= &fooupdateUniv;
  msfun->cda(thopt, &fopt, thini, &funargs);
  Rprintf("cda.               thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 4. block CDA
  msfun->fun= &foo;
  msfun->updateUniv= &fooupdateUniv;
  msfun->blockcda(thopt, &fopt, thini);
  Rprintf("blockcda.          thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 5. cdaNewton not using funargs (requires function foo)
  msfun->fun= &foo;
  msfun->gradhessUniv= &foogradhess;
  msfun->cdaNewton(thopt, &fopt, thini, 1);
  Rprintf("cdaNewton.         thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 7. cdaNewton using funargs (requires functions fooargs and fooupdate, and allocating memory for all elements in funargs)
  msfun->fun= &fooargs;
  msfun->funupdate= &fooupdate;
  msfun->gradhessUniv= &foogradhessargs;
  msfun->cdaNewton(thopt, &fopt, thini, &funargs, 1);
  Rprintf("cdaNewton.         thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 8. blockcdaNewton
  msfun->fun= &foo;
  msfun->gradhessUniv= &foogradhess;
  msfun->blockcdaNewton(thopt, &fopt, thini, 1);
  Rprintf("blockcdaNewton.    thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  free_ivector(sel, 0,nsel); free_dvector(thini, 0,nsel); free_dvector(thopt, 0,nsel);
  delete msfun;

}
*/




//*************************************************************************************
//SETTING PRIOR & MARGINALS
//*************************************************************************************


int mspriorCode(int *prCoef, int *prGroup, struct marginalPars *pars) {
  //Returns a two-digit code indicating the prior on regression coefficients. The 1st digit is the prior on individual coef; The 2nd digit the prior on groups of coefficients
  // Input
  // - prCoef: 0 for pMOM; 1 for piMOM; 2 for peMOM; 3 for Zellner; 4 for normalid; 10 for group pMOM; 13 for group Zellner
  // - prGroup: 0 for pMOM; 1 for piMOM; 2 for peMOM; 3 for Zellner; 4 for normalid; 10 for group pMOM; 11 for group iMOM; 12 for group eMOM; 13 for group Zellner
  // Output
  //   0: pMOM on all coef
  //   1: peMOM on all coef
  //   2: piMOM on all coef
  //   3: Zellner on all coef
  //   4: normalid on all coef
  //   5: group pMOM (same as pMOM, standardized by n / X'X)
  //   9: group Zellner on all coef
  //  10: pMOM + group MOM
  //  13: pMOM + group Zellner
  //  32: peMOM + group eMOM
  //  33: peMOM + group Zellner
  //  43: Zellner + group Zellner
  //  50: group MOM + group MOM
  //  53: group MOM + group Zellner
  //  63: group Zellner + group Zellner
  //  73: normalid + group Zellner
  // 100: BIC (no prior, tells marginal likelihood routines to return -0.5 BIC, the BIC approx to the marginal likelihood)
  bool hasgroups= (*((*pars).ngroups)) < (*((*pars).p));
  int ans;
  if (*prCoef==100) {
    ans= 100;  //BIC
  } else {
    if (!hasgroups) {
      if (*prCoef==0) {        //pMOM on all coef
        ans= 0;
      } else if (*prCoef==1) { //piMOM on all coef
        ans= 1;
      } else if (*prCoef==2) { //peMOM on all coef
        ans= 2;
      } else if (*prCoef==3) { //Zellner on all coef
        ans= 3;
      } else if (*prCoef==4) { //normalid on all coef
        ans= 4;
      } else if (*prCoef==10) {
        ans= 5;                //group pMOM 
      } else if (*prCoef==13) { //block Zellner on all coef
        ans= 9;
      } else {
        Rf_error("Prior specified by priorCoef not currently implemented\n");
      }
    } else {
      if ((*prCoef==0) & (*prGroup==0)) {        //pMOM on all coef
        ans= 0;
      } else if ((*prCoef==1) & (*prGroup==1)) { //piMOM on all coef
        ans= 1;
      } else if ((*prCoef==2) & (*prGroup==2)) { //peMOM on all coef
        ans= 2;
      } else if ((*prCoef==3) & (*prGroup==3)) { //Zellner on all coef
        ans= 3;
      } else if ((*prCoef==4) & (*prGroup==4)) { //normalid on all coef
        ans= 4;
      } else if ((*prCoef==0) & (*prGroup==10)) { //pMOM + group MOM
        ans= 10;
      } else if ((*prCoef==0) & (*prGroup==13)) { //pMOM + group Zellner
        ans= 13;
      } else if ((*prCoef==2) & (*prGroup==12)) { //peMOM + group eMOM
        ans= 32;
      } else if ((*prCoef==2) & (*prGroup==13)) { //peMOM + group Zellner
        ans= 33;
      } else if ((*prCoef==3) & (*prGroup==13)) { //Zellner + group Zellner
        ans= 43;
      } else if ((*prCoef==10) & (*prGroup==10)) { //group pMOM + group pMOM
        ans= 50;
      } else if ((*prCoef==10) & (*prGroup==13)) { //group pMOM + group Zellner
        ans= 53;
      } else if ((*prCoef==13) & (*prGroup==13)) { //group Zellner + group Zellner
        ans= 63;
      } else if ((*prCoef==4) & (*prGroup==13)) { //normalid + group Zellner
        ans= 73;
      } else {
        Rf_error("Prior specified by priorCoef and priorGroup not currently implemented\n");
      }
    }
  }
  return ans;
}

pt2margFun set_marginalFunction(struct marginalPars *pars) {
  /*Returns pointer to function to compute the marginal density of the data for a given model indicator

   - priorcode: prior code returned by mspriorCode
   - knownphi: 1 if residual variance phi is known, 0 otherwise. knownphi==1 currently only allowed for Normal errors
   - family: parametric family
       1 for Normal, 2 for two-piece Normal; 3 for Laplace; 4 for two-piece Laplace, 0 to infer family within 1-4
       11 for AFT survival
       21 for logistic regression
       22 for Poisson regression

  Note: if phi known, when actually calling the returned pt2margFun, phi must be set in the parameter of type struct marginalPars *
  */
  bool hasgroups= (*((*pars).ngroups)) < (*((*pars).p));
  int priorcode= *((*pars).priorcode), family= *((*pars).family), knownphi= *((*pars).knownphi);
  pt2margFun ans=NULL;
  if ((family)==1) {  //Normal errors
    if (priorcode==0) {
      if (knownphi==1) { ans= pmomMarginalKC; } else { ans= pmomMarginalUC; }
    } else if (priorcode==1) {
      if (knownphi==1) { ans= pimomMarginalKC; } else { ans= pimomMarginalUC; }
    } else if (priorcode==2) {
      if (knownphi==1) { ans= pemomMarginalKC; } else { ans= pemomMarginalUC; }
    } else if (priorcode==3) {
      if (knownphi==1) { ans= zellnerMarginalKC; } else { ans= zellnerMarginalUC; }
    } else if (priorcode==4) {
      if (knownphi==1) { ans= normalidMarginalKC; } else { ans= normalidMarginalUC; }
    } else if (priorcode==10) {
      ans= pmomgmomMarg;
    } else if (priorcode==13) {
      ans= pmomgzellMarg;
    } else if (priorcode==32) {
      ans= pemomgemomMarg;
    } else if (priorcode==33) {
      ans= pemomgzellMarg;
    } else if (priorcode==43) {
      ans= zellgzellMarg;
    } else if (priorcode==73) {
      ans= normidgzellMarg;
    } else if (priorcode==100) {
      ans= bic_lm;
    } else {
      Rf_error("The prior in (priorCoef,priorGroup) not currently implemented for linear regression");
    }
  } else if ((family)==11) {  //AFT survival with Normal errors
    if (!hasgroups) {
      if (priorcode==0) {
        ans= pmomgzellSurvMarg;
      } else if (priorcode==1) {
        ans= pemomgzellSurvMarg;
      } else if (priorcode==3) {
        Rprintf("Zellner prior not implemented, using group Zellner prior instead\n");
        ans= gzellgzellSurvMarg;
      } else if (priorcode==5) {
        ans= gmomgmomSurvMarg;
      } else if (priorcode==9) {
        ans= gzellgzellSurvMarg;
      } else {
        Rf_error("The prior in priorCoef not implemented for the specified AFT model");
      }
    } else {
      if (priorcode==10) {
        ans= pmomgmomSurvMarg;
        if (*(*pars).method !=2) Rprintf("For priorCoef=groupmomprior() + priorGroups=groupmomprior() only method='ALA' is implemented. Using ALA instead.\n");
      } else if (priorcode==13) {
        ans= pmomgzellSurvMarg;
      } else if (priorcode==32) {
        ans= pemomgemomSurvMarg;
      } else if (priorcode==33) {
        ans= pemomgzellSurvMarg;
      } else if (priorcode==43) {
        Rprintf("Zellner prior not implemented, using group Zellner prior instead\n");
        ans= gzellgzellSurvMarg;
      } else if (priorcode==50) {
        if (*(*pars).method !=2) Rprintf("For priorCoef=groupmomprior() + priorGroups=groupzellnerprior() only method='ALA' is implemented. Using ALA instead.\n");
        ans= gmomgmomSurvMarg;
      } else if (priorcode==53) {
        if (*(*pars).method !=2) Rprintf("For priorCoef=groupmomprior() + priorGroups=groupzellnerprior() only method='ALA' is implemented. Using ALA instead. You can set method='Laplace' for priorCoef=momprior() + priorGroups=groupzellnerprior()\n");
        ans= gmomgzellSurvMarg;
      } else if (priorcode==63) {
        ans= gzellgzellSurvMarg;
      } else {
        Rf_error("The prior in (priorCoef,priorGroup) not implemented for survival data\n");
      }
    }
  } else if ((family)==2) { //Two-piece Normal errors
    if (priorcode >= 10) Rf_error("Group priors not implemented for Two-piece Normal errors");
    if (priorcode==0) {
      ans= pmomMargSkewNormU;
    } else if (priorcode==1) {
      ans= pimomMargSkewNormU;
    } else if (priorcode==2) {
      ans= pemomMargSkewNormU;
    } else if (priorcode==3) {
      Rf_error("Zellner prior with two-piece Normal errors not currently implemented");
    }
  } else if ((family)==3) { //Laplace errors
    if (priorcode >= 10) Rf_error("Group priors not implemented for Laplace errors");
    if (priorcode==0) {
      ans= pmomMargLaplU;
    } else if (priorcode==1) {
      ans= pimomMargLaplU;
    } else if (priorcode==2) {
      ans= pemomMargLaplU;
    } else if (priorcode==3) {
      Rf_error("Zellner prior with Laplace errors not currently implemented");
    }
  } else if ((family)==4) { //Asymmetric Laplace errors
    if (priorcode >= 10) Rf_error("Group priors not implemented for asymmetric Laplace errors");
    if (priorcode==0) {
      ans= pmomMargAlaplU;
    } else if (priorcode==1) {
      ans= pimomMargAlaplU;
    } else if (priorcode==2) {
      ans= pemomMargAlaplU;
    } else if (priorcode==3) {
      Rf_error("Zellner prior with asymmetric Laplace errors not currently implemented");
    }
  } else if ((family)==0) { //Normal + Two-piece Normal + Laplace + Two-piece Laplace errors
    if (priorcode >= 10) Rf_error("Group priors with family='auto' not currently implemented");
    if (priorcode==0) {
      ans= pmomMargTP;
    } else if (priorcode==1) {
      ans= pimomMargTP;
    } else if (priorcode==2) {
      ans= pemomMargTP;
    } else if (priorcode==3) {
      Rf_error("Zellner prior with family='auto' not currently implemented");
    }
  } else if ((family>=21) && (family<=22)) {  //Logistic or Poisson regression
    if ((priorcode == 3) || (priorcode==43)) {
      Rprintf("Zellner prior not implemented for GLMs. Using group Zellner instead\n");
      *((*pars).priorcode)= 63;
    }
    ans= marginal_glm;
  } else {
    Rf_error("This distribution family is not available");
  }
  return ans;
}


pt2margFun set_priorFunction(int *prDelta, int *prConstr, int *family) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - prDelta: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2margFun ans=NULL;
  if (*family != 0) {
    if (*prDelta != *prConstr) Rf_error("priorConstraints must be of the same family as priorDelta (e.g. both Binomial, both Beta-Binomial, both Complexity priors)");
    if (*prDelta==0) { ans= unifPrior; } else if (*prDelta==1) { ans= binomPrior; } else if (*prDelta==2) { ans= betabinPrior; } else if (*prDelta==3) { ans= complexityPrior; }
  } else {
    if (*prDelta==0) { ans= unifPriorTP; } else if (*prDelta==1) { ans= binomPriorTP; } else if (*prDelta==2) { ans= betabinPriorTP; } else if (*prDelta==3) { ans= complexityPrior; }
  }
  return ans;
}

pt2modavgPrior set_priorFunction_modavg(int *priorModel) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - priorModel: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2modavgPrior ans=NULL;
  if (*priorModel==0) { ans= unifPrior_modavg; } else if (*priorModel==1) { ans= binomPrior_modavg; } else if (*priorModel==2) { ans= betabinPrior_modavg; } else if (*priorModel==3) { ans= complexityPrior_modavg; }
  return ans;
}


//********************************************************************************************
// GENERAL ALGEBRA
//********************************************************************************************

//multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]
//Note: A is indexed at 1. x and sel are indexed at 0. ans is indexed at 1.
void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<=fi;_i++) {
    for (_j=_i, ans[_i]=0; _j<=fi; _j++) { ans[_i]+= A[_i][_j] * x[sel[_j-1]]; }
    for (_j= 1; _j<_i; _j++) { ans[_i]+= A[_j][_i] * x[sel[_j-1]]; }
  }
}

//multiply symmetric A[1..ncolA][1..ncolA] (formatted as vector) with x[1..nsel].
//Use only selected elems in A and all elems in x[1..nsel]
void Asel_x(crossprodmat *A, int ncolA, double *x, int nsel, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<= nsel;_i++) {
    for (_j=1, ans[_i]=0; _j<= nsel; _j++) { ans[_i]+= (A->at(sel[_j]*ncolA+sel[_i])) * x[_j]; }
  }
}


//Add constant ct to diagonal elements in XtX[sel,sel]. XtX[0..p-1][0..p-1] is formatted as a vector indexed at 0, V[sel[0]..sel[nsel-1]][sel[0]..sel[nsel-1]] as a matrix indexed at 1, sel is indexed at 0
//Note: Only diagonal & upper-diagonal elements in V are set.
void addct2XtX(double *ct, crossprodmat *XtX, int *sel, int *nsel, int *p, double **V) {
  int i,j;
  for (i=1;i<=(*nsel);i++) { V[i][i]= XtX->at(sel[i-1]*(*p)+sel[i-1]) + (*ct); }
  for (i=1;i<=(*nsel);i++) {
    for (j=i+1;j<=(*nsel);j++) {
        V[i][j]= XtX->at(sel[j-1]*(*p) + sel[i-1]);
    }
  }
}


//*************************************************************************************
// MODEL AVERAGING ROUTINES
//*************************************************************************************

void set_modavgPars(struct modavgPars *pars, int *n, int *p1, int *p2, int *isbinary, int *ybinary, double *y, double *sumy2, double *x1, double *x2, crossprodmat *XtX, double *ytX, double *cholS2, double *S2inv, double *cholS2inv, double *colsumx1sq, double *alpha, double *lambda, int *priorCoef, int *r, double *tau1, double *tau2, int *priorTau1, double *atau1, double *btau1, int *priorModel, double *prModelpar) {
  (*pars).n= n;
  (*pars).p1= p1;
  (*pars).p2= p2;
  (*pars).isbinary= isbinary;
  (*pars).ybinary= ybinary;
  (*pars).y= y;
  (*pars).sumy2= sumy2;
  (*pars).x1= x1;
  (*pars).x2= x2;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).cholS2= cholS2;
  (*pars).S2inv= S2inv;
  (*pars).cholS2inv= cholS2inv;
  (*pars).colsumx1sq= colsumx1sq;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).priorCoef= priorCoef;
  (*pars).r= r;
  (*pars).tau1= tau1;
  (*pars).tau2= tau2;
  (*pars).priorTau1= priorTau1;
  (*pars).atau1= atau1;
  (*pars).btau1= btau1;
  (*pars).priorModel= priorModel;
  (*pars).prModelpar= prModelpar;
}


//MH within Gibbs scheme to sample from the joint posterior of (theta,delta) under a product MOM prior and a linear regression model
//Input
// - pars: data, pre-computed quantities and prior parameters. See definition of struct modavgPars for details
// - niter: number of MCMC iterations
// - thinning: only 1 out of each thinning iterations is kept
// - burnin: burning
// - niniModel: number of variables initially in the model
// - iniModel: initial model (vector of length p1)
// - iniCoef1: initial coefficient values for variables being selected (length p1)
// - iniCoef2: initial coefficient values for variables always in the model (length p2)
// - iniPhi: initial residual variance value
// - iniOthers: initial values for other variables. Currently this indicates initial tau value, and is only use if (*pars).priorTau1 !=0.
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
//Output
// - postModel: MCMC saves for variable inclusion indicators
// - margpp: marginal posterior prob for inclusion of each covariate (uses MH acceptance prob, which is more precise than simply averaging inclusion indicators)
// - postCoef1: MCMC saves for regression coefficients of variables being selected
// - postCoef2: MCMC saves for regression coefficients of variables which are always in the model
// - postPhi: MCMC saves for residual variance
// - postOther: MCMC saves for other parameters. Currently saves tau values (pMOM prior precision) when (*pars).priorTau1 != 0.

SEXP pmomLM_I(SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP SXtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar) {
  struct modavgPars pars;
  int mcmc2save, *postModel;
  double *margpp, *postCoef1, *postCoef2, *postPhi, *postOther, tau1copy= REAL(tau1)[0];
  crossprodmat *XtX;
  SEXP ans;

  XtX= new crossprodmat(REAL(SXtX),INTEGER(n)[0],INTEGER(p1)[0],true);

  PROTECT(ans= Rf_allocVector(VECSXP, 7));
  mcmc2save= floor((INTEGER(niter)[0] - INTEGER(burnin)[0] +.0)/(INTEGER(thinning)[0] +.0));

  SET_VECTOR_ELT(ans, 0, Rf_allocVector(INTSXP, mcmc2save * INTEGER(p1)[0]));
  postModel= INTEGER(VECTOR_ELT(ans,0));

  SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, INTEGER(p1)[0]));
  margpp= REAL(VECTOR_ELT(ans,1));

  SET_VECTOR_ELT(ans, 2, Rf_allocVector(REALSXP, mcmc2save * INTEGER(p1)[0]));
  postCoef1= REAL(VECTOR_ELT(ans,2));

  SET_VECTOR_ELT(ans, 3, Rf_allocVector(REALSXP, mcmc2save * INTEGER(p2)[0]));
  postCoef2= REAL(VECTOR_ELT(ans,3));

  SET_VECTOR_ELT(ans, 4, Rf_allocVector(REALSXP, mcmc2save));
  postPhi= REAL(VECTOR_ELT(ans,4));

  if (INTEGER(priorTau1)[0] != 0) {
    SET_VECTOR_ELT(ans, 5, Rf_allocVector(REALSXP, mcmc2save));
  } else {
    SET_VECTOR_ELT(ans, 5, Rf_allocVector(REALSXP, 1));
  }
  postOther= REAL(VECTOR_ELT(ans,5));

  set_modavgPars(&pars,INTEGER(n),INTEGER(p1),INTEGER(p2),INTEGER(isbinary),INTEGER(ybinary),REAL(y),REAL(sumy2),REAL(x1),REAL(x2),XtX,REAL(ytX),REAL(cholS2),REAL(S2inv),REAL(cholS2inv),REAL(colsumx1sq),REAL(alpha),REAL(lambda),INTEGER(priorCoef),INTEGER(r),&tau1copy,REAL(tau2),INTEGER(priorTau1),REAL(atau1),REAL(btau1),INTEGER(priorModel),REAL(prModelpar));
  pmomLM(postModel, margpp, postCoef1, postCoef2, postPhi, postOther, &pars, INTEGER(niter), INTEGER(thinning), INTEGER(burnin), INTEGER(niniModel), INTEGER(iniModel), REAL(iniCoef1), REAL(iniCoef2), REAL(iniPhi), REAL(iniOthers), INTEGER(verbose));

  delete XtX;
  UNPROTECT(1);
  return ans;
}


void pmomLM(int *postModel, double *margpp, double *postCoef1, double *postCoef2, double *postPhi, double *postOther, struct modavgPars *pars, int *niter, int *thinning, int *burnin, int *niniModel, int *iniModel, double *iniCoef1, double *iniCoef2, double *iniPhi, double *iniOthers, int *verbose) {
  int i, j, k, ilow, iupper, savecnt, niterthin, niter10, nsel= *niniModel, *curModel, newdelta, n=*(*pars).n, p1=*(*pars).p1, p2=*(*pars).p2, psn, resupdate, isbinary=*(*pars).isbinary;
  double *res, *partialres, sumres2, sumpartialres2, newcoef, *curCoef1, *curCoef2, curPhi=1.0, *linpred1, *linpred2, pinclude, *temp;
  if (*verbose) Rprintf("Running MCMC");
  niterthin= (int) floor((*niter - *burnin +.0)/(*thinning +.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  if (*burnin >0) { ilow= - *burnin; savecnt=0; iupper= *niter - *burnin +1; } else { ilow=0; savecnt=1; iupper= *niter; }
  //Initialize
  curModel= ivector(0,p1); curCoef1= dvector(0,p1); curCoef2= dvector(0,p2); linpred1= dvector(0,n); linpred2= dvector(0,n);
  res= dvector(0, n); partialres= dvector(0,n);
  for (i=0; i<p1; i++) { margpp[i]= 0; curModel[i]= postModel[i*niterthin]= iniModel[i]; curCoef1[i]= postCoef1[i*niterthin]= iniCoef1[i]; }
  for (i=0; i<p2; i++) { curCoef2[i]= postCoef2[i*niterthin]= iniCoef2[i]; }
  if (isbinary) { curPhi= postPhi[0]= 1.0; } else { curPhi= postPhi[0]= *iniPhi; }
  postOther[0]= iniOthers[0];
  Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);
  Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
  if (isbinary) sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
  for (i=0, sumres2=0; i<n; i++) { res[i]= partialres[i]= (*pars).y[i] - linpred1[i] - linpred2[i]; sumres2+= res[i]*res[i]; }
  sumpartialres2= sumres2;
  //MCMC iterations
  for (i=ilow; i< iupper; i++) {
    //Sample (curCoef1,curModel)
    for (j=0; j< *(*pars).p1; j++) {
      if (curModel[j]) {
        for (k=0, sumpartialres2=0; k<n; k++) { partialres[k]= res[k] + curCoef1[j] * ((*pars).x1[n*j+k]); sumpartialres2+= partialres[k]*partialres[k]; }
      }
      MHTheta1pmom(&newdelta, &newcoef, &pinclude, &resupdate, res, partialres, &sumres2, &sumpartialres2, j, &nsel, curModel, curCoef1, &curPhi, pars);
      if (newdelta > curModel[j]) { nsel++; } else if (newdelta < curModel[j]) { nsel--; }
      curModel[j]= newdelta; curCoef1[j]= newcoef;
      if (i>=0) margpp[j]+= pinclude;
      if (resupdate) { temp= partialres; partialres= res; res=temp; }
    }
    //Sample curCoef2
    for (k=0; k<n; k++) res[k]+= linpred2[k];
    simTheta2(curCoef2, res, &curPhi, pars);
    Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
    for (k=0, sumres2=0; k<n; k++) { res[k]-= linpred2[k]; sumres2+= res[k]*res[k]; }
    //Sample phi
    if (isbinary==0) { curPhi= simPhipmom(&nsel, curModel, curCoef1, curCoef2, &sumres2, pars); }
    //Sample tau
    if (*(*pars).priorTau1 !=0) { *(*pars).tau1= simTaupmom(&nsel, curModel, curCoef1, &curPhi, pars); }
    //Sample latent variables (only for probit model)
    if (isbinary) {
      Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);  //update linpred1 (linpred2 already updated)
      sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
    }
    //Save values
    if ((i>0) && (i % (*thinning))==0) {
      for (j=0; j<p1; j++) {
        psn= niterthin*j+savecnt;
        postModel[psn]= curModel[j];
        postCoef1[psn]= curCoef1[j];
      }
      for (j=0; j<p2; j++) postCoef2[niterthin*j+savecnt]= curCoef2[j];
      postPhi[savecnt]= curPhi;
      if (*(*pars).priorTau1 !=0) postOther[savecnt]= *(*pars).tau1;
      savecnt++;
    }
    if ((*verbose ==1) && (i%niter10)==0) Rprintf(".");
  } //end MCMC for
  if (iupper>ilow) { for (j=0; j< p1; j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*verbose ==1) Rprintf("Done.\n");
  free_ivector(curModel,0,p1); free_dvector(curCoef1,0,p1); free_dvector(curCoef2,0,p2); free_dvector(linpred1,0,n); free_dvector(linpred2,0,n);
  free_dvector(res,0,n); free_dvector(partialres,0,n);
}


//Sample from the posterior of latent variables in probit model given the regression coefficients (i.e. the linear predictor)
//Input:
// - ybinary: response variable (1: success; 0: failure)
// - linpred1: linear predictor for current regression coefficients associated to variables under selection
// - linpred2: linear predictor associated to adjustment variables
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - y: sampled values for the latent variables
// - res: sampled residuals, i.e. y - linpred1 - linpred2
// - sumres2: sum(res^2)
// - (*pars).ytX: updated t(y) %*% x1
// - (*pars).sumy2: updated sum(y^2)
void sample_latentProbit(double *y, double *res, double *sumres2, int *ybinary, double *linpred1, double *linpred2, struct modavgPars *pars) {
  int i;
  double linpred, plinpred, u;
  for (i=0, *sumres2=0, *(*pars).sumy2=0; i< *(*pars).n; i++) {
    linpred= linpred1[i] + linpred2[i];
    plinpred= pnormC(-linpred,0,1);
    if (ybinary[i]) {
      u= plinpred + (1.0-plinpred) * runif();  //u ~ Unif(plinpred,1)
    } else {
      u= plinpred * runif(); //u ~ Unif(0,plinpred)
    }
    res[i]= qnormC(u,0,1);
    (*sumres2)+= res[i]*res[i];
    y[i]= linpred + res[i];
    (*(*pars).sumy2)+= y[i]*y[i];
  }
  Atvecx((*pars).x1,y,(*pars).ytX,0,*(*pars).p1 -1,0,*(*pars).n -1); //update ytX=Xty
}

//Univariate MH update of (coef,delta)
//Input:
// - res: vector with residuals. Only used if curModel[j]==0
// - partialres: vector with partial residuals removing the variable from the model. Only used if curModel[j]==1
// - sumres2: sum(res*res)
// - sumpartialres2: sum(partialres*partialres)
// - j: index of the variable for which move is to be proposed
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - newdelta: new value
// - newcoef: new coefficient
// - pinclude: probability of including the variable in the model
// - res: on input, vector with residuals given the current coefficient value. On output, residuals given the updated coefficient value.
void MHTheta1pmom(int *newdelta, double *newcoef, double *pinclude, int *resupdate, double *res, double *partialres, double *sumres2, double *sumpartialres2, int j, int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int n= *(*pars).n, logscale=1, nsel0, nsel1, deltaprop, nu, i;
  double m1, *xj, m0, logbf, logpratio, thetaprop, m, S, propPars[6], lhood, lprior, lprop, lambda=0.0, num, den, sqrtPhi=sqrt(*curPhi);
  pt2modavgPrior priorFunction= NULL;
  *resupdate= 0;
  xj= (*pars).x1+j*n; //pointer to variable j in x1
  priorFunction= set_priorFunction_modavg((*pars).priorModel);
  //Propose delta
  if (curModel[j]) {
    m1= pmomMargKuniv(partialres, xj, sumpartialres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(partialres, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel -1; nsel1= *nsel;
  } else {
    m1= pmomMargKuniv(res, xj, sumres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(res, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel; nsel1= *nsel +1;
  }
  logbf= m0-m1;
  logpratio= priorFunction(curModel, &nsel0, pars) - priorFunction(curModel, &nsel1, pars); //we use curModel in both cases as priorFunction currently only depends on nb vars
  if ((!curModel[j]) && (((*nsel)+(*(*pars).p2)) >= n)) {
    *pinclude= 0.0;
  } else {
    *pinclude= 1.0/(1.0+exp(logbf+logpratio));
  }
  if (runif() < *pinclude) { deltaprop=1; } else { deltaprop=0; }
  //Propose coef
  nu= (int) sqrt((double) n);
  if ((curModel[j]==0) && (deltaprop==0)) {  //proposal is to keep variable out of the model
    *newdelta=0; *newcoef=0;
  } else {
    S= (*pars).colsumx1sq[j] + 1.0/(*(*pars).tau1);
    if (curModel[j]) {
      for (i=0, m=0; i<n; i++) m+= xj[i]*partialres[i];
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, partialres, xj, &m1, &nu);
    } else {
      for (i=0, m=0; i<n; i++) m+= xj[i]*res[i];
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, res, xj, &m1, &nu);
    }
    if (curModel[j] && deltaprop) {  //proposal is to keep variable in the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, lhood=0; i<n; i++) {
        partialres[i]-= thetaprop*xj[i];
        lhood+= dnormC(partialres[i],0,sqrtPhi,1) - dnormC(res[i],0,sqrtPhi,1);
      }
      lprior= dmom(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1) - dmom(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lprop= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) - dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1);
      lambda= exp(lhood+lprior+lprop);
    } else if ((curModel[j]==0) && deltaprop) { //proposal is to add variable to the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, num=0; i<n; i++) {
        partialres[i]= res[i] - thetaprop*xj[i];
        num+= dnormC(partialres[i],0,sqrtPhi,1);
      }
      num+= dmom(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      den= dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1) + m1;
      lambda= exp(num-den);
    } else {    //(curModel[j] && (deltaprop==0)), i.e. proposal is to drop variable from the model
      thetaprop=0;
      num= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) + m1;
      for (i=0, den=0; i<n; i++) { den+= dnormC(res[i],0,sqrtPhi,1); }
      den+= dmom(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lambda= exp(num-den);
    }
    if (runif()<lambda) {
      *newdelta=deltaprop; *newcoef= thetaprop;
      *resupdate= 1;  //signal that res and partialres have to be interchanged after exiting the function
      for (i=0, *sumres2=0; i< n; i++) (*sumres2)+= partialres[i]*partialres[i];
    } else {
      *newdelta=curModel[j]; *newcoef= curCoef1[j];
    }
  }
}

//Find parameters for univariate pmom proposal distribution (2 component mixture of T distributions)
//   Posterior: N(e; xj*theta; phi*I) * pmom(theta; phi, r, tau) / m1
//   Mixture: w1 * T_nu(theta;mu1,sigma21) + (1-w1) * T_nu(theta;mu2,sigma22)  (sigma21, sigma22 denote variances)
//Input:
// - m, S: posterior location & scale parameters
// - phi: residual variance
// - r: product MOM power parameter
// - tau1: product MOM prior dispersion parameter
// - e: response variable
// - xj: predictor
// - m1: normalization constant
// - nu: desired degrees of freedom
//Output: means in propPars[0:1], SD in propPars[2:3], weights in propPars[4:5]
void proposalpmom(double *propPars, double *m, double *S, double *phi, int *r, double *tau1, int *n, double *e, double *xj, double *m1, int *nu) {
  int i;
  double eps, fmode, sqrtPhi=sqrt(*phi), temp, temp2, doubler=2*(*r), ct2;
  //Find modes
  eps= sqrt((*m)*(*m) + 8.0*(*r)*(*phi)/(*S));
  propPars[0]= .5*(*m - eps); propPars[1]= .5*(*m + eps);
  //Find density at the mode
  for (i=0, fmode=0; i< *n; i++) fmode+= dnormC(e[i],propPars[1]*xj[i],sqrtPhi,1);
  fmode+= dmom(propPars[1],0,*tau1,*phi,*r,1) - *m1;
  fmode= exp(fmode);
  //Proposal variances
  temp= (*S)/(*phi);
  propPars[2]= sqrt(1.0/(temp + doubler/(propPars[0]*propPars[0]))); propPars[3]= sqrt(1.0/(temp + doubler/(propPars[1]*propPars[1])));
  temp2= .5*(*nu); temp= temp2+.5;
  //Weights
  ct2= exp(gamln(&temp) - .5*log((double) (*nu)) - gamln(&temp2) - .5*log(M_PI*propPars[3]*propPars[3]));
  propPars[4]= max_xy(0,(fmode-ct2)/(dnormC(propPars[1],propPars[0],propPars[2],0) - ct2));
  propPars[5]= 1-propPars[4];
}



//Expectation of prod_j th_j^{2*power} under multivariate Normal/T with mean m, covariance S, dimension n, degrees of freedom dof (dof=-1 for Normal)

// [[Rcpp::export]]
SEXP eprod_I(SEXP m, SEXP S, SEXP n, SEXP power, SEXP dof) {
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  *REAL(ans)= mvtexpect_vec(REAL(m), REAL(S), INTEGER(n)[0], INTEGER(power)[0], REAL(dof)[0]);
  UNPROTECT(1);
  return ans;
}



//Univariate marginal density under a product MOM prior (known variance case)
// integral N(y; x*theta, phi*I) * (theta^2/(tau*phi))^r * N(theta; 0; tau*phi) / (2r-1)!! d theta
// - y: response variable (must be a vector)
// - x: design matrix (must be a vector)
// - sumy2: sum(y*y)
// - n: length of y
// - phi: residual variance
// - tau: prior variance parameter
// - logscale: if set to 1 the log of the integral is returned
double pmomMargKuniv(double *y, double *x, double *sumy2, double *sumxsq, int *n, double *phi, double *tau, int *r, int *logscale) {
  int i; double ans, m, s, I, doubler=2.0*(*r);
  s= *sumxsq + 1.0/(*tau);
  for (i=0, m=0; i< *n; i++) { m+= y[i]*x[i]; }
  m/= s;
  I= log(mnorm(doubler,m,sqrt(*phi/s)));
  ans= I -.5*(*sumy2 - s*m*m)/(*phi) - .5*(*n)*log(2*M_PI*(*phi)) - .5*(log(s)+log(*tau)) - ldoublefact(doubler-1) - (*r)*log((*tau)*(*phi));
  if (*logscale ==0) ans=exp(ans);
  return(ans);
}

//Sample from conditional posterior of coefficients for variables not undergoing selection
//Input:
// - partialres: partial residuals obtained by removing adjustment variables from the model
// - phi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - theta2: sample from conditional posterior of theta2 given the data and all other parameters
void simTheta2(double *theta2, double *partialres, double *phi, struct modavgPars *pars) {
  int i, j; double *tmp, *m, **cholS, sqrtPhi=sqrt(*phi);
  //m= S2inv * t(x2) * partialres
  tmp= dvector(0,*(*pars).p2); m= dvector(0,*(*pars).p2); cholS= dmatrix(1,*(*pars).p2,1,*(*pars).p2);
  Atvecx((*pars).x2, partialres, tmp, 0, *(*pars).p2 -1, 0, *(*pars).n -1);
  Avecx((*pars).S2inv, tmp, m, 0, *(*pars).p2, 0, *(*pars).p2);
  //S= S2inv * phi
  for (i=0; i< *(*pars).p2; i++) { for (j=0; j< *(*pars).p2; j++) { cholS[i+1][j+1]= sqrtPhi * (*pars).cholS2inv[i+j*(*(*pars).p2)]; } }
  //Generate theta2 ~ N(m,S)
  rmvnormC(theta2-1,*(*pars).p2,m-1,cholS);
  free_dvector(tmp,0,*(*pars).p2); free_dvector(m,0,*(*pars).p2); free_dmatrix(cholS,1,*(*pars).p2,1,*(*pars).p2);
}


//Sample from conditional posterior of residual variance phi under a product MOM prior
//Input:
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curCoef2: current values of the regression coefficients for variables not undergoing selection
// - ssr: residual sum of squares for current coefficient values
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of phi given the data and all other parameteres
double simPhipmom(int *nsel, int *curModel, double *curCoef1, double *curCoef2, double *ssr, struct modavgPars *pars) {
  int i; double a, b, sumth1, sumth2;
  a= *(*pars).alpha + *(*pars).n + (2*(*(*pars).r)+1)*(*nsel) + *(*pars).p2;
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  for (i=0, sumth2=0; i< *(*pars).p2; i++) { sumth2+= curCoef2[i]*curCoef2[i]; }
  b= *(*pars).lambda + sumth1/(*(*pars).tau1) + sumth2/(*(*pars).tau2) + *ssr;
  return(1.0/rgammaC(.5*a,.5*b));
}

//Sample from conditional posterior of tau given all other parameters
//Input:
// - nsel: number of variables undergoing selection currently in the model
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value for the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of tau given the data and all other parameters
double simTaupmom(int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int i; double a, b, sumth1;
  a= *(*pars).atau1 + (2*(*(*pars).r)+1)*(*nsel);
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  b= *(*pars).btau1 + sumth1/(*curPhi);
  return(1.0/rgammaC(.5*a,.5*b));
}


//********************************************************************************************
// GENERAL MARGINAL DENSITY CALCULATION ROUTINES
//********************************************************************************************

void set_marginalPars(struct marginalPars *pars, int *family, int *n,int *nuncens,int *p,double *y,int *uncens,double *sumy2,double *sumy,double *sumlogyfact,double *x,double *colsumsx,crossprodmat *XtX,double *ytX,int *method,int *adjoverdisp,int *hesstype,int *optimMethod,int *optim_maxit,int *usethinit,double *thinit,int *B,double *alpha,double *lambda,int *knownphi,double *phi,double *tau,double *taugroup,double *taualpha, double *fixatanhalpha, int *r,double *prDeltap,double *parprDeltap, double *prConstrp,double *parprConstrp, int *logscale, double *offset, int *groups, int *isgroup, int *ngroups, int *ngroupsconstr, int *nvaringroup, int *nconstraints=NULL, int *ninvconstraints=NULL, crossprodmat *XtXuncens=NULL, double *ytXuncens=NULL) {
  (*pars).family= family;
  (*pars).n= n;
  (*pars).nuncens= nuncens;
  (*pars).p= p;
  (*pars).y= y;
  (*pars).uncens= uncens;
  (*pars).sumy2= sumy2;
  (*pars).sumy= sumy;
  (*pars).sumlogyfact= sumlogyfact;
  (*pars).x= x;
  (*pars).colsumsx= colsumsx;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).XtXuncens= XtXuncens;
  (*pars).ytXuncens= ytXuncens;
  (*pars).V0inv = new covariancemat(*p);
  (*pars).method= method;
  (*pars).adjoverdisp= adjoverdisp;
  (*pars).hesstype= hesstype;
  (*pars).optimMethod= optimMethod;
  (*pars).optim_maxit= optim_maxit;
  (*pars).usethinit= usethinit;
  (*pars).thinit= thinit;
  (*pars).B= B;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).knownphi= knownphi;
  (*pars).phi= phi;
  (*pars).tau= tau;
  (*pars).taugroup= taugroup;
  (*pars).taualpha= taualpha;
  (*pars).fixatanhalpha= fixatanhalpha;
  (*pars).r= r;
  (*pars).prDeltap= prDeltap;
  (*pars).parprDeltap= parprDeltap;
  (*pars).prConstrp= prConstrp;
  (*pars).parprConstrp= parprConstrp;
  (*pars).logscale= logscale;
  (*pars).offset= offset;
  (*pars).groups= groups;
  (*pars).isgroup= isgroup;
  (*pars).ngroups= ngroups;
  (*pars).ngroupsconstr= ngroupsconstr;
  (*pars).nvaringroup= nvaringroup;
  (*pars).nconstraints= nconstraints;
  (*pars).ninvconstraints= ninvconstraints;
}

void set_f2opt_pars(double *m, double **S, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel) {
  f2opt_pars.m= m;
  f2opt_pars.S= S;
  f2opt_pars.sumy2= sumy2;
  f2opt_pars.XtX= XtX;
  f2opt_pars.ytX= ytX;
  f2opt_pars.alpha= alpha;
  f2opt_pars.lambda= lambda;
  f2opt_pars.phi= phi;
  f2opt_pars.tau= tau;
  f2opt_pars.r= r;
  f2opt_pars.n= n;
  f2opt_pars.p= p;
  f2opt_pars.sel= sel;
  f2opt_pars.nsel= nsel;
}

void set_f2int_pars(crossprodmat *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale) {
  f2int_pars.XtX= XtX;
  f2int_pars.ytX= ytX;
  f2int_pars.tau= tau;
  f2int_pars.n= n;
  f2int_pars.p= p;
  f2int_pars.sel= sel;
  f2int_pars.nsel= nsel;
  f2int_pars.y= y;
  f2int_pars.sumy2= sumy2;
  f2int_pars.method= method;
  f2int_pars.B= B;
  f2int_pars.alpha= alpha;
  f2int_pars.lambda= lambda;
  f2int_pars.logscale= logscale;
}






//********************************************************************************************
// MODEL SELECTION ROUTINES
//********************************************************************************************

//modelSelectionEnum: model selection via enumeration in linear regression for several choices of prior distribution
//Input parameters
// - nmodels: number of models to be enumerated (rows in models)
// - models: binary matrix containing all models to be enumerated, structured as a vector
// - knownphi: is residual variance phi known?
// - priorCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
// - priorDelta: 0 for uniform, 1 for binomial, 2 for binomial with beta hyper-prior for success prob
// - niter: number of Gibbs iterations
// - ndeltaini: length of deltaini
// - deltaini: vector with indexes of covariates initially in the model (both deltaini and its indexes must be indexed at 0)
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
// - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
// - family: residual distribution (1 for Normal; 2 for two-piece Normal; 3 for Laplace; 4 for two-piece Laplace). Set family==0 to perform inference on the family
//Output
// - margpp: marginal posterior probability for inclusion of each covariate
// - postMode: model with highest posterior probability amongst those enumerated
// - postModeProb: unnormalized posterior prob of posterior mode (log scale)
// - postProb: unnormalized posterior prob of each visited model (log scale)

// [[Rcpp::export]]
SEXP modelSelectionEnumCI(SEXP Snmodels, SEXP Smodels, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose) {

  bool hasXtX= LOGICAL(ShasXtX)[0];
  int i, j, idxj, logscale=1, *postMode, mycols, mycols2, nuncens, *nconstraints, *ninvconstraints, ngroupsconstr=0, *isgroup, usethinit= INTEGER(Susethinit)[0], priorcode;
  double offset=0, *postModeProb, *postProb, *ytXuncens=NULL, *thinit=NULL;
  intptrvec constraints, invconstraints;
  crossprodmat *XtX, *XtXuncens=NULL;
  struct marginalPars pars;
  SEXP ans;

  PROTECT(ans= Rf_allocVector(VECSXP, 3));
  if (INTEGER(Sfamily)[0] !=0) { mycols= mycols2= INTEGER(Sp)[0]; } else { mycols= 2 + INTEGER(Sp)[0]; mycols2= mycols+2; }

  //Allocate memory to thinit, to safeguard against cases where subroutines use it, despite usethinit==0
  thinit= dvector(0, mycols2+1);
  if (usethinit != 3) {
    for (j=0; j<= mycols2+1; j++) { thinit[j]= 0; }
  } else {
    for (j=0; j<= INTEGER(Sp)[0]; j++) { thinit[j]= REAL(Sthinit)[j]; }
  }


  SET_VECTOR_ELT(ans, 0, Rf_allocVector(INTSXP, mycols));
  postMode= INTEGER(VECTOR_ELT(ans,0));

  SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, 1));
  postModeProb= REAL(VECTOR_ELT(ans,1));

  SET_VECTOR_ELT(ans, 2, Rf_allocVector(REALSXP, INTEGER(Snmodels)[0]));
  postProb= REAL(VECTOR_ELT(ans,2));

  isgroup= ivector(0, INTEGER(Sp)[0]);
  nconstraints= ivector(0,INTEGER(Sngroups)[0]); ninvconstraints= ivector(0,INTEGER(Sngroups)[0]);
  countConstraints(nconstraints, &constraints, ninvconstraints, &invconstraints, &ngroupsconstr, isgroup, INTEGER(Sngroups), INTEGER(Snvaringroup), Sconstraints, Sinvconstraints);

  if (hasXtX) {
    XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  } else {
    XtX= new crossprodmat(REAL(Sx),INTEGER(Sn)[0],INTEGER(Sp)[0],false);
  }

  if (LENGTH(Suncens)>0) { //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
    int n=INTEGER(Sn)[0], *uncens= INTEGER(Suncens);
    double *pty= REAL(Sy), *ptx= REAL(Sx);
    for (nuncens=0; (nuncens<n) && (uncens[nuncens]==1); nuncens++) { }
    XtXuncens= new crossprodmat(REAL(Sx), INTEGER(Sn)[0], INTEGER(Sp)[0], false, nuncens, 0);
    ytXuncens= dvector(0,INTEGER(Sp)[0]);
    for (j=0; j< INTEGER(Sp)[0]; j++) { for (i=0, ytXuncens[j]=0, idxj=j*n; i< nuncens; i++) { ytXuncens[j] += pty[i] * ptx[i + idxj]; } }
  } else { nuncens= INTEGER(Sn)[0]; }


  set_marginalPars(&pars, INTEGER(Sfamily), INTEGER(Sn), &nuncens, INTEGER(Sp), REAL(Sy), INTEGER(Suncens), REAL(Ssumy2), REAL(Ssumy), REAL(Ssumlogyfact), REAL(Sx), REAL(Scolsumsx), XtX, REAL(SytX), INTEGER(Smethod), INTEGER(Sadjoverdisp), INTEGER(Shesstype), INTEGER(SoptimMethod), INTEGER(Soptim_maxit), &usethinit, thinit, INTEGER(SB), REAL(Salpha),REAL(Slambda), INTEGER(Sknownphi), REAL(Sphi), REAL(Stau), REAL(Staugroup), REAL(Staualpha), REAL(Sfixatanhalpha), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), REAL(SprConstrp), REAL(SparprConstrp), &logscale, &offset, INTEGER(Sgroups), isgroup, INTEGER(Sngroups), &ngroupsconstr, INTEGER(Snvaringroup), nconstraints, ninvconstraints, XtXuncens, ytXuncens);

  priorcode = mspriorCode(INTEGER(SpriorCoef), INTEGER(SpriorGroup), &pars);
  pars.priorcode= &priorcode;

  modelSelectionEnum(postMode, postModeProb, postProb, INTEGER(Snmodels), INTEGER(Smodels), INTEGER(SpriorDelta), INTEGER(SpriorConstr), INTEGER(Sverbose), &pars);

  delete XtX;
  free_ivector(nconstraints, 0,INTEGER(Sngroups)[0]); free_ivector(ninvconstraints, 0,INTEGER(Sngroups)[0]);
  free_dvector(thinit, 0,mycols2+1); free_ivector(isgroup, 0, INTEGER(Sp)[0]);
  if (LENGTH(Suncens)>0) { delete XtXuncens; free_dvector(ytXuncens,0,INTEGER(Sp)[0]); }
  UNPROTECT(1);
  return ans;
}



void modelSelectionEnum(int *postMode, double *postModeProb, double *postProb, int *nmodels, int *models, int *prDelta, int *prConstr, int *verbose, struct marginalPars *pars) {

  int i, j, *sel, nsel, nselplus1, niter10, nbvars, nbfamilies=4, postModeidx, *family= (*pars).family;
  double *mfamily, *pfamily;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);
  modselIntegrals *integrals;

  marginalFunction= set_marginalFunction(pars);
  priorFunction= set_priorFunction(prDelta, prConstr, family);

  mfamily= dvector(0,nbfamilies-1); pfamily= dvector(0,nbfamilies-1);
  if ((*family)==0) {  //infer error distribution (Normal/Asymmetric Normal/Laplace/Asymmetric Laplace)
    nbvars= (*(*pars).p)+1;
    integrals= new modselIntegrals(marginalFunction, priorFunction, (*(*pars).p) +4);
  } else {
    nbvars= (*(*pars).p);
    integrals= new modselIntegrals(marginalFunction, priorFunction, (*(*pars).p));
  }

  sel= ivector(0,nbvars);
  if (*verbose ==1) Rprintf(" Computing posterior probabilities\n");
  if (*nmodels >10) { niter10= (*nmodels)/10; } else { niter10= 1; }

  postModeidx= 0;
  (*postModeProb)= R_NegInf;
  //Iterate
  for (i=0; i< *nmodels; i++) {
    nsel= 0;
    for (j=0; j< *(*pars).p; j++) { if (models[(*nmodels)*j + i]==1) { sel[nsel]= j; nsel++; } }
    if (nsel <= (*(*pars).n)) {
      if ((*family)==0) {  //inference is being done on the family
        sel[nsel]= (*(*pars).p) + models[(*nmodels)*(*(*pars).p) +i] + 2*models[(*nmodels)*nbvars +i];
        nselplus1= nsel+1;
        postProb[i]= integrals->getJoint(sel,&nselplus1,pars);
      } else {  //family is fixed
        postProb[i]= integrals->getJoint(sel,&nsel,pars);
      }
      if (postProb[i] > *postModeProb) { (*postModeProb)= postProb[i]; postModeidx= i; }
    }
    if (*verbose == 1) print_iterprogress(&i, nmodels, &niter10);
    //if ((*verbose==1) && ((i%niter10)==0)) { Rprintf("."); }
  }

  //Store posterior mode
  for (j=0; j< *(*pars).p; j++) { postMode[j]= models[(*nmodels)*j + postModeidx]; }
  if ((*family)==0) {
    for (j= *(*pars).p; j< (*(*pars).p)+2; j++) { postMode[j]= models[(*nmodels)*j + postModeidx]; }
  }

  if (*verbose == 1) { Rcout << "\r Done\n"; }
  //if (*verbose==1) Rprintf(" Done.\n");

  free_ivector(sel,0,nbvars); free_dvector(mfamily,0,nbfamilies-1); free_dvector(pfamily,0,nbfamilies-1);
  delete integrals;
}




//modelSelectionGibbs: Gibbs sampler for model selection in linear regression for several choices of prior distribution
//Input parameters
// - knownphi: is residual variance phi known?
// - priorCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
// - priorDelta: 0 for uniform, 1 for binomial, 2 for binomial with beta hyper-prior for success prob
// - niter: number of Gibbs iterations
// - ndeltaini: length of deltaini
// - deltaini: vector with indexes of covariates initially in the model (both deltaini and its indexes must be indexed at 0)
// - includevars: binary vector, elements with includevars[j]==1 indicate that variable should always be included in the model
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
// - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
// - family: residual distribution (1 for Normal; 2 for two-piece Normal; 3 for Laplace; 4 for two-piece Laplace). Set family==0 to perform inference on the family
//Output
// - postSample: matrix with niter rows and p columns with posterior samples for covariate inclusion/exclusion (formatted as a vector in column order)
// - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
// - postMode: model with highest posterior probability amongst all those visited
// - postModeProb: unnormalized posterior prob of posterior mode (log scale)
// - postProb: unnormalized posterior prob of each visited model (log scale)

// [[Rcpp::export]]
SEXP modelSelectionGibbsCI(SEXP SpostModeini, SEXP SpostModeiniProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose) {

  bool hasXtX= LOGICAL(ShasXtX)[0];
  int i, j, idxj, logscale=1, mcmc2save, *postSample, *postMode, mycols, mycols2, *nconstraints, *ninvconstraints, nuncens, ngroupsconstr=0, *isgroup, usethinit=INTEGER(Susethinit)[0], priorcode;
  double offset=0, *margpp, *postModeProb, *postProb, *ytXuncens=NULL, *thinit;
  intptrvec constraints, invconstraints;
  crossprodmat *XtX, *XtXuncens=NULL;
  struct marginalPars pars;
  SEXP ans;

  PROTECT(ans= Rf_allocVector(VECSXP, 5));
  mcmc2save= floor((INTEGER(Sniter)[0] - INTEGER(Sburnin)[0] +.0)/(INTEGER(Sthinning)[0] +.0));
  if (INTEGER(Sfamily)[0] !=0) { mycols= mycols2= INTEGER(Sp)[0]; } else { mycols= 2 + INTEGER(Sp)[0]; mycols2= mycols+2; }

  thinit= dvector(0, mycols2+1);
  if (usethinit != 3) {
    for (j=0; j<= mycols2+1; j++) { thinit[j]= 0; }
  } else {
    for (j=0; j<= INTEGER(Sp)[0]; j++) { thinit[j]= REAL(Sthinit)[j]; }
  }

  SET_VECTOR_ELT(ans, 0, Rf_allocVector(INTSXP, mcmc2save * mycols));
  postSample= INTEGER(VECTOR_ELT(ans,0));
  for (j=0; j<(mcmc2save*mycols); j++) postSample[j]= 0;

  SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, mycols2));
  margpp= REAL(VECTOR_ELT(ans,1));

  SET_VECTOR_ELT(ans, 2, Rf_allocVector(INTSXP, mycols));
  postMode= INTEGER(VECTOR_ELT(ans,2));
  for (j=0; j<mycols; j++) { postMode[j]= INTEGER(SpostModeini)[j]; }

  SET_VECTOR_ELT(ans, 3, Rf_allocVector(REALSXP, 1));
  postModeProb= REAL(VECTOR_ELT(ans,3));
  postModeProb[0]= REAL(SpostModeiniProb)[0];

  SET_VECTOR_ELT(ans, 4, Rf_allocVector(REALSXP, mcmc2save));
  postProb= REAL(VECTOR_ELT(ans,4));

  isgroup= ivector(0, INTEGER(Sp)[0]);
  nconstraints= ivector(0,INTEGER(Sngroups)[0]); ninvconstraints= ivector(0,INTEGER(Sngroups)[0]);
  countConstraints(nconstraints, &constraints, ninvconstraints, &invconstraints, &ngroupsconstr, isgroup, INTEGER(Sngroups), INTEGER(Snvaringroup), Sconstraints, Sinvconstraints);

  if (hasXtX) {
    XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  } else {
    XtX= new crossprodmat(REAL(Sx),INTEGER(Sn)[0],INTEGER(Sp)[0],false);
  }

  if (LENGTH(Suncens)>0) { //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
    int n=INTEGER(Sn)[0], *uncens= INTEGER(Suncens);
    double *pty= REAL(Sy), *ptx= REAL(Sx);
    for (nuncens=0; (nuncens<n) && (uncens[nuncens]==1); nuncens++) { }
    XtXuncens= new crossprodmat(REAL(Sx), INTEGER(Sn)[0], INTEGER(Sp)[0], false, nuncens, 0);
    ytXuncens= dvector(0,INTEGER(Sp)[0]);
    for (j=0; j< INTEGER(Sp)[0]; j++) { for (i=0, ytXuncens[j]=0, idxj=j*n; i< nuncens; i++) { ytXuncens[j] += pty[i] * ptx[i + idxj]; } }
  } else { nuncens= INTEGER(Sn)[0]; }

  set_marginalPars(&pars, INTEGER(Sfamily), INTEGER(Sn), &nuncens, INTEGER(Sp), REAL(Sy), INTEGER(Suncens), REAL(Ssumy2), REAL(Ssumy), REAL(Ssumlogyfact), REAL(Sx), REAL(Scolsumsx), XtX, REAL(SytX), INTEGER(Smethod), INTEGER(Sadjoverdisp), INTEGER(Shesstype), INTEGER(SoptimMethod), INTEGER(Soptim_maxit), &usethinit, thinit, INTEGER(SB), REAL(Salpha),REAL(Slambda), INTEGER(Sknownphi), REAL(Sphi), REAL(Stau), REAL(Staugroup), REAL(Staualpha), REAL(Sfixatanhalpha), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), REAL(SprConstrp), REAL(SparprConstrp), &logscale, &offset, INTEGER(Sgroups), isgroup, INTEGER(Sngroups), &ngroupsconstr, INTEGER(Snvaringroup), nconstraints, ninvconstraints, XtXuncens, ytXuncens);

  priorcode = mspriorCode(INTEGER(SpriorCoef), INTEGER(SpriorGroup), &pars);
  pars.priorcode= &priorcode;

  modelSelectionGibbs(postSample, margpp, postMode, postModeProb, postProb, INTEGER(SpriorDelta), INTEGER(SpriorConstr), INTEGER(Sniter), INTEGER(Sthinning), INTEGER(Sburnin), INTEGER(Sndeltaini), INTEGER(Sdeltaini), INTEGER(Sincludevars), &constraints, &invconstraints, INTEGER(Sverbose), &pars);

  free_dvector(thinit, 0,mycols2+1);   free_ivector(isgroup, 0, INTEGER(Sp)[0]);
  free_ivector(nconstraints, 0,INTEGER(Sngroups)[0]); free_ivector(ninvconstraints, 0,INTEGER(Sngroups)[0]);
  delete XtX;
  UNPROTECT(1);
  return ans;
}




void modelSelectionGibbs(int *postSample, double *margpp, int *postMode, double *postModeProb, double *postProb, int *prDelta, int *prConstr, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *includevars, intptrvec *constraints, intptrvec *invconstraints, int *verbose, struct marginalPars *pars) {

  bool copylast, validmodel, nonbinary=false;
  int i, j, jgroup, jgroup2, k, niter10, niterthin, savecnt, ilow, iupper, itercount, nbvars, nbfamilies=4, curfamily, newfamily, ngroups, *nvaringroup, *firstingroup, *nconstraints= (*pars).nconstraints, *ninvconstraints= (*pars).ninvconstraints, *family= (*pars).family;
  int *addgroups, *dropgroups, naddgroups, ndropgroups, *modelidx, *sample;
  int nsel, nselnew, nselnew2, nselnew3, *sel, *selnew, *selnew2, *selnew3, nselplus1, *selaux;
  double currentJ, *newJ, *ppnew, ppnewsum, u, *mfamily, *pfamily, sumpfamily;
  //intptrvec::iterator itlist;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);
  modselIntegrals *integrals;

  marginalFunction= set_marginalFunction(pars);
  priorFunction= set_priorFunction(prDelta, prConstr, family);

  addgroups= ivector(0,1); dropgroups= ivector(0,1); newJ= dvector(0,2); ppnew= dvector(0,3); modelidx= ivector(0,3);

  mfamily= dvector(0,nbfamilies-1); pfamily= dvector(0,nbfamilies-1);
  if ((*family)==0) {
    nbvars= (*(*pars).p)+1;
    integrals= new modselIntegrals(marginalFunction, priorFunction, (*(*pars).p) +4);
    copylast= true;
  } else {
    nbvars= (*(*pars).p);
    integrals= new modselIntegrals(marginalFunction, priorFunction, (*(*pars).p));
    copylast= false;
  }
  ngroups= *((*pars).ngroups);
  nvaringroup= (*pars).nvaringroup;
  firstingroup= ivector(0,ngroups);
  for (j=1, firstingroup[0]=0; j<ngroups; j++) { firstingroup[j]= firstingroup[j-1] + nvaringroup[j-1]; }

  //Determine if there are groups that need to be jointly sampled
  //By setting jointlysample== true the Gibbs algorithm is modified to sample jointly pairs of groups related by constraints
  //Example: group 1 is linear effect of a variable, group 2 its non-linear effect.
  //         If jointlysample== false then each group indicator is sampled separately
  //         If jointlysample== true then both group indicators are sample jointly, considering 3 possible models: (no effect, only linear, both linear & non-linear)
  //By default this option is disabled, its rarely necessary and can slow sampling (it forces considering non-linear terms at every iter, rather than only when the linear term is in)
  bool jointlysample= false;
  sample= ivector(0,ngroups-1);
  if (jointlysample) {
    for (j=0; j < ngroups; j++) sample[j]= -1;
    for (j=0; j < ngroups; j++) {
      if ((nconstraints[j]==0) && (ninvconstraints[j]==1) && (includevars[firstingroup[j]]==0)) {
        k= (*invconstraints)[j][0];
        if (nconstraints[k]==1) { sample[j]= 2; sample[k]= 0; nonbinary= true; }
      } else if (sample[j]== -1) {
        sample[j]= 1;
      }
    }
  } else {
    nonbinary= false;
    for (j=0; j < ngroups; j++) sample[j]= 1;
  }

  sel= ivector(0,nbvars); selnew= ivector(0,nbvars);
  if (nonbinary) { selnew2= ivector(0,nbvars); selnew3= ivector(0, nbvars); }

  //Initialize
  if (*verbose ==1) Rprintf(" Running Gibbs sampler");
  niterthin= (int) floor((*niter - *burnin +.0)/(*thinning+.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  for (j=0; j< (*(*pars).p); j++) { margpp[j]= 0; }
  nsel= *ndeltaini;
  for (j=0; j< nsel; j++) { sel[j]= deltaini[j]; postMode[deltaini[j]]= 1; }
  if ((*family)==0) {
    sel[nsel]= (*(*pars).p);  //initialize to baseline model
    nselplus1= nsel+1;
    currentJ= integrals->getJoint(sel,&nselplus1,pars);
    margpp[(*(*pars).p)]= margpp[(*(*pars).p)+1]= margpp[(*(*pars).p)+2]= margpp[(*(*pars).p)+3]= 0;
  } else {
    currentJ= integrals->getJoint(sel,&nsel,pars);
  }
  postProb[0]= *postModeProb= currentJ;
  if (*burnin >0) {
    ilow=-(*burnin); savecnt=0; iupper= *niter - *burnin +1;
  } else {
    for (j=0; j<nsel; j++) { postSample[sel[j]*niterthin]= 1; }
    ilow=1; savecnt=1; iupper= *niter;
  } //if no burnin, start at i==1 & save initial value


  //Iterate
  for (i=ilow, itercount=0; i< iupper; i++, itercount++) {
    j= jgroup= 0;
    while (j< *(*pars).p) {

      if (sample[jgroup]>0) { //if jgroup should be sampled

        sel2selnew(jgroup,sel,&nsel,selnew,&nselnew,copylast,&ngroups,nvaringroup,firstingroup); //copy sel into selnew, adding/removing jth group
        if (nsel > nselnew) { naddgroups= 0; ndropgroups=1; dropgroups[0]= jgroup; } else { naddgroups=1 ; ndropgroups=0; addgroups[0]= jgroup; }
        validmodel= checkConstraints(addgroups,&naddgroups,dropgroups,&ndropgroups,constraints,nconstraints,invconstraints,ninvconstraints,(*pars).groups,nvaringroup,sel,&nsel);
        if (nselnew > (*(*pars).n)) validmodel= false;

        ppnew[0]= ppnewsum= 1;
        if (includevars[j]==0 && validmodel) { //if proposed model is valid
            if ((*family)==0) { nselplus1= nselnew+1; newJ[0]= integrals->getJoint(selnew,&nselplus1,pars); } else { newJ[0]= integrals->getJoint(selnew,&nselnew,pars); }
            if (newJ[0] > *postModeProb) { *postModeProb= newJ[0];  update_postMode(postMode, nselnew, selnew, *(*pars).p, *family); } //update posterior mode
            ppnew[1]= exp(newJ[0]-currentJ);
            ppnewsum+= ppnew[1];
        } else {
          ppnew[1]= -0.1;
        }

        if ((nonbinary) && sample[jgroup]==2) { //non-binary update

          jgroup2= (*invconstraints)[jgroup][0];
          sel2selnew(jgroup2,sel,&nsel,selnew2,&nselnew2,copylast,&ngroups,nvaringroup,firstingroup); //Model changing inclusion/exclusion of jgroup2

          sel2selnew(jgroup2,selnew,&nselnew,selnew3,&nselnew3,copylast,&ngroups,nvaringroup,firstingroup); //Model changing inclusion/exclusion of (jgroup,jgroup2)
          if ((nsel > nselnew) && (nselnew > nselnew3)) {         //nselnew3 dropped jgroup, dropped jgroup2
            modelidx[0]= 2; modelidx[1]= -1; modelidx[2]= 1; modelidx[3]= 0;
          } else if ((nsel > nselnew) && (nselnew < nselnew3)) {  //nselnew3 dropped jgroup, added jgroup2
            modelidx[0]= 1; modelidx[1]= 0; modelidx[2]= 2; modelidx[3]= -1;
            //  } else if ((nsel < nselnew) && (nselnew > nselnew3)) {  //nselnew3 added jgroup, dropped jgroup2 (this case cannot happen, it would mean that current model was invalid)
      //    modelidx[0]= -1; modelidx[1]= 2; modelidx[2]= 0; modelidx[3]= 1;
          } else {                                                //nselnew3 added jgroup, added jgroup2
            modelidx[0]= 0; modelidx[1]= 1; modelidx[2]= -1; modelidx[3]= 2;
          }

          if ((modelidx[2]!= -1) && (nselnew2 < (*(*pars).n))) {
            if ((*family)==0) { nselplus1= nselnew2+1; newJ[1]= integrals->getJoint(selnew2,&nselplus1,pars); } else { newJ[1]= integrals->getJoint(selnew2,&nselnew2,pars); }
            if (newJ[1] > *postModeProb) { *postModeProb= newJ[1];  update_postMode(postMode, nselnew2, selnew2, *(*pars).p, *family); }
            ppnew[2]= exp(newJ[1]-currentJ);
            ppnewsum+= ppnew[2];
          } else {
            ppnew[2]= -0.1;
          }

          if ((modelidx[3]!= -1) && (nselnew3 < (*(*pars).n))) {
            if ((*family)==0) { nselplus1= nselnew3+1; newJ[2]= integrals->getJoint(selnew3,&nselplus1,pars); } else { newJ[2]= integrals->getJoint(selnew3,&nselnew3,pars); }
            if (newJ[2] > *postModeProb) { *postModeProb= newJ[2];  update_postMode(postMode, nselnew3, selnew3, *(*pars).p, *family); }
            ppnew[3]= exp(newJ[2]-currentJ);
            ppnewsum+= ppnew[3];
          } else {
            ppnew[3]= -0.1;
          }

          //update model
          u= runif();
          ppnew[0] /= ppnewsum; ppnew[1] /= ppnewsum; ppnew[2] /= ppnewsum; ppnew[3] /= ppnewsum;
          if (u< ppnew[1]) {
            selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; currentJ= newJ[0];
          } else {
            ppnewsum= max_xy(ppnew[1],0.0);
            if (u < ppnewsum + ppnew[2]) {
              selaux= sel; sel=selnew2; selnew2=selaux; nsel=nselnew2; currentJ= newJ[1];
            } else {
              ppnewsum += max_xy(ppnew[2],0.0);
              if (u < ppnewsum + ppnew[3]) { selaux= sel; sel=selnew3; selnew3=selaux; nsel=nselnew3; currentJ= newJ[2]; }
            }
          }

          //update Rao-Blackwellized marginal inclusion probabilities
          for (k=0; k<=3; k++) {
            if (modelidx[k]==1) {
              margpp[firstingroup[jgroup]] += ppnew[k];
            } else if (modelidx[k]==2) {
              margpp[firstingroup[jgroup]] += ppnew[k];
              margpp[firstingroup[jgroup2]] += ppnew[k];
            }
          }

        } else {  //binary update

          if (includevars[j]==0 && validmodel) { //if proposed model is valid

            ppnew[1] /= ppnewsum;
            if (i>=0) { if (nselnew>nsel) { margpp[j]+= ppnew[1]; } else { margpp[j]+= (1-ppnew[1]); } } //update Rao-Blackwellized inclusion probabilities
            u= runif();
            if (u < ppnew[1]) {  selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; currentJ= newJ[0]; } //update model

          } else {

            if ((i>=0) && (ndropgroups>0) && (!validmodel)) { margpp[j] += 1.0 ; } //conditional marginal inclusion prob=1

          }

        }

      } //end if jgroup should be sampled

      j += nvaringroup[jgroup];
      jgroup++;

    }  //end j for

    if ((*family)==0) {  //update family
      nselplus1= nsel+1;
      curfamily= sel[nsel] - (*(*pars).p);
      sumpfamily= 0;
      for (j=0; j<nbfamilies; j++) {
        if (j==curfamily) {
          mfamily[j]= currentJ;
          pfamily[j]= 1.0;
          sumpfamily += 1.0;
        } else {
          sel[nsel]= (*(*pars).p) + j;
          mfamily[j]= integrals->getJoint(sel,&nselplus1,pars);
          pfamily[j]= exp(mfamily[j] - currentJ);
          sumpfamily += pfamily[j];
        }
      }
      for (j=0; j<nbfamilies; j++) { pfamily[j] /= sumpfamily; }
      rmultinomial(1, nbfamilies, pfamily, &newfamily);
      sel[nsel]= (*(*pars).p) + newfamily;
      if (i>=0) {
        for (j=0; j<nbfamilies; j++) margpp[(*(*pars).p) +j]+= pfamily[j];
      }
      currentJ= mfamily[newfamily];
    }

    if ((i>0) && ((i%(*thinning))==0)) {
      for (j=0; j<nsel; j++) { postSample[sel[j]*niterthin+savecnt]= 1; }
      if ((*family)==0) {
        if (sel[nsel]== (*(*pars).p)) { //Normal residuals
          postSample[(*(*pars).p)*niterthin + savecnt]= 0;
          postSample[((*(*pars).p +1))*niterthin + savecnt]= 0;
        } else if (sel[nsel]== (*(*pars).p) +1) { //Asymmetric Normal residuals
          postSample[(*(*pars).p)*niterthin + savecnt]= 1;
          postSample[((*(*pars).p +1))*niterthin + savecnt]= 0;
        } else if (sel[nsel]== (*(*pars).p) +2) { //Laplace residuals
          postSample[(*(*pars).p)*niterthin + savecnt]= 0;
          postSample[((*(*pars).p +1))*niterthin + savecnt]= 1;
        } else { //Asymmetric Laplace residuals
          postSample[(*(*pars).p)*niterthin + savecnt]= 1;
          postSample[((*(*pars).p +1))*niterthin + savecnt]= 1;
        }
      }
      postProb[savecnt]= currentJ;
      savecnt++;
    }
    if (*verbose == 1) print_iterprogress(&itercount, niter, &niter10);
    //if ((*verbose==1) && ((i%niter10)==0)) { Rprintf("."); }
  }
  if (iupper>ilow) { for (j=0; j< (*(*pars).p); j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*family ==0) {
    margpp[(*(*pars).p)] /= (iupper-imax_xy(0,ilow)+.0); margpp[(*(*pars).p)+1] /= (iupper-imax_xy(0,ilow)+.0); margpp[(*(*pars).p)+2] /= (iupper-imax_xy(0,ilow)+.0); margpp[(*(*pars).p)+3] /= (iupper-imax_xy(0,ilow)+.0);
  }
  if (*verbose) { Rcout << "\r Done\n"; }
  //if (*verbose==1) Rprintf(" Done.\n");

  //Copy marginal probabilities from 1st variable in each group to all other variables in the group
  for (jgroup=0; jgroup<ngroups; jgroup++) {
    for (j=1; j<nvaringroup[jgroup]; j++) {
      margpp[firstingroup[jgroup]+j]= margpp[firstingroup[jgroup]];
    }
  }

  free_ivector(addgroups, 0,1); free_ivector(dropgroups, 0,1); free_dvector(newJ, 0,2); free_dvector(ppnew, 0,3); free_ivector(modelidx, 0,3);
  free_ivector(firstingroup,0,ngroups);
  free_ivector(sel,0,nbvars); free_ivector(selnew,0,nbvars); free_ivector(sample, 0,ngroups-1);
  if (nonbinary) { free_ivector(selnew2, 0,nbvars); free_ivector(selnew3, 0,nbvars); }
  free_dvector(mfamily,0,nbfamilies-1); free_dvector(pfamily,0,nbfamilies-1);
  delete integrals;
}





//greedyVarSelC: greedy search for posterior mode in variable selection.
//               Similar to Gibbs sampling, except that deterministic updates are made iff there is an increase in post model prob
//               The scheme proceeds until no variable is included/excluded or niter iterations are reached
// Input arguments: same as in modelSelectionC.

// [[Rcpp::export]]
SEXP greedyVarSelCI(SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP ShasXtX, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP SpriorConstr, SEXP SprConstrp, SEXP SparprConstrp, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Sverbose) {

  bool hasXtX= LOGICAL(ShasXtX)[0];
  int i, j, idxj, logscale=1, mycols, *postMode, *nconstraints, *ninvconstraints, nuncens, ngroupsconstr=0, *isgroup, usethinit=INTEGER(Susethinit)[0], priorcode;
  double offset=0, *postModeProb, *ytXuncens=NULL, *thinit;
  intptrvec constraints, invconstraints;
  crossprodmat *XtX, *XtXuncens=NULL;
  struct marginalPars pars;
  SEXP ans;

  mycols= INTEGER(Sp)[0];

  thinit= dvector(0,mycols+1);
  if (usethinit != 3) {
    for (j=0; j<= mycols+1; j++) { thinit[j]= 0; }
  } else {
    for (j=0; j<= INTEGER(Sp)[0]; j++) { thinit[j]= REAL(Sthinit)[j]; }
  }


  PROTECT(ans= Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ans, 0, Rf_allocVector(INTSXP, mycols));
  postMode= INTEGER(VECTOR_ELT(ans,0));
  for (j=0; j<mycols; j++) postMode[j]= 0;

  SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, 1));
  postModeProb= REAL(VECTOR_ELT(ans,1));

  isgroup= ivector(0, INTEGER(Sp)[0]);
  nconstraints= ivector(0,INTEGER(Sngroups)[0]); ninvconstraints= ivector(0,INTEGER(Sngroups)[0]);
  countConstraints(nconstraints, &constraints, ninvconstraints, &invconstraints, &ngroupsconstr, isgroup, INTEGER(Sngroups), INTEGER(Snvaringroup), Sconstraints, Sinvconstraints);

  if (hasXtX) {
    XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  } else {
    XtX= new crossprodmat(REAL(Sx),INTEGER(Sn)[0],INTEGER(Sp)[0],false);
  }

  if (LENGTH(Suncens)>0) { //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
    int n=INTEGER(Sn)[0], *uncens= INTEGER(Suncens);
    double *pty= REAL(Sy), *ptx= REAL(Sx);
    for (nuncens=0; (nuncens<n) && (uncens[nuncens]==1); nuncens++) { }
    XtXuncens= new crossprodmat(REAL(Sx), INTEGER(Sn)[0], INTEGER(Sp)[0], false, nuncens, 0);
    ytXuncens= dvector(0,INTEGER(Sp)[0]);
    for (j=0; j< INTEGER(Sp)[0]; j++) { for (i=0, ytXuncens[j]=0, idxj=j*n; i< nuncens; i++) { ytXuncens[j] += pty[i] * ptx[i + idxj]; } }
  } else { nuncens= INTEGER(Sn)[0]; }

  set_marginalPars(&pars, INTEGER(Sfamily), INTEGER(Sn), &nuncens, INTEGER(Sp), REAL(Sy), INTEGER(Suncens), REAL(Ssumy2), REAL(Ssumy), REAL(Ssumlogyfact), REAL(Sx), REAL(Scolsumsx), XtX, REAL(SytX), INTEGER(Smethod), INTEGER(Sadjoverdisp), INTEGER(Shesstype), INTEGER(SoptimMethod), INTEGER(Soptim_maxit), &usethinit, thinit, INTEGER(SB), REAL(Salpha), REAL(Slambda), INTEGER(Sknownphi), REAL(Sphi), REAL(Stau), REAL(Staugroup), REAL(Staualpha), REAL(Sfixatanhalpha), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), REAL(SprConstrp), REAL(SparprConstrp), &logscale, &offset, INTEGER(Sgroups), isgroup, INTEGER(Sngroups), &ngroupsconstr, INTEGER(Snvaringroup), nconstraints, ninvconstraints, XtXuncens, ytXuncens);

  priorcode = mspriorCode(INTEGER(SpriorCoef), INTEGER(SpriorGroup), &pars);
  pars.priorcode= &priorcode;

  greedyVarSelC(postMode,postModeProb,INTEGER(SpriorDelta),INTEGER(SpriorConstr),INTEGER(Sniter),INTEGER(Sndeltaini),INTEGER(Sdeltaini),INTEGER(Sincludevars),&constraints,&invconstraints,INTEGER(Sverbose),&pars);

  free_dvector(thinit, 0,mycols+1); free_ivector(isgroup, 0, INTEGER(Sp)[0]);
  free_ivector(nconstraints, 0,INTEGER(Sngroups)[0]); free_ivector(ninvconstraints, 0,INTEGER(Sngroups)[0]);
  delete XtX;
  UNPROTECT(1);
  return ans;

}

void greedyVarSelC(int *postMode, double *postModeProb, int *prDelta, int *prConstr, int *niter, int *ndeltaini, int *deltaini, int *includevars, intptrvec *constraints, intptrvec *invconstraints, int *verbose, struct marginalPars *pars) {

  int i, j, k, l, jgroup, *sel, *selnew, *selaux, nsel, nselnew, nchanges, ngroups, *nvaringroup, *firstingroup, *nconstraints= (*pars).nconstraints, *family= (*pars).family;
  double newJ;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(pars);
  priorFunction= set_priorFunction(prDelta, prConstr, family);
  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  ngroups= *((*pars).ngroups);
  nvaringroup= (*pars).nvaringroup;
  firstingroup= ivector(0,ngroups);
  for (j=1, firstingroup[0]=0; j<ngroups; j++) { firstingroup[j]= firstingroup[j-1] + nvaringroup[j-1]; }

  //Initialize
  if (*verbose==1) Rprintf(" Greedy searching posterior mode... ");
  for (j=0, nsel=*ndeltaini; j< nsel; j++) { sel[j]= deltaini[j]; }
  *postModeProb= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);

  //Iterate
  for (i=0, nchanges=1; (i< *niter) && (nchanges>0); i++) {
    j= jgroup= 0; nchanges= 0;

    while(j< *(*pars).p) {
      sel2selnew(jgroup,sel,&nsel,selnew,&nselnew,false,&ngroups,nvaringroup,firstingroup); //copy sel into selnew, adding/removing jth group
      //if (nsel > nselnew) { naddgroups= 0; ndropgroups=1; dropgroups= jgroup; } else { naddgroups=1 ; ndropgroups=0; addgroups= jgroup; }
      //validmodel= checkConstraints(&addgroups,&naddgroups,&dropgroups,&ndropgroups,constraints,nconstraints,invconstraints,ninvconstraints,(*pars).groups,nvaringroup,sel,&nsel);

      if (includevars[j]==0) {
      //if (includevars[j]==0 && validmodel) {
        newJ= marginalFunction(selnew,&nselnew,pars) + priorFunction(selnew,&nselnew,pars);
        if (newJ > *postModeProb) {
          *postModeProb= newJ;  //update post mode prob
          selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; //update model indicator
          nchanges++;
        }
      }
      j += nvaringroup[jgroup];
      jgroup++;
    } //end j for
  }
  for (j=0; j<nsel; j++) { postMode[sel[j]]= 1; }

  //Enforce hierarchical constraints
  nchanges= 1; nselnew= nsel;
  while (nchanges > 0) {
    nchanges= 0;
    for (j=0; j < ngroups; j++) {
      if ((postMode[firstingroup[j]]==1) && (nconstraints[j]>0)) {

        for (i=0; i < nconstraints[j]; i++) {
          k= (*constraints)[j][i]; //include all variables in group k
          if (postMode[firstingroup[k]]==0) {
            nchanges++;
            for (l=0; l < nvaringroup[k]; l++) { postMode[firstingroup[k]+l]= 1; nselnew += nvaringroup[k]; }
          }
        }

      }
    }
  }

  if (nselnew > nsel) {
    for (j=0, i=0; (j < *(*pars).p) && (i < nselnew); j++) { if (postMode[j]==1) { selnew[i]=j; i++; } } //copy postMode into selnew
    (*postModeProb)= marginalFunction(selnew,&nselnew,pars) + priorFunction(selnew,&nselnew,pars);
  }


  if (*verbose==1) Rprintf("Done.\n");

  free_ivector(firstingroup,0,ngroups);free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}



//Update posterior mode by copying selnew into postMode
void update_postMode(int *postMode, int nselnew, int *selnew, int p, int family) {
  int k;
  for (k=0; k< p; k++) { postMode[k]= 0; }
  for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; }
  if (family==0) {
    if (selnew[nselnew]== p) { //Normal residuals
      postMode[p]= 0;
      postMode[p +1]= 0;
    } else if (selnew[nselnew]== p+1) { //Asymmetric Normal residuals
      postMode[p]= 1;
      postMode[p +1]= 0;
    } else if (selnew[nselnew]== (p) +2) { //Laplace residuals
      postMode[p]= 0;
      postMode[p +1]= 1;
    } else { //Asymmetric Laplace residuals
      postMode[p]= 1;
      postMode[p +1]= 1;
    }
  }
}


bool checkConstraints(int *addgroups, int *naddgroups, int *dropgroups, int *ndropgroups, intptrvec *constraints, int *nconstraints, intptrvec *invconstraints, int *ninvconstraints, int *groups, int *nvaringroup, int *sel, int *nsel) {
  /* Check if adding variables in addgroups and dropping variables from dropgroups to current model (sel) gives a valid model satisfying all constraints
   Input

   - addgroups: id of groups to be added to sel (vector of length naddgroups)
   - naddgroups:  number of groups to be added
   - dropgroups: id of groups to be dropped from sel (vector of length ndropgroups)
   - ndropgroups: number of groups to be dropped
   - constraints: constraints[j] is a vector indicating all groups required by group j. Assumed to be ordered increasingly
   - nconstraints: nconstraints[j] is the length of constraints[j] (number of constraints required by group j)
   - invconstraints: invconstraints[j] is a vector indicating all inverse constraints of group j (ids of other groups requiring group j)
   - ninvconstraints: ninvconstraints[j] is the length of invconstraints[j] (number of inverse constraints of group j)
   - firstingroup: index of the first variable in all groups
   - sel: index of variables currently in the model. Assumed to be ordered increasingly
   - nsel: number of variables currently in the model
  */

  bool valid= true;
  int j, k, l, *curconstraints, curgroup, nvalid;

  //For any group we want to add to sel, check that its constraints are also in sel
  for (k=0; (k< *naddgroups) && valid; k++) {
    nvalid= j= l= 0;
    curgroup= addgroups[k];
    curconstraints= (*constraints)[curgroup];
    while ((l < nconstraints[curgroup]) && (j < *nsel) && valid) {
      if (groups[sel[j]] > curconstraints[l]) {
        valid= false;
      } else if (groups[sel[j]] == curconstraints[l]) {  //add groups and nvaringroup as parameters
        l++; nvalid++;
      } else {
        j += nvaringroup[sel[j]];
      }
    }
    if (nvalid < nconstraints[curgroup]) { valid= false; }
  }

  //For any group we want to drop from sel, check that its inverse constraints are not in sel
  for (k=0; (k< *ndropgroups) && valid; k++) {
    j= l= 0;
    curgroup= dropgroups[k];
    curconstraints= (*invconstraints)[curgroup];
    while ((l < ninvconstraints[curgroup]) && (j < *nsel) && valid) {
      if (groups[sel[j]] > curconstraints[l]) {
        l++;
      } else if (groups[sel[j]] == curconstraints[l]) {
        valid= false;
      } else {
        j+= nvaringroup[groups[sel[j]]];
      }
    }
  }

  return valid;
}


/* OLD VERSION, IT DID NOT WORK CORRECTLY
bool checkConstraints(int *constraints, int *nconstraints, int *firstingroup, int *sel, int *nsel) {
  // Check if current model (sel) contains all the variable groups indicated in constraints
  // Input
  // - constraints: list of groups that should be included in sel, e.g. group 1,3. Assumed to be ordered increasingly
  // - nconstraints: number of groups (elements in constraints)
  // - firstingroup: index of the first variable in all groups
  // - sel: index of variables currently in the model. Assumed to be ordered increasingly
  // - nsel: number of variables currently in the model
  bool valid= true;
  int j, jgroup, firstvar;
  jgroup= 0;
  firstvar= firstingroup[constraints[jgroup]];
  for (j=0; (j<(*nsel)) && valid && (jgroup< *nconstraints); j++) {
    if (sel[j]>firstvar) {
      valid= false;
    } else if (sel[j]==firstvar) {
      jgroup++;
      if (jgroup< *nconstraints) { firstvar= firstingroup[constraints[jgroup]]; }
    }
  }
  return valid;
}
*/


void sel2selnew(int newgroup, int *sel, int *nsel, int *selnew, int *nselnew, bool copylast, int *ngroups, int *nvaringroup, int *firstingroup) {
//Copy sel into selnew. Elements are always kept ordered so that sel[0] < sel[1] < sel[2] ...
// - If newgroup in sel, don't copy it in selnew and set nselnew=nsel-1.
// - If newgroup not in sel, add it to selnew and set nselnew=nsel+1.
// - If copylast==true, copy last element sel[nsel] into selnew[nselnew]
  bool found;
  int i, ii, iii;
  for (i=0, found=false; (i< *nsel) && (sel[i]<=firstingroup[newgroup]) && (!found); i++) {
    selnew[i]= sel[i]; found= (sel[i]==firstingroup[newgroup]);
  }
  if (!found) { //add new group
    for (ii=0; ii<nvaringroup[newgroup]; ii++) { selnew[i+ii]= firstingroup[newgroup]+ii; }
    for (iii=0; (i+iii)<*nsel ; iii++) { selnew[i+ii+iii]= sel[i+iii]; }
    (*nselnew)= (*nsel) + nvaringroup[newgroup];
  } else {  //remove new elem
    for (ii=i-1+nvaringroup[newgroup]; ii< *nsel; ii++) { selnew[ii-nvaringroup[newgroup]]= sel[ii]; }
    (*nselnew)= (*nsel) - nvaringroup[newgroup];
  }
  if (copylast) selnew[*nselnew]= sel[*nsel];
}

/* OLD VERSION. ELEMENTS IN sel WERE NOT ORDERED, ALSO VARIABLE GROUPS WERE NOT SUPPORTED
void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew, bool copylast) {
//Copy sel into selnew.
// - If j in sel, don't copy it in selnew and set nselnew=nsel-1.
// - If j not in sel, add it to selnew and set nselnew=nsel+1.
// - If copylast==true, copy last element sel[nsel] into selnew[nselnew]
  int i, ii, found;
  for (i=0, found=0; (i< *nsel) && (found==0); i++) { selnew[i]= sel[i]; found= (newelem==sel[i]); }
  if (found==0) { //add newelem
    selnew[*nsel]= newelem; (*nselnew)= (*nsel)+1;
  } else {  //remove new elem
    for (ii=i; ii< *nsel; ii++) { selnew[ii-1]= sel[ii]; }
    (*nselnew)= (*nsel)-1;
  }
  if (copylast) selnew[*nselnew]= sel[*nsel];
  }*/


void findselgroups(double *nvarinselgroups, double *firstingroup, double *nselgroups, double *selgroups, int *sel, int *nsel, int *nvaringroup, int *ngroups) {
  //Return sub-vector of nvaringroup corresponding to selected groups (as indicated by sel), and the total number of selected groups
  // Input
  // - sel: vector indicating selected variables. IMPORTANT: these must be given in increasing order, i.e. sel[0] < sel[1] < sel[2] ...
  // - nsel: number of selected variables (length of sel)
  // - nvaringroup: number of variables in each group
  // - ngroup: number of groups (length of nvaringroup)
  // Output
  // - nvarinselgroups: subset of elements in nvaringroup containing only groups with variables in sel
  // - firstingroup: index of first variable in each group, e.g if nvarinselgroups=(2,1,1) then firstingroup=(0,2,3)
  // - nselgroups: number of selected groups (length of nvarinselgroups)
  // - selgroups: group id for each variable in sel (starting from 0,1,...)
  //
  // Example:  Suppose that nvaringroup= (1,1,2), indicating that there's 3 groups with (1,1,2) variables in each.
  //           Hence variable 0 is in group 1, variable 1 in group 2, variables 2-3 in group 3
  //
  //           Further suppose that sel[0]=0, sel[1]=2, sel[2]=3, indicating that variables (0,2,3) are selected
  //
  //           Then we selected 2 groups with (1,2) variables respectively, hence nselgroups=2 and nvarinselgroups= (1,2)
  int j, l, groupid, nselgroupsi, firstidx;
  (*nselgroups) = 0;
  groupid= firstidx= j= 0;
  while (j< *nsel) {
    while ((groupid < *ngroups) && (firstidx < sel[j])) {
      firstidx += nvaringroup[groupid];
      groupid++;
    }
    nselgroupsi= (int) (*nselgroups +.1);
    nvarinselgroups[nselgroupsi]= (double) nvaringroup[groupid];
    if (nselgroupsi==0) { firstingroup[0]= 0; } else  { firstingroup[nselgroupsi]= firstingroup[nselgroupsi-1] + nvarinselgroups[nselgroupsi-1]; }
    for (l=0; l< nvaringroup[groupid]; l++) selgroups[j+l]= (double) (*nselgroups);
    j += nvaringroup[groupid];
    (*nselgroups)++;
  }
}


void nselConstraints(int *ngroups0, int *ngroups1, int *sel, int *nsel, int *group, int *nconstraints, int *nvaringroup) {
  /* Return number of selected groups that have hierarchical constraints

     Input
     - sel: indexes of selected variables
     - nsel: number of selected variables (length of sel)
     - group: group[sel[j]] indicates the group that variable j belongs to. Variables are assumed to be ordered by groups (group 1, group 2 etc)
     - nconstraints: nconstraints[l] is the number of constraints for group l (nconstraints[l]==0 indicates no constraints)

     Output:
     - ngroups0: number of groups in sel that do not have hierarchical constraints
     - ngroups1: number of groups in sel that have hierarchical constraints
  */
  int j=0, g;
  (*ngroups0)= (*ngroups1)= 0;
  while (j< *nsel) {
    g= group[sel[j]];
    if (nconstraints[g]==0) { (*ngroups0)++; } else { (*ngroups1)++; }
    j += nvaringroup[g];
  }
}


void countConstraints(int *nconstraints, intptrvec *constraints, int *ninvconstraints, intptrvec *invconstraints, int *ngroupsconstr, int *isgroup, int *ngroups, int *nvaringroup, SEXP Sconstraints, SEXP Sinvconstraints) {
  /*Count number of constraints, number of groups with constraints, determine which variables are in a group

  Input:
  - ngroups: number of groups
  - nvaringroup: nvaringroup[j] is the number of variables in group j=0,1...,ngroups-1
  - Scontraints: list with length ngroups indicating the constraints for each group
  - Sinvcontraints: list with length ngroups indicating the inverse constraints for each group (i.e. if groups (2,3) require group 1, then the inverse constraints for group 1 are (2,3))

  Output:
  - nconstraints: nconstraints[j] is the number of constraints in group j (length of Sconstraints[[j]])
  - constraints: element j is a vector of length nconstraints[j] indicating the constraints in Sconstraints[[j]]
  - ninvconstraints: ninvconstraints[j] is the number of inverse constraints in group j (length of Sinvconstraints[[j]])
  - invconstraints: element j is a vector of length ninvconstraints[j] indicating the inverse constraints in Sinvconstraints[[j]]
  - ngroupsconstr: number of groups with constraints (that is, with nconstraints[j]>0)
  - isgroup: isgroup[i]==1 indicates that variable i=0,...,sum(nvaringroup)-1 belongs to a group (isgroup[i]==0 that it doesn't belong to a group)
  */
  int i, j, jj;
  for (j=0,jj=0; j< *ngroups; j++) {
    nconstraints[j]= LENGTH(VECTOR_ELT(Sconstraints, j));
    ninvconstraints[j]= LENGTH(VECTOR_ELT(Sinvconstraints, j));
    constraints->push_back(INTEGER(VECTOR_ELT(Sconstraints, j)));
    invconstraints->push_back(INTEGER(VECTOR_ELT(Sinvconstraints, j)));
    if (nconstraints[j]>0) (*ngroupsconstr)++;
    isgroup[jj]= ((int) (nvaringroup[j] +.1)) >1;
    jj++;
    for (i=1; i< nvaringroup[j]; i++, jj++) isgroup[jj]= isgroup[jj-1];
  }
}




/* ALA orthogonal approximation to posterior expectation of pMOM penalty on individual coefficients (grouped coefficients are not penalized, for that see gmompenalty_approx below)

   The MOM prior across individual coefficients g is theta ~ prod_g theta_g^2 / tau N(theta; 0, tau)

   The ALA approximation is    prod_g (thopt_g^2 + 1/Hinv_g) / tau, where thopt is the posterior mode the underlying normal prior, and Hinv the inverse of the posterior hessian at thopt

   INPUT
   - thopt: posterior mode under the underlying Normal prior
   - Hinv: inverse of posterior hessian at thopt
   - thlength: length of thopt
   - nvaringroups: number of variables in each selected group, as returned by function cholSini_indexes
   - firstingroup: index of first variable in each group, as returned by function cholSini_indexes
*/

double pmompenalty_approx(double *thopt, double **Hinv, double *tau, int thlength, double *nvaringroups, double *firstingroup) {
  int i, idx;
  double ans= 0, logtau= log(*tau);
  for (i=0; i< thlength; i++) {
    if (nvaringroups[i] < 1.1) {
      idx= (int) (firstingroup[i]+.1);
      ans += log(thopt[idx]*thopt[idx] + Hinv[idx+1][idx+1]) - logtau;
    }
  }
  return ans;
}

/* ALA orthogonal approximation to posterior expectation of group MOM penalty. 

   The MOM prior across individual/grouped coefficients g is theta ~ prod_g (theta_g' Sinv theta_g) / phi N(theta; 0, phi Sinv^{-1})

   The ALA approximation is    prod_g trace(Sinv_g %*% Hinv_g) + thopt' Sinv thopt / phi, where thopt is the posterior mode the underlying normal prior, and Hinv the inverse of the posterior hessian at thopt

   INPUT
   - momsingle: if true, the MOM penalty is computed across individual coefficients
   - momgroup: if true, the MOM penalty is computed across grouped coefficients
   - thopt: posterior mode under the underlying Normal prior
   - Hinv: inverse of posterior hessian at thopt
   - phi: prior dispersion parameter, assumed known
   - thlength: length of thopt
   - nsel: number of parameters in theta
   - nselgroups: number of groups in theta
   - nvaringroups: number of variables in each selected group, as returned by function cholSini_indexes
   - firstingroup: index of first variable in each group, as returned by function cholSini_indexes
   - cholSini: Cholesky decomposition of S, as returned by function cholSini_indexes
*/

double gmompenalty_approx(bool momsingle, bool momgroup, double *thopt, double **Hinv, double *Sinv, double phi, int thlength, int nsel, int nselgroups, double *nvaringroups, double *firstingroup, double *cholSini) {

  int i, j, groupid, ningroup, idx;
  double ans=0, trgroup, th2group, Sij;

    //For each group, compute trace(Sinv %*% Hinv) + thopt' Sinv thopt / phi
    for (groupid=0; groupid< nselgroups; groupid++) {

      ningroup= (int) (nvaringroups[groupid] +.1);

      if (((ningroup==1) && momsingle) || ((ningroup>1) && momgroup)) {  //if MOM prior was set on the current group
        idx= (int) (firstingroup[groupid]+.1);
         
        trgroup= th2group= 0;
        for (i=1; i<=ningroup; i++) {
          Sij= getelem_Sinv(groupid, i, i, Sinv, cholSini, ningroup);
          trgroup += Sij * Hinv[idx+i][idx+i];
          th2group += thopt[idx + i -1] * thopt[idx + i -1] * Sij;
          for (j=i+1; j<=ningroup; j++) {
            Sij= getelem_Sinv(groupid, i, j, Sinv, cholSini, ningroup);
            trgroup += 2.0 * Sij * Hinv[idx+j][idx+i];  //Sinv[i][j] * Hinv[j][i]
            th2group += 2.0 * thopt[idx + i -1] * thopt[idx + j -1] * Sij;
          }
        }
         
        th2group= th2group / phi;
         
        ans += log( (trgroup + th2group)/(ningroup +.0) );
      }

    } //end for each group

  return ans;

}




//********************************************************************************************
// PRIORS ON MODEL SPACE (always return on log scale)
//********************************************************************************************

double unifPrior(int *sel, int *nsel, struct marginalPars *pars) { return -((*(*pars).ngroups) +.0) * log(2.0); }
double unifPriorTP(int *sel, int *nsel, struct marginalPars *pars) { return -((*(*pars).ngroups) +2.0) * log(2.0); }
double unifPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) { return 0.0; }

double vectBinom(int *sel, int *nsel, int len_prDeltap, int len_prConstrp, struct marginalPars *pars) {
  int i, sel_i=0, delta_i=0, constr_i=0, ngroups=*(*pars).ngroups, *groups=(*pars).groups, *nconstraints=(*pars).nconstraints, *nvaringroup=(*pars).nvaringroup;
  double ans=0, *prDeltap=(*pars).prDeltap, *prConstrp=(*pars).prConstrp;

  if (*nsel == 0) {
    for (i=0; i<len_prDeltap; i++) ans += log(1-prDeltap[(len_prDeltap > 1) ? i : 0]);
    if (*(*pars).ngroupsconstr > 0) {
      for (i=0; i<len_prConstrp; i++) ans += log(1-prConstrp[(len_prConstrp > 1) ? i : 0]);
    }
  } else {
    for (i = 0; i<ngroups; i++) {
      if (nconstraints[i] == 0) {
        if (i == groups[sel[sel_i]]) {
          ans += log(prDeltap[delta_i]);
          if (sel_i < *nsel-1) sel_i += nvaringroup[groups[i]];
        } else {
          ans += log(1 - prDeltap[delta_i]);
        }
        if (len_prDeltap > 1) delta_i++;
      } else { // constrained, use prConstrp
        if (i == groups[sel[sel_i]]) {
          ans += log(prConstrp[constr_i]);
          if (sel_i < *nsel-1) sel_i += nvaringroup[groups[i]];
        } else {
          ans += log(1 - prConstrp[constr_i]);
        }
        if (len_prConstrp > 1) constr_i++;
      }
    }
  }
  return ans;
}

//nsel ~ Binom(p,prDeltap)
double binomPrior(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1, n_notconstr, n_constr=*(*pars).ngroupsconstr,
      len_prDeltap=(int) *(*pars).parprDeltap, len_prConstrp=(int) *(*pars).parprConstrp;
  double ans, *prDeltap=(*pars).prDeltap, *prConstrp=(*pars).prConstrp;
  nselConstraints(&ngroups0, &ngroups1, sel, nsel, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  n_notconstr= *(*pars).ngroups - n_constr;
  if ((len_prDeltap == 1) & (len_prConstrp == 1)) {
    ans= (ngroups0+.0)*log(*prDeltap) + (n_notconstr-ngroups0+.0)*log(1-*prDeltap);
    if ((n_constr) > 0) {
      ans += (ngroups1+.0)*log(*prConstrp) + (n_constr-ngroups1+.0)*log(1-*prConstrp);
    }
  } else {
    ans = vectBinom(sel, nsel, len_prDeltap, len_prConstrp, pars);
  }
  return ans;
  //return dbinomial(*nsel,*(*pars).p,*(*pars).prDeltap,1);
}

double binomPriorTP(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1, nselminus= *nsel -1;
  double ans;
  nselConstraints(&ngroups0, &ngroups1, sel, &nselminus, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  ans = dbinomial(ngroups0,*(*pars).ngroups - *(*pars).ngroupsconstr,*(*pars).prDeltap,1);
  if ((*(*pars).ngroupsconstr) >0) ans += dbinomial(ngroups1,*(*pars).ngroupsconstr,*(*pars).prDeltap,1);
  return ans - 2.0*log(2.0);
  //return dbinomial(*nsel -1,*(*pars).p,*(*pars).prDeltap,1) - 2.0*log(2.0);
}

double binomPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return dbinomial(*nsel,*(*pars).p1,(*pars).prModelpar[0],1);
}

//nsel ~ Beta-Binomial(prModelpar[0],prModelPar[1])
double betabinPrior(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1;
  double ans;
  nselConstraints(&ngroups0, &ngroups1, sel, nsel, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  ans = bbPrior(ngroups0,*(*pars).ngroups - *(*pars).ngroupsconstr, (*pars).parprDeltap[0], (*pars).parprDeltap[1],1);
  if ((*(*pars).ngroupsconstr) >0) ans += bbPrior(ngroups1,*(*pars).ngroupsconstr,(*pars).parprConstrp[0], (*pars).parprConstrp[1],1);
  return ans;
  //return bbPrior(*nsel, *(*pars).p, (*pars).parprDeltap[0], (*pars).parprDeltap[1],1);
}

double betabinPriorTP(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1, nselminus= *nsel -1;
  double ans;
  nselConstraints(&ngroups0, &ngroups1, sel, &nselminus, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  ans = bbPrior(ngroups0,*(*pars).ngroups - *(*pars).ngroupsconstr, (*pars).parprDeltap[0], (*pars).parprDeltap[1],1);
  if ((*(*pars).ngroupsconstr) >0) ans += bbPrior(ngroups1,*(*pars).ngroupsconstr,(*pars).parprConstrp[0], (*pars).parprConstrp[1],1);
  return ans - 2.0*log(2.0);
  //return bbPrior(*nsel -1, *(*pars).p, (*pars).parprDeltap[0], (*pars).parprDeltap[1], 1) - 2.0*log(2.0);
}

double betabinPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return bbPrior(*nsel,*(*pars).p1,(*pars).prModelpar[0],(*pars).prModelpar[1],1);
}

double complexityPrior(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1;
  double ans;
  nselConstraints(&ngroups0, &ngroups1, sel, nsel, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  ans = complexPrior(ngroups0,*(*pars).ngroups - *(*pars).ngroupsconstr, (*pars).prDeltap[0],1);
  if ((*(*pars).ngroupsconstr) >0) ans += complexPrior(ngroups1,*(*pars).ngroupsconstr,(*pars).prDeltap[0],1);
  return ans;
  //return complexPrior(*nsel,*(*pars).p,(*pars).prDeltap[0],1);
}

double complexityPriorTP(int *sel, int *nsel, struct marginalPars *pars) {
  int ngroups0, ngroups1, nselminus= *nsel -1;
  double ans;
  nselConstraints(&ngroups0, &ngroups1, sel, &nselminus, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
  ans = complexPrior(ngroups0,*(*pars).ngroups - *(*pars).ngroupsconstr, (*pars).prDeltap[0],1);
  if ((*(*pars).ngroupsconstr) >0) ans += complexPrior(ngroups1,*(*pars).ngroupsconstr,(*pars).prDeltap[0],1);
  return ans - 2.0*log(2.0);
  //return complexPrior(*nsel -1,*(*pars).p,(*pars).prDeltap[0],1) - 2.0*log(2.0);
}

double complexityPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return complexPrior(*nsel,*(*pars).p1,(*pars).prModelpar[0],1);
}



//*************************************************************************************
// PRIOR DENSITIES ON PARAMETERS
//*************************************************************************************


void dmomgzell(double *ans, double *th, double *tau, double *nvaringroup, double *ngroups, double *ldetSinv, double *cholSinv, double *cholSini, bool logscale) {
  /*Evaluate pMOM(tau) + block Zellner(taugroup) prior density at th

     prod_j pMOM(beta_j; tau)  prod_j N(delta_j; 0, S_j^{-1})

     where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
     corresponding to coefficients for individual variables and grouped variables (respectively)

    Input
     - th: point at which to evaluate the density is th[0],...,th[*nsel]
     - tau: prior pMOM dispersion parameter (drives prior on individual coefficients, i.e. groups of size 1)
     - nvaringroup: number of variables in each group
     - ngroups: number of selected groups (length of nvaringroup)
     - ldetSinv: log-determinant of S_j^{-1} for each group j=1,...,ngroups
     - cholSinv: Cholesky decomp of S_j^{-1} for all groups stored as a single vector
     - cholSini: Chol decomp of S_j^{-1} for group j starts at cholSinv[cholSini[j]]
     - logscale: if true return log-density, else return the density

    Output
     - ans: log-density (if logscale==true) or density (if logscale==false)
  */
  int i, firstingroup=0, ningroup, ngroupsi= (int) (*ngroups +.1);

  (*ans)= 0;
  for (i=0; i< ngroupsi; i++) {
    ningroup= (int) (nvaringroup[i] +.1);
    if (ningroup== 1) {
      (*ans) += dmom(th[firstingroup], 0, *tau, 1, 1, true);
    } else {
      (*ans) += dmvnorm0(th+firstingroup-1, ningroup, cholSinv + (int) (cholSini[i]+.1), ldetSinv[i], true, true);
    }
    firstingroup += ningroup;
  }
  if (!logscale) (*ans)= exp(*ans);
}


void demomgzell(double *ans, double *th, double *tau, double *nvaringroup, double *ngroups, double *ldetSinv, double *cholSinv, double *cholSini, bool logscale) {
  /*Evaluate peMOM(tau) + block Zellner(taugroup) prior density at th

     prod_j peMOM(beta_j; tau)  prod_j N(delta_j; 0, S_j^{-1})

     where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
     corresponding to coefficients for individual variables and grouped variables (respectively)

    Input
     - th: point at which to evaluate the density is th[0],...,th[*nsel]
     - tau: prior pMOM dispersion parameter (drives prior on individual coefficients, i.e. groups of size 1)
     - nvaringroup: number of variables in each group
     - ngroups: number of selected groups (length of nvaringroup)
     - ldetSinv: log-determinant of S_j^{-1} for each group j=1,...,ngroups
     - cholSinv: Cholesky decomp of S_j^{-1} for all groups stored as a single vector
     - cholSini: Chol decomp of S_j^{-1} for group j starts at cholSinv[cholSini[j]]
     - logscale: if true return log-density, else return the density

    Output
     - ans: log-density (if logscale==true) or density (if logscale==false)
  */
  int i, firstingroup=0, ningroup, ngroupsi= (int) (*ngroups +.1);

  (*ans)= 0;
  for (i=0; i< ngroupsi; i++) {
    ningroup= (int) (nvaringroup[i] +.1);
    if (ningroup== 1) {
      (*ans) += demom(th[firstingroup], *tau, 1, true);
    } else {
      (*ans) += dmvnorm0(th+firstingroup-1, ningroup, cholSinv + (int) (cholSini[i]+.1), ldetSinv[i], true, true);
    }
    firstingroup += ningroup;
  }
  if (!logscale) (*ans)= exp(*ans);
}


void dgzellgzell(double *ans, double *th, double *nvaringroup, double *ngroups, double *ldetSinv, double *cholSinv, double *cholSini, bool logscale) {
  /*Evaluate Zellner + block Zellner prior density at th

     prod_j N(beta_j; 0; S_j^{-1})  prod_j N(delta_j; 0, S_j^{-1})

     where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
     corresponding to coefficients for individual variables and grouped variables (respectively)

    Input
     - th: point at which to evaluate the density is th[0],...,th[*nsel]
     - nvaringroup: number of variables in each group
     - ngroups: number of selected groups (length of nvaringroup)
     - ldetSinv: log-determinant of S_j^{-1} for each group j=1,...,ngroups
     - cholSinv: Cholesky decomp of S_j^{-1} for all groups stored as a single vector
     - cholSini: Chol decomp of S_j^{-1} for group j starts at cholSinv[cholSini[j]]
     - logscale: if true return log-density, else return the density

    Output
     - ans: log-density (if logscale==true) or density (if logscale==false)
  */
  int i, firstingroup=0, ningroup, ngroupsi= (int) (*ngroups +.1);

  (*ans)= 0;
  for (i=0; i< ngroupsi; i++) {
    ningroup= (int) (nvaringroup[i] +.1);
    (*ans) += dmvnorm0(th+firstingroup-1, ningroup, cholSinv + (int) (cholSini[i]+.1), ldetSinv[i], true, true);
    firstingroup += ningroup;
  }
  if (!logscale) (*ans)= exp(*ans);
}



/*Compute inverse covariance matrix in group Zellner's prior, its Cholesky decomposition and determinant
  INPUT
  - ngroups: total number of groups
  - nvaringroups: nvaringroups[g] is the number of variables in group g=0,...,ngroups-1
  - sel: sel[j] is the index of the jth selected variable
  - cholSini: cholSini[g] is the index at which the Sinv and cholSinv for group g should start
  - XtX: XtX[sel[j],sel[l]] is the entry in the Gram matrix corresponding to variables (sel[j],sel[l])
  - tau: prior dispersion parameter for groups with 1 variable
  - taugroup: prior dispersion parameter for groups with >1 variables

  OUTPUT. Sinv, cholSinv and ldetSinv store the group Zellner's prior precision matrix, its Cholesky decomposition and log-determinants.

  If A is the submatrix of XtX corresponding to group j, then the precision matrix for group g is

  (1 / tau) A, if nvaringroups ==1 and orthoapprox==false

  (1 / tau), if nvaringroups ==1 and orthoapprox==true

  (nvaringroups / taugroup) A, if nvaringroups >1

  and it is stored in Sinv[cholSini[g]...cholSini[g+1]-1], its Cholesky decomp in cholSinv[cholSini[g]...cholSini[g+1]-1], its log-determinant is stored in ldetSinv[g]
*/
void gzell_Sinv(double *Sinv, double *cholSinv, double *ldetSinv, int *ngroups, double *nvaringroups, int *sel, double *cholSini, crossprodmat *XtX, double *tau, double *taugroup, bool orthoapprox) {
  bool posdef;
  int i, j, l, k, firstingroup, groupsize, idxini, Sidxini;
  double ct, sqrtct, tauinv=1, sqrttauinv=1, neglogtau= 0;
  if (orthoapprox) {
    tauinv= 1.0/(*tau);
    sqrttauinv= sqrt(tauinv);
    neglogtau= -log(*tau);
  }

  for (i=0, firstingroup=0; i< *ngroups; i++) {
    groupsize= (int) (nvaringroups[i] + .1);
    idxini= sel[firstingroup]; Sidxini= (int) (cholSini[i]+.1);
    XtX->choldc(idxini, idxini+groupsize-1, cholSinv + Sidxini, ldetSinv+i, &posdef);
    if (groupsize==1) {
      ct= 1.0 / (*tau);
    } else {
      ct= nvaringroups[i] / (*taugroup);
    }
    sqrtct= sqrt(ct);
    if ((!orthoapprox) || (groupsize>1)) {
      for (j=0, l=0, k=0; j< groupsize*(groupsize+1)/2; j++) {
        (*(cholSinv + Sidxini + j)) *= sqrtct;
        Sinv[Sidxini + j]= ct * XtX->at(sel[firstingroup+l],sel[firstingroup+k]);
        if (l< groupsize-1) { l++; } else { k++; l=k; }
      }
      ldetSinv[i] = log(ldetSinv[i]) + 2.0 * nvaringroups[i] * log(sqrtct);
    } else {
      cholSinv[Sidxini]= sqrttauinv;
      Sinv[Sidxini]= tauinv;
      ldetSinv[i] = neglogtau;
    }

    firstingroup += groupsize;
  }
}


//Same as gzell_Sinv, adjusting the prior dispersion for pMOM and groupMOM
//
//  priorcode (options currently implemented here)
//     0: pMOM on all coef
//     5: group pMOM (same as pMOM, standardized by n / X'X)
//     9: group Zellner on all coef
//    10: pMOM + group pMOM
//    13: pMOM + group Zellner
//    43: group Zellner + group Zellner
//    50: group pMOM + group pMOM
//    53: group pMOM + group Zellner
//    63: group Zellner + group Zellner
void gzell_Sinv_byprior(double *Sinv, double *cholSinv, double *ldetSinv, int *ngroups, double *nvaringroups, int *sel, double *cholSini, crossprodmat *XtX, int *n, double *tau, double *taugroup, int *priorcode) {
  bool posdef, idcov, groupmom, groupzell;
  int i, j, l, k, firstingroup, groupsize, idxini, Sidxini;
  double ct, sqrtct, ct1=1, tauinv=1, sqrttauinv=1, neglogtau=0;

  idcov= ((*priorcode==0) || (*priorcode == 10) || (*priorcode == 13));  //is the covariance for individual coef tau * 1
  groupmom= ((*priorcode== 10) || (*priorcode==50));                     //was a group MOM prior set on groups
  groupzell= ((*priorcode==9) || (*priorcode==13) || (*priorcode== 43) || (*priorcode== 53) || (*priorcode==63)); //was a group Zellner set on groups

  if (idcov)  {
    tauinv= 1.0/(*tau);
    sqrttauinv= sqrt(tauinv);
    neglogtau= -log(*tau);
  }

  //Pre-compute scaling factor for individual parameters
  if ((*priorcode==9) || (*priorcode==43) || (*priorcode==63)) ct1= 1.0 / (*tau);           //group Zellner on individual coef
  if ((*priorcode==5) || (*priorcode==50) || (*priorcode==53)) ct1= 3.0 / ((*n) * (*tau));  //group pMOM on individual coef

  //Iterate over groups
  for (i=0, firstingroup=0; i< *ngroups; i++) {
    groupsize= (int) (nvaringroups[i] + .1);
    idxini= sel[firstingroup]; Sidxini= (int) (cholSini[i]+.1);

    if ((groupsize==1) && idcov) {

      cholSinv[Sidxini]= sqrttauinv;
      Sinv[Sidxini]= tauinv;
      ldetSinv[i] = neglogtau;

    } else { 

      XtX->choldc(idxini, idxini+groupsize-1, cholSinv + Sidxini, ldetSinv+i, &posdef);

      if (groupsize==1) {
        ct= ct1;
      } else {
        if (groupmom) { //group MOM
          ct= (nvaringroups[i] + 2.0) / ((*n) * (*taugroup));
        } else if (groupzell) { //group Zellner
          ct= nvaringroups[i] / (*taugroup);
        } else {
          Rf_error("prior not implemented in gzell_Sinv\n");
        }
      }
      sqrtct= sqrt(ct);

      for (j=0, l=0, k=0; j< groupsize*(groupsize+1)/2; j++) {
        (*(cholSinv + Sidxini + j)) *= sqrtct;
        Sinv[Sidxini + j]= ct * XtX->at(sel[firstingroup+l],sel[firstingroup+k]);
        if (l< groupsize-1) { l++; } else { k++; l=k; }
      }
      ldetSinv[i] = log(ldetSinv[i]) + 2.0 * nvaringroups[i] * log(sqrtct);

    }

    firstingroup += groupsize;
  }
}

/* Return element (k,l) in Sinv of groupid, as stored by gzell_Sinv_byprior */
double getelem_Sinv(int groupid, int k, int l, double *Sinv, double *cholSini, int ningroup) {
  int k0, l0, kk, idxini;

  idxini= (int) (cholSini[groupid]+.1);

  if (l >= k) { k0= k-1; l0= l-1; } else { k0= l-1; l0= k-1; } //ensure that l>=k

  kk= idxini + k0*ningroup - k0*(k0-1)/2;
  return *(Sinv + kk + l0-k0);  //elem (k,l) in Sinv of groupid, works for l>=k

}

/*Return indexes at which the Cholesky decomposition starts for each group when storing the lower-triangular Cholesky decomp as vectors
  Input
   - ngroups: number of groups
   - nvaringroups: number of variables in each group
  Output
   - cholSini: indexes at which the Cholesky decomp starts
   - cholSsize: total size

  Example: suppose there's three groups with nvaringroups=(2,1,1).
           That is, the 1st group has 2 variables (requiring 3 entries in the Cholesky decomp), the 2nd and 3rd groups have 1 variable.
           Then cholSini[0]=0, cholSini[1]=3, cholSini[2]=4, and cholSsize=5
*/
void cholSini_indexes(double *cholSini, int *cholSsize, int ngroups, double *nvaringroups) {
  int i, groupsize, nelem;
  cholSini[0]= (*cholSsize)= 0;
  for (i=0; i< ngroups; i++) {
    groupsize= (int) (nvaringroups[i] + .1);
    nelem= groupsize * (groupsize+1) / 2;
    (*cholSsize) += nelem;
    if (i< ngroups-1) cholSini[i+1]= cholSini[i] + nelem;
  }
}




//*************************************************************************************
// PRIORS, GRADIENTS AND HESSIANS IN STANDARDIZED FORMAT
//*************************************************************************************


//Evaluate negative log-likelihood + log-prior (pMOM + group Zellner) and initialize funargs
void pmomgzell_log (double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  dmomgzell(f, th, (*pars).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  (*f) = -(*f);

}

//Gradient wrt th[j] of log pMOM + group Zellner prior
void pmomgzell_grad (double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau);
      (*priorgrad)= -2.0/th[j] + th[j]/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    }

}

//Gradient and hessian wrt th[j] of negative log pMOM + group Zellner prior
void pmomgzell_gradhess (double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau);;
      (*priorgrad)= -2.0/th[j] + th[j]/tau;
      (*priorhess)= 2.0/(th[j]*th[j]) + 1.0/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      (*priorhess)= Sinv[jj];
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    }

}


void pmomgzell_hess (double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  /* NEGATIVE HESSIAN FOR pMOM - GROUP ZELLNER PRIOR IS ADDED TO hess

     Input/Output

     - hess: the input value is modified by adding the negative hessian. For example, if on input hess contains the negative hessian of the log-likelihood, then on output hess contains the negative hessian of the log-joint (log-likelihood + log prior)

*/
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    pmomgzell_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] += priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}


//Evaluate negative log-likelihood + log-prior (peMOM + group Zellner) and initialize funargs
void pemomgzell_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  demomgzell(f, th, (*pars).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  (*f)= -(*f);

}

void pemomgzell_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
/* NEGATIVE GRADIENT AND HESSIAN WRT th[j] OF THE peMOM - GROUP ZELLNER PRIOR */

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau), th2= th[j]*th[j];
      (*priorgrad)= -2.0 * tau/(th2 * th[j]) + th[j]/tau;
      (*priorhess)= 6.0*tau/(th2*th2) + 1.0/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      (*priorhess)= Sinv[jj];
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    }

}

void pemomgzell_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
/* NEGATIVE GRADIENT WRT th[j] OF THE peMOM - GROUP ZELLNER PRIOR */

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau), th2= th[j]*th[j];
      (*priorgrad)= -2.0 * tau/(th2 * th[j]) + th[j]/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    }

}

void pemomgzell_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  /* NEGATIVE HESSIAN FOR peMOM - GROUP ZELLNER PRIOR IS ADDED TO hess

     Input/Output

     - hess: the input value is modified by adding the negative hessian. For example, if on input hess contains the negative hessian of the log-likelihood, then on output hess contains the negative hessian of the log-joint (log-likelihood + log prior)

*/
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    pemomgzell_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] += priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}




//Evaluate negative log-likelihood + log-prior (group Zellner + group Zellner) and initialize funargs
void gzellgzell_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  dgzellgzell(f, th, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  (*f)= -(*f);

}

void gzellgzell_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
/* NEGATIVE GRADIENT AND HESSIAN WRT th[j] OF THE GROUP ZELLNER - GROUP ZELLNER PRIOR */

    int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
    double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
    ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
    firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
    idxini= (int) (cholSini[groupid] +.1);
    jingroup= j - firstingroup;
    jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
    (*priorhess)= Sinv[jj];
    for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
    for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }

}

void gzellgzell_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
/* NEGATIVE GRADIENT WRT th[j] OF THE GROUP ZELLNER - GROUP ZELLNER PRIOR */
    int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
    double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
    ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
    firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
    idxini= (int) (cholSini[groupid] +.1);
    jingroup= j - firstingroup;
    jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
    for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
    for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }

}

void gzellgzell_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  /* NEGATIVE HESSIAN FOR GROUP ZELLNER - GROUP ZELLNER PRIOR IS ADDED TO hess

     Input/Output

     - hess: the input value is modified by adding the negative hessian. For example, if on input hess contains the negative hessian of the log-likelihood, then on output hess contains the negative hessian of the log-joint (log-likelihood + log prior)

*/
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    gzellgzell_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] += priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}



//Evaluate negative log-likelihood + log-prior (pMOM + group Zellner + inverse gamma) and initialize funargs
void pmomgzellig_log (double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  pmomgzell_log(f, th, sel, thlength, pars, funargs);
  (*f) -= dinvgammaC(exp(-2.0*th[*thlength -1]), *((*pars).alpha)/2.0, *((*pars).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];

}

//Gradient and hessian wrt th[j] of log pMOM + group Zellner prior
void pmomgzellig_gradhess (double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau);
      (*priorgrad)= 2.0/th[j] - th[j]/tau;
      (*priorhess)= -2.0/(th[j]*th[j]) - 1.0/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      (*priorhess)= - Sinv[jj];
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
      (*priorgrad) = -(*priorgrad);
    }

  } else { //if exp(th[j]) is the residual precision
    (*priorhess)= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * (*priorhess) + (*((*pars).alpha));
  }
}


//Gradient wrt th[j] of log pMOM + group Zellner prior
void pmomgzellig_grad (double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau);
      (*priorgrad)= 2.0/th[j] - th[j]/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
      (*priorgrad) = -(*priorgrad);
    }

  } else { //if exp(th[j]) is the residual precision
    double priorhess= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * priorhess + (*((*pars).alpha));
  }
}


void pmomgzellig_hess (double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    pmomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] -= priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}


//Evaluate negative log-likelihood + log-prior (peMOM + group Zellner + inverse gamma) and initialize funargs
void pemomgzellig_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  pemomgzell_log(f, th, sel, thlength, pars, funargs);
  (*f) -= dinvgammaC(exp(-2.0*th[*thlength -1]), *((*pars).alpha)/2.0, *((*pars).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];

}

void pemomgzellig_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau), th2= th[j]*th[j];
      (*priorgrad)= 2.0 * tau/(th2 * th[j]) - th[j]/tau;
      (*priorhess)= -6.0*tau/(th2*th2) - 1.0/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      (*priorhess)= - Sinv[jj];
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
      (*priorgrad) = -(*priorgrad);
    }

  } else { //if exp(th[j]) is the residual precision
    (*priorhess)= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * (*priorhess) + (*((*pars).alpha));
  }
}

void pemomgzellig_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    if (((*pars).isgroup)[sel[j]] == 0) {
      double tau= *((*pars).tau), th2= th[j]*th[j];
      (*priorgrad)= 2.0 * tau/(th2 * th[j]) - th[j]/tau;
    } else {
      int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
      double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
      ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
      firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
      idxini= (int) (cholSini[groupid] +.1);
      jingroup= j - firstingroup;
      jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
      for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
      for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
      (*priorgrad) = -(*priorgrad);
    }

  } else { //if exp(th[j]) is the residual precision
    double priorhess= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * priorhess + (*((*pars).alpha));
  }
}


void pemomgzellig_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    pemomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] -= priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}



//Evaluate negative log-likelihood + log-prior (group Zellner + group Zellner + inverse gamma) and initialize funargs
void gzellgzellig_log(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {

  gzellgzell_log(f, th, sel, thlength, pars, funargs);
  (*f) -= dinvgammaC(exp(-2.0*th[*thlength -1]), *((*pars).alpha)/2.0, *((*pars).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];

}

void gzellgzellig_gradhess(double *priorgrad, double *priorhess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
    double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
    ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
    firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
    idxini= (int) (cholSini[groupid] +.1);
    jingroup= j - firstingroup;
    jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
    (*priorhess)= - Sinv[jj];
    for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
    for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    (*priorgrad) = -(*priorgrad);

  } else { //if exp(th[j]) is the residual precision
    (*priorhess)= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * (*priorhess) + (*((*pars).alpha));
  }
}


void gzellgzellig_grad(double *priorgrad, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {

  if (j < *thlength -1) { //if th[j] is a regression coefficient

    int jingroup, jj, l, ll, idxini, ningroup, firstingroup, groupid= ((*funargs)["selgroups"])[j];
    double *Sinv= (*funargs)["Sinv"], *cholSini= (*funargs)["cholSini"];
    ningroup= (int) (((*funargs)["nvarinselgroups"])[groupid] +.1);
    firstingroup= (int) (((*funargs)["firstingroup"])[groupid] +.1);
    idxini= (int) (cholSini[groupid] +.1);
    jingroup= j - firstingroup;
    jj= idxini + jingroup*ningroup - jingroup*(jingroup-1)/2;
    for (l=0, (*priorgrad)=0; l< jingroup; l++) { ll= idxini + l*ningroup - l*(l-1)/2; (*priorgrad) +=  Sinv[ll+jingroup-l] * th[firstingroup+l]; } //Sinv[jingroup,] %*% th
    for (l=jingroup; l< ningroup; l++) { (*priorgrad) +=  Sinv[jj+l-jingroup] * th[firstingroup+l]; }
    (*priorgrad) = -(*priorgrad);

  } else { //if exp(th[j]) is the residual precision
    double priorhess= -2.0 * (*((*pars).lambda) * exp(2.0 * th[j]));
    (*priorgrad)= 0.5 * priorhess + (*((*pars).alpha));
  }
}


void gzellgzellig_hess(double **hess, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  /* HESSIAN FOR AFT LOG-LIKELIHOOD Xobs + LOG-LIKELIHOOD Xcens evaluated at any th   */
  int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  for (j=0; j< *thlength; j++) {
    gzellgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, pars, funargs);
    hess[j+1][j+1] -= priorhess;
  }
  //Add group Zellner's prior contribution to off-diagonal hessian elements
  ngroups= (int) (*((*funargs)["nselgroups"]) +.1);
  firstingroup=0;
  for (j=0; j< ngroups; j++) {
    ningroup= (int) (nvaringroup[j] +.1);
    if (ningroup>1) {
      idxini= (int) (cholSini[j]+.1);
      for (k=0; k< ningroup; k++) {
        kk= idxini + k*ningroup - k*(k-1)/2;
        for (l=k+1; l< ningroup; l++) {
          hess[firstingroup+k+1][firstingroup+l+1] +=  *(Sinv + kk + l-k);  //elem (k,l) in Sinv of group j
        }
      }
    }
    firstingroup += ningroup;
  }

}






//*************************************************************************************
// LEAST SQUARES
//*************************************************************************************


void leastsquares(double *theta, double *phi, double *ypred, double *y, double *x, crossprodmat *XtX, double *ytX, int *n, int *p, int *sel, int *nsel) {
  //Least squares estimate for y= x[,sel] %*% theta + e, where e ~ N(0,phi)   (phi is the variance)
  //Input
  // - y: observed response
  // - x: predictors in vector format
  // - n: length(y)
  // - sel: variables in x to be included in the model are given by sel[0], sel[1], etc.
  // - nsel: length(sel)
  //Output
  // - theta[1...*nsel]: least squares estimate for theta
  // - phi: MLE for residual variance (i.e. SSR/n). If <1.0e-10 then phi=1.0e-10 is returned
  // - ypred[0..n-1]: predicted y, i.e. x[,sel] %*% theta where theta is the least squares estimate
  int i;
  double zero=0, **S, **Sinv, detS, e;

  (*phi)= 0;

  if ((*nsel)>0) {
    //Least squares
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,XtX,sel,nsel,p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,ytX,sel,theta);
    free_dmatrix(S, 1,*nsel,1,*nsel); free_dmatrix(Sinv, 1,*nsel,1,*nsel);

    //MLE for residual variance
    Aselvecx(x, theta+1, ypred, 0, (*n) -1, sel, nsel);
    for (i=0; i< (*n); i++) { e= y[i]-ypred[i]; (*phi) += e*e; }

  } else {

    for (i=0; i< (*n); i++) { (*phi) += y[i]*y[i]; }

  }

  (*phi)= (*phi) / (*n);
  if ((*phi) < 1.0e-10) { (*phi)= 1.0e-10; }
}



//*************************************************************************************
// MARGINAL LIKELIHOOD UNDER NORMAL ERRORS
//*************************************************************************************

// pMOM on individual coef, group Zellner on groups
double pmomgzellMarg(int *sel, int *nsel, struct marginalPars *pars) {
  int i_var, i, j_var, j, p_i, varcount, groupcount, *groupsel, singlevarcount=0, *isgroup=pars->isgroup, nu;
  double num, den, ans=0.0, aux, term1, *m, **S, **Sinv, **Vinv, **Vinv_chol, detS, detVinvtau, logdetVinv, tau= *(*pars).tau, tauinv=1/tau, taugroup=*(*pars).taugroup, taugroupinv=1/taugroup, logtaus, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0, *nvarinselgroups, *firstingroup, nselgroups, *selgroups;
  covariancemat *V0inv=(*pars).V0inv;
  bool posdef;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    nvarinselgroups= dvector(0, min_xy(*nsel, *((*pars).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*pars).ngroups))); selgroups= dvector(0, *nsel -1);
    findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*pars).nvaringroup, (*pars).ngroups); //copy subset of nvaringroup into nvarinselgroups
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    Vinv = dmatrix(1, *nsel, 1, *nsel); Vinv_chol = dmatrix(1, *nsel, 1, *nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    for (i=1; i<=*nsel; i++) {
      for (j=i; j<=*nsel; j++) {
        Vinv[i][j] = 0;
      }
    }
    for (varcount=1, groupcount=0; varcount <= *nsel; groupcount++) {
      p_i = (int) nvarinselgroups[groupcount];
      if (p_i==1) {
        Vinv[varcount][varcount] = tauinv;
        S[varcount][varcount]+=tauinv;
        varcount++;
        singlevarcount++;
      } else {
        groupsel = ivector(0,p_i);
        for (i=0; i<p_i; i++) {  groupsel[i] = sel[varcount-1+i]; }
        for (i_var=varcount, i=0; i<p_i; i_var++, i++) {
          for (j_var=i_var, j=i; j<p_i; j_var++, j++) {
            if (V0inv->computed_at(groupsel[i], groupsel[j])) {
              Vinv[i_var][j_var]= V0inv->at(groupsel[i], groupsel[j]);
              S[i_var][j_var]+=Vinv[i_var][j_var];
            } else {
              aux = (*pars).XtX->at(groupsel[i], groupsel[j])*taugroupinv*p_i;
              V0inv->set(groupsel[i], groupsel[j], aux);
              Vinv[i_var][j_var]= aux;
              S[i_var][j_var]+=aux;
            }
          }
        }
        varcount = varcount + p_i;
        free_ivector(groupsel,0,p_i);
      }
    }
    choldc(Vinv, *nsel, Vinv_chol, &posdef);
    detVinvtau = choldc_det(Vinv_chol, *nsel);
    logtaus = singlevarcount*log(tau) + (*nsel - singlevarcount) * log(taugroup);
    logdetVinv = log(detVinvtau) + logtaus;
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    /* nuhalf= .5*(*(*pars).n + *(*pars).alpha); */
    nuhalf= (*(*pars).r)*singlevarcount + .5*(*(*pars).n + *(*pars).alpha);
    nu= (int) (2.0*nuhalf);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI + log(detS) - logdetVinv) + .5 * logtaus + gamln(&alphahalf);
    ans= num - den;

    // Orthogonal approx
    term1= ss / ((double) (nu-2));
    for (i=1; i<=(*nsel); i++) {
      if (isgroup[sel[i-1]] == 0) {
        ans+= log(pow(m[i],2.0) + Sinv[i][i] * term1);
      }
    }

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
    free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vinv_chol,1,*nsel,1,*nsel);

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}

// pMOM on individual coef, groups MOM on groups
double pmomgmomMarg(int *sel, int *nsel, struct marginalPars *pars) {
  int i_var, i, j_var, j, p_i, nu, varcount, groupcount, *groupsel, singlevarcount=0, n=*(*pars).n;
  double num, den, ans=0.0, aux, trSV, term1, *m, *mj, **S, **Sinv, **Vinv, **Vinv_chol, **Vinvj,
         detS, detVinvtau, logdetVinv, tau= *(*pars).tau, tauinv=1/tau, taugroup=*(*pars).taugroup*((double) n),
         taugroupinv=1/taugroup, logtaus, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda),
         ss, zero=0, *nvarinselgroups, *firstingroup, nselgroups, *selgroups;
  covariancemat *V0inv=(*pars).V0inv;
  bool posdef;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    nvarinselgroups= dvector(0, min_xy(*nsel, *((*pars).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*pars).ngroups))); selgroups= dvector(0, *nsel -1);
    findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*pars).nvaringroup, (*pars).ngroups); //copy subset of nvaringroup into nvarinselgroups
    free_dvector(firstingroup, 0, min_xy(*nsel, *((*pars).ngroups)));
    free_dvector(selgroups, 0, *nsel -1);
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    Vinv = dmatrix(1, *nsel, 1, *nsel); Vinv_chol = dmatrix(1, *nsel, 1, *nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    for (i=1; i<=*nsel; i++) {
      for (j=i; j<=*nsel; j++) {
        Vinv[i][j] = 0;
      }
    }
    for (varcount=1, groupcount=0; varcount <= *nsel; groupcount++) {
      p_i = (int) nvarinselgroups[groupcount];
      if (p_i==1) {
        Vinv[varcount][varcount] = tauinv;
        S[varcount][varcount]+= tauinv;
        varcount++;
        singlevarcount++;
      } else {
        groupsel = ivector(0,p_i-1);
        for (i=0; i<p_i; i++) {  groupsel[i] = sel[varcount-1+i]; }
        for (i_var=varcount, i=0; i<p_i; i_var++, i++) {
          for (j_var=i_var, j=i; j<p_i; j_var++, j++) {
            if (V0inv->computed_at(groupsel[i], groupsel[j])) {
              Vinv[i_var][j_var]= V0inv->at(groupsel[i], groupsel[j]);
              S[i_var][j_var]+=Vinv[i_var][j_var];
            } else {
              aux = (*pars).XtX->at(groupsel[i], groupsel[j])*taugroupinv*(p_i+2);
              V0inv->set(groupsel[i], groupsel[j], aux);
              Vinv[i_var][j_var]= aux;
              S[i_var][j_var]+=aux;
            }
          }
        }
        varcount = varcount + p_i;
        free_ivector(groupsel,0,p_i-1);
      }
    }
    choldc(Vinv, *nsel, Vinv_chol, &posdef);
    detVinvtau = choldc_det(Vinv_chol, *nsel);
    logtaus = singlevarcount*log(tau) + (*nsel - singlevarcount) * log(taugroup);
    logdetVinv = log(detVinvtau) + logtaus;
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    /* nuhalf= .5*(*(*pars).n + *(*pars).alpha); */
    nuhalf= (*(*pars).r)*singlevarcount + .5*(*(*pars).n + *(*pars).alpha);
    nu= (int) (2.0*nuhalf);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI + log(detS) - logdetVinv) + .5 * logtaus + gamln(&alphahalf);
    ans= num - den;

    // Orthogonal approx
    term1= ss / ((double) (nu-2));
    for (varcount=1, groupcount=0; varcount <= *nsel; groupcount++) {
      p_i = (int) nvarinselgroups[groupcount];
      if (p_i==1) {
        ans += log(pow(m[varcount],2.0) + Sinv[varcount][varcount] * term1);
        varcount++;
      } else {
        trSV = 0;
        Vinvj = dmatrix(1, p_i, 1, p_i);
        mj = dvector(1, p_i);
        for (i=1, i_var=varcount; i<=p_i; i++, i_var++) {
          Vinvj[i][i] = Vinv[i_var][i_var] / p_i;
          trSV += Vinvj[i][i] * Sinv[i_var][i_var];
          mj[i] = m[i_var];
          for (j=i+1, j_var=i_var+1; j<=p_i; j++, j_var++) {
            Vinvj[i][j] = Vinv[i_var][j_var] / p_i;
            trSV += 2 * Vinvj[i][j] * Sinv[i_var][j_var];
          }
        }
        aux = quadratic_xtAx(mj, Vinvj, 1, p_i);
        ans += log(trSV + aux / term1);
        varcount = varcount + p_i;
        free_dmatrix(Vinvj,1,p_i,1,p_i);
        free_dvector(mj,1,p_i);
      }
    }
    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
    free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vinv_chol,1,*nsel,1,*nsel);
    free_dvector(nvarinselgroups, 0, min_xy(*nsel, *((*pars).ngroups)));
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}

// peMOM on individual coef, group eMOM on groups
double pemomgemomMarg(int *sel, int *nsel, struct marginalPars *pars) {
  Rf_error("peMOM + group eMOM not currently implemented for linear regression");
}


// peMOM on individual coef, block Zellner on groups
double pemomgzellMarg(int *sel, int *nsel, struct marginalPars *pars) {
  Rf_error("peMOM + block Zellner not currently implemented for linear regression");
}

// Zellner on individual coef, block Zellner on groups
double zellgzellMarg (int *sel, int *nsel, struct marginalPars *pars) {
  int i_var, i, j_group, j_var, j, p_i, varcount, groupcount, *groupsel, singlevarcount=0;
  double num, den, ans=0.0, aux, term1, *m, **S, **Sinv, **Vinv, **Vinv_chol, detS, detVinvtau, logdetVinv, tau= *(*pars).tau, taugroup=*(*pars).taugroup, logtaus, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0, *nvarinselgroups, *firstingroup, nselgroups, *selgroups;
  covariancemat *V0inv=(*pars).V0inv;
  bool posdef;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    nvarinselgroups= dvector(0, min_xy(*nsel, *((*pars).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*pars).ngroups))); selgroups= dvector(0, *nsel -1);
    findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*pars).nvaringroup, (*pars).ngroups); //copy subset of nvaringroup into nvarinselgroups
    free_dvector(selgroups, 0, *nsel -1);
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    Vinv = dmatrix(1, *nsel, 1, *nsel); Vinv_chol = dmatrix(1, *nsel, 1, *nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    for (i=1; i<=*nsel; i++) {
      for (j=i; j<=*nsel; j++) {
        Vinv[i][j] = 0;
      }
    }
    for (varcount=1, groupcount=0; varcount <= *nsel; groupcount++) {
      p_i = (int) nvarinselgroups[groupcount];
      if (p_i==1) {
        for (j_group=groupcount; j_group< nselgroups; j_group++) {
          if (nvarinselgroups[j_group] == 1) {
            j_var = firstingroup[j_group] + 1;
            if (V0inv->computed_at(sel[varcount-1], sel[j_var-1])) {
              Vinv[varcount][j_var] = V0inv->at(sel[varcount-1], sel[j_var-1]);
              S[varcount][j_var]+=Vinv[varcount][j_var];
            } else {
              aux = (*pars).XtX->at(sel[varcount-1], sel[j_var-1])/tau;
              V0inv->set(sel[varcount-1], sel[j_var-1], aux);
              Vinv[varcount][j_var]= aux;
              S[varcount][j_var]+=aux;
            }
          }
        }
        varcount++;
        singlevarcount++;
      } else {
        groupsel = ivector(0,p_i);
        for (i=0; i<p_i; i++) {  groupsel[i] = sel[varcount-1+i]; }
        for (i_var=varcount, i=0; i<p_i; i_var++, i++) {
          for (j_var=i_var, j=i; j<p_i; j_var++, j++) {
            if (V0inv->computed_at(groupsel[i], groupsel[j])) {
              Vinv[i_var][j_var]= V0inv->at(groupsel[i], groupsel[j]);
              S[i_var][j_var]+=Vinv[i_var][j_var];
            } else {
              aux = (*pars).XtX->at(groupsel[i], groupsel[j])/taugroup*p_i;
              V0inv->set(groupsel[i], groupsel[j], aux);
              Vinv[i_var][j_var]= aux;
              S[i_var][j_var]+=aux;
            }
          }
        }
        varcount = varcount + p_i;
        free_ivector(groupsel,0,p_i);
      }
    }
    choldc(Vinv, *nsel, Vinv_chol, &posdef);
    detVinvtau = choldc_det(Vinv_chol, *nsel);
    logtaus = singlevarcount*log(tau) + (*nsel - singlevarcount) * log(taugroup);
    logdetVinv = log(detVinvtau) + logtaus;
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= .5*(*(*pars).n + *(*pars).alpha);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI + log(detS) - logdetVinv) + .5 * logtaus + gamln(&alphahalf);
    ans= num - den;

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
    free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vinv_chol,1,*nsel,1,*nsel);
    free_dvector(nvarinselgroups, 0, min_xy(*nsel, *((*pars).ngroups)));
    free_dvector(firstingroup, 0, min_xy(*nsel, *((*pars).ngroups)));

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}

// Zellner on individual coef, normalid on groups
double normidgzellMarg (int *sel, int *nsel, struct marginalPars *pars) {
  int i_var, i, j_var, j, p_i, varcount, groupcount, *groupsel, singlevarcount=0;
  double num, den, ans=0.0, aux, term1, *m, **S, **Sinv, **Vinv, **Vinv_chol, detS, detVinvtau, logdetVinv, tau= *(*pars).tau, tauinv=1/tau, taugroup=*(*pars).taugroup, taugroupinv=1/taugroup, logtaus, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0, *nvarinselgroups, *firstingroup, nselgroups, *selgroups;
  covariancemat *V0inv=(*pars).V0inv;
  bool posdef;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    nvarinselgroups= dvector(0, min_xy(*nsel, *((*pars).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*pars).ngroups))); selgroups= dvector(0, *nsel -1);
    findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*pars).nvaringroup, (*pars).ngroups); //copy subset of nvaringroup into nvarinselgroups
    free_dvector(firstingroup, 0, min_xy(*nsel, *((*pars).ngroups)));
    free_dvector(selgroups, 0, *nsel -1);
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    Vinv = dmatrix(1, *nsel, 1, *nsel); Vinv_chol = dmatrix(1, *nsel, 1, *nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    for (i=1; i<=*nsel; i++) {
      for (j=i; j<=*nsel; j++) {
        Vinv[i][j] = 0;
      }
    }
    for (varcount=1, groupcount=0; varcount <= *nsel; groupcount++) {
      p_i = (int) nvarinselgroups[groupcount];
      if (p_i==1) {
        Vinv[varcount][varcount] = tauinv;
        S[varcount][varcount]+=tauinv;
        varcount++;
        singlevarcount++;
      } else {
        groupsel = ivector(0,p_i);
        for (i=0; i<p_i; i++) {  groupsel[i] = sel[varcount-1+i]; }
        for (i_var=varcount, i=0; i<p_i; i_var++, i++) {
          for (j_var=i_var, j=i; j<p_i; j_var++, j++) {
            if (V0inv->computed_at(groupsel[i], groupsel[j])) {
              Vinv[i_var][j_var]= V0inv->at(groupsel[i], groupsel[j]);
              S[i_var][j_var]+=Vinv[i_var][j_var];
            } else {
              aux = (*pars).XtX->at(groupsel[i], groupsel[j])*taugroupinv*p_i;
              V0inv->set(groupsel[i], groupsel[j], aux);
              Vinv[i_var][j_var]= aux;
              S[i_var][j_var]+=aux;
            }
          }
        }
        varcount = varcount + p_i;
        free_ivector(groupsel,0,p_i);
      }
    }
    choldc(Vinv, *nsel, Vinv_chol, &posdef);
    detVinvtau = choldc_det(Vinv_chol, *nsel);
    logtaus = singlevarcount*log(tau) + (*nsel - singlevarcount) * log(taugroup);
    logdetVinv = log(detVinvtau) + logtaus;
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= .5*(*(*pars).n + *(*pars).alpha);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI + log(detS) - logdetVinv) + .5 * logtaus + gamln(&alphahalf);
    ans= num - den;

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
    free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vinv_chol,1,*nsel,1,*nsel);
    free_dvector(nvarinselgroups, 0, min_xy(*nsel, *((*pars).ngroups)));

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}



//BIC for the Gaussian linear model
double bic_lm (int *sel, int *nsel, struct marginalPars *pars) {
  int n= *((*pars).n);
  double ans=0.0, ndouble= (double) n, *m, penalty, phi, sumy2= *(*pars).sumy2, *ypred;

  if (*((*pars).tau) >= 0) {
    penalty= *((*pars).tau); //if tau>=0, it stores the information criteria's penalty on dimension
  } else {
    penalty= log(ndouble);   //if tau<0, use the BIC penalty
  }

  if (*nsel ==0) {

    phi= sumy2 / ((double) n);
    ans= -0.5 * (((double) n) + ((double) n) * log(2 * M_PI * phi)) + penalty;  //Model dimension is 1 (variance parameter)

  } else {

    m= dvector(1,*nsel); ypred= dvector(0, n-1);

    leastsquares(m, &phi, ypred, (*pars).y, (*pars).x, (*pars).XtX, (*pars).ytX, (*pars).n, (*pars).p, sel, nsel);

    ans= -0.5* (ndouble + ndouble * log(2 * M_PI * phi) + ((double) *nsel + 1) * penalty); //Model dim is nsel + 1 (variance parameter)

    free_dvector(m, 1, *nsel); free_dvector(ypred, 0, n-1);

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }

  return ans;
}







//*************************************************************************************
// MARGINAL LIKELIHOOD UNDER NORMAL / TWO-PIECE NORMAL / LAPLACE / TWO-PIECE LAPLACE RESIDUALS
//*************************************************************************************

//The first nsel elements in sel indicate variables in/out of the model
//Element nsel+1 indicates single-piece if equal to p+1, two-piece if equal to p+2
// - Example1: If nsel=2,p=10 and sel=(2,5,11) then variables 2 and 5 are in the model, and residuals are Normal
// - Example2: If nsel=2,p=10 and sel=(2,5,12) then variables 2 and 5 are in the model, and residuals are two-piece Normal
// - Example3: If nsel=1,p=10 and sel=(0,11) then variables 0 is in the model and residuals are Normal


double pmomMargTP(int *sel, int *nsel, struct marginalPars *pars) {
  int p= (*((*pars).p)), nvars= *nsel -1;
  double ans;

  if (sel[nvars] == p) { //Normal residuals

    ans= pmomMarginalUC(sel, &nvars, pars);

  } else if (sel[nvars]== (p+1)) { //Two-piece Normal residuals

    int prior=1, symmetric=0;
    ans= nlpMargSkewNorm(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+2)) { //Laplace residuals

    int prior=1, symmetric=1;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+3)) { //Asymmetric Laplace residuals

    int prior=1, symmetric=0;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else {
    Rf_error("Invalid residual distribution\n");
  }

  return ans;
}


double pimomMargTP(int *sel, int *nsel, struct marginalPars *pars) {
  int p= (*((*pars).p)), nvars= *nsel -1;
  double ans;

  if (sel[nvars] == p) { //Normal residuals

    ans= pimomMarginalUC(sel, &nvars, pars);

  } else if (sel[nvars]== (p+1)) { //Two-piece Normal residuals

    int prior=2, symmetric=0;
    ans= nlpMargSkewNorm(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+2)) { //Laplace residuals

    int prior=2, symmetric=1;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+3)) { //Asymmetric Laplace residuals

    int prior=2, symmetric=0;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else {
    Rf_error("Invalid residual distribution\n");
  }

  return ans;
}


double pemomMargTP(int *sel, int *nsel, struct marginalPars *pars) {
  int p= (*((*pars).p)), nvars= *nsel -1;
  double ans;

  if (sel[nvars] == p) { //Normal residuals

    ans= pemomMarginalUC(sel, &nvars, pars);

  } else if (sel[nvars]== (p+1)) { //Two-piece Normal residuals

    int prior=3, symmetric=0;
    ans= nlpMargSkewNorm(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+2)) { //Laplace residuals

    int prior=3, symmetric=1;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else if (sel[nvars]== (p+3)) { //Asymmetric Laplace residuals

    int prior=3, symmetric=0;
    ans= nlpMargAlapl(sel, &nvars, pars, &prior, &symmetric);

  } else {
    Rf_error("Invalid residual distribution\n");
  }

  return ans;
}



//*************************************************************************************
// TWO-PIECE LAPLACE ROUTINES
//*************************************************************************************
double pmomMargLaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=1, symmetric=1;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}

double pimomMargLaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=2, symmetric=1;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}

double pemomMargLaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=3, symmetric=1;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}

double pmomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=1, symmetric=0;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}

double pimomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=2, symmetric=0;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}

double pemomMargAlaplU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=3, symmetric=0;
  return nlpMargAlapl(sel, nsel, pars, &prior, &symmetric);
}


double nlpMargAlapl(int *sel, int *nsel, struct marginalPars *pars, int *prior, int *symmetric) {
//Integrated likelihood for linear regression model y= Xtheta + e where e ~ asymmetric Laplace(0,vartheta,alpha) (vartheta prop to variance, alpha gives asymmetry)
// and priors theta ~ pMOM/piMOM/peMOM(0,tau*vartheta), vartheta ~ IG(alpha/2,lambda/2), atanh(alpha) ~ pMOM(0,taualpha)
// Input
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - pars: parameters needed to compute marginal likelihood
// - prior: prior==1 for pMOM, prior==2 for piMOM, prior==3 for peMOM
// - symmetric: symmetric==1 for Laplace residuals, symmetric==0 for asymmetric Laplace residuals
// Output: integrated likelihood
// IMPORTANT: it is assumed that prior dispersion tau was elicited on theta/sqrt(2*vartheta), but lower-level functions operate on theta/sqrt(vartheta), hence we set taulapl= 2*tau. Similarly for vartheta we set lambdalapl= 2*lambda

  bool posdef;
  int maxit= 100, p, n= (*((*pars).n)), *hesstype= ((*pars).hesstype), fixedalpha;
  double ans, *thmode, fmode, **hess, **cholhess, det, *ypred, taulapl, lambdalapl, ftol=0.001, thtol=0.0001;

  taulapl= 2.0 * (*(*pars).tau);
  lambdalapl= 2.0 * (*(*pars).lambda);
  if (*((*pars).fixatanhalpha) > -9999) { fixedalpha= 1; } else { fixedalpha= 0; }
  if ((*symmetric ==0) & (!fixedalpha)) { p= *nsel +2; } else { p= *nsel +1; }
  thmode= dvector(1,p+fixedalpha); hess= dmatrix(1, p+fixedalpha, 1, p+fixedalpha); ypred=dvector(0,n-1);

  postmodeAlaplCDA(thmode, &fmode, hess, sel, nsel, (*pars).n, (*pars).p, (*pars).y, (*pars).x, (*pars).XtX, (*pars).ytX, &maxit, &ftol, &thtol, &taulapl, (*pars).taualpha, (*pars).fixatanhalpha, (*pars).alpha, &lambdalapl, prior, hesstype, symmetric);

  int method= *((*pars).method);
  if ((method!=0) & (method!=1)) method= 0; //If unrecognized method, set to Laplace

  cholhess= dmatrix(1,p,1,p);
  choldc(hess,p,cholhess,&posdef);

  if (!posdef) {
    int i;
    double lmin=0, *vals;
    vals= dvector(1,p);
    eigenvals(hess,p,vals);
    for (i=1; i<=p; i++) if (vals[i]<lmin) lmin= vals[i];
    lmin = -lmin + .01;
    for (i=1; i<=p; i++) hess[i][i] += lmin;
    choldc(hess,p,cholhess,&posdef);
    free_dvector(vals,1,p);
  }
  det= choldc_det(cholhess, p);

  if (method==0) { //Laplace

    ans= -fmode + 0.5 * p * LOG_M_2PI - 0.5*log(det);

  } else if (method==1) { //Monte Carlo

    int i, j, nu=3;
    double *thsim, **cholV, **cholVinv, ctnu= sqrt((nu-2.0)/(nu+.0)), detVinv, term1, term2;

    thsim= dvector(1, p+fixedalpha); cholV= dmatrix(1,p,1,p); cholVinv= dmatrix(1,p,1,p);

    thmode[*nsel +1]= log(thmode[*nsel +1]);
    //if (*symmetric ==0) { thmode[p]= atanh(thmode[p]); }
    if ((*symmetric ==0) & (fixedalpha==0)) { thmode[p]= atanh(thmode[p]); } else if ((*symmetric ==0) & (fixedalpha==1)) { thmode[p+1]= *((*pars).fixatanhalpha); }
    cholS_inv(cholhess, p, cholV);
    for (i=1; i<=p; i++) {
      for (j=1; j<=i; j++) {
        cholV[i][j]= cholV[i][j] * ctnu;
        cholVinv[i][j]= cholhess[i][j] / ctnu;
      }
    }
    detVinv= exp(log(det) - 2*p*log(ctnu));

    ans= 0;
    for (i=1; i<= (*(*pars).B); i++) {
      rmvtC(thsim, p, thmode, cholV, nu);
      if ((*symmetric ==0) & (fixedalpha==1)) { thsim[p+1]= *((*pars).fixatanhalpha); }
      fnegAlapl(&term1,ypred,thsim,sel,nsel,(*pars).n,(*pars).y,(*pars).x,&taulapl,(*pars).taualpha,(*pars).alpha,&lambdalapl,prior,true,symmetric,fixedalpha);
      term1 -= thsim[*nsel +1];
      term2= -dmvtC(thsim, p, thmode, cholVinv, detVinv, nu, 1);
      ans += exp(-term1 + fmode + term2);
    }
    ans= log(ans / ((*(*pars).B)+.0)) - fmode;

    free_dvector(thsim, 1,p+fixedalpha); free_dmatrix(cholV, 1,p,1,p); free_dmatrix(cholVinv, 1,p,1,p);
  }

  free_dmatrix(cholhess, 1,p,1,p);

  if (*((*pars).logscale) == 0) ans= exp(ans);

  free_dvector(thmode, 1,p+fixedalpha); free_dmatrix(hess, 1,p+fixedalpha,1,p+fixedalpha); free_dvector(ypred,0,n-1);
  return(ans);

}


void postmodeAlaplCDA(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *ftol, double *thtol, double *tau, double *taualpha, double *fixatanhalpha, double *alphaphi, double *lambdaphi, int *prior, int *hesstype, int *symmetric) {

  bool useinit= false;
  int i, j, jj, it, p, maxitmle=20, fixedalpha;
  double err, ferr, g, H, delta, fnew, *thnew, *ypred, *fudgeh;

  if (*fixatanhalpha > -9999) { fixedalpha= 1; } else { fixedalpha= 0; }
  if ((*symmetric ==0) & (fixedalpha==0)) { p= *nsel +2; } else { p= *nsel +1; }
  ypred= dvector(0,*n -1); thnew= dvector(1,p+fixedalpha); fudgeh= dvector(1,p);
  for (j=1; j<=p; j++) fudgeh[j]= 1.0;

  //Initialize at MLE
  mleAlaplCDA(thmode,fmode,ypred,sel,nsel,n,pvar,y,x,XtX,ytX,&maxitmle,useinit,symmetric,fixatanhalpha);

  for (i=1; i<=(*nsel); i++) { thnew[i]= thmode[i]; }
  thnew[*nsel +1]= thmode[*nsel +1]; //phi
  if ((*symmetric ==0) & (fixedalpha==0)) {   //alpha
    //(silly avoiding 0 posterior at alpha=0)
    //if (fabs(thmode[p])>0.01) {
    //  thnew[p]= thmode[p];
    //} else { if (thmode[p]<=0) { thmode[p]= thnew[p]= -0.01; } else { thmode[p]= thnew[p]= 0.01; } }
    //Improve initial guess for alpha
    loglnegGradHessAlaplUniv(p-1,&g,&H,thmode,nsel,sel,n,pvar,y,ypred,x,XtX,symmetric);
    if (*prior ==1) {  //pMOM prior
      double s= (1.0 + 1.0/(H*(*taualpha)));
      double ss= sqrt(thmode[p]*thmode[p] + 8.0*(1.0/H)*s);
      if (thmode[p]>0) { thmode[p]= thnew[p]= 0.5*(thmode[p] + ss)/s; } else { thmode[p]= thnew[p]= 0.5*(thmode[p] - ss)/s; }
    } else {  //piMOM prior (also used for peMOM, although no longer exact)
      bool found;
      int root_count;
      double *coef, *real_vector, *imag_vector;
      Polynomial poly;
      PolynomialRootFinder::RootStatus_T status;

      coef= dvector(0,4); real_vector= dvector(0,4); imag_vector= dvector(0,4);
      coef[0]= 2.0*(*taualpha); coef[1]= 0; coef[2]= -2.0; coef[3]= H*thmode[p]; coef[4]= -H;
      poly.SetCoefficients(coef, 4);
      status= poly.FindRoots(real_vector,imag_vector,&root_count);

      if (status == PolynomialRootFinder::SUCCESS) {
        j=0; found= false;
        while ((!found) & (j<=4)) {
          if (fabs(imag_vector[j])<1.0e-5) {
            if (((real_vector[j]>0) & (thmode[p]>0)) | ((real_vector[j]<=0) & (thmode[p]<=0))) {
              thmode[p]= thnew[p]= real_vector[j];
              found= true;
            }
          }
          j++;
        }
      }
      free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
    }
  }

  it=1; err= ferr= 1;
  fnegAlapl(fmode,ypred,thmode,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,true,symmetric,fixedalpha);
  (*fmode) -= thmode[*nsel +1];

  while ((err> *thtol) & (it<(*maxit)) & (ferr> *ftol)) {

    err= ferr= 0;
    for (j=1; j<=p; j++) {

      fpnegAlaplUniv(j,&g,&H,thmode,ypred,sel,nsel,n,pvar,y,x,XtX,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //gradient and hessian
      if (j== *nsel +1) g-= 1.0;
      delta= g/H;
      thnew[j]= thmode[j] - fudgeh[j]*delta;
      fnegAlapl(&fnew,ypred,thnew,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,true,symmetric,fixedalpha);
      fnew -= thnew[*nsel +1];

      if ((fnew< *fmode) & (fudgeh[j]<1)) fudgeh[j]*= 2;
      jj=1;
      while ((fnew> *fmode) && (jj<5)) {
        fudgeh[j]= fudgeh[j]/2;
        thnew[j]= thmode[j]- fudgeh[j]*delta;
        fnegAlapl(&fnew,ypred,thnew,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,true,symmetric,fixedalpha);
        fnew -= thnew[*nsel +1];
        jj++;
      }

      //If new value improves target function, update thmode, fmode
      if (fnew<(*fmode)) {
        err= max_xy(err,fabs(thmode[j]-thnew[j]));
        ferr+= *fmode - fnew;
        thmode[j]= thnew[j];
        (*fmode)= fnew;
      } else {
        Aselvecx(x, thmode+1, ypred, 0, (*n) -1, sel, nsel);
        thnew[j]= thmode[j];
      }

    }

    it++;

  }

  fppnegAlapl(hess,thmode,ypred,sel,nsel,n,pvar,y,x,XtX,tau,taualpha,alphaphi,lambdaphi,prior,symmetric,hesstype); //Hessian

  thmode[*nsel +1]= exp(thmode[*nsel +1]);
  if ((*symmetric== 0) & (fixedalpha==0)) { thmode[p]= tanh(thmode[p]); } else if ((*symmetric ==0) & (fixedalpha==1)) { thmode[p+1]= tanh(*fixatanhalpha); }  //Note: tanh(z)= -1 + 2/(1+exp(-2*z))

  free_dvector(ypred, 0,*n -1); free_dvector(thnew,1,p+fixedalpha); free_dvector(fudgeh,1,p);

}


void fppnegAlapl(double **H, double *th, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric, int *hesstype) {
  int i, j, one=1, nselplus1= (*nsel)+1;
  double **Hprior, *hprioralpha, zero=0;

  Hprior= dmatrix(1,nselplus1,1,nselplus1);
  hprioralpha= dvector(1,1);

  loglnegHessAlapl(H,th,nsel,sel,n,p,y,ypred,x,XtX,symmetric,hesstype);

  if ((*prior)==1) {

    dmomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
        H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    if (*symmetric ==0) {
      dmomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else if ((*prior)==2) {

    dimomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
        H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    if (*symmetric ==0) {
      dimomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else if ((*prior)==3) {

    demomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) {
        H[i][j]= H[j][i]= H[i][j] - Hprior[i][j];
      }
    }

    if (*symmetric ==0) {
      demomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  free_dmatrix(Hprior,1,nselplus1,1,nselplus1);
  free_dvector(hprioralpha,1,1);

}


void mleAlaplCDA(double *thmode, double *fmode, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, bool useinit, int *symmetric, double *fixatanhalpha) {
  //MLE for linear regression with asymmetric Laplace errors using a Coordinate Descent Algorithm
  //Input
  // - useinit: if true then thmode is used as initial value (ypred should contain linear predictor for thmode), else it is initializes at least squares estimator
  //Output
  // - thmode: MLE
  // - fmode: log-likelihood evaluated at the MLE
  // - ypred[0.. n-1] contains linear predictor at MLE
  int i, ii, j, jj, fixedalpha;
  double *thnew, fnew, err, scale, alpha, *fudgeh, g, H, s1, s2;

  if (*fixatanhalpha > -9999) { fixedalpha= 1; } else { fixedalpha= 0; }
  //if ((*symmetric ==0) & (fixedalpha==0)) { p= *nsel +2; } else { p= *nsel +1; }

  //Initialize
  thnew= dvector(1,*nsel +2); fudgeh= dvector(1,*nsel +2);

  if ((*symmetric ==0) & (fixedalpha==0)) {   //if alpha must be estimated, init regression coef to median regression after 5 iter

    int fiveiter=1, issymmetric= 1;
    mleAlaplCDA(thmode,fmode,ypred,sel,nsel,n,p,y,x,XtX,ytX,&fiveiter,false,&issymmetric,fixatanhalpha);
    thmode[*nsel +2]= thnew[*nsel +2]= 0;

  } else {

    if ((*nsel)>0) {
      if (!useinit) { leastsquares(thmode,thmode+(*nsel)+1,ypred,y,x,XtX,ytX,n,p,sel,nsel); }
      for (j=1; j<= *nsel; j++) thnew[j]= thmode[j];
    } else {
      for (i=0; i< *n; i++) ypred[i]=0 ;
    }
    thmode[*nsel +1]= thnew[*nsel +1]= 0;

    if ((*symmetric ==0) & (fixedalpha==1)) { thmode[*nsel +2]= thnew[*nsel +2]= *fixatanhalpha; }
  }

  if ((*symmetric ==0) & (fixedalpha==0)) { thmode[*nsel +2]= thnew[*nsel +2]= 0; } else if ((*symmetric ==0) & (fixedalpha==1)) { thmode[*nsel +2]= thnew[*nsel +2]= *fixatanhalpha; }

  scale= exp(thmode[*nsel +1]);
  if ((*symmetric ==0) & (fixedalpha==0)) {
    alpha= tanh(thmode[*nsel +2]);
  } else if ((*symmetric ==0) & (fixedalpha==1)) {
    alpha= tanh(*fixatanhalpha);
  } else { alpha= 0; }
  loglAlapl(fmode,ypred,thmode,nsel,sel,n,&scale,&alpha,y,x,symmetric);

  //Coordinate descent
  ii=0; err= 1;
  for (j=1; j<= *nsel +2; j++) fudgeh[j]= 1.0;
  while ((err>0.0001) && (ii<(*maxit))) {

    ii++; err= 0;
    //Update theta
    if (*nsel >0) {
      for (j=1; j<= *nsel; j++) {
        loglnegGradHessAlaplUniv(j-1,&g,&H,thmode,nsel,sel,n,p,y,ypred,x,XtX,symmetric);
        thnew[j]= thmode[j]-fudgeh[j]*g/H;
        loglAlapl(&fnew,ypred,thnew,nsel,sel,n,&scale,&alpha,y,x,symmetric);
        jj=1;
        while ((fnew< *fmode) && (jj<5)) {
          fudgeh[j]= fudgeh[j]/2;
          thnew[j]= thmode[j]- fudgeh[j]*g/H;
          loglAlapl(&fnew,ypred,thnew,nsel,sel,n,&scale,&alpha,y,x,symmetric);
          jj++;
        }
        if (fnew > *fmode) {
          err= max_xy(err,fabs(thnew[j]-thmode[j]));
          thmode[j]= thnew[j];
          (*fmode)= fnew;
        } else {
          Aselvecx(x, thmode+1, ypred, 0, (*n) -1, sel, nsel);
          thnew[j]= thmode[j];
        }
      }
    }

    //Update vartheta and alpha
    for (i=0, s1=s2=0; i< *n; i++) { if (y[i]<ypred[i]) { s1+= ypred[i]-y[i]; } else { s2+= y[i]-ypred[i]; } }

    if ((*symmetric ==0) & (fixedalpha==0)) {
      thnew[*nsel +2]= atanh((sqrt(s1) - sqrt(s2))/(sqrt(s1) + sqrt(s2))); //alpha
      thnew[*nsel +1]= log(0.25) - 2.0*log(*n +.0) + 4*log(sqrt(s1) + sqrt(s2)); //vartheta
      err= max_xy(err,max_xy(fabs(thnew[*nsel +1]-thmode[*nsel +1]), fabs(thnew[*nsel +2]-thmode[*nsel +2])));
      thmode[*nsel +2]= thnew[*nsel +2]; thmode[*nsel +1]= thnew[*nsel +1];
      alpha= tanh(thmode[*nsel +2]);
    } else if ((*symmetric ==0) & (fixedalpha==1)) {
      thnew[*nsel +1]= log(s1/(1+alpha)  + s2/(1-alpha)) - log(*n + .0); //vartheta
      err= max_xy(err,fabs(thnew[*nsel +1]-thmode[*nsel +1]));
      thmode[*nsel +1]= thnew[*nsel +1];
    } else {
      thnew[*nsel +1]= 2.0*log(s1+s2) - 2.0*log(*n +.0);
      err= max_xy(err,fabs(thnew[*nsel +1]-thmode[*nsel +1]));
      thmode[*nsel +1]= thnew[*nsel +1];
      alpha= 0;
    }
    scale= exp(thmode[*nsel +1]);
    loglAlapl(fmode,ypred,thmode,nsel,sel,n,&scale,&alpha,y,x,symmetric);

  }

  free_dvector(thnew, 1, *nsel +2); free_dvector(fudgeh, 1, *nsel +2);

}


void fnegAlapl(double *ans, double *ypred, double *th, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, bool logscale, int *symmetric, int fixedalpha) {
//Negative log-joint for two-piece Laplace under MOM/eMOM/iMOM prior on coef and IG on variance
//Note: log-joint evaluated for vartheta, if log-joint for log(vartheta) is desired you need to substract -th[nsel+1] to consider the Jacobian term
// Input
// - th[1..nsel+2]: (theta, log(vartheta), atanh(alpha)) where theta=regression coef, vartheta \propto variance and alpha=asymmetry parameter in [-1,1]
// - Other parameters as in postmodeSkewNorm
// Output:
// - ans: value of minus the log-joint evaluated at th
// - ypred: linear predictor x %*% th
  double scale, alpha;

  scale= exp(th[*nsel +1]);
  if (*symmetric ==0) { alpha= tanh(th[*nsel +2]); } else { alpha= 0; }
  loglAlapl(ans, ypred, th, nsel, sel, n, &scale, &alpha, y, x, symmetric);
  (*ans)= -(*ans);

  if ((*prior)==1) {

    if ((*symmetric ==0) & (fixedalpha==0)) {
      if ((*nsel)>0) {
        (*ans) += -dmomvec(th+1,*nsel,0.0,*tau,scale,1,1) - dmom(th[*nsel +2],0.0,*taualpha,1.0,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -dmom(th[*nsel +2],0.0,*taualpha,1.0,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    } else {
      if ((*nsel)>0) {
        (*ans) += -dmomvec(th+1,*nsel,0.0,*tau,scale,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    }

  } else if ((*prior)==2) {

    if ((*symmetric ==0) & (fixedalpha==0)) {
      if ((*nsel)>0) {
        (*ans) += -dimomvec(th+1,*nsel,0.0,*tau,scale,1) - dimom(th[*nsel +2],0.0,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -dimom(th[*nsel +2],0.0,*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    } else {
      if ((*nsel)>0) {
        (*ans) += -dimomvec(th+1,*nsel,0.0,*tau,scale,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    }

  } else if ((*prior)==3) {

    if ((*symmetric ==0) & (fixedalpha==0)) {
      if ((*nsel)>0) {
        (*ans) += -demomvec(th+1,*nsel,*tau,scale,1) - demom(th[*nsel +2],*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -demom(th[*nsel +2],*taualpha,1.0,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    } else {
      if ((*nsel)>0) {
        (*ans) += -demomvec(th+1,*nsel,*tau,scale,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      } else {
        (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
      }
    }

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  if (!logscale) { (*ans)= exp(*ans); }
}


void fpnegAlaplUniv(int j, double *g, double *H, double *th, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric) {
  //Univariate gradient and hessian of fnegAlapl wrt j^th element in th
  //Note: the -1.0 term when j=nsel+1 coming from the jacobian of tvartheta=log(vartheta) is not included, it must be added separately after calling this function
  int i;
  double gprior, hprior, zero=0, sumth2, suminvth2;

  loglnegGradHessAlaplUniv(j-1,g,H,th,nsel,sel,n,p,y,ypred,x,XtX,symmetric);

  if ((*prior)==1) {   //MOM prior

    if (j <= (*nsel)) {
      gprior= dmomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0; i<=(*nsel); i++) { sumth2 += th[i]*th[i]; }
      gprior= -1.5*(*nsel) - 0.5*(*alphaphi) -1.0 + 0.5*(sumth2/(*tau) + *lambdaphi) * exp(-th[*nsel +1]);
    } else {
      gprior= dmomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }
    (*g) -= gprior;

    if (j <= (*nsel)) {
      hprior= dmomhessuniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0;i<=(*nsel);i++) sumth2 += pow(th[i],2.0);
      hprior= -0.5 * exp(-th[(*nsel)+1]) * (sumth2/(*tau)+(*lambdaphi));
    } else {
      hprior= dmomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else if ((*prior)==2) {   //iMOM prior

    if (j <= (*nsel)) {
      gprior= dimomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, suminvth2=0; i<=(*nsel); i++) { suminvth2 += 1.0/(th[i]*th[i]); }
      gprior= 0.5*(*nsel) - 0.5*(*alphaphi) -1.0 + 0.5*(*lambdaphi)*exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      gprior= dimomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }
    (*g) -= gprior;

    if (j <= (*nsel)) {
      hprior= dimomhessuniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, suminvth2=0; i<=(*nsel); i++) suminvth2 += pow(1.0/th[i],2.0);
      hprior= -0.5*exp(-th[(*nsel)+1])*(*lambdaphi) - (*tau) * exp(th[(*nsel)+1]) * suminvth2;
    } else {
      hprior= dimomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else if ((*prior)==3) {   //eMOM prior

    if (j <= (*nsel)) {
      gprior= demomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0, suminvth2=0; i<=(*nsel); i++) { sumth2 += th[i]*th[i]; suminvth2 += 1.0/(th[i]*th[i]); }
      gprior= -0.5*(*nsel) - 0.5*(*alphaphi) - 1.0 +0.5*(sumth2/(*tau) + *lambdaphi) * exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      gprior= demomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }

    (*g) -= gprior;

    if (j<=(*nsel)) {
      hprior= demomhessuniv(th+j, th+(*nsel)+1, tau);

    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0, suminvth2=0; i<=(*nsel); i++) { sumth2+= pow(th[i],2.0); suminvth2 += pow(1.0/th[i],2.0); }
      hprior= -0.5*(*nsel) - 0.5*(*alphaphi) -1.0 + 0.5*(sumth2/(*tau) + (*lambdaphi)) * exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      hprior= demomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

}




void loglAlapl(double *ans, double *ypred, double *th, int *nsel, int *sel, int *n, double *scale, double *alpha, double *y, double *x, int *symmetric) {
  //Log-likelihood function of a linear model with two-piece Laplace errors evaluated at th=(theta,scale,alpha)
  //Output
  // - ans: value of the log-likelihood evaluated at th
  // - ypred: linear predictor x %*% th
  int i;
  double w1, w2;

  (*ans)= 0;
  if (*symmetric ==0) {

    w1= 1.0 / ((1.0 + (*alpha)) * sqrt(*scale));
    w2= 1.0 / ((1.0 - (*alpha)) * sqrt(*scale));
    if ((*nsel)>0) {
      Aselvecx(x, th+1, ypred, 0, (*n) -1, sel, nsel); //ypred= x %*% th
      for (i=0; i<(*n); i++) {
        if (y[i]<ypred[i]) { (*ans) -= w1 * (ypred[i]-y[i]); } else { (*ans) -= w2 * (y[i]-ypred[i]); }
      }
    } else {
      for (i=0; i<(*n); i++) {
        if (y[i]<0) { (*ans) -= w1 * fabs(y[i]); } else { (*ans) -= w2 * fabs(y[i]); }
      }
    }

  } else {

    if ((*nsel)>0) {
      Aselvecx(x, th+1, ypred, 0, (*n) -1, sel, nsel); //ypred= x %*% th
      for (i=0; i<(*n); i++) { (*ans) -= fabs(y[i]-ypred[i]); }
    } else {
      for (i=0; i<(*n); i++) { (*ans) -= fabs(y[i]); }
    }
    (*ans)= (*ans) / sqrt(*scale);

  }

  (*ans)+= -(*n +.0)*log(2.0) -0.5*(*n +.0)*log(*scale);
}


void loglnegGradHessAlaplUniv(int j, double *g, double *H, double *th, int *nsel, int *sel, int *n, int *p, double *y, double *ypred, double *x, crossprodmat *XtX, int *symmetric) {
  //Gradient for th[j] of minus the log-likelihood function of a linear model with two-piece Laplace errors
  //Note: the -1.0 term when j=nsel+1 coming from the jacobian of tvartheta=log(vartheta) is not included, it must be added separately after calling this function
  int i, colidx;
  double scale, sqscale, alpha, alphasq, w1, w2, ws1, ws2, wsbar1, wsbar2, tmp;

  scale= exp(th[*nsel +1]);
  sqscale= sqrt(scale);
  (*g)= (*H)= 0;

  if (*symmetric ==0) {
    alpha= tanh(th[*nsel +2]);
    alphasq= alpha*alpha;
    w1= 1.0 / (1.0 + alpha);
    w2= 1.0 / (1.0 - alpha);

    if (j< *nsel) {  //derivative wrt theta

      colidx= sel[j]*(*n);
      for (i=0; i< *n; i++) {
        if (y[i]<ypred[i]) { (*g)+= w1 * x[colidx+i]; } else { (*g)-= w2 * x[colidx+i]; }
      }
      (*g)= (*g)/sqscale;
      (*H)= (XtX->at(sel[j]*(*p)+sel[j])) / (scale*(1-alphasq));

    } else if (j== *nsel) {  //derivative wrt vartheta

      for (i=0; i< *n; i++) { if (y[i]<ypred[i]) { (*g)+= w1 * (ypred[i]-y[i]); } else { (*g)+= w2 * (y[i]-ypred[i]); } }
      (*H)= 0.25 * (*g) / sqscale;
      (*g)= (*n - (*g)/sqscale) * 0.5;

    } else {  //derivative wrt alpha

      ws1= exp(-2.0 * th[*nsel +2]); ws2= exp(2.0 * th[*nsel +2]);
      wsbar1= ws1; wsbar2= -ws2;
      for (i=0; i< *n; i++) {
        tmp= y[i]-ypred[i];
        if (tmp<0) {
          (*g)-= wsbar1*tmp;
          (*H)-= ws1*tmp;
        } else {
          (*g)+= wsbar2*tmp;
          (*H)+= ws2*tmp;
        }
      }
      (*g)= -(*g)/sqscale;
      (*H)*= 2.0/sqscale;

    }

  } else {  // symmetric Laplace residuals

    if (j< *nsel) {  //derivative wrt theta

      colidx= sel[j]*(*n);
      for (i=0; i< *n; i++) { if (y[i]<ypred[i]) { (*g)+= x[colidx+i]; } else { (*g)-= x[colidx+i]; } }
      (*g)= (*g)/sqscale;
      (*H)= (XtX->at(sel[j]*(*p)+sel[j])) / scale;

    } else {  //derivative wrt vartheta

      for (i=0; i< *n; i++) { if (y[i]<ypred[i]) { (*g)+= (ypred[i]-y[i]); } else { (*g)+= (y[i]-ypred[i]); } }
      (*H)= 0.25 * (*g) / sqscale;
      (*g)= (*n - (*g)/sqscale) * 0.5;

    }

  }

}


void loglnegHessAlapl(double **H, double *th, int *nsel, int *sel, int *n, int *p, double *y, double *ypred, double *x, crossprodmat *XtX, int *symmetric, int *hesstype) {
  int i, j, k, npar;
  double alpha, alphasq, scale, sqscale, wy0, sumwsy0, sumwsbary0, w1=1, w2=1, ws1, ws2, wsbar1, wsbar2, *Xtwbar, *Xtws, *hdiag, *D, *y0;
  scale= exp(th[*nsel +1]);
  sqscale= sqrt(scale);

  Xtwbar= dvector(0,*nsel);
  y0= dvector(0,*n);

  if (*symmetric ==0) { //asymmetric Laplace errors

    Xtws= dvector(0,*nsel);
    alpha= tanh(th[*nsel +2]);
    alphasq= alpha*alpha;
    w1= 1.0 / (1.0 + alpha);
    w2= 1.0 / (1.0 - alpha);

    for (j=0; j< *nsel; j++) {  //hessian wrt theta
      H[j+1][j+1]= (XtX->at(sel[j]*(*p)+sel[j])) / (scale*(1-alphasq));
      for (k=0; k<j; k++) { H[j+1][k+1]= H[k+1][j+1]= (XtX->at(sel[k]*(*p)+sel[j])) / (scale*(1-alphasq)); }
      Xtwbar[j]= Xtws[j]= 0;
    }

    ws1= exp(-2.0 * th[*nsel +2]); ws2= exp(2.0 * th[*nsel +2]);
    wsbar1= ws1; wsbar2= -ws2;
    wy0= sumwsy0= sumwsbary0=  0;
    for (i=0; i< *n; i++) {
      y0[i]= y[i]-ypred[i];
      if (y[i]<ypred[i]) {
        wy0-= w1 * y0[i];
        sumwsy0-= ws1 * y0[i];
        sumwsbary0-= wsbar1 * y0[i];
        for (j=0; j< *nsel; j++) { Xtwbar[j]+= w1* x[sel[j]*(*n) +i]; Xtws[j]+= ws1* x[sel[j]*(*n) +i]; }
      } else {
        wy0+= w2 * y0[i];
        sumwsy0+= ws2 * y0[i];
        sumwsbary0+= wsbar2 * y0[i];
        for (j=0; j< *nsel; j++) { Xtwbar[j]-= w2* x[sel[j]*(*n) +i]; Xtws[j]+= ws2* x[sel[j]*(*n) +i]; }
      }
    }

    H[*nsel +1][*nsel +1]= 0.25*wy0/sqscale; //hessian wrt vartheta
    H[*nsel +2][*nsel +2]= 2*sumwsy0/sqscale; //hessian wrt alpha
    H[*nsel +1][*nsel +2]= H[*nsel +2][*nsel +1]= 0.5*sumwsbary0/sqscale; //hessian wrt vartheta, alpha

    for (j=0; j< *nsel; j++) {
      H[j+1][*nsel +1]= H[*nsel +1][j+1]= -0.5*Xtwbar[j]/sqscale; //hessian wrt theta, vartheta
      H[j+1][*nsel +2]= H[*nsel +2][j+1]= -Xtws[j]/sqscale; //hessian wrt theta, alpha
    }

    free_dvector(Xtws,0,*nsel);

  } else {  //symmetric Laplace errors

    for (j=0; j< *nsel; j++) {  //hessian wrt theta
      H[j+1][j+1]= (XtX->at(sel[j]*(*p)+sel[j])) / scale;
      for (k=0; k<j; k++) { H[j+1][k+1]= H[k+1][j+1]= (XtX->at(sel[k]*(*p)+sel[j])) / scale; }
      Xtwbar[j]= 0;
    }

    wy0= 0;
    for (i=0; i< *n; i++) {
      y0[i]= y[i]-ypred[i];
      if (y[i]<ypred[i]) {
        wy0-= y0[i];
        for (j=0; j< *nsel; j++) { Xtwbar[j]+= x[sel[j]*(*n) +i]; }
      } else {
        wy0+= y0[i];
        for (j=0; j< *nsel; j++) { Xtwbar[j]-= x[sel[j]*(*n) +i]; }
      }
    }

    H[*nsel +1][*nsel +1]= 0.25*wy0/sqscale; //hessian wrt vartheta

    for (j=0; j< *nsel; j++) {
      H[j+1][*nsel +1]= H[*nsel +1][j+1]= -0.5*Xtwbar[j]/sqscale; //hessian wrt theta, vartheta
    }

  }

  free_dvector(Xtwbar,0,*nsel);

  if (*hesstype ==2) {

    if (*symmetric ==0) { npar= *nsel +2; } else { npar= *nsel +1; }
    hdiag= dvector(1,npar); D= dvector(1,npar);

    quadapproxALaplace(hdiag, H, nsel, sel, n, y0, x, th, &scale, &alpha, &wy0, symmetric, &w1, &w2);
    for (i=1; i<= *nsel; i++) D[i]= sqrt(hdiag[i]/H[i][i]);
    hdiag[*nsel +1]= H[*nsel +1][*nsel +1]; D[*nsel +1]= 1;
    if (*symmetric ==0) { hdiag[*nsel +2]= H[*nsel +2][*nsel +2]; D[*nsel +2]= 1; }

    for (i=1; i<=npar; i++) {
      H[i][i]= hdiag[i];
      for (j=1; j< i; j++) H[i][j]= H[j][i]= H[i][j] * D[i] * D[j];
    }
    free_dvector(hdiag,1,npar); free_dvector(D,1,npar);

  }

  free_dvector(y0,0,*n);
}


void quadapproxALaplace(double *hdiag, double **H, int *nsel, int *sel, int *n, double *y0, double *x, double *th, double *vartheta, double *alpha, double *wy0, int *symmetric, double *w1, double *w2) {
  //Diagonal elements of the hessian in a quadratic approximation to asymmetric Laplace log-likelihood
  // Input
  // - H: asymptotic hessian
  // - y0: residuals y - X th where th is the MLE
  // - x: matrix with predictors
  // - th: vector containing MLE for theta
  // - vartheta: MLE for vartheta parameter (must be >0)
  // - alpha: MLE for asymmetry parameter (must be in (-1,1))
  // - wy0: sum of weighted absolute errors at MLE (i.e. log-likelihood ignoring terms depending on vartheta)
  // Output
  // - hdiag: diagonal terms of the hessian
  int i, j, k, colidx;
  double *l, *fl, f0, *e, l2, suml2, suml4, ct;

  l= dvector(1,2);
  fl= dvector(1,2);
  e= dvector(0,*n -1);

  f0= *wy0;
  ct= 2.0/sqrt(*vartheta);
  for (j=1; j<= *nsel; j++) {
    l[2]= 1.96/sqrt(H[j][j]);
    l[1]= -l[2];
    colidx= sel[j-1]* (*n);
    suml2= suml4= 0;
    for (k=1; k<=2; k++) {
      fl[k]= 0;
      if (*symmetric ==0) { //asymmetric Laplace errors
        for (i=0; i< *n; i++) {
          e[i]= y0[i] - l[k] * x[i + colidx];
          if (e[i]<0) { fl[k]-= (*w1) * e[i]; } else { fl[k]+= (*w2) * e[i]; }
        }
      } else { //symmetric Laplace errors
        for (i=0; i< *n; i++) {
          e[i]= y0[i] - l[k] * x[i + colidx];
          if (e[i]<0) { fl[k]-= e[i]; } else { fl[k]+= e[i]; }
        }
      }
      l2= l[k]*l[k];
      suml2+= l2 * (fl[k]-f0);
      suml4+= l2*l2;
      hdiag[j]= ct * suml2 / suml4;
    }
  }

  free_dvector(l,1,2);
  free_dvector(fl,1,2);
  free_dvector(e,0,*n -1);

}



//*************************************************************************************
// TWO-PIECE NORMAL ROUTINES
//*************************************************************************************
double pmomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=1, symmetric=0;
  return nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);
}

double pimomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=2, symmetric=0;
  return nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);
}

double pemomMargSkewNormU(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=3, symmetric=0;
  return nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);
}


double nlpMargSkewNorm(int *sel, int *nsel, struct marginalPars *pars, int *prior, int *symmetric) {
//Integrated likelihood for linear regression model y= Xtheta + e where e ~ two-piece Normal(0,vartheta,alpha) (vartheta prop to variance, alpha gives asymmetry)
// and priors theta ~ pMOM/piMOM/peMOM(0,tau*vartheta), vartheta ~ IG(alpha/2,lambda/2), atanh(alpha) ~ pMOM(0,taualpha)
// Input
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - pars: parameters needed to compute marginal likelihood
// - prior: prior==1 for pMOM, prior==2 for piMOM, prior==3 for peMOM
// - symmetric: symmetric==1 for Normal residuals, symmetric==0 for two-piece normal residuals
// Output: integrated likelihood

  bool initmle=true, posdef;
  int maxit= 50, p, n= (*((*pars).n));
  double ans, *thmode, fmode, **hess, **cholhess, det, *ypred, ftol=0.001, thtol=0.0001;

  if (*symmetric ==0) { p= *nsel +2; } else { p= *nsel +1; }
  thmode= dvector(1,p); hess= dmatrix(1, p, 1, p); ypred=dvector(0,n-1);

  if ((*symmetric ==1) | (*((*pars).optimMethod) != 1)) { //Coordinate Descent Algorithm
    postmodeSkewNormCDA(thmode,&fmode,hess,sel,nsel,(*pars).n,(*pars).p,(*pars).y,(*pars).x,(*pars).XtX,(*pars).ytX,&maxit,&ftol,&thtol,(*pars).tau,(*pars).taualpha,(*pars).alpha,(*pars).lambda,prior,symmetric);
  } else {  //LMA (modified Newton-Raphson)
    postmodeSkewNorm(thmode,&fmode,hess,sel,nsel,(*pars).n,(*pars).p,(*pars).y,(*pars).x,(*pars).XtX,(*pars).ytX,&maxit,(*pars).tau,(*pars).taualpha,(*pars).alpha,(*pars).lambda,&initmle, prior);
  }

  int method= *((*pars).method);
  if ((method!=0) & (method!=1)) method= 0; //If unrecognized method, set to Laplace

  cholhess= dmatrix(1,p,1,p);
  choldc(hess,p,cholhess,&posdef);
  if (!posdef) {
    int i;
    double lmin=0, *vals;
    vals= dvector(1,p);
    eigenvals(hess,p,vals);
    for (i=1; i<=p; i++) if (vals[i]<lmin) lmin= vals[i];
    lmin = -lmin + .01;
    for (i=1; i<=p; i++) hess[i][i] += lmin;
    choldc(hess,p,cholhess,&posdef);
    free_dvector(vals,1,p);
  }
  det= choldc_det(cholhess, p);

  if (method==0) { //Laplace

    ans= -fmode + 0.5 * p * LOG_M_2PI - 0.5*log(det);

  } else if (method==1) { //Monte Carlo

    int i, j, nu=3;
    double *thsim, **cholV, **cholVinv, ctnu= sqrt((nu-2.0)/(nu+.0)), detVinv, term1, term2;

    thsim= dvector(1, p); cholV= dmatrix(1,p,1,p); cholVinv= dmatrix(1,p,1,p);

    thmode[*nsel +1]= log(thmode[*nsel +1]);
    if (*symmetric ==0) thmode[p]= atanh(thmode[p]);
    cholS_inv(cholhess, p, cholV);
    for (i=1; i<=p; i++) {
      for (j=1; j<=i; j++) {
        cholV[i][j]= cholV[i][j] * ctnu;
        cholVinv[i][j]= cholhess[i][j] / ctnu;
      }
    }
    detVinv= exp(log(det) - 2*p*log(ctnu));

    ans= 0;
    for (i=1; i<= (*(*pars).B); i++) {
      rmvtC(thsim, p, thmode, cholV, nu);
      fnegSkewnorm(&term1,ypred,thsim,sel,nsel,(*pars).n,(*pars).y,(*pars).x,(*pars).XtX,(*pars).tau,(*pars).taualpha,(*pars).alpha,(*pars).lambda,prior,true,symmetric);
      term1 -= thsim[*nsel +1];
      term2= -dmvtC(thsim, p, thmode, cholVinv, detVinv, nu, 1);
      ans += exp(-term1 + fmode + term2);
    }
    ans= log(ans / ((*(*pars).B)+.0)) - fmode;

    free_dvector(thsim, 1,p); free_dmatrix(cholV, 1,p,1,p); free_dmatrix(cholVinv, 1,p,1,p);
  }

  free_dmatrix(cholhess, 1,p,1,p);

  if (*((*pars).logscale) == 0) ans= exp(ans);

  free_dvector(thmode, 1,p); free_dmatrix(hess, 1,p,1,p); free_dvector(ypred,0,n-1);
  return(ans);

}



void postmodeSkewNorm(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *tau, double *taualpha, double *alpha, double *lambda, bool *initmle, int *prior) {
//Posterior mode for two-piece normal under pMOM, piMOM or peMOM prior on (theta,atanh(alpha)) and vartheta ~ IG(alpha/2,lambda/2)
//
//Maximization is done via modified Newton-Raphson (LMA)
//
  // Input
  // - y: observed response
  // - x: design matrix
  // - maxit: maximum number of iterations
  // - tau: dispersion parameter for MOM prior on theta
  // - tau.alpha: dispersion parameter for MOM prior on atanh(alpha)
  // - alpha, lambda: prior on vartheta ~ IG(alpha/2,lambda/2)
  // - init: init=='mle' to initialize at MLE, else initialize at least squares
  // Ouput
  // - thmode[1..nsel]: posterior mode for (theta,vartheta,alpha) (i.e. in original parameterization)
  // - fmode: minus log-joint evaluated at thmode
  // - hess: hessian evaluated at thmode

  bool posdef;
  int i, ii, j, p=(*nsel)+2, itermle=10, symmetric=0;
  double err, ferr, damp, *g, **H, **Hinv, *delta, lmin, *vals, fnew, *thnew, *ypred;

  ypred= dvector(0,*n -1);

  if (*initmle) {  //Initialize at MLE

    mleSkewnorm(thmode, ypred, sel, nsel, n, pvar, y, x, XtX, ytX, &itermle, false);

  } else {  //Initialize at least-squares for theta; set (phi,alpha)= argmax likelihood for given theta

    double s1=0, s2=0, pows1, pows2;

    leastsquares(thmode, thmode+(*nsel)+1, ypred, y, x, XtX, ytX, n, pvar, sel, nsel);

    for (i=0; i<(*n); i++) {
      if (y[i]<=ypred[i]) { s1+= pow(y[i]-ypred[i], 2.0); } else { s2+= pow(y[i]-ypred[i], 2.0); }
    }

    pows1= pow(s1, 1.0/3.0); pows2= pow(s2, 1.0/3.0);
    thmode[p+2]= (pows1 - pows2)/(pows1 + pows2);  //estimate for alpha
    thmode[p+1]= (0.25/((*n)+.0)) * pow(pows1 + pows2, 3.0);  //estimate for phi

  }

  thmode[p-1]= log(thmode[p-1]); //phi
  thmode[p]= atanh(thmode[p]);   //alpha (Note: atanh(z)= 0.5*(log(1+z)-log(1-z)))

  g= dvector(1,p); delta= dvector(1,p); thnew= dvector(1,p);
  H= dmatrix(1,p,1,p); Hinv= dmatrix(1,p,1,p);

  i=ii=1; err=ferr=1; damp=2.0;

  fnegSkewnorm(fmode,ypred,thmode,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true,&symmetric);
  (*fmode) -= thmode[p-1];

  while ((err>0.001) & (i<(*maxit)) & (ferr>0.001)) {

    fpnegSkewnorm(g,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alpha,lambda,prior); //gradient
    g[p-1]-= 1.0; //jacobian term from tvartheta=log(vartheta)
    fppnegSkewnorm(H,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alpha,lambda,prior,&symmetric); //Hessian

    inv_posdef(H,p,Hinv,&posdef);

    if (posdef) {
      Ax(Hinv,g,delta,1,p,1,p);
    } else {
      //Ensure H is posdef
      vals= dvector(1,p);
      eigenvals(H,p,vals);
      lmin= vals[1];
      for (j=2; j<=p; j++) if (vals[j]<lmin) lmin= vals[j];
      lmin = -lmin + .01;
      for (j=1; j<=p; j++) H[j][j] += lmin;
      choldc_inv(H,p,Hinv,&posdef);
      Ax(Hinv,g,delta,1,p,1,p);
      free_dvector(vals,1,p);
    }

    for (j=1; j<=p; j++) { thnew[j]= thmode[j] - delta[j]; }
    fnegSkewnorm(&fnew,ypred,thnew,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true,&symmetric);
    fnew -= thnew[p-1];

    //If Newton update fails, use Levenberg-Marquardt (LMA)
    ii= 1;
    if (fnew > (*fmode)) {
      while ((fnew > (*fmode)) & (ii<5)) {
        for (j=1; j<=p; j++) H[j][j] *= damp;
        inv_posdef(H,p,Hinv,&posdef);
        Ax(Hinv,g,delta,1,p,1,p);
        for (j=1; j<=p; j++) { thnew[j]= thmode[j] - delta[j]; }
        fnegSkewnorm(&fnew,ypred,thnew,sel,nsel,n,y,x,XtX,tau,taualpha,alpha,lambda,prior,true,&symmetric);
        fnew -= thnew[p-1];
        ii++;
      }
    }

    //If new value improves target function, update thmode, fmode
    if (fnew<(*fmode)) {
      ferr= *fmode - fnew;
      err= 0;
      for (j=1; j<=p; j++) { err= max_xy(err,fabs(delta[j])); thmode[j]= thnew[j]; }
      (*fmode)= fnew;
      i++;
    } else {
      i= (*maxit);
      for (j=1; j<=p; j++) { thnew[j]= thmode[j]; }
    }
  }

  thmode[p-1]= exp(thmode[p-1]);
  thmode[p]= tanh(thmode[p]); //Note: tanh(z)= -1 + 2/(1+exp(-2*z))

  if (ii==1) { //LMA update not used in last iteration
    for (i=1; i<=p; i++) {
      hess[i][i]= H[i][i];
      for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= H[i][j]; }
    }
  } else {  //LMA update used, which modifies H
    damp= pow(damp,ii - 1.0);
    for (i=1; i<=p; i++) {
      hess[i][i]= H[i][i] / damp;
      for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= H[i][j]; }
    }
  }

  free_dvector(ypred, 0,*n -1); free_dvector(g,1,p); free_dvector(delta,1,p); free_dvector(thnew,1,p);
  free_dmatrix(H,1,p,1,p); free_dmatrix(Hinv,1,p,1,p);
}



void postmodeSkewNormCDA(double *thmode, double *fmode, double **hess, int *sel, int *nsel, int *n, int *pvar, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, double *ftol, double *thtol, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric) {
//Posterior mode for two-piece normal under pMOM, piMOM or peMOM prior on (theta,atanh(alpha)) and vartheta ~ IG(alpha/2,lambda/2)
//
//Maximization is done via Coordinate Descent Algorithm (i.e. iterative univariate minimization)
//
  // Input
  // - y: observed response
  // - x: design matrix
  // - maxit: maximum number of iterations
  // - tau: dispersion parameter for MOM prior on theta
  // - tau.alpha: dispersion parameter for MOM prior on atanh(alpha)
  // - alphaphi, lambdaphi: prior on vartheta ~ IG(alphaphi/2,lambdaphi/2)
  // - init: init=='mle' to initialize at MLE, else initialize at least squares
  // - symmetric: set symmetric==0 for two-piece normal residuals, otherwise normal residuals are assumed
  // Ouput
  // - thmode[1..nsel]: posterior mode for (theta,vartheta,alpha) (i.e. in original parameterization)
  // - fmode: minus log-joint evaluated at thmode
  // - hess: hessian evaluated at thmode

  int i, j, it, p;
  double err, ferr, g, H, delta, fnew, *thnew, *ypred, s1=0, s2=0, pows1, pows2, acur, aa, bb, sumth2=0;

  if (*symmetric ==0) { p= *nsel +2; } else { p= *nsel +1; }
  ypred= dvector(0,*n -1); thnew= dvector(1,p);

  //Initialize at least squares estimate
  leastsquares(thmode, thmode+(*nsel)+1, ypred, y, x, XtX, ytX, n, pvar, sel, nsel);
  for (i=0; i<(*n); i++) { if (y[i]<=ypred[i]) { s1+= pow(y[i]-ypred[i], 2.0); } else { s2+= pow(y[i]-ypred[i], 2.0); } }
  if (*symmetric ==0) {
    pows1= pow(s1, 1.0/3.0); pows2= pow(s2, 1.0/3.0);
    thmode[p]= (pows1 - pows2)/(pows1 + pows2);
    thmode[*nsel +1]= (0.25/((*n)+.0)) * pow(pows1 + pows2, 3.0);
  } else {
    thmode[*nsel +1]= (s1+s2)/(*n +.0);
  }

  for (i=1; i<=(*nsel); i++) { thnew[i]= thmode[i]; }
  thnew[*nsel +1]= thmode[*nsel +1]= log(thmode[*nsel +1]);
  if (*symmetric ==0) {
    thnew[p]= thmode[p]= atanh(thmode[p]); //Note: atanh(z)= 0.5*(log(1+z)-log(1-z))
    //Improve initial guess for alpha
    loglnegGradSkewNormUniv(p,&g,thmode,nsel,sel,n,y,ypred,x,symmetric); //gradient
    loglnegHessSkewNormUniv(p,&H,thmode,nsel,sel,n,y,ypred,x,symmetric);  //hessian
    if (*prior ==1) {  //pMOM prior
      double s= (1.0 + 1.0/(H*(*taualpha)));
      double ss= sqrt(thmode[p]*thmode[p] + 8.0*(1.0/H)*s);
      if (thmode[p]>0) { thmode[p]= thnew[p]= 0.5*(thmode[p] + ss)/s; } else { thmode[p]= thnew[p]= 0.5*(thmode[p] - ss)/s; }
    } else {  //piMOM prior (also used for peMOM, although no longer exact)
      bool found;
      int root_count;
      double *coef, *real_vector, *imag_vector;
      Polynomial poly;
      PolynomialRootFinder::RootStatus_T status;

      coef= dvector(0,4); real_vector= dvector(0,4); imag_vector= dvector(0,4);
      coef[0]= 2.0*(*taualpha); coef[1]= 0; coef[2]= -2.0; coef[3]= H*thmode[p]; coef[4]= -H;
      poly.SetCoefficients(coef, 4);
      status= poly.FindRoots(real_vector,imag_vector,&root_count);

      if (status == PolynomialRootFinder::SUCCESS) {
        j=0; found= false;
        while ((!found) & (j<=4)) {
          if (fabs(imag_vector[j])<1.0e-5) {
            if (((real_vector[j]>0) & (thmode[p]>0)) | ((real_vector[j]<=0) & (thmode[p]<=0))) {
              thmode[p]= thnew[p]= real_vector[j];
              found= true;
            }
          }
          j++;
        }
      }
      free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
    }
  }

  it=1; err= ferr= 1;

  fnegSkewnorm(fmode,ypred,thmode,sel,nsel,n,y,x,XtX,tau,taualpha,alphaphi,lambdaphi,prior,true,symmetric);
  (*fmode) -= thmode[*nsel +1];

  while ((err> *thtol) & (it<(*maxit)) & (ferr> *ftol)) {

    err= ferr= 0; sumth2= 0;
    for (j=1; j<=p; j++) {

      if (j== *nsel +1) {  //update for phi
        if (*prior ==1) { //under MOM use exact max for phi
          for (i=0, s1=0, s2=0; i<(*n); i++) { if (y[i]<=ypred[i]) { s1+= pow(y[i]-ypred[i], 2.0); } else { s2+= pow(y[i]-ypred[i], 2.0); } }
          for (i=1, sumth2=0; i<=(*nsel); i++) { sumth2 += thnew[i]*thnew[i]; }
          if (*symmetric ==0) { acur= tanh(thnew[p]); } else { acur= 0; }
          aa= (*n + 3*(*nsel) + (*alphaphi));
          bb= (s1/pow(1+acur,2.0) + s2/pow(1-acur,2.0) + sumth2/(*tau) + *lambdaphi);
          thnew[j]= log(bb/aa);
        } else {
          fpnegSkewnormUniv(j,&g,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //gradient
          g -= 1.0; //jacobian term from tvartheta=log(tvartheta)
          fppnegSkewnormUniv(j,&H,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //Hessian
          delta= g/H;
          thnew[j]= thmode[j] - delta;
        }
      } else {  //update for theta, alpha
        fpnegSkewnormUniv(j,&g,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //gradient
        fppnegSkewnormUniv(j,&H,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //Hessian
        delta= g/H;
        thnew[j]= thmode[j] - delta;
      }

      fnegSkewnorm(&fnew,ypred,thnew,sel,nsel,n,y,x,XtX,tau,taualpha,alphaphi,lambdaphi,prior,true,symmetric);
      fnew -= thnew[*nsel +1];
      //If new value improves target function, update thmode, fmode
      if (fnew<(*fmode)) {
        err= max_xy(err,fabs(thmode[j]-thnew[j]));
        thmode[j]= thnew[j];
        ferr+= *fmode - fnew;
        (*fmode)= fnew;
      } else {
        thnew[j]= thmode[j];
      }

    }

    it++;

  }

  fppnegSkewnorm(hess,thmode,ypred,sel,nsel,n,y,x,tau,taualpha,alphaphi,lambdaphi,prior,symmetric); //Hessian

  thmode[*nsel +1]= exp(thmode[*nsel +1]);
  if (*symmetric ==0) { thmode[p]= tanh(thmode[p]); } //Note: tanh(z)= -1 + 2/(1+exp(-2*z))
  //Rprintf("--- niter=%d, Posterior mode= %f %f %f \n",it,thmode[1],thmode[2],thmode[p]);

  free_dvector(ypred, 0,*n -1); free_dvector(thnew,1,p);
}



void fnegSkewnorm(double *ans, double *ypred, double *th, int *sel, int *nsel, int *n, double *y, double *x, crossprodmat *XtX, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, bool logscale, int *symmetric) {
//Negative log-joint for two-piece Normal under MOM/eMOM/iMOM prior on coef and IG on variance
//Note: log-joint evaluated for vartheta, if log-joint for log(vartheta) is desired you need to substract -th[nsel+1] to consider the Jacobian term
// Input
// - th[1..nsel+2]: (theta, log(vartheta), atanh(alpha)) where theta=regression coef, vartheta \propto variance and alpha=asymmetry parameter in [-1,1]
// - Other parameters as in postmodeSkewNorm
// Output:
// - ans: value of minus the log-joint evaluated at th
// - ypred: linear predictor x %*% th
  double scale, alpha;

  scale= exp(th[*nsel +1]);
  if (*symmetric ==0) { alpha= tanh(th[*nsel +2]); } else { alpha= 0; }
  loglSkewnorm(ans, ypred, th, nsel, sel, n, &scale, &alpha, y, x, XtX);
  (*ans)= -(*ans);

  if ((*prior)==1) {

    if ((*nsel)>0) {
      (*ans) += -dmomvec(th+1,*nsel,0.0,*tau,scale,1,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }
    if (*symmetric ==0) (*ans)-= dmom(th[*nsel +2],0.0,*taualpha,1.0,1,1);

  } else if ((*prior)==2) {

    if ((*nsel)>0) {
      (*ans) += -dimomvec(th+1,*nsel,0.0,*tau,scale,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }
    if (*symmetric ==0) (*ans)-= dimom(th[*nsel +2],0.0,*taualpha,1.0,1);

  } else if ((*prior)==3) {

    if ((*nsel)>0) {
      (*ans) += -demomvec(th+1,*nsel,*tau,scale,1) - dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    } else {
      (*ans) += -dinvgammaC(scale,0.5*(*alphaphi),0.5*(*lambdaphi),1);
    }
    if (*symmetric ==0) (*ans)-= demom(th[*nsel +2],*taualpha,1.0,1);

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  if (!logscale) { (*ans)= exp(*ans); }
}



void fpnegSkewnorm(double *g, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior) {
 //Gradient of fnegSkewnorm
 //Note: g[*nsel +1] does not include the -1.0 term coming from the jacobian of tvartheta=log(vartheta)
  int i, one=1, nselplus1= (*nsel)+1;
  double *gprior, zero=0;

  gprior= dvector(1,(*nsel)+2);

  loglnegGradSkewNorm(g,th,nsel,sel,n,y,ypred,x);

  if ((*prior)==1) {

    dmomiggrad(gprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    dmomgrad(gprior+(*nsel)+1,&one,th+(*nsel)+1,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else if ((*prior)==2) {

    dimomiggrad(gprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    dimomgrad(gprior+(*nsel)+1,&one,th+(*nsel)+1,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else if ((*prior)==3) {

    demomiggrad(gprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) { g[i] -= gprior[i]; }

    demomgrad(gprior+(*nsel)+1,&one,th+(*nsel)+1,&zero,taualpha);
    g[(*nsel)+2] -= gprior[(*nsel)+2];

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  free_dvector(gprior,1,(*nsel)+2);
}


void fpnegSkewnormUniv(int j, double *g, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric) {
  //Gradient of fnegSkewnorm
  //Note: the -1.0 term when j=nsel+1 coming from the jacobian of tvartheta=log(vartheta) is not included, it must be added separately after calling this function
  int i;
  double gprior, zero=0, sumth2, suminvth2;

  loglnegGradSkewNormUniv(j,g,th,nsel,sel,n,y,ypred,x,symmetric);

  if ((*prior)==1) {

    if (j <= (*nsel)) {
      gprior= dmomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0; i<=(*nsel); i++) { sumth2 += th[i]*th[i]; }
      gprior= -1.5*(*nsel) - 0.5*(*alphaphi) -1 + 0.5*(sumth2/(*tau) + *lambdaphi) * exp(-th[*nsel +1]);
    } else {
      gprior= dmomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }
    (*g) -= gprior;

  } else if ((*prior)==2) {

    if (j <= (*nsel)) {
      gprior= dimomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, suminvth2=0; i<=(*nsel); i++) { suminvth2 += 1.0/(th[i]*th[i]); }
      gprior= 0.5*(*nsel) - 0.5*(*alphaphi) -1 + 0.5*(*lambdaphi)*exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      gprior= dimomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }
    (*g) -= gprior;

  } else if ((*prior)==3) {

    if (j <= (*nsel)) {
      gprior= demomgraduniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0, suminvth2=0; i<=(*nsel); i++) { sumth2 += th[i]*th[i]; suminvth2 += 1.0/(th[i]*th[i]); }
      gprior= -0.5*(*nsel) - 0.5*(*alphaphi) - 1 +0.5*(sumth2/(*tau) + *lambdaphi) * exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      gprior= demomgraduniv(th+(*nsel)+2, &zero, taualpha);
    }
    (*g) -= gprior;

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

}


void fppnegSkewnorm(double **H, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric) {
 //Hessian of fnegSkewnorm

  int i, j, one=1, nselplus1= (*nsel)+1;
  double **Hprior, *hprioralpha, zero=0;

  Hprior= dmatrix(1,nselplus1,1,nselplus1);
  hprioralpha= dvector(1,1);

  loglnegHessSkewNorm(H,th,nsel,sel,n,y,ypred,x,symmetric);

  if ((*prior)==1) {

    dmomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) { H[i][j]= H[j][i]= H[i][j] - Hprior[i][j]; }
    }

    if (*symmetric ==0) {
      dmomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else if ((*prior)==2) {

    dimomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) { H[i][j]= H[j][i]= H[i][j] - Hprior[i][j]; }
    }

    if (*symmetric ==0) {
      dimomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else if ((*prior)==3) {

    demomighess(Hprior,&nselplus1,th,th+(*nsel)+1,tau,alphaphi,lambdaphi);
    for (i=1; i<= (*nsel)+1; i++) {
      H[i][i] -= Hprior[i][i];
      for (j=1; j<i; j++) { H[i][j]= H[j][i]= H[i][j] - Hprior[i][j]; }
    }

    if (*symmetric ==0) {
      demomhess(hprioralpha,&one,th+(*nsel)+1,&zero,taualpha);
      H[(*nsel)+2][(*nsel)+2] -= hprioralpha[1];
    }

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

  free_dmatrix(Hprior,1,nselplus1,1,nselplus1);
  free_dvector(hprioralpha,1,1);
}


void fppnegSkewnormUniv(int j, double *H, double *th, double *ypred, int *sel, int *nsel, int *n, double *y, double *x, double *tau, double *taualpha, double *alphaphi, double *lambdaphi, int *prior, int *symmetric) {
 //Hessian of fnegSkewnorm

  int i;
  double hprior, zero=0, sumth2, suminvth2;

  loglnegHessSkewNormUniv(j,H,th,nsel,sel,n,y,ypred,x,symmetric);

  if ((*prior)==1) {

    if (j<=(*nsel)) {
      hprior= dmomhessuniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0;i<=(*nsel);i++) sumth2 += pow(th[i],2.0);
      hprior= -0.5 * exp(-th[(*nsel)+1]) * (sumth2/(*tau)+(*lambdaphi));
    } else {
      hprior= dmomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else if ((*prior)==2) {

    if (j<=(*nsel)) {
      hprior= dimomhessuniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, suminvth2=0; i<=(*nsel); i++) suminvth2 += pow(1.0/th[i],2.0);
      hprior= -0.5*exp(-th[(*nsel)+1])*(*lambdaphi) - (*tau) * exp(th[(*nsel)+1]) * suminvth2;
    } else {
      hprior= dimomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else if ((*prior)==3) {

    if (j<=(*nsel)) {
      hprior= demomhessuniv(th+j, th+(*nsel)+1, tau);
    } else if (j== (*nsel)+1) {
      for (i=1, sumth2=0, suminvth2=0; i<=(*nsel); i++) { sumth2+= pow(th[i],2.0); suminvth2 += pow(1.0/th[i],2.0); }
      hprior= -0.5*(*nsel) - 0.5*(*alphaphi) -1.0 + 0.5*(sumth2/(*tau) + (*lambdaphi)) * exp(-th[*nsel +1]) - exp(th[*nsel +1])*(*tau)*suminvth2;
    } else {
      hprior= demomhessuniv(th+(*nsel)+2,&zero,taualpha);
    }
    (*H) -= hprior;

  } else {

    Rf_error("prior must be 'mom', 'imom' or 'emom'");

  }

}



void loglSkewnorm(double *ans, double *ypred, double *th, int *nsel, int *sel, int *n, double *scale, double *alpha, double *y, double *x, crossprodmat *XtX) {
  //Log-likelihood function of a linear model with two-piece normal errors evaluated at th=(theta,scale,alpha)
  //Output
  // - ans: value of the log-likelihood evaluated at th
  // - ypred: linear predictor x %*% th
  int i;
  double w1, w2;

  w1= 0.5 / (pow(1.0 + (*alpha),2) * (*scale));
  w2= 0.5 / (pow(1.0 - (*alpha),2) * (*scale));
  (*ans)= -0.5*(*n)*(LOG_M_2PI + log(*scale));

  if ((*nsel)>0) {

    Aselvecx(x, th+1, ypred, 0, (*n) -1, sel, nsel); //ypred= x %*% th

    for (i=0; i<(*n); i++) {

      if (y[i]<ypred[i]) { (*ans) -= w1 * pow(y[i]-ypred[i],2); } else { (*ans) -= w2 * pow(y[i]-ypred[i],2); }

    }

  } else {

    for (i=0; i<(*n); i++) {

      if (y[i]<0) { (*ans) -= w1 * pow(y[i],2); } else { (*ans) -= w2 * pow(y[i],2); }

    }

  }

}


void loglnegGradSkewNorm(double *g, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x) {
  //Gradient of minus the log-likelihood function of a linear model with two-piece Normal errors
  int i;
  double sigma, alpha, alphat, w1, w2, ws1, ws2, *y0, *Wy0, y0Wy0=0, y0Wsy0=0;

  Wy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  alphat= tanh(th[*nsel +2]);
  alpha= th[*nsel+2];

  w1= 1.0 / (pow(1.0 + alphat,2));
  w2= 1.0 / (pow(1.0 - alphat,2));
  ws1= -2.0 / (pow(cosh(alpha),2) * pow(1.0+alphat,3));
  ws2=  2.0 / (pow(cosh(alpha),2) * pow(1.0-alphat,3));

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);

      if (y[i]<ypred[i]) {
        Wy0[i]= w1 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws1;
      } else {
        Wy0[i]= w2 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws2;
      }
      y0Wy0+= y0[i] * Wy0[i];

    }

    Atselvecx(x, Wy0, g+1, 0, (*n)-1, sel, nsel);  //g[1:nsel]= t(x) %*% Wy0
    for (i=1; i<=(*nsel); i++) g[i]= -g[i]/sigma;
    free_dvector(y0,0,*n -1);

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
        Wy0[i]= w1 * y[i];
        y0Wsy0+= pow(y[i],2)*ws1;
      } else {
        Wy0[i]= w2 * y[i];
        y0Wsy0+= pow(y[i],2)*ws2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

  }

  g[*nsel +1]= 0.5*(*n) - 0.5*y0Wy0/sigma;

  g[*nsel +2]= 0.5*y0Wsy0/sigma;

  free_dvector(Wy0,0,*n -1);
}


void loglnegGradSkewNormUniv(int j, double *g, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric) {
  //Gradient for th[j] of minus the log-likelihood function of a linear model with two-piece Normal errors
  int i;
  double sigma, alpha, alphat, w1, w2, ws1, ws2, *y0, *Wy0, y0Wy0=0, y0Wsy0=0;

  Wy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  if (*symmetric ==0) { alphat= tanh(th[*nsel +2]); alpha= th[*nsel+2]; } else { alphat= alpha= 0; }

  w1= 1.0 / (pow(1.0 + alphat,2));
  w2= 1.0 / (pow(1.0 - alphat,2));
  ws1= -2.0 / (pow(cosh(alpha),2) * pow(1.0+alphat,3));
  ws2=  2.0 / (pow(cosh(alpha),2) * pow(1.0-alphat,3));

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);

      if (y[i]<ypred[i]) {
        Wy0[i]= w1 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws1;
      } else {
        Wy0[i]= w2 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws2;
      }
      y0Wy0+= y0[i] * Wy0[i];

    }

    if (j<=(*nsel)) {
      int selj= sel[j-1], nselj= 1;
      Atselvecx(x, Wy0, g, 0, (*n)-1, &selj, &nselj);  //j^th elem in t(x) %*% Wy0
      (*g)= -(*g)/sigma;
    }

    free_dvector(y0,0,*n -1);

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
        Wy0[i]= w1 * y[i];
        y0Wsy0+= pow(y[i],2)*ws1;
      } else {
        Wy0[i]= w2 * y[i];
        y0Wsy0+= pow(y[i],2)*ws2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

  }

  if (j== *nsel +1) {
    (*g)= 0.5*(*n) - 0.5*y0Wy0/sigma;
  } else if (j== *nsel +2) {
    (*g)= 0.5*y0Wsy0/sigma;
  }

  free_dvector(Wy0,0,*n -1);
}



void loglnegHessSkewNorm(double **H, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric) {
  //Hessian of minus the log-likelihood function of a linear model with two-piece Normal errors
  //NOTE: log.lik.hess function in R, but we need to change the sign of the output
  int i, j, k, idxi, idxj;
  double sigma, alphat, alpha, w, w1, w2, ws1, ws2, wss1, wss2, *y0, *Wy0, *Wsy0, y0Wy0=0, y0Wsy0=0, y0Wssy0=0, *tX0Wy0;

  Wy0= dvector(0,*n -1); Wsy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  if (*symmetric ==0) { alphat= tanh(th[*nsel +2]); alpha= th[*nsel+2]; } else { alphat= alpha= 0; }

  w1= 1.0 / (pow(1.0 + alphat,2.0));
  w2= 1.0 / (pow(1.0 - alphat,2.0));
  ws1= -2.0 / (pow(cosh(alpha),2.0) * pow(1.0+alphat,3.0));
  ws2=  2.0 / (pow(cosh(alpha),2.0) * pow(1.0-alphat,3.0));
  wss1= 2.0 * exp(-2.0*alpha) + 4.0*exp(-4.0*alpha);
  wss2= 2.0 * exp(2.0*alpha) + 4.0*exp(4.0*alpha);

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);

      if (y[i]<ypred[i]) {
        Wy0[i]= w1 * y0[i];
        Wsy0[i]= ws1 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws1;
        y0Wssy0+= pow(y0[i],2)*wss1;
      } else {
        Wy0[i]= w2 * y0[i];
        Wsy0[i]= ws2 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws2;
        y0Wssy0+= pow(y0[i],2)*wss2;
      }
      y0Wy0+= y0[i] * Wy0[i];

    }

    free_dvector(y0,0,*n -1);

    //Compute H[1:*nsel,1:*nsel] <- t(X0)%*%W%*%X0/sigma
    for (i=1; i<=(*nsel); i++) {
      idxi= (*n)*sel[i-1];
      for (j=i; j<=(*nsel); j++) {
        idxj= (*n)*sel[j-1];
        H[i][j]= 0;
        for (k=0; k<(*n); k++) {
          if (y[k]<ypred[k]) { w= w1; } else { w= w2; }
          H[i][j] += x[k+ idxi] * x[k +idxj] * w;  //x[k][i] * x[k][j] * w
        }
        H[i][j] /= sigma;
        H[j][i]= H[i][j];
      }
    }

    //Compute H[1:*nsel,*nsel+1] <- t(X0)%*%Wy0/sigma
    tX0Wy0= dvector(1,*nsel);
    Atselvecx(x, Wy0, tX0Wy0+1, 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +1]= tX0Wy0[i]/sigma; H[*nsel +1][i]= H[i][*nsel +1]; }

    //Compute H[1:*nsel,*nsel+2] <- -t(X0)%*%Wsy0/sigma
    if (*symmetric ==0) {
      Atselvecx(x, Wsy0, tX0Wy0+1, 0, (*n)-1, sel, nsel);
      for (i=1; i<=(*nsel); i++) { H[i][*nsel +2]= -tX0Wy0[i]/sigma; H[*nsel +2][i]= H[i][*nsel +2]; }
    }
    free_dvector(tX0Wy0, 1,*nsel);

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
        Wy0[i]= w1 * y[i];
        Wsy0[i]= ws1 * y[i];
        y0Wsy0+= pow(y[i],2)*ws1;
        y0Wssy0+= pow(y[i],2)*wss1;
      } else {
        Wy0[i]= w2 * y[i];
        Wsy0[i]= ws2 * y[i];
        y0Wsy0+= pow(y[i],2)*ws2;
        y0Wssy0+= pow(y[i],2)*wss2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

    //Compute H[1:*nsel,*nsel+1] <- t(X0)%*%Wy0/sigma
    Atselvecx(x, Wy0, (double *)H +(*nsel)*(*nsel), 0, (*n)-1, sel, nsel);
    for (i=1; i<=(*nsel); i++) { H[i][*nsel +1]= H[i][*nsel +1]/sigma; H[*nsel +1][i]= H[i][*nsel +1]; }

    //Compute H[1:*nsel,*nsel+2] <- -t(X0)%*%Wsy0/sigma
    if (*symmetric ==0) {
      Atselvecx(x, Wsy0, (double *)H +(*nsel)*(*nsel +1), 0, (*n)-1, sel, nsel);
      for (i=1; i<=(*nsel); i++) { H[i][*nsel +2]= -H[i][*nsel +2]/sigma; H[*nsel +2][i]= H[i][*nsel +2]; }
    }

  }

  H[*nsel +1][*nsel +1]= 0.5 * y0Wy0 / sigma;
  if (*symmetric ==0) {
    H[*nsel +2][*nsel +2]= 0.5 * y0Wssy0 / sigma;
    H[*nsel +1][*nsel +2]= H[*nsel +2][*nsel +1] = -0.5 * y0Wsy0 / sigma;
  }

  free_dvector(Wy0, 0,*n -1); free_dvector(Wsy0, 0,*n -1);
}


void loglnegHessSkewNormUniv(int jj, double *H, double *th, int *nsel, int *sel, int *n, double *y, double *ypred, double *x, int *symmetric) {
  //Hessian for th[jj] of minus the log-likelihood function of a linear model with two-piece Normal errors
  //NOTE: log.lik.hess function in R, but we need to change the sign of the output
  int i, k, idxi;
  double sigma, alphat, alpha, w, w1, w2, ws1, ws2, wss1, wss2, *y0, *Wy0, *Wsy0, y0Wy0=0, y0Wsy0=0, y0Wssy0=0;

  Wy0= dvector(0,*n -1); Wsy0= dvector(0,*n -1);
  sigma= exp(th[*nsel +1]);
  if (*symmetric ==0) { alphat= tanh(th[*nsel +2]); alpha= th[*nsel+2]; } else { alphat= alpha= 0; }

  w1= 1.0 / (pow(1.0 + alphat,2.0));
  w2= 1.0 / (pow(1.0 - alphat,2.0));
  ws1= -2.0 / (pow(cosh(alpha),2.0) * pow(1.0+alphat,3.0));
  ws2=  2.0 / (pow(cosh(alpha),2.0) * pow(1.0-alphat,3.0));
  wss1= 2.0 * exp(-2.0*alpha) + 4.0*exp(-4.0*alpha);
  wss2= 2.0 * exp(2.0*alpha) + 4.0*exp(4.0*alpha);

  if ((*nsel)>0) {

    y0= dvector(0,*n -1);

    for (i=0; i<(*n); i++) {
      y0[i]= (y[i]-ypred[i]);
      if (y[i]<ypred[i]) {
        Wy0[i]= w1 * y0[i];
        Wsy0[i]= ws1 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws1;
        y0Wssy0+= pow(y0[i],2)*wss1;
      } else {
        Wy0[i]= w2 * y0[i];
        Wsy0[i]= ws2 * y0[i];
        y0Wsy0+= pow(y0[i],2)*ws2;
        y0Wssy0+= pow(y0[i],2)*wss2;
      }
      y0Wy0+= y0[i] * Wy0[i];
    }

    free_dvector(y0,0,*n -1);

    if (jj<= *nsel) {
      idxi= (*n)*sel[jj-1];
      (*H)= 0;
      for (k=0; k<(*n); k++) {
        if (y[k]<ypred[k]) { w= w1; } else { w= w2; }
        (*H) += pow(x[k+ idxi],2.0) * w;
      }
      (*H) /= sigma;
    }

  } else {

    for (i=0; i<(*n); i++) {
      if (y[i]<0) {
        Wy0[i]= w1 * y[i];
        Wsy0[i]= ws1 * y[i];
        y0Wsy0+= pow(y[i],2)*ws1;
        y0Wssy0+= pow(y[i],2)*wss1;
      } else {
        Wy0[i]= w2 * y[i];
        Wsy0[i]= ws2 * y[i];
        y0Wsy0+= pow(y[i],2)*ws2;
        y0Wssy0+= pow(y[i],2)*wss2;
      }
      y0Wy0+= y[i] * Wy0[i];
    }

  }

  if (jj== *nsel +1) {
    (*H)= 0.5 * y0Wy0 / sigma;
  } else if (jj== *nsel +2) {
    (*H)= 0.5 * y0Wssy0 / sigma;
  }

  free_dvector(Wy0, 0,*n -1); free_dvector(Wsy0, 0,*n -1);
}


void mleSkewnorm(double *thmode, double *ypred, int *sel, int *nsel, int *n, int *p, double *y, double *x, crossprodmat *XtX, double *ytX, int *maxit, bool useinit) {
  //Find MLE for linear regression with skew-normal residuals
  //Output:
  // - thmode[1..nsel+2] contains MLE for regression parameters, residual dispersion and asymmetry
  // - ypred[0.. n-1] contains linear predictor at MLE
  //If useinit==false thmode is initialized at least squares, else thmode is used
  bool posdef;
  int i, ii=1, j, k, idxi, idxj, nseluniv=1, seluniv=sel[0], maxituniv=100;
  double *Xtwy, **XtwX, **XtwXinv, *thnew, err=1.0, s1=0, s2=0, s1pow=1, s2pow=1, w, w1, w2, *e, *epred, *thuniv, difth1;

  if ((*nsel)>0) {  //There are covariates

    Xtwy= dvector(1,*nsel);
    XtwX= dmatrix(1,*nsel,1,*nsel);
    XtwXinv= dmatrix(1,*nsel,1,*nsel);
    thnew= dvector(1,*nsel);

    if (!useinit) {

      leastsquares(thmode,thmode+(*nsel)+1,ypred,y,x,XtX,ytX,n,p,sel,nsel);

      //Refine estimate for intercept (assumed to be sel[0]) via several univariate updates
      e= dvector(0,*n -1); epred= dvector(0,*n -1); thuniv= dvector(1,3);
      idxj= (*n)*sel[0];
      for (i=0; i< (*n); i++) { epred[i]= thmode[1]*x[i +idxj]; e[i]= y[i] - ypred[i] + epred[i]; }
      thuniv[1]= thmode[1]; thuniv[2]= thmode[*nsel +1]; thuniv[3]= thmode[*nsel +2];
      mleSkewnorm(thuniv,epred,&seluniv,&nseluniv,n,p,e,x,XtX,ytX,&maxituniv,true);
      difth1= thuniv[1]-thmode[1];
      thmode[1]= thuniv[1]; thmode[*nsel +1]= thuniv[2]; thmode[*nsel +2]= thuniv[3];
      for (i=0; i< (*n); i++) { ypred[i]+= difth1*x[i +idxj]; }  //update linear predictor
      if ((*nsel)==1) { //only 1 covariate, we're done
        err= -1;
        for (i=0, s1=0, s2=0; i< *n; i++) { if (y[i]<=ypred[i]) { s1 += pow(y[i]-ypred[i],2.0); } else { s2 += pow(y[i]-ypred[i],2.0); } }
        s1pow= pow(s1,1.0/3.0); s2pow= pow(s2,1.0/3.0);
        thmode[*nsel +2]= (s1pow-s2pow)/(s1pow+s2pow); //alpha
      }
      free_dvector(e, 0,*n -1); free_dvector(epred, 0,*n -1); free_dvector(thuniv,1,3);
    }

    while ((err>0.0001) && (ii<(*maxit))) {

      for (i=0, s1=0, s2=0; i< *n; i++) { if (y[i]<=ypred[i]) { s1 += pow(y[i]-ypred[i],2.0); } else { s2 += pow(y[i]-ypred[i],2.0); } }

      s1pow= pow(s1,1.0/3.0); s2pow= pow(s2,1.0/3.0);
      thmode[*nsel +2]= (s1pow-s2pow)/(s1pow+s2pow); //alpha

      w1= 1.0/pow(1.0+thmode[*nsel +2],2.0); w2= 1.0/pow(1.0-thmode[*nsel +2],2.0);

      //Compute t(X) %*% W %*% y, t(X) %*% W %*% X
      for (i=1; i<=(*nsel); i++) {
        idxi= (*n)*sel[i-1];
        //Find t(X) %*% W %*% y
        Xtwy[i]= 0;
        for (k=0; k<(*n); k++) {
          if (y[k]<ypred[k]) { w= w1; } else { w= w2; }
          Xtwy[i] += x[k+ idxi] * y[k] * w;  //x[k][i] * y[k] * w
        }
        //Find t(X) %*% W %*% X
        for (j=i; j<=(*nsel); j++) {
          idxj= (*n)*sel[j-1];
          XtwX[i][j]= 0;
          for (k=0; k<(*n); k++) {
            if (y[k]<ypred[k]) { w= w1; } else { w= w2; }
            XtwX[i][j] += x[k+ idxi] * x[k +idxj] * w;  //x[k][i] * x[k][j] * w
          }
        }
        for (j=1; j<i; j++) { XtwX[i][j]= XtwX[j][i]; }
      }

      inv_posdef(XtwX, *nsel, XtwXinv, &posdef);
      Ax(XtwXinv, Xtwy,thnew,1,*nsel,1,*nsel);  //thmode[1..nsel]= (t(X) * W * X)^{-1} t(X) * W * y

      err= fabs(thnew[1]-thmode[1]);
      thmode[1]= thnew[1];
      for (j=2; j<= (*nsel); j++) {
        err= max_xy(err,fabs(thnew[j]-thmode[j]));
        thmode[j]= thnew[j];
      }
      ii++;
      if ((err>0.0001) && (ii<(*maxit))) {
        Aselvecx(x, thmode+1, ypred, 0, (*n) -1, sel, nsel);
      }
    }

    free_dvector(Xtwy,1,*nsel);
    free_dmatrix(XtwX,1,*nsel,1,*nsel);
    free_dmatrix(XtwXinv,1,*nsel,1,*nsel);
    free_dvector(thnew,1,*nsel);

  } else {  //No covariates

    for (i=0; i< *n; i++) {
      if (y[i]<=0) { s1 += pow(y[i],2.0); } else { s2 += pow(y[i],2.0); }
    }

    s1pow= pow(s1,1.0/3.0); s2pow= pow(s2,1.0/3.0);
    thmode[*nsel +2]= (s1pow-s2pow)/(s1pow+s2pow); //alpha

  }

  if (err>=0) thmode[(*nsel)+1]= (0.25/((*n)+.0)) * pow(s1pow + s2pow, 3);  //MLE for residual dispersion

}



//*************************************************************************************
// PRODUCT MOM ROUTINES
//*************************************************************************************

double f2opt_mom(double *th) {
  return fmomNegC_non0(th+1,f2opt_pars.m+1,f2opt_pars.S,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.r,f2opt_pars.n,f2opt_pars.nsel);
}

//Note: th and m are assumed to be indexed at 0; S indexed at 1
double fmomNegC_non0(double *th, double *m, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i;
  double ans, sumlogth, *z;
  z= dvector(0,*nsel);
  for (i=0, sumlogth=0; i<(*nsel); i++) { sumlogth+= log(th[i]*th[i]); z[i]= th[i]-m[i]; }
  ans= .5*quadratic_xtAx(z-1,S,1,*nsel)/(*phi) - (*r +.0)*sumlogth;
  free_dvector(z,0,*nsel);
  return ans;
}

void fppmomNegC_non0(double **ans, double *th, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i, j;
  for (i=1; i<=(*nsel); i++) { ans[i][i]= S[i][i]/(*phi) + 2.0*(*r)/(th[i]*th[i]); }
  for (i=1; i<=(*nsel); i++) { for (j=i+1; j<=(*nsel); j++) { ans[i][j]= ans[j][i]= S[i][j]/(*phi); } }
}

void momIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *n, int *nsel, double *m, double **S, double *detS, double *phi, double *tau, int *r, int *logscale) {
  int i, emptyint, iter, maxit=100;
  double emptydouble, ftol= 1.0e-5, **Vopt, detVopt, **dirth;
  crossprodmat *XtX;

  XtX= new crossprodmat(&emptydouble, 0, 0, true);
  Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  set_f2opt_pars(m,S,&emptydouble,XtX,&emptydouble,&emptydouble,&emptydouble,phi,tau,r,n,nsel,&emptyint,nsel);

  //Minimization
  for (i=1; i<=(*nsel); i++) { thopt[i]= m[i]; }  //init
  ddiag(dirth,1,*nsel);
  minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_mom, maxit);

  //Laplace approx
  fppmomNegC_non0(Vopt,thopt,S,phi,tau,r,n,nsel);
  invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);

  (*ILaplace)= -(*fopt) + .5*(log(*detS)-log(detVopt)- (*nsel)*log(*phi)) ;

  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
  delete XtX; free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
}


//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ N(m,Sinv)
double MC_mom_normal(double *m,double **Sinv,int *r,int *nsel, int *B) {
  bool posdef;
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv,&posdef); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvnormC(thsim,*nsel,m,cholSinv);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac);
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}

//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ T_nu(m,Sinv)
double MC_mom_T(double *m,double **Sinv,int *nu,int *r,int *nsel, int *B) {
  bool posdef;
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv,&posdef); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvtC(thsim,*nsel,m,cholSinv,*nu);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac);
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}


// PRODUCT MOMENT MARGINALS
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - r: MOM power parameter
// - method==0 for Laplace; method==1 for Monte Carlo; method==2 for plug-in; method== -1 for exact calculation if p<20
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

// [[Rcpp::export]]
SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Sngroups, SEXP Snvaringroup) {
  struct marginalPars pars;
  int SoptimMethod= 1, Soptim_maxit= -1, emptyint=1, usethinit=0, knownphi=1;
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  crossprodmat *XtX;
  SEXP ans;

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);

  set_marginalPars(&pars,&emptyint,INTEGER(Sn),INTEGER(Sn),INTEGER(Sp),REAL(Sy),&emptyint,REAL(Ssumy2),&emptydouble,&emptydouble,&emptydouble,&emptydouble,XtX,REAL(SytX),INTEGER(Smethod),&emptyint,&emptyint,&SoptimMethod,&Soptim_maxit,&usethinit,&emptydouble,INTEGER(SB),&emptydouble,&emptydouble,&knownphi,REAL(Sphi),REAL(Stau),&emptydouble,taualpha,taualpha,INTEGER(Sr),NULL,NULL,NULL,NULL,INTEGER(Slogscale),&offset,NULL,NULL,INTEGER(Sngroups),NULL,INTEGER(Snvaringroup));
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalKC(INTEGER(Ssel),INTEGER(Snsel),&pars);

  delete XtX;
  UNPROTECT(1);
  return ans;
}


// Function to compute r * sum(log(th^2))
double rsumlogsq(double *th, int *r, int *nsel) {
  int i; double ans;
  for (i=1, ans=0; i<=(*nsel); i++) { ans+= log(th[i]*th[i]); }
  ans*= (*r);
  return(ans);
}

double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int i,j;
  double *m, s, **S, **Sinv, detS, num, den, logtau= log(*(*pars).tau), tauinv= 1.0/(*(*pars).tau), logphi= log(*(*pars).phi), ans=0.0, *thopt, **Voptinv, fopt;

  if (*nsel ==0) {
    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);
  } else {

    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*logtau) + (*nsel)*(*(*pars).r)*(logtau+logphi+ldoublefact(2*(*(*pars).r)-1));

    if ((*(*pars).method ==0) | ((*(*pars).method == -1) & ((*nsel)>10)))  { //Laplace

      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,(*pars).phi,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    } else if (*(*pars).method ==1) { //MC

      for (i=1; i<=(*nsel); i++) { Sinv[i][i]= (*(*pars).phi)*Sinv[i][i]; for (j=i+1; j<=(*nsel); j++) { Sinv[i][j]=Sinv[j][i]= (*(*pars).phi)*Sinv[i][j]; } }
      ans= MC_mom_normal(m,Sinv,(*pars).r,nsel,(*pars).B);

    } else if (*(*pars).method ==2) { //Plug-in

      ans= rsumlogsq(m,(*pars).r,nsel);

    } else if ((*(*pars).method == -1) & ((*nsel)<=10)) { //Exact

      Voptinv= dmatrix(1,*nsel,1,*nsel);
      for (i=1; i<= *nsel; i++) for (j=i; j<= *nsel; j++) Voptinv[i][j]= Voptinv[j][i]= Sinv[i][j] * (*(*pars).phi);
      ans= log(mvtexpect(m, Voptinv, *nsel, 2, -1));
      free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    }
    ans+= num - den;
    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


// [[Rcpp::export]]
SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP Sngroups, SEXP Snvaringroup) {
  int SoptimMethod= 1, Soptim_maxit= -1, emptyint=1, usethinit=0, knownphi=0;
  double *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  crossprodmat *XtX;
  SEXP ans;

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);

  set_marginalPars(&pars,&emptyint,INTEGER(Sn),INTEGER(Sn),INTEGER(Sp),REAL(Sy),&emptyint,REAL(Ssumy2),&emptydouble,&emptydouble,REAL(Sx),&emptydouble,XtX,REAL(SytX),INTEGER(Smethod),&emptyint,&emptyint,&SoptimMethod,&Soptim_maxit,&usethinit,&emptydouble,INTEGER(SB),REAL(Salpha),REAL(Slambda),&knownphi,&emptydouble,REAL(Stau),&emptydouble,taualpha,taualpha,INTEGER(Sr),NULL,NULL,NULL,NULL,INTEGER(Slogscale),&offset,NULL,NULL,INTEGER(Sngroups),NULL,INTEGER(Snvaringroup));
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalUC(INTEGER(Ssel), INTEGER(Snsel), &pars);

  delete XtX;
  UNPROTECT(1);
  return ans;
}


double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, nu, *isgroup=pars->isgroup, varingroup, singlevarcount=0, maxvarsexact= 3;
  double num, den, ans=0.0, term1, *m, **S, **Sinv, **Voptinv, detS, tau=*(*pars).tau, tauinv= 1.0/tau, taugroup=*(*pars).taugroup, taugroupinv=1.0/taugroup, logtaus, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0;
  bool hasgroups= (*((*pars).ngroups)) < (*((*pars).p));

  if (*nsel == 0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    if (*(*pars).method ==0)  { //Laplace
      if ((abs(tau - taugroup) > 1e-7) and hasgroups) {
        Rprintf("Laplace method does not support different tau values between priorCoef and priorGroup\nUsing only tau from priorCoef");
      }

      int prior=1, symmetric=1;
      ans= nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);

    } else {

      m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
      addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);
      for (i = 1; i <= *nsel; i++ ) {
        varingroup = isgroup[sel[i-1]];
        if (varingroup==0) {
          S[i][i] += tauinv;
          singlevarcount++;
        } else {
          S[i][i] += taugroupinv;
        }
      }
      invdet_posdef(S,*nsel,Sinv,&detS);
      Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
      nuhalf= (*(*pars).r)*(*nsel) + .5*(*(*pars).n + *(*pars).alpha);
      nu= (int) (2.0*nuhalf);

      ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
      num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
      logtaus = (.5 + *(*pars).r)*(singlevarcount*log(tau) + (*nsel - singlevarcount)*log(taugroup));
      den= (*nsel)*ldoublefact(2*(*(*pars).r)-1.0) + .5*(*(*pars).n * LOG_M_2PI + log(detS)) + logtaus + gamln(&alphahalf);

      if (*(*pars).method ==1) {  //MC

        term1= (*(*pars).lambda + *(*pars).sumy2 - quadratic_xseltAxsel((*pars).ytX,Sinv,1,nsel,sel))/(nu+.0);
        for (i=1; i<= *nsel; i++) { for (j=i; j<= *nsel; j++) { Sinv[i][j]= Sinv[j][i]= Sinv[i][j]*term1; } } //Vinv matrix
        ans= MC_mom_T(m,Sinv,&nu,(*pars).r,nsel,(*pars).B);

      } else if ((*(*pars).method ==2) | ((*(*pars).method == -1) & ((*nsel)>maxvarsexact)))  { //Orthogonal approx

        //ans= rsumlogsq(m,(*pars).r,nsel); //old version
        term1= ss / ((double) (nu-2)); // (ss/nu) * nu / (nu-2)
        for (i=1, ans=0; i<=(*nsel); i++) { ans+= log(pow(m[i],2.0) + Sinv[i][i] * term1); }

      } else if ((*(*pars).method == -1) & ((*nsel)<=maxvarsexact)) { //Exact

        Voptinv= dmatrix(1,*nsel,1,*nsel);
        for (i=1; i<= *nsel; i++) for (j=i; j<= *nsel; j++) Voptinv[i][j]= Voptinv[j][i]= Sinv[i][j] * ss / (nu+.0);
        ans= log(mvtexpect(m, Voptinv, *nsel, 2, nu));
        free_dmatrix(Voptinv,1,*nsel,1,*nsel);

      }
      ans+= num - den;
      free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
    }
  }

  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}



double pmomMarginalUC_old(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, nu;
  double num, den, ans=0.0, term1, *m, **S, **Sinv, detS, *thopt, **Voptinv, fopt, phiadj, tauinv= 1.0/(*(*pars).tau), nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss;
  if (*nsel == 0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
  } else {
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= (*(*pars).r)*(*nsel) + .5*(*(*pars).n + *(*pars).alpha);
    nu= (int) (2.0*nuhalf);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= (*nsel)*ldoublefact(2*(*(*pars).r)-1.0) + .5*(*(*pars).n * LOG_M_2PI + log(detS)) + (*nsel)*(.5 + *(*pars).r)*log(*(*pars).tau) + gamln(&alphahalf);

    //Rprintf("%d variables, method=%d \n", *nsel, *((*pars).method));
    if ((*(*pars).method ==0) | ((*(*pars).method == -1) & ((*nsel)>10)))  { //Laplace

      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      phiadj= (nu+.0)/(nu-2.0);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,&phiadj,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    } else if (*(*pars).method ==1) {  //MC

      term1= (*(*pars).lambda + *(*pars).sumy2 - quadratic_xseltAxsel((*pars).ytX,Sinv,1,nsel,sel))/(nu+.0);
      for (i=1; i<= *nsel; i++) { for (j=i; j<= *nsel; j++) { Sinv[i][j]= Sinv[j][i]= Sinv[i][j]*term1; } } //Vinv matrix
      ans= MC_mom_T(m,Sinv,&nu,(*pars).r,nsel,(*pars).B);

    } else if (*(*pars).method ==2) {  //Plug-in

      ans= rsumlogsq(m,(*pars).r,nsel);

    } else if ((*(*pars).method == -1) & ((*nsel)<=10)) { //Exact

      Voptinv= dmatrix(1,*nsel,1,*nsel);
      for (i=1; i<= *nsel; i++) for (j=i; j<= *nsel; j++) Voptinv[i][j]= Voptinv[j][i]= Sinv[i][j] * ss / (nu+.0);
      ans= log(mvtexpect(m, Voptinv, *nsel, 2, nu));
      free_dmatrix(Voptinv,1,*nsel,1,*nsel);

    }
    ans+= num - den;
    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


//********************************************************************************************
// PRODUCT IMOM ROUTINES
//********************************************************************************************

//fimomNeg: minus log integrand of the function needed to compute product iMOM marginal density of the data conditional under a given linear model
//
// fimomNegC
// Input
// - th: theta value at which to evaluate the function (includes coef for all covariates, even those excluded in the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - n: sample size (length of y)
// - p: number of columns in XtX
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// Output: scalar evaluating the function for a single value of theta
double fimomNegC(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[sel[i]];
    th2= th[sel[i]] * th[sel[i]];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xseltAselxsel(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


double f2opt_imom(double *th) {
  double ans;
  ans= fimomNegC_non0(th+1,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomNegC_non0(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
//same as fimomNegC but loops over all elem in th (i.e. th has length *nsel and contains non-zero elements only). th is indexed at 0.
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xtAselx(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


//Hessian of fimomNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)][1:(*nsel)])
// - th: th[1:(*nsel)] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomNegC_non0(double **ans, double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2;

  for (i=1; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= (XtX->at(sel[i-1]*(*p)+sel[i-1])) / (*phi) + 6.0*(*tau)*(*phi)/(th2*th2) - 2.0/th2;
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= (XtX->at(sel[i-1]*(*p)+sel[j-1])) /(*phi);
    }
  }
}


void imomModeK(double *th, PolynomialRootFinder::RootStatus_T *status, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *sel, int *nsel, int *p) {
  //piMOM mode when phi is known using gradient algorithm
  // - th: contains initial estimate at input and mode at output
  // - status: indicates if root finding has been successful
  bool found=false;
  int i, j, niter=0, root_count;
  double err= 1.0, *coef, *real_vector, *imag_vector;
  Polynomial poly;

  coef= dvector(0,4);
  real_vector= dvector(0,4);
  imag_vector= dvector(0,4);

  coef[0]= 2.0 * (*tau) * (*phi);
  coef[1]= 0.0;
  coef[2]= -2;
  while ((err > 1.0e-5) & (niter<50)) {
    err= 0;
    for (i=1; i<=(*nsel); i++) {
      coef[3]= ytX[sel[i-1]];
      for (j=1; j<i; j++) { coef[3]-= (XtX->at(sel[i-1]*(*p)+sel[j-1])) * th[j]; }
      for (j=i+1; j<=(*nsel); j++) { coef[3]-= (XtX->at(sel[i-1]*(*p)+sel[j-1])) * th[j]; }
      coef[3]= coef[3]/(*phi);
      coef[4]= -(XtX->at(sel[i-1]*(*p)+sel[i-1])) / (*phi);
      poly.SetCoefficients(coef, 4);
      (*status)= poly.FindRoots(real_vector,imag_vector,&root_count);

      j=0; found= false;
      while ((!found) & (j<=4)) {
        if (fabs(imag_vector[j])<1.0e-5) {
          if (((real_vector[j]>0) & (th[i]>0)) | ((real_vector[j]<0) & (th[i]<0))) {
            err= max_xy(err,fabs(th[i] - real_vector[j]));
            th[i]= real_vector[j];
            found= true;
          }
        }
        j++;
      }

    }
    niter++;
  }

  free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
}


void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *logscale, int *hessian) {
  bool posdef;
  int iter, maxit=100, emptyint;
  double **V, **Vinv, ftol= 1.0e-5, **dirth, **Vopt, detVopt, emptydouble=0, **emptymatrix;
  PolynomialRootFinder::RootStatus_T status;

  V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel); Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to XtX diagonal, store in V
  inv_posdef_upper(V,*nsel,Vinv,&posdef);
  Asym_xsel(Vinv,*nsel,ytX,sel,thopt);  //product Vinv * selected elements in ytX
  //Minimization
  imomModeK(thopt,&status,XtX,ytX,phi,tau,sel,nsel,p);
  set_f2opt_pars(&emptydouble,emptymatrix,&emptydouble,XtX,ytX,&emptydouble,&emptydouble,phi,tau,&emptyint,n,p,sel,nsel);
  if (status == PolynomialRootFinder::SUCCESS) {
    (*fopt)= f2opt_imom(thopt);
  } else {
    ddiag(dirth,1,*nsel);
    minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_imom, maxit);
  }

  if (*hessian == 1) {
    //Laplace approx
    fppimomNegC_non0(Vopt,thopt,XtX,ytX,phi,tau,n,p,sel,nsel);
    invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);
    (*ILaplace)= -(*fopt) - 0.5*log(detVopt);
  } else {
    (*ILaplace)= -(*fopt) - 0.5*(*nsel)*log(*n +.0);  //BIC-type approximation
  }

  free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
  free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
}


//Product iMOM marginal density for known phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - method: method to approximate the marginal. method==0 for Laplace approximation (may underestimate the true value), method==1 for Importance Sampling Monte Carlo
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

// [[Rcpp::export]]
SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Sngroups, SEXP Snvaringroup) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), SoptimMethod=1, Soptim_maxit= -1, *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1, emptyint=1, usethinit=0, knownphi=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *ytX=REAL(SytX), *phi=REAL(Sphi), *tau=REAL(Stau), *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  crossprodmat *XtX;
  SEXP ans;

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  set_marginalPars(&pars,&emptyint,n,n,p,y,&emptyint,sumy2,&emptydouble,&emptydouble,&emptydouble,&emptydouble,XtX,ytX,method,&emptyint,&emptyint,&SoptimMethod,&Soptim_maxit,&usethinit,&emptydouble,B,&emptydouble,&emptydouble,&knownphi,phi,tau,&emptydouble,taualpha,taualpha,&r,NULL,NULL,NULL,NULL,logscale,&offset,NULL,NULL,INTEGER(Sngroups),NULL,INTEGER(Snvaringroup));
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalKC(sel, nsel, &pars);
  delete XtX;
  UNPROTECT(1);
  return ans;
}


double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int one=1, hessian;
  double k, ans, m, s, ILaplace, *thopt, **Voptinv, fopt;
  thopt= dvector(1,*nsel);
  Voptinv= dmatrix(1,*nsel,1,*nsel);
  if ((*nsel)==0) {
    m= 0;
    s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m,s,1);
  } else {
    if (*(*pars).method == 2) { hessian=0; } else { hessian=1; }
    imomIntegralApproxC(&ILaplace,thopt,Voptinv,&fopt,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,&one,&hessian);
    k= .5*((*nsel)*log(*(*pars).tau) - (*(*pars).sumy2)/(*(*pars).phi) - (*(*pars).n +.0) *LOG_M_2PI - (*(*pars).n - *nsel)*log(*(*pars).phi) - (*nsel)*LOG_M_PI);
    if (((*(*pars).method)==0) || ((*(*pars).method)==2)) {
      ans= k + ILaplace;
    } else {
      ans= k + IS_imom(thopt,Voptinv,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,(*pars).B);
    }
  }
  if ((*(*pars).logscale)!=1) { ans= exp(ans); }
  free_dvector(thopt,1,*nsel);
  free_dmatrix(Voptinv,1,*nsel,1,*nsel);
  return(ans);
}


//Evaluation of iMOM integral via Importance Sampling (result is returned in log-scale)
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *B) {
  bool posdef;
  int i,j;
  double *sdprop, **Vprop, *sopt, **cholVprop, **cholVpropinv, detVpropinv, *mprop, *thsim, *logr, maxlogr, ans;

  sdprop= dvector(1,*nsel); sopt= dvector(1,*nsel);
  mprop= dvector(1,*nsel); thsim= dvector(1, *nsel);
  logr= dvector(0,999);
  Vprop= dmatrix(1,*nsel,1,*nsel); cholVprop= dmatrix(1,*nsel,1,*nsel); cholVpropinv= dmatrix(1,*nsel,1,*nsel);

  for (i=1; i<=(*nsel); i++) {
    mprop[i]= 0;
    sopt[i]= sqrt(Voptinv[i][i]);
    sdprop[i]= .5*fabs(thopt[i] + 2*dsign(thopt[i])*sopt[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i; j<=(*nsel); j++) {
      Vprop[i][j]= Vprop[j][i]= sdprop[i]*sdprop[j]*Voptinv[i][j]/(sopt[i]*sopt[j]);
    }
  }
  choldc(Vprop,*nsel,cholVprop,&posdef);
  choldc_inv(Vprop,*nsel,cholVpropinv,&posdef);
  detVpropinv= choldc_det(cholVpropinv, *nsel);
  rmvtC(thsim, *nsel, mprop, cholVprop, 1);
  maxlogr= logr[0]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
  for (i=1;i<1000;i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    logr[i]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
    if (logr[i]>maxlogr) { maxlogr= logr[i]; }
  }
  for (i=0, ans=0; i<1000; i++) { ans+= exp(logr[i]-maxlogr+500); }
  for (i=1000;i<(*B);i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    ans+= exp(-fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1) -maxlogr+500);
  }
  ans= log(ans/(.0+ (*B))) + maxlogr-500;

  free_dvector(sdprop,1,*nsel); free_dvector(sopt,1,*nsel);
  free_dvector(mprop, 1,*nsel); free_dvector(thsim, 1, *nsel);
  free_dvector(logr,0,999);
  free_dmatrix(Vprop,1,*nsel,1,*nsel); free_dmatrix(cholVprop,1,*nsel,1,*nsel); free_dmatrix(cholVpropinv,1,*nsel,1,*nsel);
  return(ans);
}




double f2opt_imomU(double *th) {
  //last element in th corresponds to eta=log(tau), i.e. log residual variance
  double ans;
  ans= fimomUNegC_non0(th+1,f2opt_pars.sumy2,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.alpha,f2opt_pars.lambda,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomUNegC_non0(double *th, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
//loops over all elem in th (i.e. th has length *nsel+1 and contains non-zero elements only). th is indexed at 0.
//Note: last element in th corresponds to eta=log(tau), i.e. log residual variance
  int i;
  double ans, ytXth, sumlogth, suminvth, th2, eta, phi;
  eta= th[*nsel]; phi= exp(eta);
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth + sumlogth + .5*eta*(*n - *nsel + *alpha);
  return ans;
}


//Hessian of fimomUNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)+1][1:(*nsel)+1])
// - th: th[1:(*nsel)+1] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomUNegC_non0(double **ans, double *th, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2, eta, phi, suminvth, ytXth, *XtXth;

  XtXth= dvector(1,*nsel);
  eta= th[*nsel +1]; phi= exp(eta);
  Asel_x(XtX,*p,th,*nsel,sel-1,XtXth);
  for (i=1, ytXth=0, suminvth=0; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= (XtX->at(sel[i-1]*(*p)+sel[i-1])) / phi + 6.0*(*tau)*phi/(th2*th2) - 2.0/th2;
    ans[i][*nsel+1]= ans[*nsel+1][i]= -2.0*(*tau)*phi/(th2*th[i]) - (XtXth[i]-ytX[sel[i-1]])/phi;
    ytXth+= ytX[sel[i-1]] * th[i];
    suminvth+= 1/(th[i]*th[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= (XtX->at(sel[i-1]*(*p)+sel[j-1]))/phi;
    }
  }
  ans[*nsel+1][*nsel+1]= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th+1, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth;
  free_dvector(XtXth,1,*nsel);
}

void imomModeU(double *th, PolynomialRootFinder::RootStatus_T *status, double *sumy2, crossprodmat *XtX, double *ytX, double *tau, double *alpha, double *lambda, int *sel, int *nsel, int *n, int *p) {
  //piMOM mode when phi is unknown using gradient algorithm
  // - th: contains initial estimate (theta,log(phi)) at input and mode at output
  // - status: indicates if root finding has been successful
  bool found=false;
  int i, j, niter=0, root_count;
  double err= 1.0, *coef, *real_vector, *imag_vector, phi, phinew, a, b, b2, c, d, suminvth2, *XtXth;
  Polynomial poly;

  phi= exp(th[*nsel +1]);
  b= (*n -(*nsel) + *alpha)/2.0;
  b2= b*b;

  coef= dvector(0,4);
  real_vector= dvector(0,4);
  imag_vector= dvector(0,4);
  XtXth= dvector(1,*nsel);

  coef[1]= 0.0;
  coef[2]= -2;
  while ((err > 1.0e-5) & (niter<50)) {
    coef[0]= 2.0 * (*tau) * phi;
    suminvth2= 0.0;
    err= 0;
    //Update th
    for (i=1; i<=(*nsel); i++) {
      coef[3]= ytX[sel[i-1]];
      for (j=1; j<i; j++) { coef[3]-= (XtX->at(sel[i-1]*(*p)+sel[j-1])) * th[j]; }
      for (j=i+1; j<=(*nsel); j++) { coef[3]-= (XtX->at(sel[i-1]*(*p)+sel[j-1])) * th[j]; }
      coef[3]= coef[3]/phi;
      coef[4]= -(XtX->at(sel[i-1]*(*p)+sel[i-1]))/phi;
      poly.SetCoefficients(coef, 4);
      (*status)= poly.FindRoots(real_vector,imag_vector,&root_count);

      j=0; found= false;
      while ((!found) & (j<=4)) {
        if (fabs(imag_vector[j])<1.0e-5) {
          if (((real_vector[j]>0) & (th[i]>0)) | ((real_vector[j]<0) & (th[i]<0))) {
            err= max_xy(err, fabs(th[i] - real_vector[j]));
            th[i]= real_vector[j];
            suminvth2 += 1.0/(th[i]*th[i]);
            found= true;
          }
        }
        j++;
      }
    }

    //Update phi
    a= (*tau) * suminvth2;
    c= 0;
    Asel_x(XtX,*p,th,*nsel,sel-1,XtXth);
    for (i=1; i<=(*nsel); i++) { c += -2.0*ytX[sel[i-1]]*th[i] + th[i]*XtXth[i]; }
    c= -.5*(*lambda + *sumy2 + c);
    d= sqrt(b2 - 4.0*a*c);

    if (-b > d) { phinew= (-b-d)/(2.0*a); } else { phinew= (-b+d)/(2.0*a); }
    err= max_xy(err,fabs(phi-phinew));
    phi= phinew;

    niter++;
  }

  th[*nsel +1]= log(phi);

  free_dvector(coef,0,4); free_dvector(real_vector,0,4); free_dvector(imag_vector,0,4);
  free_dvector(XtXth,1,*nsel);
}


void imomUIntegralApproxC(double *ILaplace, double *thopt, int *sel, int *nsel, int *n, int *p, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *logscale, int *hessian) {
  int iter, maxit=100, emptyint;
  double ftol= 1.0e-10, **dirth, **Vopt, **Voptinv, detVopt, emptydouble=0, **emptymatrix, fopt;
  PolynomialRootFinder::RootStatus_T status;

  Vopt= dmatrix(1,*nsel +1,1,*nsel +1); Voptinv= dmatrix(1,*nsel +1,1,*nsel +1);
  dirth= dmatrix(1,*nsel +1,1,*nsel +1);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  set_f2opt_pars(&emptydouble,emptymatrix,sumy2,XtX,ytX,alpha,lambda,&emptydouble,tau,&emptyint,n,p,sel,nsel);

  //Minimization
  imomModeU(thopt,&status,sumy2,XtX,ytX,tau,alpha,lambda,sel,nsel,n,p);
  set_f2opt_pars(&emptydouble,emptymatrix,sumy2,XtX,ytX,alpha,lambda,&emptydouble,tau,&emptyint,n,p,sel,nsel);
  if (status == PolynomialRootFinder::SUCCESS) {
    fopt= f2opt_imomU(thopt);
  } else {
    ddiag(dirth,1,*nsel +1);
    minimize(thopt, dirth, *nsel +1, ftol, &iter, &fopt, f2opt_imomU, maxit);
  }

  if (*hessian ==1) {
    //Laplace approx
    fppimomUNegC_non0(Vopt,thopt,sumy2,XtX,ytX,alpha,lambda,tau,n,p,sel,nsel);
    invdet_posdef(Vopt,*nsel +1,Voptinv,&detVopt);
    (*ILaplace)= -fopt - 0.5*log(detVopt) + .5*(*nsel)*log(2.0*(*tau));
  } else {
    (*ILaplace)= -fopt - 0.5*(*nsel)*log(*n +.0) + .5*(*nsel)*log(2.0*(*tau));  //BIC-type approximation
  }

  free_dmatrix(Vopt,1,*nsel +1,1,*nsel +1); free_dmatrix(Voptinv,1,*nsel +1,1,*nsel +1);
  free_dmatrix(dirth,1,*nsel +1,1,*nsel+1); free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
}



//Product iMOM marginal density for unknown phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - x: design matrix (includes all covariates, even those excluded under the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - tau: prior dispersion parameter
// - method: method to approximate the integral for known phi. 0 for Laplace approx which may underestimate true value, 1 for exact evaluation which can be very computationally expensive. 2 ('Hybrid') integrates wrt phi numerically (Romberg) and wrt theta via Laplace approx (Laplace error is adjusted based on exact evaluation for single phi value)
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale
// - alpha, lambda: prior for phi (residual variance) is Inverse Gamma (.5*alpha,.5*lambda)
double f2int_imom(double phi) {
  int one=1, *inputlog= f2int_pars.logscale;
  double ans, *inputphi= f2int_pars.phi;
  f2int_pars.phi= &phi;
  f2int_pars.logscale= &one;
  ans= exp(pimomMarginalKC(f2int_pars.sel,f2int_pars.nsel,&f2int_pars) + dinvgammaC(phi,.5*(*f2int_pars.alpha),.5*(*f2int_pars.lambda),1) - *(f2int_pars.offset));
  f2int_pars.phi= inputphi;
  f2int_pars.logscale= inputlog;
  return(ans);
}


// [[Rcpp::export]]
SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP Sngroups, SEXP Snvaringroup) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1, SoptimMethod=1, Soptim_maxit= -1, emptyint=1, usethinit=0, knownphi=0;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *x=REAL(Sx), *ytX=REAL(SytX), *tau=REAL(Stau), *alpha=REAL(Salpha), *lambda=REAL(Slambda), *rans, emptydouble=0, offset=0, *taualpha=NULL;
  struct marginalPars pars;
  crossprodmat *XtX;
  SEXP ans;

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  set_marginalPars(&pars,&emptyint,n,n,p,y,&emptyint,sumy2,&emptydouble,&emptydouble,x,&emptydouble,XtX,ytX,method,&emptyint,&emptyint,&SoptimMethod,&Soptim_maxit,&usethinit,&emptydouble,B,alpha,lambda,&knownphi,&emptydouble,tau,&emptydouble,taualpha,taualpha,&r,NULL,NULL,NULL,NULL,logscale,&offset,NULL,NULL,INTEGER(Sngroups),NULL,INTEGER(Snvaringroup));
  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalUC(sel, nsel, &pars);
  delete XtX;
  UNPROTECT(1);
  return ans;
}


double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  bool posdef;
  int i, j, zero=0, one=1, *inputlog, hessian;
  double ans=0, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, *inputphi, num, den, term1, alphahalf=.5*(*(*pars).alpha);

  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*LOG_M_PI + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    if ((*(*pars).method)==0) {  //Laplace
      int prior=2, symmetric=1;
      ans= nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);
    } else {

      V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel); thest= dvector(1,*nsel+1);

      addct2XtX((*pars).tau,(*pars).XtX,sel,nsel,(*pars).p,V); //add tau to diagonal elem of XtX
      inv_posdef_upper(V,*nsel,Vinv,&posdef);
      Asym_xsel(Vinv,*nsel,(*pars).ytX,sel,thest);
      for (i=0, sumer2=0; i<(*(*pars).n); i++) {
        for (j=1, ypred=0; j<=(*nsel); j++) { ypred += (*pars).x[i + (*(*pars).n)*sel[j-1]] * thest[j]; }
        er= (*pars).y[i] - ypred;
        sumer2+= er*er;
      }
      phiest= (sumer2 + (*(*pars).lambda))/(*(*pars).alpha + *(*pars).n);
      if ((*(*pars).method)==2) {  //Plug-in

        hessian=0;
        thest[*nsel +1]= log(phiest);
        imomUIntegralApproxC(&ans,thest,sel,nsel,(*pars).n,(*pars).p,(*pars).sumy2,(*pars).XtX,(*pars).ytX,(*pars).alpha,(*pars).lambda,(*pars).tau,&one,&hessian);
        ans= ans + alphahalf*log(.5*(*(*pars).lambda)) - .5*(*(*pars).n)*LOG_M_2PI - gamln(&alphahalf);

      } else if ((*(*pars).method)==1) {  //MC for each fixed phi + univariate integration
        set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
        inputphi= (*pars).phi; (*pars).phi= &phiest;
        (*(*pars).method)= 0; inputlog= (*pars).logscale; (*pars).logscale= &one; //Laplace approx for phi=phiest
        intlapl= pimomMarginalKC(sel, nsel, pars);
        (*pars).phi= inputphi; (*(*pars).method)= 1; (*pars).logscale= inputlog;  //reset input values for phi, method
        f2int_pars.offset= &intlapl; //f2int_imom returns result divided by exp(intlapl) to avoid numerical overflow
        ans= intlapl + log(qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));

      } else if ((*(*pars).method)==3) {  //Hybrid Laplace - MC - Univariate integration
        set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
        inputphi= (*pars).phi; (*pars).phi= &phiest;
        (*(*pars).method)= 1; //IS evaluation of marginal for phi=phiest
        intmc= pimomMarginalKC(sel, nsel, pars);
        (*(*pars).method)= 0; //Laplace approx for phi=phiest
        intlapl= pimomMarginalKC(sel, nsel, pars);
        (*pars).phi= inputphi; (*(*pars).method)= 2;  //reset input values for phi, method
        if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
        f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx
        f2int_pars.offset= &intlapl; //f2int_imom returns result divided by exp(intlapl) to avoid numerical overflow
        ans= intmc + log(qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf)); //adjusment is intmc - intlapl, but intlapl is the offset so needs to added back in
      }
      free_dmatrix(V,1,*nsel,1,*nsel);
      free_dmatrix(Vinv,1,*nsel,1,*nsel);
      free_dvector(thest,1,*nsel+1);
    }
  }
  if ((*(*pars).logscale)==0) ans= exp(ans);
  return(ans);
}



//*************************************************************************************
// Product eMOM routines
//*************************************************************************************
double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}

double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int prior=3, symmetric=1;
  double ans;
  ans= nlpMargSkewNorm(sel, nsel, pars, &prior, &symmetric);
  return ans;
}



//*************************************************************************************
// Zellner's prior routines
//*************************************************************************************

// Marginal likelihood for linear models under Zellner's prior
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - logscale: if set to 1 result is returned in log scale
double zellnerMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int i,j;
  double *m, s, **S, **Sinv, detS, num, den, adj, tau= *(*pars).tau, logphi= log(*(*pars).phi), ans=0.0, zero=0;

  if (*nsel ==0) {

    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);

  } else {

    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX into S
    adj= (tau+1)/tau;
    for (i=1; i<=(*nsel); i++) {
      S[i][i]= S[i][i] * adj;
      for (j=i+1; j<=(*nsel); j++) { S[i][j]= S[i][j] * adj; }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + (*nsel)*log(tau+1.0));
    //den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*logtau);
    ans= num - den;

    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}

double zellnerMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j;
  double num, den, ans=0.0, term1, *m, **S, **Sinv, detS, adj, tau= *(*pars).tau, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    adj= (tau+1)/tau;
    for (i=1; i<=(*nsel); i++) {
      S[i][i]= S[i][i] * adj;
      for (j=i+1; j<=(*nsel); j++) { S[i][j]= S[i][j] * adj; }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= .5*(*(*pars).n + *(*pars).alpha);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI) + .5 * (*nsel) *log((*(*pars).tau)+1.0) + gamln(&alphahalf);
    //den= .5*(*(*pars).n * LOG_M_2PI + log(detS)) + .5 * (*nsel) *log(*(*pars).tau) + gamln(&alphahalf);
    ans= num - den;

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


//*************************************************************************************
// Normal's prior routines
//*************************************************************************************

// Marginal likelihood for linear models under Normal's prior. For now, the covariance matrix is assumed to be the identity
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - logscale: if set to 1 result is returned in log scale
double normalidMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  double *m, s, **S, **Sinv, detS, num, den, tau= *(*pars).tau, tauinv=1/tau, taugroup=*(*pars).taugroup, taugroupinv=1/taugroup, logphi= log(*(*pars).phi), ans=0.0, zero=0;
  int i, *isgroup=(*pars).isgroup, varingroup, singlevarcount=0;

  if (*nsel ==0) {

    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);

  } else {

    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX into S
    for (i = 1; i <= *nsel; i++ ) {
      varingroup = isgroup[sel[i-1]];
      if (varingroup==0) {
        S[i][i] += tauinv;
        singlevarcount++;
      } else {
        S[i][i] += taugroupinv;
      }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*log(tau));
    ans= num - den;

    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}

double normalidMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  double num, den, ans=0.0, term1, *m, **S, **Sinv, detS, tau= *(*pars).tau, tauinv=1/tau, taugroup=*(*pars).taugroup, taugroupinv = 1/taugroup, nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda), ss, zero=0;
  int i, singlevarcount=0, *isgroup=(*pars).isgroup, varingroup;
  if (*nsel ==0) {

    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*(LOG_M_PI) + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);

  } else {

    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&zero,(*pars).XtX,sel,nsel,(*pars).p,S);  //copy XtX onto S
    for (i = 1; i <= *nsel; i++ ) {
      varingroup = isgroup[sel[i-1]];
      if (varingroup==0) {
        S[i][i] += tauinv;
        singlevarcount++;
      } else {
        S[i][i] += taugroupinv;
      }
    }
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= .5*(*(*pars).n + *(*pars).alpha);

    ss= *(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel);
    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2.0) - log(ss));
    den= .5*(*(*pars).n * LOG_M_2PI + log(detS)) + .5 * (singlevarcount * log(tau) + (*nsel - singlevarcount) * log(taugroup)) + gamln(&alphahalf);
    ans= num - den;

    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);

  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


//*************************************************************************************
// ADDITIVE MODELS
//*************************************************************************************

// [[Rcpp::export]]
SEXP bsplineCI(SEXP x, SEXP degree, SEXP knots) {
  int nknots=LENGTH(knots), nx=LENGTH(x);
  SEXP ans;
  PROTECT(ans= Rf_allocVector(REALSXP, (nknots-INTEGER(degree)[0]-1) * nx));

  bspline_vec(REAL(ans), REAL(x), &nx, INTEGER(degree), REAL(knots), &nknots);

  UNPROTECT(1);
  return ans;
}
