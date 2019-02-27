#ifndef MODSELFUNCTION
#define MODSELFUNCTION 1

#include <RcppArmadillo.h>
#include <map>
#include <string>
#include "modelSel.h"
#include "cstat.h"
using namespace std;

//*************************************************************************************************************
// TYPEDEF: Pointers to functions returning minus log-integrand (e.g. -log(likelihood) - log(prior)), its gradients and hessians
//*************************************************************************************************************

typedef void (*pt2updateUniv)(double *thnew, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);

typedef void (*pt2fun)(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);  //objective function
typedef void (*pt2funupdate)(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs);

typedef void (*pt2gradUniv)(double *grad, int j, double *th, int *sel, int *thlength, struct marginalPars *, std::map<string, double*> *funargs);
typedef void (*pt2gradhessUniv)(double *grad, double *hess, int j, double *th, int *sel, int *thlength, struct marginalPars *, std::map<string, double*> *funargs);
typedef void (*pt2hess)(double **H, double *th, int *sel, int *thlength, struct marginalPars *, std::map<string, double*> *funargs);


//*************************************************************************************************************
//
// FOR A FULL EXAMPLE OF USAGE SEE BELOW
//
// USER-PROVIDED INPUT WHEN CREATING AN OBJECT OF THIS CLASS
//
// - fun: evaluates the objective function (optionally computing auxiliary arguments "funargs")
// - funupdate: updates the objective function value after changing th[j] to thjnew, using "funargs", and also updates "funargs"
// - updateUniv: function giving a rule to update th to thnew
// - gradUniv, gradhessUniv: functions computing only gradient or both (gradient,hessian) wrt to th[j]
// - hess: function computing the full hessian wrt th[1]^2, th[1]th[2] etc.
//
// NOTES
//
// - No need to specify all these inputs, just those required by the algorithms one intends to use (e.g. cda only requires updateUniv)
// - In all functions there's thlength parameters th[0], ..., th[thlength-1]. The corresponding variable indexes are sel[0],...,sel[thlength-1]
// - "funargs" is an optional argument to share computations between fun, funupdate and gradhessUniv. The idea is that fun and funupdate compute not only the function at the provided th but also "funargs", and then funargs can be re-used when evaluating the function, gradient and/or hessian. For example suppose that
//
// fun= sum_i (y[i] - ypred[i])^2 where ypred[i]= sum_l x[i][l] th[l]. Its gradient wrt th[j] is
//
// grad[j]= sum_i -2(y[i] - ypred[i]) x[i][j] th[j]
//
// When first evaluating fun we can store funargs[1]= ypred.
// If we update th[j] to thnew it's wasteful to call fun again, instead we use funupdate to update ypred[i] --> ypred[i] - x[i][j] (thnew - th[j]), avoiding the sum over l != j.
// Then to compute the gradient at thnew we can pass "funargs" to gradhess so it uses the readily available ypred
//
//
// METHODS PROVIDED BY THE CLASS
//
// - cda, blockcda: Coordinate Descent Algorithms using updateUniv
// - cdaNewton, blockcdaNewton: Newton-based CDA using gradhess





//*************************************************************************************************************
// CLASS DEFINITION
//*************************************************************************************************************

class modselFunction {

public:

  //Constructor and destructor
  modselFunction(int *sel, int thlength, struct marginalPars *pars, pt2fun fun);
  ~modselFunction();

  //PARAMETERS THAT CAN BE ACCESSED/MODIFIED BY THE USER
  int maxiter; //Maximum number of iterations in optimization algorithms (1 iter corresponds to updating all parameters)
  double ftol;  //Tolerance for objective function. Optimization stops when improvement is < ftol
  double thtol; //Tolerance for parameter values. Optim stops when largest change in th < thtol

  //FUNCTIONS THAT SHOULD BE SET BY THE USER (TYPICALLY ONLY A SUBSET OF THESE IS NEEDED, SEE EXAMPLES)
  pt2updateUniv updateUniv;  //rule to update th[j] into thjnew

  pt2fun fun;  //evaluate objective function
  pt2funupdate funupdate; //evaluate objective function by updating its value at the previous th  (optional, typically much faster than fun)
  pt2gradUniv gradUniv; //evaluate gradient
  pt2gradhessUniv gradhessUniv; //evaluate gradient/hessian wrt th[j]
  pt2hess hess; //evaluate full hessian matrix of -fun (rather than fun). Important: returned H should have indexes [1..thlength][1..thlength] (rather than 0-indexed th used in other functions)

  //PUBLIC METHODS PROVIDED BY THE CLASS
  void evalfun(double *f, double *th, std::map<string, double *> *funargs); //Evaluate fun at th (and optionally return the value of funargs)
  void evalfunupdate(double *fnew, double *thjnew, int j, double *f, double *th, std::map<string, double *> *funargs); //Eval fun at thjnew by updating its value f at th[j]

  void cda(double *thopt, double *fopt, double *thini);  //Coordinate Descent Algorithm (uses updateUniv)
  void cda(double *thopt, double *thini);  //same but does not evaluate objective function (stopping depends only on change in thopt)
  void cda(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs);
  void blockcda(double *thopt, double *fopt, double *thini);  //Block CDA jointly updating all parameters (uses updateUniv)
  void cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //CDA with approx updates given by Newton's method (uses gradhess)
  void cdaNewton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps);
  void blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Block CDA with Newton method updates (uses gradhess)

  double laplaceapprox(double *thopt, double *fopt, double **H); //Laplace approximation to int exp(-fun(th)) dth
  double laplaceapprox(double *thopt, double *fopt, std::map<string, double *> *funargs);
  double laplaceapprox(double *thopt, std::map<string, double *> *funargs);


private:

  int thlength;  //total number of parameters
  int *sel;   //sel[0], ..., sel[thlength -1] contain the indexes of active variables
  struct marginalPars *pars;

};

#endif



//*************************************************************************************************************
// EXAMPLE OF USAGE
//*************************************************************************************************************

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
void foo(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {
  int k, l;
  for (k=0, (*f)=0; k< *thlength; k++) {
    (*f) += (double)(sel[k]+1) * th[k] * th[k];
    for (l=k+1; l< *thlength; l++) { (*f) += th[k] * th[l]; }
  }
}

//Compute gradient and hessian wrt th[j], not using funargs
void foogradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int l;
  (*hess)= 2.0 * (double)(sel[j]+1);
  (*grad)= (*hess) * th[j];
  for (l=0; l< j; l++) { (*grad)+= th[l]; }
  for (l=j+1; l< *thlength; l++) { (*grad)+= th[l]; }
}

//Return univariate optimum for th[j], that is thnew= -0.5/(sel[j]+1) * sum_{l \neq j} th[l]. Not using funargs
void fooupdateUniv(double *thnew, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int l;
  *thnew= 0;
  for (l=0; l< j; l++) (*thnew)-= th[l];
  for (l=j+1; l< *thlength; l++) (*thnew)-= th[l];
  (*thnew) *= 0.5/((double)(sel[j]+1));
}


//Evaluate function and funargs
void fooargs(double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {
  int k, l;
  double sumth=0, sumth2=0, sumcrossprod=0;
  for (k=0, (*f)=0; k< *thlength; k++) {
    sumth += th[k];
    sumth2 += (double)(sel[k]+1) * th[k] * th[k];
    for (l=k+1; l< *thlength; l++) { sumcrossprod += th[k] * th[l]; }
  }
  (*f)= sumth2 + sumth;
  *(*funargs)["sumth"]= sumth;
  *(*funargs)["sumth2"]= sumth2;
  *(*funargs)["sumcrossprod"]= sumcrossprod;
}

//Update function and funargs from changing th[j] to thjnew
void fooupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double *> *funargs) {
  double thdif= *thjnew - th[j];
  *(*funargs)["sumth"] += thdif;
  *(*funargs)["sumth2"] += (double)(sel[j]+1) * (pow(*thjnew,2) - pow(th[j],2));
  *(*funargs)["sumcrossprod"] += (*thjnew - th[j]) * (*(*funargs)["sumth"] - *thjnew);
  (*fnew)= *(*funargs)["sumth2"] + *(*funargs)["sumcrossprod"];
}

//Compute gradient and hessian wrt th[j], using funargs
void foogradhessargs(double *grad, double *hess, int j, double *th, int *sel, int *thlength, struct marginalPars *pars, std::map<string, double*> *funargs) {
  (*hess)= 2.0 * (double)(sel[j]+1);
  (*grad)= (*hess) * th[j] + (*(*funargs)["sumth"]) - th[j];
}



void testfunction() {

  int thlength=2, *sel;
  double *thini, *thopt, fopt;
  struct marginalPars *pars= NULL;
  modselFunction *msfun;
  std::map<string, double *> funargs;

  //Alloc memory for elements in funargs. For vector arguments use dvector
  double sumth= 0, sumth2= 0, sumcrossprod= 0;
  funargs["sumth"]= &sumth; funargs["sumth2"]= &sumth2; funargs["sumcrossprod"]= &sumcrossprod;

  sel= ivector(0,thlength); thini= dvector(0,thlength); thopt= dvector(0,thlength);
  sel[0]= 0; sel[1]= 2;
  thini[0]= 1; thini[1]= 1;
  msfun= new modselFunction(sel, thlength, pars, NULL);

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

  free_ivector(sel, 0,thlength); free_dvector(thini, 0,thlength); free_dvector(thopt, 0,thlength);
  delete msfun;

}


*/
