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

typedef void (*pt2updateUniv)(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs);

typedef void (*pt2fun)(double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs);  //objective function
typedef void (*pt2funupdate)(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs);

typedef void (*pt2gradUniv)(double *grad, int j, double *th, int *sel, int *nsel, struct marginalPars *, std::map<string, double*> *funargs);
typedef void (*pt2gradhessUniv)(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *, std::map<string, double*> *funargs);
typedef void (*pt2hess)(double **hess, double *th, int *sel, int *nsel, struct marginalPars *);


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
// - In all functions there's nsel parameters th[0], ..., th[nsel-1]. The corresponding variable indexes are sel[0],...,sel[nsel-1]
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

  modselFunction(int *sel, int *nsel, struct marginalPars *pars, pt2fun fun);
  ~modselFunction();

  int maxiter; //Maximum number of iterations in optimization algorithms (1 iter corresponds to updating all parameters)
  double ftol;  //Tolerance for objective function. Optimization stops when improvement is < ftol
  double thtol; //Tolerance for parameter values. Optim stops when largest change in th < thtol

  void evalfun(double *f, double *th, std::map<string, double *> *funargs); //Evaluate fun at th (and optionally return the value of funargs)
  void evalfunupdate(double *fnew, double *thjnew, int j, double *f, double *th, std::map<string, double *> *funargs); //Eval fun at thjnew by updating its value f at th[j]

  //Optimization algorithms
  void cda(double *thopt, double *fopt, double *thini);  //Coordinate Descent Algorithm (uses updateUniv)
  void cda(double *thopt, double *thini);  //same but does not evaluate objective function (stopping depends only on change in thopt)
  void cda(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs);
  void blockcda(double *thopt, double *fopt, double *thini);  //Block CDA jointly updating all parameters (uses updateUniv)
  void cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //CDA with approx updates given by Newton's method (uses gradhess)
  void cdaNewton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps);
  void blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Block CDA with Newton method updates (uses gradhess)

  //pointers to functions computing gradient, hessian and updates
  pt2updateUniv updateUniv;

  pt2fun fun;  //evaluate objective function
  pt2funupdate funupdate; //evaluate objective function by updating its value at the previous th  (optional, typically much faster than fun)
  pt2gradUniv gradUniv; //evaluate gradient
  pt2gradhessUniv gradhessUniv; //evaluate gradient/hessian
  pt2hess hess;

private:

  int *nsel;  //total number of parameters
  int *sel;   //sel[0], ..., sel[*nsel -1] contain the indexes of active variables
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
void foo(double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs) {
  int k, l;
  for (k=0, (*f)=0; k< *nsel; k++) {
    (*f) += (double)(sel[k]+1) * th[k] * th[k];
    for (l=k+1; l< *nsel; l++) { (*f) += th[k] * th[l]; }
  }
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

//Update value of function and funargs when changing th[j] to thjnew
void fooupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double *> *funargs) {
  double thdif= *thjnew - th[j];
  *(*funargs)["sumth"] += thdif;
  *(*funargs)["sumth2"] += (double)(sel[j]+1) * (pow(*thjnew,2) - pow(th[j],2));
  *(*funargs)["sumcrossprod"] += (*thjnew - th[j]) * (*(*funargs)["sumth"] - *thjnew);
  (*fnew)= *(*funargs)["sumth2"] + *(*funargs)["sumcrossprod"];
}

//Evaluate gradient and hessian
void foogradhess(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<string, double*> *funargs) {
  int l;
  (*hess)= 2.0 * (double)(sel[j]+1);
  (*grad)= (*hess) * th[j];
  for (l=0; l< j; l++) { (*grad)+= th[l]; }
  for (l=j+1; l< *nsel; l++) { (*grad)+= th[l]; }
}


void testfunction() {

  int nsel=2, *sel;
  double *thini, *thopt, fopt;
  struct marginalPars *pars= NULL;
  modselFunction *msfun;

  sel= ivector(0,nsel); thini= dvector(0,nsel); thopt= dvector(0,nsel);
  sel[0]= 0; sel[1]= 2;
  thini[0]= 1; thini[1]= 1;
  msfun= new modselFunction(sel, &nsel, pars, NULL);

  //Option 1. CDA without providing foo
  msfun->updateUniv= &fooupdate;
  msfun->cda(thopt, thini);
  Rprintf("cda. thopt= %f %f\n", thopt[0], thopt[1]);

  //Option 2. CDA providing foo (stops when foo doesn't improve further)
  msfun->fun= &foo;
  msfun->updateUniv= &fooupdate;
  msfun->cda(thopt, &fopt, thini);
  Rprintf("cda. thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 3. block CDA
  msfun->fun= &foo;
  msfun->updateUniv= &fooupdate;
  msfun->blockcda(thopt, &fopt, thini);
  Rprintf("blockcda. thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 4. cdaNewton not using funargs (requires function foo)
  msfun->fun= &foo;
  msfun->gradhessUniv= &foogradhess;
  msfun->cdaNewton(thopt, &fopt, thini, 1);
  Rprintf("cdaNewton. thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 5. cdaNewton using funargs (requires functions fooargs and fooupdate, and allocating memory for all elements in funargs)
  std::map<string, double *> funargs;
  double sumth= 0, sumth2= 0, sumcrossprod= 0; //if funargs has vector arguments use dvector to allocate them
  funargs["sumth"]= &sumth; funargs["sumth2"]= &sumth2; funargs["sumcrossprod"]= &sumcrossprod;
  msfun->fun= &fooargs;
  msfun->funupdate= &fooupdate;
  msfun->gradhessUniv= &foogradhess;
  msfun->cdaNewton(thopt, &fopt, thini, &funargs, 1);
  Rprintf("cdaNewton. thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  //Option 6. blockcdaNewton
  msfun->fun= &foo;
  msfun->gradhessUniv= &foogradhess;
  msfun->blockcdaNewton(thopt, &fopt, thini, 1);
  Rprintf("blockcdaNewton. thopt= %f %f; fopt=%f\n", thopt[0], thopt[1], fopt);

  free_ivector(sel, 0,nsel); free_dvector(thini, 0,nsel); free_dvector(thopt, 0,nsel);
  delete msfun;

}


*/
