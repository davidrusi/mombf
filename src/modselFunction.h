#ifndef MODSELFUNCTION
#define MODSELFUNCTION 1

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <map>
#include <string>
#include "modelSel.h"
#include "cstat.h"
using namespace std;

//*************************************************************************************************************
// TYPEDEF: Pointers to functions returning minus log-integrand (e.g. -log(likelihood) - log(prior)), its gradients and hessians
//
//- Any function requiring gradient and hessian wrt th[j] uses gradhessUniv if available, assumed to be faster than calling gradUniv and hessUniv separately
//- No need to specify all these when creating an object of class modselFunction, just those used by the algorithms one intends to use (e.g. cda only requires updateUniv)
//- In all functions there's nsel parameters th[0], ..., th[nsel-1]. The corresponding variable indexes are sel[0],...,sel[nsel-1].
//*************************************************************************************************************

typedef void (*pt2updateUniv)(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars); //function updating th[j] to thnew

typedef void (*pt2fun)(double *f, double *th, int *sel, int *nsel, struct marginalPars *pars, std::map<char, double *> *funargs);  //log-function to be maximized/integrated (e.g. log-likelihood + log-prior)
typedef void (*pt2gradUniv)(double *grad, int j, double *th, int *sel, int *nsel, struct marginalPars *);  //gradient wrt th[j]. Parameters (grad,j,th,sel,nsel,pars)
typedef void (*pt2gradhessUniv)(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *);  //gradient wrt th[j] and hessian wrt th[j]^2
typedef void (*pt2hess)(double **hess, double *th, int *sel, int *nsel, struct marginalPars *);  //full hessian matrix wrt th[1]^2, th[1]th[2] etc.


//Optionally you can provide a way to compute the (function,grad,hessian) due from updating th[j] to thnew[j] by updating their previous values
//Since all th's other than th[j] remain unchanged, these can be vastly faster than *pt2fun, pt2gradhessUniv
typedef void (*pt2funupdate)(double *fnew, double *thjnew, double *f, double *th, int *j, int *sel, int *nsel, struct marginalPars *pars, std::map<char, double *> *funargs);
typedef void (*pt2gradhessupdate)(double *gradnew, double *hessnew, double *thjnew, double *th, double *f, double *grad, double *hess, std::map<char, double *> *funargs);




//*************************************************************************************************************
// EXAMPLE OF USAGE
//*************************************************************************************************************

/*
  double foo(double *th, int *sel, int *nsel, struct marginalPars *pars) {  //returns sum_j th[sel[j]]^2 + sum_{l \neq j} th[sel[j]] th[sel[l]]
    int j, k; double ans;
    for (j=0, ans=0; j< *nsel; j++) { for (k=0; k< *nsel; k++) { (*thnew) += th[sel[j]] * th[sel[k]]; } }
    return(ans);
  }

  void foogradhess(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *pars) {  //returns grad= 2 th[sel[j]] + sum_{l \neq j} th[sel[l]]; hess= 2.0;
    int l;
    (*grad)= th[sel[j]];
    for (l=0; l< *nsel; l++) { (*grad)+= th[sel[j]]; }
    (*hess)= 2.0;
  }

  void fooupdate(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars) { //returns univariate optimum 0.5 * sum_{l \neq j} th[sel[j]]
    int l;
    for (l=0; l< j; l++) (*thnew)= th[sel[l]];
    for (l=j+1; l< *nsel; l++) (*thnew)= th[sel[l]];
    (*thnew) *= 0.5;
  }

  int nsel=2, *sel;
  double *thini, *thopt, fopt;
  modselFunction *msfun;

  sel= ivector(0,nsel); thini= dvector(0,nsel); thopt= dvector(0,nsel);
  sel[0]= 0; sel[1]= 1;
  thini[0]= 1; thini[1]= 1;
  msfun= new modselFunction(sel, &nsel, pars);

  //Option 1. Use fooupdate to run CDA
  msfun->updateUniv= &fooupdate;
  msfun->cda(thopt, thini);

  //Option 2. Use foograd and foohess to run cdaNewton
  msfun->fun= &foo;
  msfun->gradhessUniv= &foogradhess;
  msfun->cdaNewton(thopt, &fopt, thini);

  free_ivector(sel, 0,nsel); free_dvector(thini, 0,nsel); free_dvector(thopt, 0,nsel);
  delete msfun;
*/

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

  void evalfun(double *f, double *th); //Evaluate fun at th
  void evalfun(double *f, double *th, std::map<char, double *> *funargs); //Evaluate fun at th and return the value of funargs
  void evalfunupdate(double *fnew, double *thjnew, double *f, double *th, int *j, std::map<char, double *> *funargs); //Eval fun at thjnew by updating its value f at th[j]

  //Optimization algorithms
  void cda(double *thopt, double *fopt, double *thini);  //Coordinate Descent Algorithm (uses updateUniv)
  void cda(double *thopt, double *thini);  //same but does not evaluate objective function (faster but stopping depends only on change in thopt)
  void blockcda(double *thopt, double *fopt, double *thini);  //Block CDA jointly updating all parameters (uses updateUniv)
  void cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Classical CDA with approx updates given by Newton's method (uses gradhess)
  void blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Block CDA with Newton method updates (uses gradhess)

  //pointers to functions computing gradient, hessian and updates
  pt2updateUniv updateUniv;

  pt2fun fun;  //evaluate objective function
  pt2funupdate funupdate; //evaluate objective function by updating its value at the previous th  (optional, typically much faster than fun)
  pt2gradUniv gradUniv; //evaluate gradient
  pt2gradhessUniv gradhessUniv; //evaluate gradient/hessian
  pt2gradhessupdate gradhessupdate; //evaluate grad/hessian by updating their values at the previous th  (optional, typically much faster than gradhessUniv)
  pt2hess hess;

private:

  int *nsel;  //total number of parameters
  int *sel;   //sel[0], ..., sel[*nsel -1] contain the indexes of active variables
  struct marginalPars *pars;

};

#endif

