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

typedef double (*pt2jointFun)(double *th, int *sel, int *nsel, struct marginalPars *pars);  //log-function to be maximized/integrated (e.g. log-likelihood + log-prior)
typedef void (*pt2gradUniv)(double *grad, int j, double *th, int *sel, int *nsel, struct marginalPars *);  //gradient wrt th[j]. Parameters (grad,j,th,sel,nsel,pars)
typedef void (*pt2gradhessUniv)(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *);  //gradient wrt th[j] and hessian wrt th[j]^2
typedef void (*pt2hess)(double **hess, double *th, int *sel, int *nsel, struct marginalPars *);  //full hessian matrix wrt th[1]^2, th[1]th[2] etc.
typedef void (*pt2updateUniv)(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars); //function updating th[j] to thnew
typedef void (*pt2updateBlock)(double *thnew, double *th, int *sel, int *nsel, struct marginalPars *pars); //function updating th[j] to thnew



//*************************************************************************************************************
// EXAMPLE OF USAGE
//*************************************************************************************************************

/*
  double foo(double *th, int *sel, int *nsel, struct marginalPars *pars) {  //returns sum_j th[sel[j]]^2 + sum_{l \neq j} th[sel[j]] th[sel[l]]
    int j, k; double ans;
    for (j=0, ans=0; j< *nsel; j++) { for (k=0; k< *nsel; k++) { (*thnew) += th[sel[j]] * th[sel[k]]; } }
    return(ans);
  }

  void foogradhess(double *grad, double *hess, int j, double *th, int *sel, int *nsel) {  //returns grad= 2 th[sel[j]] + sum_{l \neq j} th[sel[l]]; hess= 2.0;
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

  modselFunction(int *sel, int *nsel, struct marginalPars *pars, pt2jointFun fun);
  modselFunction(int *sel, int *nsel, struct marginalPars *pars) : modselFunction(sel, nsel, pars, NULL) {};
  ~modselFunction();

  int maxiter; //Maximum number of iterations in optimization algorithms (1 iter corresponds to updating all parameters)
  double ftol;  //Tolerance for objective function. Optimization stops when improvement is < ftol
  double thtol; //Tolerance for parameter values. Optim stops when largest change in th < thtol

  double evalfun(double *th); //Evaluate fun at th

  //Optimization algorithms
  void cda(double *thopt, double *fopt, double *thini);  //Classical Coordinate Descent Algorithm sequentially updating each parameter (uses updateUniv)
  void cda(double *thopt, double *thini);  //same but does not evaluate objective function (faster but stopping depends only on change in thopt)
  void blockcda(double *thopt, double *fopt, double *thini);  //Block CDA jointly updating all parameters (uses updateBlock)
  void cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Classical CDA with approx updates given by Newton's method (uses gradhess)
  void blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps); //Block CDA with Newton method updates (uses gradhess)

  //pointers to functions computing gradient, hessian and updates
  pt2jointFun fun;
  pt2updateUniv updateUniv;
  pt2updateBlock updateBlock;
  pt2gradUniv gradUniv;
  pt2gradhessUniv gradhessUniv;
  pt2hess hess;
  //double (*fun)(double *th, int *sel, int *nsel, struct marginalPars *pars); //Pointer to function returns log-likelihood + log-prior
  //void (*updateUniv)(double *thnew, int j, double *th, int *sel, int *nsel, struct marginalPars *pars); //Pointer to function returning univariate update for th[j] given other elements in th
  //void (*gradhessUniv)(double *grad, double *hess, int j, double *th, int *sel, int *nsel, struct marginalPars *pars);  //Pointer to function returning gradient wrt th[j]
  //void (*gradUniv)(double *grad, int j, double *th, int *sel, int *nsel, struct marginalPars *pars);  //Pointer to function returning gradient wrt th[j]


private:

  int *nsel;  //total number of parameters
  int *sel;   //sel[0], ..., sel[*nsel -1] contain the indexes of active variables
  struct marginalPars *pars;

};

#endif

