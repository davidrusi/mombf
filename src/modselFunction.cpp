#include "modselFunction.h"
using namespace std;

modselFunction::modselFunction(int *sel, int *nsel, struct marginalPars *pars, pt2jointFun fun=NULL) {

  this->nsel= nsel;
  this->sel= sel;
  this->pars= pars;
  this->maxiter= 50;
  this->ftol= 0.001;
  this->thtol= 0.0001;
  this->fun= fun;

  this->updateUniv= NULL;

  this->gradhessUniv= NULL;
  this->gradUniv= NULL;
  this->hess= NULL;

  this->funupdate= NULL;
  this->gradhessupdate= NULL;

}


modselFunction::~modselFunction() {

}


//Evaluate fun at th
double modselFunction::evalfun(double *th) {

  return (fun(th, this->sel, this->nsel, this->pars));

}


//*******************************************************************************************************************
//OPTIMIZATION ALGORITHMS
//*******************************************************************************************************************

//Classical Coordinate Descent Algorithm sequentially updating each parameter (uses updateUniv)
// Input
// - thini: initial parameter value
// Output
// - thopt: final parameter value
// - fopt: value of objective function (fun) at thopt
void modselFunction::cda(double *thopt, double *fopt, double *thini) {

  int j, iter=0;
  double therr=1, ferr=1, thnew, fnew;

  if ((this->fun)==NULL) Rf_error("To run CDA you need to specify evalfun");
  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");

  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) thopt[j]= thini[j];
  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
    fnew= this->evalfun(thopt);
    ferr= (*fopt) - fnew;
    (*fopt)= fnew;
    iter++;
  }
}


//Same but does not evaluate objective function (faster but stopping depends only on change in thopt)
void modselFunction::cda(double *thopt, double *thini) {

  int j, iter=0;
  double therr=1, thnew;

  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");

  for (j=0; j< *nsel; j++) thopt[j]= thini[j];
  while ((iter< this->maxiter) & (therr > this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
    iter++;
  }
}


//BLOCK CDA JOINTLY UPDATING ALL PARAMETERS
//In contrast to cda here th[j] is updated without updating first th[0], ..., th[j-1].
//Hence even if CDA were guaranteed to converge blockcda may not, as it cannot be interpreted as a sequence of univariate optimizations
void modselFunction::blockcda(double *thopt, double *fopt, double *thini) {

  int j, iter=0;
  double *thnew, fnew, therr=1, ferr=1;

  if ((this->fun)==NULL) Rf_error("To run blockcda you need to specify evalfun");
  thnew= dvector(0,*nsel);

  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0; j< *nsel; j++) { (*(this->updateUniv))(thnew+j, j, thopt, this->sel, this->nsel, this->pars); }

    fnew= this->evalfun(thnew);
    ferr= (*fopt) - fnew;
    if (ferr>0) {
      (*fopt)= fnew;
      for (j=0,therr=0; j< *nsel; j++) {
	therr= max_xy(therr, fabs(thnew[j] - thopt[j]));
	thopt[j]= thnew[j];
      }
    }
    iter++;

  }

  free_dvector(thnew,0,*nsel);

}





//Classical CDA with approx updates given by Newton's method (uses gradhess)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double thcur, therr=1, ferr=1, fnew, delta, g, H;

  if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify evalfun");
  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");

  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0, therr=ferr=0; j< *nsel; j++) {

      gradhessUniv(&g, &H, j, thopt, sel, nsel, pars);
      delta= g/H;

      nsteps= 1; found= false;
      thcur= thopt[j];
      while (!found & (nsteps<=maxsteps)) {

	thopt[j] -= delta;
	fnew= this->evalfun(thopt);

	if (fnew < *fopt) {
	  found= true;
	  ferr+= *fopt - fnew;
	  (*fopt)= fnew;
	  therr= max_xy(therr, fabs(delta));
	} else {
	  thopt[j]= thcur;
	  delta /= 2.0;
	  nsteps++;
	}

      } //end while !found

    } //end for j
    iter++;

  } //end while iter


}



//BLOCK CDA WITH NEWTON METHOD UPDATES (USES gradhess)
//Each parameter is updated to th[j] - 0.5^k grad[j]/hess[j] as in cdaNewton, but here grad[j] and hess[j] are evaluated at the current th prior to updating any th
//In contrast, in cdaNewton grad[j] and hess[j] are evaluated after updating th[0], ..., th[j-1]
void modselFunction::blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double *thcur, therr=1, ferr=1, fnew, *delta, *g, *H;

  if ((this->fun)==NULL) Rf_error("To run blockcdaNewton you need to specify evalfun");
  if ((this->gradhessUniv)==NULL) Rf_error("To run blockcdaNewton you need to specify either gradhessUniv");
  thcur= dvector(0,*nsel); delta= dvector(0,*nsel); g= dvector(0,*nsel); H= dvector(0,*nsel);

  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thcur[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    therr= ferr= 0;
    for (j=0; j< *nsel; j++) {
      gradhessUniv(g+j, H+j, j, thopt, sel, nsel, pars);
      delta[j]= g[j]/H[j];
    }

    nsteps= 1; found= false;
    while (!found & (nsteps<=maxsteps)) {

      for (j=0; j< *nsel; j++) { thopt[j] -= delta[j]; therr= max_xy(therr, fabs(delta[j])); }
      fnew= this->evalfun(thopt);

      if (fnew < *fopt) {
	found= true;
	ferr+= *fopt - fnew;
	(*fopt)= fnew;
      } else {
	for (j=0; j< *nsel; j++) { thopt[j]= thcur[j]; delta[j] /= 2.0; }
	therr= 0;
	nsteps++;
      }

    } //end while !found
    iter++;

  } //end while iter

  free_dvector(thcur, 0,*nsel); free_dvector(delta, 0,*nsel); free_dvector(g, 0,*nsel); free_dvector(H, 0,*nsel);

}


