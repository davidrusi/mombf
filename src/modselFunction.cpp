#include "modselFunction.h"
using namespace std;

modselFunction::modselFunction(int *sel, int *nsel, struct marginalPars *pars, pt2jointFun fun) {

  this->nsel= nsel;
  this->sel= sel;
  this->pars= pars;
  this->maxiter= 50;
  this->ftol= 0.001;
  this->thtol= 0.0001;
  this->fun= fun;

  this->updateUniv= NULL;
  this->updateBlock= NULL;
  this->gradhessUniv= NULL;
  this->gradUniv= NULL;
  this->hess= NULL;

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

  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");

  (*fopt)= this->evalfun(thini);
  while ((iter< this->maxiter) & (ferr < this->ftol) & (therr < this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
    fnew= this->evalfun(thopt);
    ferr= (*fopt) - fnew;
    (*fopt)= fnew;
  }
}


//Same but does not evaluate objective function (faster but stopping depends only on change in thopt)
void modselFunction::cda(double *thopt, double *thini) {
  
  int j, iter=0;
  double therr=1, thnew;

  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");
  
  while ((iter< this->maxiter) & (therr < this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
  }
}


//Block CDA jointly updating all parameters (uses updateBlock)
void modselFunction::blockcda(double *thopt, double *fopt, double *thini) {

  int j, iter=0;
  double *thnew, fnew, therr=1, ferr=1;

  if ((this->updateBlock)==NULL) Rf_error("To run CDA you need to specify updateBlock");
  thnew= dvector(0,*nsel);

  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }
  
  while ((iter< this->maxiter) & (ferr < this->ftol) & (therr < this->thtol)) {
    
    (*(this->updateBlock))(thnew, thopt, this->sel, this->nsel, this->pars);
    fnew= this->evalfun(thnew);
    ferr= (*fopt) - fnew;
    if (ferr>0) {
      (*fopt)= fnew;
      for (j=0,therr=0; j< *nsel; j++) {
	therr= max_xy(therr, fabs(thnew[j] - thopt[j]));
	thopt[j]= thnew[j];
      }
    }

  }

  free_dvector(thnew,0,*nsel);

}





//Classical CDA with approx updates given by Newton's method (uses gradhess)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double thcur, therr=1, ferr=1, fnew, delta, g, H;

  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");


  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr < this->ftol) & (therr < this->thtol)) {

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

  } //end while iter

  
}



//Block CDA with Newton method updates (uses gradhess)
void modselFunction::blockcdaNewton(double *thopt, double *fopt, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double *thcur, therr=1, ferr=1, fnew, *delta, *g, *H;
  
  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");
  thcur= dvector(0,*nsel); delta= dvector(0,*nsel); g= dvector(0,*nsel); H= dvector(0,*nsel);
  
  (*fopt)= this->evalfun(thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thcur[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr < this->ftol) & (therr < this->thtol)) {

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

  } //end while iter

  free_dvector(thcur, 0,*nsel); free_dvector(delta, 0,*nsel); free_dvector(g, 0,*nsel); free_dvector(H, 0,*nsel);
  
}


