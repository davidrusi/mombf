#include "modselFunction.h"
using namespace std;

modselFunction::modselFunction(int *sel, int *nsel, struct marginalPars *pars, pt2fun fun=NULL) {

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
  //this->gradhessupdate= NULL;

}


modselFunction::~modselFunction() {

}



//Evaluate fun at th and return the value of funargs
void modselFunction::evalfun(double *f, double *th, std::map<string, double *> *funargs= NULL) {

  fun(f, th, this->sel, this->nsel, this->pars, funargs);

}

//Evaluate fun at new value thjnew by updating its old value f at th[j]. Also update the value of funargs
// Input
//  - thjnew: new value for th[j]
//  - f: value of fun at th
//  - th: current values for th[0], ..., th[*nsel -1]
//  - j: index of parameter th[j] being updated
// Output
//  - fnew: value of fun at thnew= th[0],...,th[j-1],thjnew,th[j+1],...,th[*nsel -1]
// Input/Output
//  - funargs: on input these are arguments needed to evaluate fun at th, at output arguments needed to evaluate fun at thnew
void modselFunction::evalfunupdate(double *fnew, double *thjnew, int j, double *f, double *th, std::map<string, double *> *funargs) {

  funupdate(fnew, thjnew, j, f, th, this->sel, this->nsel, this->pars, funargs);

}




//*******************************************************************************************************************
// OPTIMIZATION ALGORITHMS
//*******************************************************************************************************************

//Coordinate Descent Algorithm (uses updateUniv)
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

  this->evalfun(fopt, thini);
  for (j=0; j< *nsel; j++) thopt[j]= thini[j];

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars, NULL);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
    this->evalfun(&fnew, thopt);
    ferr= (*fopt) - fnew;
    (*fopt)= fnew;
    iter++;
  }
}



//Same but does not evaluate objective function (stopping depends only on change in thopt)
void modselFunction::cda(double *thopt, double *thini) {

  int j, iter=0;
  double therr=1, thnew;

  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");

  for (j=0; j< *nsel; j++) thopt[j]= thini[j];
  while ((iter< this->maxiter) & (therr > this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars, NULL);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      thopt[j]= thnew;
    }
    iter++;
  }
}


void modselFunction::cda(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs) {
  int j, iter=0;
  double therr=1, ferr=1, thnew, fnew;

  if ((this->fun)==NULL) Rf_error("To run CDA you need to specify evalfun");
  if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");

  this->evalfun(fopt, thini, funargs);
  for (j=0; j< *nsel; j++) thopt[j]= thini[j];

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
    for (j=0, therr=0; j< *nsel; j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, this->nsel, this->pars, funargs);
      therr= max_xy(therr, fabs(thnew - thopt[j]));
      evalfunupdate(&fnew,&thnew,j,fopt,thopt,funargs); //Eval fun at thjnew, update funargs
      thopt[j]= thnew;
    }
    ferr= (*fopt) - fnew;
    (*fopt)= fnew;
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

  this->evalfun(fopt,thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0; j< *nsel; j++) { (*(this->updateUniv))(thnew+j, j, thopt, this->sel, this->nsel, this->pars, NULL); }

    this->evalfun(&fnew,thnew);
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





//CDA with approx updates given by Newton's method (uses gradhess and funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double thjnew, thjcur, therr=1, ferr=1, fnew, delta, g, H;

  if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify fun");
  if ((this->funupdate)==NULL) Rf_error("To run cdaNewton you need to specify funupdate");
  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");

  this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0, therr=ferr=0; j< *nsel; j++) {

      gradhessUniv(&g, &H, j, thopt, this->sel, this->nsel, this->pars, funargs);
      delta= g/H;

      nsteps= 1; found= false;
      while (!found & (nsteps<=maxsteps)) {

	thjnew= thopt[j] - delta;
	evalfunupdate(&fnew,&thjnew,j,fopt,thopt,funargs); //Eval fun at thjnew, update funargs

	if (fnew < *fopt) {
	  found= true;
	  ferr+= *fopt - fnew;
	  (*fopt)= fnew;
	  therr= max_xy(therr, fabs(delta));
	  thopt[j]= thjnew;
	} else {
	  delta /= 2.0;
	  nsteps++;
	  thjcur= thopt[j]; thopt[j]= thjnew;
  	  evalfunupdate(fopt,&thjcur,j,&fnew,thopt,funargs); //revert funargs to earlier th
	  thopt[j]= thjcur;
	}

      } //end while !found

    } //end for j
    iter++;

  } //end while iter

}


//CDA with approx updates given by Newton's method (uses gradhess but not funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double thcur, therr=1, ferr=1, fnew, delta, g, H;

  if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify evalfun");
  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");

  this->evalfun(fopt,thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0, therr=ferr=0; j< *nsel; j++) {

      gradhessUniv(&g, &H, j, thopt, this->sel, this->nsel, this->pars, NULL);
      delta= g/H;

      nsteps= 1; found= false;
      thcur= thopt[j];
      while (!found & (nsteps<=maxsteps)) {

	thopt[j] -= delta;
	this->evalfun(&fnew,thopt);

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

  this->evalfun(fopt,thini);
  for (j=0; j< *nsel; j++) { thopt[j]= thcur[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    therr= ferr= 0;
    for (j=0; j< *nsel; j++) {
      gradhessUniv(g+j, H+j, j, thopt, sel, nsel, pars, NULL);
      delta[j]= g[j]/H[j];
    }

    nsteps= 1; found= false;
    while (!found & (nsteps<=maxsteps)) {

      for (j=0; j< *nsel; j++) { thopt[j] -= delta[j]; therr= max_xy(therr, fabs(delta[j])); }
      this->evalfun(&fnew,thopt);

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



