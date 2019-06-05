#include "modselFunction.h"
using namespace std;

modselFunction::modselFunction(int *sel, int thlength, struct marginalPars *pars, pt2fun fun=NULL) {

  this->thlength= thlength;
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

  fun(f, th, this->sel, &(this->thlength), this->pars, funargs);

}

//Evaluate fun at new value thjnew by updating its old value f at th[j]. Also update the value of funargs
// Input
//  - thjnew: new value for th[j]
//  - f: value of fun at th
//  - th: current values for th[0], ..., th[thlength -1]
//  - j: index of parameter th[j] being updated
// Output
//  - fnew: value of fun at thnew= th[0],...,th[j-1],thjnew,th[j+1],...,th[thlength -1]
// Input/Output
//  - funargs: on input these are arguments needed to evaluate fun at th, at output arguments needed to evaluate fun at thnew
void modselFunction::evalfunupdate(double *fnew, double *thjnew, int j, double *f, double *th, std::map<string, double *> *funargs) {

  funupdate(fnew, thjnew, j, f, th, this->sel, &(this->thlength), this->pars, funargs);

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
  for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
    for (j=0, therr=0; j< (this->thlength); j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->pars, NULL);
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

  for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];
  while ((iter< this->maxiter) & (therr > this->thtol)) {
    for (j=0, therr=0; j< (this->thlength); j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->pars, NULL);
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
  for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
    for (j=0, therr=0; j< (this->thlength); j++) {
      (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->pars, funargs);
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
  thnew= dvector(0,this->thlength);

  this->evalfun(fopt,thini);
  for (j=0; j< (this->thlength); j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0; j< this->thlength; j++) { (*(this->updateUniv))(thnew+j, j, thopt, this->sel, &(this->thlength), this->pars, NULL); }

    this->evalfun(&fnew,thnew);
    ferr= (*fopt) - fnew;
    if (ferr>0) {
      (*fopt)= fnew;
      for (j=0,therr=0; j< this->thlength; j++) {
	therr= max_xy(therr, fabs(thnew[j] - thopt[j]));
	thopt[j]= thnew[j];
      }
    }
    iter++;

  }

  free_dvector(thnew,0,this->thlength);

}





//CDA with approx updates given by Newton's method (uses gradhess and funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps=5) {

  bool found;
  int j, iter=0, nsteps;
  double thjnew, thjcur, therr=1, ferr=1, fnew, delta, g, H;

  if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify fun");
  if ((this->funupdate)==NULL) Rf_error("To run cdaNewton you need to specify funupdate");
  if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");

  this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs
  for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0, therr=ferr=0; j< this->thlength; j++) {

      gradhessUniv(&g, &H, j, thopt, this->sel, &(this->thlength), this->pars, funargs);
      if (H>0) { delta= g/H; } else { delta= g/max_xy(-H,.001); }  //if H<0 then target is -def, fix to ensure step is in the direction of -gradient

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

  //Rprintf("nparam= %d, niter=%d\n",this->thlength, iter); //debug
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
  for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    for (j=0, therr=ferr=0; j< this->thlength; j++) {

      gradhessUniv(&g, &H, j, thopt, this->sel, &(this->thlength), this->pars, NULL);
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
void modselFunction::blockcdaNewton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double therr=1, ferr=1, fnew, *delta, *g, *H;

  if ((this->fun)==NULL) Rf_error("To run blockcdaNewton you need to specify evalfun");
  if ((this->gradhessUniv)==NULL) Rf_error("To run blockcdaNewton you need to specify either gradhessUniv");
  delta= dvector(0,this->thlength); g= dvector(0,this->thlength); H= dvector(0,this->thlength);

  this->evalfun(fopt,thini,funargs);
  for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    therr= ferr= 0;
    for (j=0; j< this->thlength; j++) {
      gradhessUniv(g+j, H+j, j, thopt, sel, &(this->thlength), this->pars, funargs);
      delta[j]= g[j]/H[j];
    }

    nsteps= 1; found= false;
    for (j=0; j< this->thlength; j++) { thopt[j] -= delta[j]; therr= max_xy(therr, fabs(delta[j])); }
    while (!found & (nsteps<=maxsteps)) {

      this->evalfun(&fnew,thopt,funargs);

      if (fnew < *fopt) {
	found= true;
	ferr= *fopt - fnew;
	(*fopt)= fnew;
      } else {
	for (j=0; j< this->thlength; j++) { delta[j] /= 2.0; thopt[j] += delta[j]; }
	ferr= 0;
	nsteps++;
      }

    } //end while !found
    iter++;

  } //end while iter

  free_dvector(delta, 0,this->thlength); free_dvector(g, 0,this->thlength); free_dvector(H, 0,this->thlength);

}




/* 

Newton-Raphson optimization (modifying hessian to be +def when needed)

*/
void modselFunction::Newton(double *thopt, double *fopt, double *thini, std::map<string, double *> *funargs, int maxsteps=5) {

  bool posdef;
  int j, iter=0;
  double *thnew, therr=1, ferr=1, fnew, *delta, *g, **H, **Hinv;

  if ((this->fun)==NULL) Rf_error("To run Newton you need to specify fun");
  if ((this->hess)==NULL) Rf_error("To run Newton you need to specify hess");
  if ((this->gradUniv)==NULL) Rf_error("To run Newton you need to specify gradUniv");

  thnew= dvector(1,this->thlength); delta= dvector(1,this->thlength); g= dvector(1,this->thlength);
  H= dmatrix(1,this->thlength,1,this->thlength); Hinv= dmatrix(1,this->thlength,1,this->thlength);
  
  this->evalfun(fopt, thini, funargs); //call evalfun and initialize funargs
  for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }

  while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {

    this->hess(H, thopt, this->sel, &(this->thlength), this->pars, funargs);
    inv_posdef(H, this->thlength, Hinv, &posdef);
    if (!posdef) {
      int i;
      double lmin=0, *vals;
      vals= dvector(1,this->thlength);
      eigenvals(H,this->thlength,vals);
      for (i=1; i<= this->thlength; i++) if (vals[i]<lmin) lmin= vals[i];
      lmin = -lmin + .01;
      for (i=1; i<= this->thlength; i++) H[i][i] += lmin;
      free_dvector(vals,1,this->thlength);
    }
    
    for (j=0; j< this->thlength; j++) { this->gradUniv(g+1+j, j, thopt, this->sel, &(this->thlength), this->pars, funargs); }
    Ax(Hinv,g,delta-1,1,this->thlength,1,this->thlength);

    for (j=0; j< this->thlength; j++) { thnew[j]= thopt[j] - delta[j]; }
    
    this->evalfun(&fnew, thnew, funargs); //call evalfun and initialize funargs

    if (fnew < *fopt) {
      
      for (j=0; j< this->thlength; j++) { therr= max_xy(therr, fabs(thopt[j]-thopt[j])); thopt[j]= thnew[j]; }
      ferr= *fopt - fnew;
      (*fopt)= fnew;

    } else {

      ferr= 0; //causes exit

    }

    iter++;

  }

  free_dvector(thnew, 1,this->thlength); free_dvector(delta,1,this->thlength); free_dvector(g,1,this->thlength);
  free_dmatrix(H, 1,this->thlength,1,this->thlength); free_dmatrix(Hinv, 1,this->thlength,1,this->thlength);
  
}



/*Laplace approximation to int exp(-fun(th)) dth

Input

- thopt: argmin_th fun(th)
- H: hessian of -fun at th=thopt. If not positive definite then +.01 - lmin is added to the diagonal of H, where lmin is the smallest eigenvalue of H

Ouput: logarithm of Laplace approximation

  -fun(thopt) + 0.5 * dim(th) * log(2 pi) - 0.5*log(det(H));

If returnH==false, it is assumed that H contains the pre-computed hessian matrix at the mode
If returnH==true, then H is computed and returned, and its Cholesky decomp cholH is also returned (provided the pointer cholH is non-null)

 */
double modselFunction::laplaceapprox(double *thopt, double *fopt, double **H, double **cholH= NULL, bool returnH=false, std::map<string, double *> *funargs=NULL) {
  bool posdef;
  double ans, logdetH, **mycholH;

  if (returnH) this->hess(H, thopt, this->sel, &(this->thlength), this->pars, funargs);

  if (cholH == NULL) {
    mycholH= dmatrix(1,this->thlength,1,this->thlength);
  } else {
    mycholH = cholH;
  }

  choldc(H,this->thlength,mycholH,&posdef);
  if (!posdef) {
    make_posdef(H,this->thlength);
    choldc(H,this->thlength,mycholH,&posdef);
  }
  
  logdetH= logcholdc_det(mycholH, this->thlength);
  ans= - (*fopt) + 0.5 * (this->thlength) * LOG_M_2PI - 0.5*logdetH;

  if (cholH== NULL) free_dmatrix(mycholH, 1,this->thlength,1,this->thlength);
  return ans;
}


double modselFunction::laplaceapprox(double *thopt, double *fopt, std::map<string, double *> *funargs=NULL) {
  double ans, **H;

  if ((this->hess)==NULL) Rf_error("To run laplaceapprox you need to specify hess");
  H= dmatrix(1,this->thlength,1,this->thlength);

  this->hess(H, thopt, this->sel, &(this->thlength), this->pars, funargs);

  ans= this->laplaceapprox(thopt, fopt, H);

  free_dmatrix(H, 1,this->thlength,1,this->thlength);
  return ans;
}


double modselFunction::laplaceapprox(double *thopt, std::map<string, double *> *funargs=NULL) {
  double ans, fopt;

  if ((this->hess)==NULL) Rf_error("To run laplaceapprox you need to specify hess");

  if (funargs==NULL) {
    this->evalfun(&fopt, thopt);
    ans= this->laplaceapprox(thopt, &fopt);
  } else {
    this->evalfun(&fopt, thopt, funargs);
    ans= this->laplaceapprox(thopt, &fopt, funargs);
  }

  return ans;
}
