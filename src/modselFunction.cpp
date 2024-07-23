#include "modselFunction.h"
using namespace std;




/* EVALUATE LOG-JOINT (LOG-LIKELIHOOD + LOG-PRIOR) AND INITIALIZE funargs

   INPUT

   - th: 0-indexed vector with values at which to evaluate the log-likelihood
   - thlength: length of th
   - sel: sel[0], sel[1] etc indicate the selected variables
   - lm: further parameters needed to evaluate the likelihood and prior, see lmObject in cstat.h

   OUTPUT

   - f: value of -loglikelihood(th) - logprior(th)

   INPUT / OUTPUT

   - funargs: optional input arguments storing model-specific calculations needed to evaluate the log-likelihood or prior, such as determinants of prior covariance matrices, etc. funargs can also have output parameters to store calculations done by logl, such as the linear predictor at th, to speed up subsequent log-likelihood evaluations by fjoint_update

*/
void fjoint(pt2fun logl, pt2fun logprior, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens=0;

  logl(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs

  logprior(&priordens, th, sel, thlength, lm, funargs); //evaluate -log(prior)

  (*f) += priordens;

}


//Update log-joint and funargs due to changing th[j] into thjnew
void fjoint_update(pt2funupdate logl_update, pt2fun logprior, double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens= 0, thtmp;

  logl_update(fnew, thjnew, j, f, th, sel, thlength, lm, funargs); //-loglikelihood at thnew and update funargs["residuals"]

  thtmp= th[j];
  th[j]= *thjnew;
  logprior(&priordens, th, sel, thlength, lm, funargs); //evaluate -log(prior)
  th[j]= thtmp;

  (*fnew) += priordens;

}


/*Minus log-joint gradient and hessian wrt th[j]. The function computes 

      grad= logl_grad + loglprior_grad

      hess= logl_hess + loglprior_hess

   where

   logl_grad: gradient of the negative log-likelihood
   logl_hess: hessian of the negative log-likelihood

   logprior_grad: gradient of the negative log-prior
   logprior_hess: hessian of the negative log-prior
*/
void fjoint_gradhess(pt2gradhessUniv logl_gradhess, pt2gradhessUniv logprior_gradhess, double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
   double gradprior=0, hessprior=0;

   logl_gradhess(grad, hess, j, th, sel, thlength, lm, funargs);  //store loglikelihood gradient and hessian wrt th[j]

   logprior_gradhess(&gradprior, &hessprior, j, th, sel, thlength, lm, funargs); //compute logprior gradient and hessian wrt th[j]

   (*grad) += gradprior;
   (*hess) += hessprior;

}

/* Minus log-joint gradient and hessian wrt th[j]

   The function computes logl_grad + loglprior_grad, where

   logl_grad: gradient of the negative log-likelihood
   logprior_grad: gradient of the negative log-prior
*/
void fjoint_grad(pt2gradUniv logl_grad, pt2gradUniv logprior_grad, double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) { 
   double gradprior=0;

   logl_grad(grad, j, th, sel, thlength, lm, funargs);  //store loglikelihood gradient and hessian wrt th[j]

   logprior_grad(&gradprior, j, th, sel, thlength, lm, funargs); //compute logprior gradient and hessian wrt th[j]

   (*grad) += gradprior;

}


//Hessian matrix of the minus log-joint (- loglikelihood - logprior)
// INPUT
// - logl_hess: function to compute the minus log-likelihood hessian
// - logprior_hess: function that adds to the output of logl_hess the hessian of the minus log-prior
void fjoint_hess(pt2hess logl_hess, pt2hess logprior_hess, double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

   logl_hess(hess, th, sel, thlength, lm, funargs);  //store loglikelihood hessian in hess

   logprior_hess(hess, th, sel, thlength, lm, funargs); //add log-prior to hess

}




//*************************************************************************************
// CLASS modselFunction
//*************************************************************************************


modselFunction::modselFunction(int *sel, int thlength, lmObject *lm, pt2fun fun=NULL) {

  this->thlength= thlength;
  this->sel= sel;
  this->lm= lm;
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

  fun(f, th, this->sel, &(this->thlength), this->lm, funargs);

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

  funupdate(fnew, thjnew, j, f, th, this->sel, &(this->thlength), this->lm, funargs);

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
// - converged: true if convergence criteria satisfied, false otherwise

void modselFunction::cda(double *thopt, double *fopt, bool *converged, double *thini) {

  int j, iter=0;
  double therr=1, ferr=1, thnew, fnew;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run CDA you need to specify evalfun");
    if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");
     
    (*converged)= false;
    this->evalfun(fopt, thini);
    for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
      for (j=0, therr=0; j< (this->thlength); j++) {
        (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->lm, NULL);
        therr= max_xy(therr, fabs(thnew - thopt[j]));
        thopt[j]= thnew;
      }
      this->evalfun(&fnew, thopt);
      ferr= (*fopt) - fnew;
      (*fopt)= fnew;
      iter++;
    }
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;

  } else {

    (*converged)= true;
    this->evalfun(fopt, thini);

  }

}



//Same but does not evaluate objective function (stopping depends only on change in thopt)
void modselFunction::cda(double *thopt, bool *converged, double *thini) {

  int j, iter=0;
  double therr=1, thnew;

  if (this->thlength > 0) {

    if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");
     
    (*converged)= false;
    for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];
    while ((iter< this->maxiter) & (therr > this->thtol)) {
      for (j=0, therr=0; j< (this->thlength); j++) {
        (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->lm, NULL);
        therr= max_xy(therr, fabs(thnew - thopt[j]));
        thopt[j]= thnew;
      }
      iter++;
    }
     
    if (therr < this->thtol) (*converged)= true;

  } else {

    (*converged)= true;

  }
}


void modselFunction::cda(double *thopt, double *fopt, bool *converged, double *thini, std::map<string, double *> *funargs) {
  int j, iter=0;
  double therr=1, ferr=1, thnew, fnew;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run CDA you need to specify evalfun");
    if ((this->updateUniv)==NULL) Rf_error("To run CDA you need to specify updateUniv");
     
    (*converged)= false;
    this->evalfun(fopt, thini, funargs);
    for (j=0; j< (this->thlength); j++) thopt[j]= thini[j];
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
      for (j=0, therr=0; j< (this->thlength); j++) {
        (*(this->updateUniv))(&thnew, j, thopt, this->sel, &(this->thlength), this->lm, funargs);
        therr= max_xy(therr, fabs(thnew - thopt[j]));
        evalfunupdate(&fnew,&thnew,j,fopt,thopt,funargs); //Eval fun at thjnew, update funargs
        thopt[j]= thnew;
      }
      ferr= (*fopt) - fnew;
      (*fopt)= fnew;
      iter++;
    }
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;

  } else {

    (*converged)= true;
    this->evalfun(fopt, thini, funargs);

  }
}


//BLOCK CDA JOINTLY UPDATING ALL PARAMETERS
//In contrast to cda here th[j] is updated without updating first th[0], ..., th[j-1].
//Hence even if CDA were guaranteed to converge blockcda may not, as it cannot be interpreted as a sequence of univariate optimizations
void modselFunction::blockcda(double *thopt, double *fopt, bool *converged, double *thini) {

  int j, iter=0;
  double *thnew, fnew, therr=1, ferr=1;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run blockcda you need to specify evalfun");
    (*converged)= false;
    thnew= dvector(0,this->thlength);
     
    this->evalfun(fopt,thini);
    for (j=0; j< (this->thlength); j++) { thopt[j]= thini[j]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
      for (j=0; j< this->thlength; j++) { (*(this->updateUniv))(thnew+j, j, thopt, this->sel, &(this->thlength), this->lm, NULL); }
     
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
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;
     
    free_dvector(thnew,0,this->thlength);

  } else {

    (*converged)= true;
    this->evalfun(fopt,thini);

  }

}





//CDA with approx updates given by Newton's method (uses gradhess and funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, bool *converged, double *thini, std::map<string, double *> *funargs, int maxsteps=5) {

  bool found;
  int j, iter=0, nsteps;
  double thjnew, thjcur, therr=1, ferr=1, fnew, delta, g, H;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify fun");
    if ((this->funupdate)==NULL) Rf_error("To run cdaNewton you need to specify funupdate");
    if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");
     
    (*converged)= false;
    this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs
    for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
      for (j=0, ferr=therr=0; j< this->thlength; j++) {
     
        gradhessUniv(&g, &H, j, thopt, this->sel, &(this->thlength), this->lm, funargs);
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
            ferr= therr= 1;
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
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;

  } else {

    (*converged)= true;
    this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs

  }

}


//CDA with approx updates given by Newton's method (uses gradhess but not funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
void modselFunction::cdaNewton(double *thopt, double *fopt, bool *converged, double *thini, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double thcur, therr=1, ferr=1, fnew, delta, g, H;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run cdaNewton you need to specify evalfun");
    if ((this->gradhessUniv)==NULL) Rf_error("To run cdaNewton you need to specify either gradhessUniv");
     
    (*converged)= false;
    this->evalfun(fopt,thini);
    for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
      for (j=0, therr=ferr=0; j< this->thlength; j++) {
     
        gradhessUniv(&g, &H, j, thopt, this->sel, &(this->thlength), this->lm, NULL);
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
            ferr= therr= 1;
     	  thopt[j]= thcur;
     	  delta /= 2.0;
     	  nsteps++;
     	}
     
        } //end while !found
     
      } //end for j
      iter++;
     
    } //end while iter
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;

  } else {

    (*converged)= true;
    this->evalfun(fopt,thini);

  }

}



//BLOCK CDA WITH NEWTON METHOD UPDATES (USES gradhess)
//Each parameter is updated to th[j] - 0.5^k grad[j]/hess[j] as in cdaNewton, but here grad[j] and hess[j] are evaluated at the current th prior to updating any th
//In contrast, in cdaNewton grad[j] and hess[j] are evaluated after updating th[0], ..., th[j-1]
void modselFunction::blockcdaNewton(double *thopt, double *fopt, bool *converged, double *thini, std::map<string, double *> *funargs, int maxsteps=1) {

  bool found;
  int j, iter=0, nsteps;
  double therr=1, ferr=1, fnew, *delta, *g, *H;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run blockcdaNewton you need to specify evalfun");
    if ((this->gradhessUniv)==NULL) Rf_error("To run blockcdaNewton you need to specify either gradhessUniv");
    delta= dvector(0,this->thlength); g= dvector(0,this->thlength); H= dvector(0,this->thlength);
     
    (*converged)= false;
    this->evalfun(fopt,thini,funargs);
    for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
      therr= ferr= 0;
      for (j=0; j< this->thlength; j++) {
        gradhessUniv(g+j, H+j, j, thopt, sel, &(this->thlength), this->lm, funargs);
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
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;
     
    free_dvector(delta, 0,this->thlength); free_dvector(g, 0,this->thlength); free_dvector(H, 0,this->thlength);

  } else {

    (*converged)= true;
    this->evalfun(fopt,thini,funargs);

  }

}




/* 

Newton-Raphson optimization (modifying hessian to be +def when needed)

*/
void modselFunction::Newton(double *thopt, double *fopt, bool *converged, double *thini, std::map<string, double *> *funargs, int maxsteps=5) {

  bool posdef, found;
  int j, iter=0, nsteps;
  double *thnew, therr=1, ferr=1, fnew, *delta, *g, **H, **Hinv;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run Newton you need to specify fun");
    if ((this->hess)==NULL) Rf_error("To run Newton you need to specify hess");
    if ((this->gradUniv)==NULL) Rf_error("To run Newton you need to specify gradUniv");
     
    (*converged)= false;
    thnew= dvector(0,this->thlength -1); delta= dvector(1,this->thlength); g= dvector(1,this->thlength);
    H= dmatrix(1,this->thlength,1,this->thlength); Hinv= dmatrix(1,this->thlength,1,this->thlength);
     
    this->evalfun(fopt, thini, funargs); //call evalfun and initialize funargs
    for (j=0; j< this->thlength; j++) { thopt[j]= thini[j]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
      this->hess(H, thopt, this->sel, &(this->thlength), this->lm, funargs);
      inv_posdef(H, this->thlength, Hinv, &posdef);
      if (!posdef) { //if not positive definite, make +def by adding smallest eigenvalue to diagonal
        int i;
        double lmin=0, *vals;
        vals= dvector(1,this->thlength);
        eigenvals(H,this->thlength,vals);
        for (i=1; i<= this->thlength; i++) if (vals[i]<lmin) lmin= vals[i];
        lmin = -lmin + .01;
        for (i=1; i<= this->thlength; i++) H[i][i] += lmin;
        free_dvector(vals,1,this->thlength);
      }
      
      for (j=0; j< this->thlength; j++) { this->gradUniv(g+1+j, j, thopt, this->sel, &(this->thlength), this->lm, funargs); }
      Ax(Hinv,g,delta,1,this->thlength,1,this->thlength);
     
      nsteps= 1; found= false;
      while (!found & (nsteps<=maxsteps)) {
     
     
        for (j=0; j< this->thlength; j++) { thnew[j]= thopt[j] - delta[j+1]; }
         
        this->evalfun(&fnew, thnew, funargs); //call evalfun and initialize funargs
         
        if (fnew < *fopt) {
          
          found= true;
          for (j=0, therr=0; j< this->thlength; j++) { therr= max_xy(therr, fabs(thopt[j]-thnew[j])); thopt[j]= thnew[j]; }
          ferr= *fopt - fnew;
          (*fopt)= fnew;
         
        } else {
         
     	for (j=0; j< this->thlength; j++) { delta[j+1] /= 2.0; }
     	nsteps++;
         
        }
         
        iter++;
     
      }
     
    }
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;
     
    free_dvector(thnew, 0,this->thlength -1); free_dvector(delta,1,this->thlength); free_dvector(g,1,this->thlength);
    free_dmatrix(H, 1,this->thlength,1,this->thlength); free_dmatrix(Hinv, 1,this->thlength,1,this->thlength);

  } else {

    (*converged)= true;
    this->evalfun(fopt, thini, funargs); //call evalfun and initialize funargs

  }
 
}



//Univariate Newton-Raphson on th[j] (uses gradhess and funupdate)
// Each th[j] is updated to th[j] - 0.5^k g[j]/H[j]; where k in {1,...,maxsteps} is the smallest value improving the objective function
//
// Output: 
// - thj: optimal value of th[j], for all other parameters evaluated at thini
// - fopt: value of objective function at thj

void modselFunction::Newtonuniv(double *thj, int j, double *fopt, bool *converged, double *thini, std::map<string, double *> *funargs, int maxsteps=5) {

  bool found;
  int i, iter=0, nsteps;
  double thjnew, thjcur, therr=1, ferr=1, fnew, delta, g, H, *thopt;

  if (this->thlength > 0) {

    if ((this->fun)==NULL) Rf_error("To run Newtonuniv you need to specify fun");
    if ((this->funupdate)==NULL) Rf_error("To run Newtonuniv you need to specify funupdate");
    if ((this->gradhessUniv)==NULL) Rf_error("To run Newtonuniv you need to specify gradhessUniv");
     
    thopt= dvector(0, this->thlength);
     
    (*converged)= false;
    this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs
    for (i=0; i< this->thlength; i++) { thopt[i]= thini[i]; }
     
    while ((iter< this->maxiter) & (ferr > this->ftol) & (therr > this->thtol)) {
     
        gradhessUniv(&g, &H, j, thopt, this->sel, &(this->thlength), this->lm, funargs);
        if (H>0) { delta= g/H; } else { delta= g/max_xy(-H,.001); }  //if H<0 then target is -def, fix to ensure step is in the direction of -gradient
     
        nsteps= 1; found= false;
        while (!found & (nsteps<=maxsteps)) {
     
     	thjnew= thopt[j] - delta;
     	evalfunupdate(&fnew,&thjnew,j,fopt,thopt,funargs); //Eval fun at thjnew, update funargs
     
     	if (fnew < *fopt) {
     	  found= true;
     	  ferr= *fopt - fnew;
     	  (*fopt)= fnew;
     	  therr= fabs(delta);
     	  thopt[j]= thjnew;
     	} else {
     	  delta /= 2.0;
     	  nsteps++;
     	  thjcur= thopt[j]; thopt[j]= thjnew;
     	  evalfunupdate(fopt,&thjcur,j,&fnew,thopt,funargs); //revert funargs to earlier th
     	  thopt[j]= thjcur;
     	}
     
        } //end while !found
     
      iter++;
     
    } //end while iter
     
    (*thj)= thopt[j];
     
    if ((ferr < this->ftol) | (therr < this->thtol)) (*converged)= true;
     
    free_dvector(thopt, 0, this->thlength);

  } else {

    (*converged)= true;
    this->evalfun(fopt,thini,funargs); //eval fun at thini, initialize funargs

  }

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

  if (this->thlength >0) {

    if (returnH) this->hess(H, thopt, this->sel, &(this->thlength), this->lm, funargs);
     
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

  } else {

    logdetH= 0;

  }

  ans= - (*fopt) + 0.5 * (this->thlength) * LOG_M_2PI - 0.5*logdetH;

  if (cholH== NULL) free_dmatrix(mycholH, 1,this->thlength,1,this->thlength);
  return ans;
}


double modselFunction::laplaceapprox(double *thopt, double *fopt, std::map<string, double *> *funargs=NULL) {
  double ans, **H;

  if ((this->hess)==NULL) Rf_error("To run laplaceapprox you need to specify hess");
  H= dmatrix(1,this->thlength,1,this->thlength);

  this->hess(H, thopt, this->sel, &(this->thlength), this->lm, funargs);

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



double modselFunction::BIC(double *fopt, int *n) {
  double ans;

  ans= - (*fopt) - 0.5 * (this->thlength) * log(*n);

  return ans;
}



/* APPROXIMATE LAPLACE APPROXIMATION */


/*Laplace approximation to int exp(-fun(th)) dth, based on quadratic Taylor expansion at th0

Input

- th0: initial parameter value
- f0: value of fun at th=th0
- g0: gradient of fun at th=th0. Indexed g0[1], g0[2],... 
- H0: hessian of fun at th=th0. Indexed H0[1][1], H0[1][2], ... If not positive definite then +.01 - lmin is added to the diagonal of H0, where lmin is its smallest eigenvalue
- returng0: if returng0=true then g0 is computed, else g0 is assumed to be a pre-computed input parameter
- returnH0: if returnH0=true then H0 is computed, else H0 is assumed to be a pre-computed input parameter

Ouput: logarithm of approximate Laplace approximation

  -fun(th0) + 0.5 * dim(th) * log(2 pi) - 0.5*log(det(H0)) + 0.5 g0' H0 g0;

If returnH==false, it is assumed that H contains the pre-computed hessian matrix at the mode

If returnH==true, then H is computed and returned, and its Cholesky decomp cholH and inverse H0inv are also returned (provided the pointers to cholH0 and H0inv are non-null)

 */


double modselFunction::ALA(double *th0, double *f0, double *g0, double **H0, double **cholH0= NULL, double **H0inv=NULL, bool returng0= false, bool returnH0=false, double adjfactor=1, std::map<string, double *> *funargs= NULL) { 

  bool posdef;
  int j;
  double ans, logdetH0, **mycholH0, **myH0inv, g0norm;

  if (this->thlength >0) {
    if (returng0) {
      if ((this->gradUniv) != NULL) {
        for (j=0; j< this->thlength; j++) { this->gradUniv(g0+1+j, j, th0, this->sel, &(this->thlength), this->lm, funargs); }
      } else { 
        double h;
        for (j=0; j< this->thlength; j++) { this->gradhessUniv(g0+1+j, &h, j, th0, this->sel, &(this->thlength), this->lm, funargs); }
      }
    }
     
    if (returnH0) this->hess(H0, th0, this->sel, &(this->thlength), this->lm, funargs);
     
    if (cholH0 == NULL) {
      mycholH0= dmatrix(1,this->thlength,1,this->thlength);
    } else {
      mycholH0 = cholH0;
    }
     
    if (H0inv == NULL) {
      myH0inv= dmatrix(1,this->thlength,1,this->thlength);
    } else {
      myH0inv = H0inv;
    }
     
    choldc(H0,this->thlength,mycholH0,&posdef);
    if (!posdef) {
      make_posdef(H0,this->thlength);
      choldc(H0,this->thlength,mycholH0,&posdef);
    }
     
    logdetH0= logcholdc_det(mycholH0, this->thlength);
    inv_posdef(H0, this->thlength, myH0inv, &posdef, mycholH0, NULL);
     
    g0norm= quadratic_xtAx(g0, myH0inv, 1, this->thlength);
     
    ans= - (*f0) + 0.5 * ((this->thlength) * (LOG_M_2PI - log(adjfactor)) - logdetH0 + g0norm / adjfactor);
     
    if (cholH0== NULL) free_dmatrix(mycholH0, 1,this->thlength,1,this->thlength);
    if (H0inv== NULL) free_dmatrix(myH0inv, 1,this->thlength,1,this->thlength);

  } else {

    ans= - (*f0);

  }

  return ans;
}


double modselFunction::ALA(double *th0, double *f0, double adjfactor=1, std::map<string, double *> *funargs=NULL) {
  int j;
  double ans, *g0, **H0;

  if (((this->gradUniv)==NULL) && ((this->gradhessUniv)==NULL)) Rf_error("To run ALA you need to specify gradUniv or gradhessUniv");
  if ((this->hess)==NULL) Rf_error("To run ALA you need to specify hess");

  g0= dvector(1,this->thlength); H0= dmatrix(1,this->thlength,1,this->thlength);

  if ((this->gradUniv) != NULL) {
    for (j=0; j< this->thlength; j++) { this->gradUniv(g0+1+j, j, th0, this->sel, &(this->thlength), this->lm, funargs); }
  } else {
    double h;
    for (j=0; j< this->thlength; j++) { this->gradhessUniv(g0+1+j, &h, j, th0, this->sel, &(this->thlength), this->lm, funargs); }
  }

  this->hess(H0, th0, this->sel, &(this->thlength), this->lm, funargs);

  ans= this->ALA(th0, f0, g0, H0, NULL, NULL, false, false, adjfactor);

  free_dvector(g0, 1,this->thlength); free_dmatrix(H0, 1,this->thlength,1,this->thlength);
  return ans;
}


double modselFunction::ALA(double *th0, double adjfactor=1, std::map<string, double *> *funargs=NULL) {
  double ans, f0;

  if (funargs==NULL) {
    this->evalfun(&f0, th0);
    ans= this->ALA(th0, &f0, adjfactor);
  } else {
    this->evalfun(&f0, th0, funargs);
    ans= this->ALA(th0, &f0, adjfactor, funargs);
  }

  return ans;
}
