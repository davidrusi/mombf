#include "glm.h"


double marginal_glm(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= marginal_glm(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return(ans);
}

double marginal_glm(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood for GLMs

   */

  std::map<string, double *> funargs;
  bool orthoapprox=false, nonlocal, momsingle, momgroup, converged;
  int i, nselgroupsint, cholSsize, thlength= *nsel, n= *(lm->n), priorcode= lm->priorcode, optimMethod= *(lm->optimMethod), family= *(lm->family);
  double ans, *linpred, *ytlinpred, *ypred, nselgroups, *nvarinselgroups, *firstingroup, *selgroups, *ldetSinv, *cholSini, *cholSinv, *Sinv, *thini, *thopt, fini, fopt, *y, **H, **Hinv, **cholH;
  modselFunction *msfun;
  pt2fun fjoint=NULL, fjoint0=NULL;
  pt2funupdate fjoint_update=NULL; 
  pt2gradUniv fjoint_grad=NULL;
  pt2gradhessUniv fjoint_gradhess=NULL, fjoint_gradhess0=NULL;
  pt2hess fjoint_hess=NULL, fjoint_hess0=NULL;

  set_logjoint_glm(&fjoint, &fjoint_update, &fjoint_grad, &fjoint_gradhess, &fjoint_hess, &fjoint0, &fjoint_gradhess0, &fjoint_hess0, &orthoapprox, &nonlocal, &momsingle, &momgroup, &family, priorcode, lm->method);

  //Create object of class modselFunction
  msfun= new modselFunction(sel, thlength, lm, NULL);
  msfun->ftol= 0.001; msfun->thtol= 0.001;
  if (*(lm->optim_maxit) >= 0) msfun->maxiter= *(lm->optim_maxit);

  //Allocate memory
  H= dmatrix(1,thlength,1,thlength); Hinv= dmatrix(1,thlength,1,thlength); cholH= dmatrix(1,thlength,1,thlength);
  thopt= dvector(0, *nsel); thini= dvector(0, *nsel); 
  y= lm->y;

  //Initialize static elements in funargs (not changed by msfun)
  // Pre-compute prior covariances (determinant, Cholesky decomp) and book-keeping (number of selected groups, number of variables in each group etc)

  nvarinselgroups= dvector(0, min_xy(*nsel, *(lm->ngroups))); firstingroup= dvector(0, min_xy(*nsel, *(lm->ngroups))); selgroups= dvector(0, *nsel -1);
  findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, lm->nvaringroup, lm->ngroups); //copy subset of nvaringroup into nvarinselgroups
  funargs["nvarinselgroups"]= nvarinselgroups;
  funargs["firstingroup"]= firstingroup;
  funargs["nselgroups"]= &nselgroups;
  funargs["selgroups"]= selgroups;
  nselgroupsint= (int) (nselgroups +.1);
     
  ldetSinv= dvector(0, nselgroupsint); cholSini= dvector(0, nselgroupsint);   //Obtain Cholesky decomp and determinant of prior scale covariances for each group
  cholSini_indexes(cholSini, &cholSsize, nselgroupsint, nvarinselgroups);
  cholSinv= dvector(0, cholSsize); Sinv= dvector(0, cholSsize);
     
  funargs["cholSini"]= cholSini; //cholSini[j] is the index in cholSinv at which Sinv_j starts
  if (priorcode != 100) gzell_Sinv_byprior(Sinv, cholSinv, ldetSinv, &nselgroupsint, nvarinselgroups, sel, cholSini, lm->XtX, lm->n, lm->tau, lm->taugroup, priorcode); //priorcode==100 is the BIC, no prior hessians to compute
  funargs["ldetSinv"]= ldetSinv; funargs["cholSinv"]= cholSinv; funargs["Sinv"]= Sinv;


  //Initialize dynamic elements in funargs (changed by msfun)
  linpred= dvector(0, n); ypred= dvector(0, n); ytlinpred= dvector(0, thlength);
  funargs["linpred"]= linpred; funargs["ypred"]= ypred; funargs["ytlinpred"]= ytlinpred;

  //Initialize posterior mode and funargs
  if (*(lm->usethinit) == 3) {
    for (i=0; i< thlength; i++) thini[i]= (lm->thinit)[sel[i]];
  } else {
    for (i=0; i< thlength; i++) thini[i]= 0;
  }

  
  if (*(lm->method) == 2) {
    // APPROXIMATE LAPLACE APPROXIMATION AT theta=0

    msfun->fun= fjoint0; msfun->gradhessUniv= fjoint_gradhess0; msfun->hess= fjoint_hess0;

    double adjfactor=1, ybar= *(lm->sumy) / ((double) n);

    if ((family == 21) || (family == 22)) {
      if (*(lm->adjoverdisp) == 1) {  //estimate over-dispersion from intercept-only model

        adjfactor= *(lm->sumy2) / ((double) n)  - pow(ybar, 2);
        if (family == 21) { 
          adjfactor /= ybar * (1 - ybar);
        } else if (family == 22) {
          adjfactor /= ybar;
        }

      } else if (*(lm->adjoverdisp) == 2) { //estimate over-dispersion from residuals

        Rf_error("This over-dispersion adjustment method is not implemented yet\n");
        double ss=0;

        if (*(lm->usethinit) != 3) {
          if (family==21) {
            get_thini_glm(thopt, thini, H, Hinv, negloglgradhess00_logreg, negloglhess00_logreg, sel, &thlength, nonlocal, orthoapprox, &funargs, lm);
          } else if (family==22) {
            get_thini_glm(thopt, thini, H, Hinv, negloglgradhess00_poisson, negloglhess00_poisson, sel, &thlength, nonlocal, orthoapprox, &funargs, lm);
          }
        }
        
        Aselvecx(lm->x, thopt, ypred, 0, n-1, sel, &thlength); //Returns ypred=x[,sel] %*% thini

        if (family == 21) { 
          for (i=0; i<n; i++) { ss += pow( y[i] - 1.0/(1.0 + exp(-ypred[i])), 2); } 
          adjfactor= 4; //1/0.25
        } else if (family == 22) {
          for (i=0; i<n; i++) { ss += pow( y[i] - exp(ypred[i]), 2); } 
          adjfactor= 1;
        }
        adjfactor *= ss / ((double) n);
      }
    }

    ans= msfun->ALA(thini, adjfactor, &funargs);


  } else {
    // LAPLACE APPROXIMATION or APPROXIMATE LAPLACE APPROX AT theta != 0

    msfun->fun = fjoint; msfun->funupdate = fjoint_update; msfun->gradUniv = fjoint_grad; msfun->gradhessUniv = fjoint_gradhess; msfun->hess = fjoint_hess;

    msfun->evalfun(&fini, thini, &funargs); //initialize funargs

    if (*(lm->usethinit) != 3) get_thini_glm(thini, thini, H, Hinv, fjoint_gradhess0, fjoint_hess0, sel, &thlength, nonlocal, orthoapprox, &funargs, lm);
   
    if (nonlocal && !orthoapprox) {   //if it's a non-local prior, avoid exact zeroes (0 prior density)
     
      for (i=0; i< *nsel; i++) {
        if (fabs(thini[i]) < 1.0e-5) {
          double fminus, fplus;
          thini[i]= -1.0e-5; msfun->evalfun(&fminus, thini, &funargs);
          thini[i]=  1.0e-5; msfun->evalfun(&fplus, thini, &funargs);;
          if (fminus<=fplus) { thini[i]= -1.0e-5; } else { thini[i]= 1.0e-5; }
        }
      }
     
    }
     
    //Obtain posterior mode and Laplace approx
    if (((optimMethod==0) && (*nsel >=15)) || (optimMethod==2)) {
      msfun->cdaNewton(thopt, &fopt, &converged, thini, &funargs, 10);
    } else {
      msfun->Newton(thopt, &fopt, &converged, thini, &funargs, 10);
    }

    //Deprecated: Laplace approx equivalent to ALA, upon convergence. Use ALA in case convergence failed (then ALA is often more accurate than Laplace)
    //ans= msfun->laplaceapprox(thopt, &fopt, H, cholH, true, &funargs); //Laplace approx (also returns H and cholH)

    if (priorcode != 100) {

      double *g;
      g= dvector(1, thlength);
      ans= msfun->ALA(thopt, &fopt, g, H, cholH, Hinv, true, true, 1, &funargs); //ALA (also returns g, H, cholH and Hinv)
      free_dvector(g, 1, thlength);

    } else {

      ans= msfun->BIC(&fopt, lm->n);

    }

  }

  //For MOM priors, add log-penalty to log-marginal likelihood
  if ((momsingle | momgroup) & orthoapprox) {
    double pen;

    if (momgroup) {
      pen = gmompenalty_approx(momsingle, momgroup, thopt, Hinv, Sinv, exp(thopt[*sel]), thlength, *nsel, nselgroupsint, nvarinselgroups, firstingroup, cholSini);
    } else {
      pen = pmompenalty_approx(thopt, Hinv, lm->tau, nselgroupsint, nvarinselgroups, firstingroup);
    }
    ans += pen;

  }

  //Free memory and return output
  delete msfun;
  free_dmatrix(H,1,thlength,1,thlength); free_dmatrix(Hinv,1,thlength,1,thlength); free_dmatrix(cholH,1,thlength,1,thlength);
  free_dvector(thopt,0,*nsel); free_dvector(thini,0,*nsel);

  free_dvector(nvarinselgroups, 0, min_xy(*nsel, *(lm->ngroups))); free_dvector(firstingroup, 0, min_xy(*nsel, *(lm->ngroups))); free_dvector(selgroups,0,*nsel -1);
  free_dvector(ldetSinv, 0, nselgroupsint); free_dvector(cholSini, 0, nselgroupsint);
  free_dvector(cholSinv, 0, cholSsize); free_dvector(Sinv, 0, cholSsize);

  free_dvector(linpred, 0, n);free_dvector(ypred, 0, n); free_dvector(ytlinpred, 0, thlength);

  return ans;

}



/* Initialize parameter value in GLMs to thini= H0^{-1} g0, where g0 and H0 are the gradient and hessian at 0 of the log-joint (for local priors) or the log-likelihood (for non-local priors)

  - thini: initial parameter value

  OUPUT

  - thhat: thini - Hinv g. You can set thhat=thini, then both get updated on output
  - H: hessian at 0
  - Hinv: inverse of H

*/
void get_thini_glm(double *thhat, double *thini, double **H, double **Hinv, pt2gradhessUniv fjoint_gradhess0, pt2hess fjoint_hess0, int *sel, int *thlength, bool nonlocal, bool orthoapprox, std::map<string, double *> *funargs, lmObject *lm) {
  bool posdef;
  int i;
  double *g, *h;

  g= dvector(1,*thlength); h= dvector(1,*thlength); 

  if (!nonlocal || orthoapprox) { //if it's a local prior, evaluate log-posterior gradient and hessian
     
    fjoint_hess0(H, thini, sel, thlength, lm, funargs);
    for (i=0; i< *thlength; i++) { fjoint_gradhess0(g+1+i, h+1+i, i, thini, sel, thlength, lm, funargs); g[i+1]= -g[i+1]; }
     
  } else { //if it's a non-local prior, evaluate log-likelihood gradient and hessian
     
    pt2fun logl=NULL, logl0=NULL;
    pt2funupdate logl_update=NULL; 
    pt2gradUniv logl_grad=NULL;
    pt2gradhessUniv logl_gradhess=NULL, logl_gradhess0=NULL;
    pt2hess logl_hess=NULL, logl_hess0=NULL;
     
    set_logl_glm(&logl, &logl_update, &logl_grad, &logl_gradhess, &logl_hess, &logl0, &logl_gradhess0, &logl_hess0, lm->family);
    fjoint_hess0(H, thini, sel, thlength, lm, funargs);
    for (i=0; i< *thlength; i++) { fjoint_gradhess0(g+1+i, h+1+i, i, thini, sel, thlength, lm, funargs); g[i+1]= -g[i+1]; }
     
  }
     
  inv_posdef(H, *thlength, Hinv, &posdef);   //initialize posterior mode with single Newton step
  Ax(Hinv,g,thhat-1,1,*thlength,1,*thlength);

  free_dvector(g,1,*thlength); free_dvector(h,1,*thlength); 

}


//*************************************************************************************
//  Generic function returning a pointer to log-joint (or log-likelihood for BIC), gradient and hessians for a given priorcode
//*************************************************************************************

/*Given a likelihood family and a priorcode (as returned by mspriorCode) and whether an ALA approximation is desired for non-local priors (orthoapprox)

  Return 
  - fjoint, fjoint_update, fjoint_gradhess, fjoint_hess: pointers to functions evaluating the negative log-joint density, updating its value, gradient and hessian
  - nonlocal: boolean indicating whether the prior has any non-local component (either in single coefficients or in grouped coefficients0
  - momsingle: boolean variable indicating whether a pMOM or groupMOM prior was set single coefficients 
  - momgroup: boolean variable indicating whether a groupMOM prior was set on groups of coefficients
  - orthoapprox: set to true if non-local contribution to the integrated likelihood is to be approximated via ALA, hence the prior must be set to the underlying Normal prior

  Priorcodes
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
  // 100: BIC (rather than computing the integrated likelihood, return -0.5 BIC)

*/
void set_logjoint_glm(pt2fun *fjoint, pt2funupdate *fjoint_update, pt2gradUniv *fjoint_grad, pt2gradhessUniv *fjoint_gradhess, pt2hess *fjoint_hess, pt2fun *fjoint0, pt2gradhessUniv *fjoint_gradhess0, pt2hess *fjoint_hess0, bool *orthoapprox, bool *nonlocal, bool *momsingle, bool *momgroup, int *family, int priorcode, int *method) {


  if (priorcode == 100) {  // If no prior was specified (e.g. BIC), return pointers to log-likelihood and derivatives

    set_logl_glm(fjoint, fjoint_update, fjoint_grad, fjoint_gradhess, fjoint_hess, fjoint0, fjoint_gradhess0, fjoint_hess0, family);
    (*momsingle)= (*momgroup)= (*orthoapprox)= (*nonlocal)= false;

  } else {  //If a prior was specified, return pointers to log-joint and derivatives

    //Record if pMOM or groupMOM on single coef was set
    (*momsingle)= ((priorcode==0) || (priorcode==5) || (priorcode==10) || (priorcode==13) || (priorcode==50) || (priorcode==53)); 
    //Record if groupMOM on groups of coef was set
    (*momgroup)= ((priorcode==10) || (priorcode)==50); 
     
    //Record if the prior has any non-local component
    (*nonlocal)= (priorcode == 0) || (priorcode == 1) || (priorcode == 2) || (priorcode == 5) || (priorcode == 10) || (priorcode == 13) || (priorcode == 32) || (priorcode == 33) || (priorcode == 50) || (priorcode == 53);
     
    //Record if orthogonal approx to posterior expected MOM penalty should be used
    if (*momsingle || *momgroup) {
      if ((*method ==2) || (*method == -1))  { //For MOM priors approx the mean of products via the product of means
        (*orthoapprox)= true; 
      } else if (*method ==0) { //If Laplace approx requested, do not approx the mean of product
        (*orthoapprox)= false;
      } else {
        Rf_error("For GLMs and MOM priors only method='auto', 'Laplace' and 'ALA' implemented\n");
      }
    }
     
     
    if (*family ==21) {  //logistic regression
     
      if ((*momsingle) && !(*momgroup) && !(*orthoapprox)) { //MOM prior on indiv coef, orthogonal approx not wanted
     
        (*fjoint)= &fjoint_logreg_pmomgzell; (*fjoint_update)= &fjointu_logreg_pmomgzell; 
        (*fjoint_grad)= &fjointg_logreg_pmomgzell; (*fjoint_gradhess)= &fjointgh_logreg_pmomgzell; (*fjoint_hess)= &fjointh_logreg_pmomgzell;
        (*fjoint0)= &neglogl0_logreg; (*fjoint_gradhess0)= &negloglgradhess0_logreg; (*fjoint_hess0)= &negloglhess0_logreg;
     
      } else   if ((priorcode == 1) || (priorcode == 33)) { //eMOM-eMOM or eMOM-gZellner priors
     
        (*fjoint)= &fjoint_logreg_pemomgzell; (*fjoint_update)= &fjointu_logreg_pemomgzell; 
        (*fjoint_grad)= &fjointg_logreg_pemomgzell; (*fjoint_gradhess)= &fjointgh_logreg_pemomgzell; (*fjoint_hess)= &fjointh_logreg_pemomgzell;
        (*fjoint0)= &neglogl0_logreg; (*fjoint_gradhess0)= &negloglgradhess0_logreg; (*fjoint_hess0)= &negloglhess0_logreg;
     
      } else   if ((priorcode == 9) || (priorcode== 63) || ( (*momsingle || *momgroup) && *orthoapprox)) { //Zellner or MOM priors with orthogonal approx
     
        (*fjoint)= &fjoint_logreg_gzellgzell; (*fjoint_update)= &fjointu_logreg_gzellgzell; 
        (*fjoint_grad)= &fjointg_logreg_gzellgzell; (*fjoint_gradhess)= &fjointgh_logreg_gzellgzell; (*fjoint_hess)= &fjointh_logreg_gzellgzell;
        (*fjoint0)= &fjoint0_logreg_gzellgzell; (*fjoint_gradhess0)= &fjointgh0_logreg_gzellgzell; (*fjoint_hess0)= &fjointh0_logreg_gzellgzell;
     
      } else {
        Rf_error("The specified method to obtain the integrated likelihood is not implemented in GLMs for this prior");
      }
     
    } else if (*family ==22) {  //Poisson regression
     
      if ((*momsingle) && !(*momgroup) && !(*orthoapprox)) { //MOM prior on indiv coef, orthogonal approx not wanted
     
        (*fjoint)= &fjoint_poisson_pmomgzell; (*fjoint_update)= &fjointu_poisson_pmomgzell; 
        (*fjoint_grad)= &fjointg_poisson_pmomgzell; (*fjoint_gradhess)= &fjointgh_poisson_pmomgzell; (*fjoint_hess)= &fjointh_poisson_pmomgzell;
        (*fjoint0)= &neglogl0_poisson; (*fjoint_gradhess0)= &negloglgradhess0_poisson; (*fjoint_hess0)= &negloglhess0_poisson;
     
      } else   if ((priorcode == 1) || (priorcode == 33)) { //eMOM-eMOM or eMOM-gZellner priors
     
        (*fjoint)= &fjoint_poisson_pemomgzell; (*fjoint_update)= &fjointu_poisson_pemomgzell; 
        (*fjoint_grad)= &fjointg_poisson_pemomgzell; (*fjoint_gradhess)= &fjointgh_poisson_pemomgzell; (*fjoint_hess)= &fjointh_poisson_pemomgzell;
        (*fjoint0)= &neglogl0_poisson; (*fjoint_gradhess0)= &negloglgradhess0_poisson; (*fjoint_hess0)= &negloglhess0_poisson;
     
      } else   if ((priorcode == 9) || (priorcode== 63) || ( (*momsingle || *momgroup) && *orthoapprox)) { //Zellner or MOM priors with orthogonal approx
     
        (*fjoint)= &fjoint_poisson_gzellgzell; (*fjoint_update)= &fjointu_poisson_gzellgzell; 
        (*fjoint_grad)= &fjointg_poisson_gzellgzell; (*fjoint_gradhess)= &fjointgh_poisson_gzellgzell; (*fjoint_hess)= &fjointh_poisson_gzellgzell;
        (*fjoint0)= &fjoint0_poisson_gzellgzell; (*fjoint_gradhess0)= &fjointgh0_poisson_gzellgzell; (*fjoint_hess0)= &fjointh0_poisson_gzellgzell;
     
      } else {
        Rf_error("The specified method to obtain the integrated likelihood is not implemented in GLMs for this prior");
      }
     
    } else {
     
      Rf_error("This likelihood family is not implemented");
     
    }

  }

}



void set_logl_glm(pt2fun *logl, pt2funupdate *logl_update, pt2gradUniv *logl_grad, pt2gradhessUniv *logl_gradhess, pt2hess *logl_hess, pt2fun *logl0, pt2gradhessUniv *logl_gradhess0, pt2hess *logl_hess0, int *family) {

  if (*family ==21) {  //logistic regression

      (*logl)= &neglogl_logreg; (*logl_update)= &negloglupdate_logreg; 
      (*logl_grad)= &negloglgrad_logreg; (*logl_gradhess)= &negloglgradhess_logreg; (*logl_hess)= &negloglhess_logreg;
      (*logl0)= &neglogl0_logreg; (*logl_gradhess0)= &negloglgradhess0_logreg; (*logl_hess0)= &negloglhess0_logreg;

  } else if (*family ==22) {  //Poisson regression

      //Uncomment after adding Poisson log-likelihood and derivatives
      (*logl)= &neglogl_poisson; (*logl_update)= &negloglupdate_poisson; 
      (*logl_grad)= &negloglgrad_poisson; (*logl_gradhess)= &negloglgradhess_poisson; (*logl_hess)= &negloglhess_poisson;
      (*logl0)= &neglogl0_poisson; (*logl_gradhess0)= &negloglgradhess0_poisson; (*logl_hess0)= &negloglhess0_poisson;

  } else {

    Rf_error("This likelihood family is not implemented");

  }


}




//*************************************************************************************
// LOGISTIC REGRESSION
//*************************************************************************************


//Function to evaluate the negative loglikelihood, and update funargs
void neglogl_logreg(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 
  int i, n= *(lm->n), nvars= *thlength;
  double *ypred, *linpred, *ytlinpred, *ytX= lm->ytX, sumlog=0;

  ypred= (*funargs)["ypred"];
  linpred= (*funargs)["linpred"];
  ytlinpred= (*funargs)["ytlinpred"];
  (*ytlinpred)= 0;

  if (*thlength >0) {

    for (i=0; i< nvars; i++) { (*ytlinpred) += ytX[sel[i]] * th[i]; }
    Aselvecx(lm->x, th, linpred, 0, n-1, sel, &nvars); //Returns linpred= x[,sel] %*% th
    for (i=0; i< n; i++) {
      sumlog += log(1.0 + exp(linpred[i])); 
      ypred[i]= 1.0 / (1.0 + exp(-linpred[i]));
    }

    (*f) = -(*ytlinpred) + sumlog;

  } else {

    (*ytlinpred)= 0;
    for (i=0; i< n; i++) { linpred[i]= 0; ypred[i]= 0.5; }

    (*f) = - n * log(.5);
    //neglogl0_logreg(f, th, sel, thlength, lm, funargs);

  }

}


//Update the negative loglikelihood due to changing th[j] into thjnew, and update funargs
void negloglupdate_logreg(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  int i, n= *(lm->n), idxj;
  double *linpred, *ypred, *ytlinpred, *ytX= lm->ytX, *x= lm->x, sumlog=0, thdif;

  linpred= (*funargs)["linpred"];
  ypred= (*funargs)["ypred"];
  ytlinpred= (*funargs)["ytlinpred"];

  if (*thlength >0) {

    //Update inner product of ytX and linear predictor
    idxj= n * sel[j];
    thdif= *thjnew - th[j];
    (*ytlinpred) += ytX[sel[j]] * thdif; 

    //Update linear predictor
    for (i=0; i< n; i++) { 
      linpred[i] += x[i + idxj] * thdif;
      ypred[i]= 1.0 / (1.0 + exp(-linpred[i]));
      sumlog += log(1.0 + exp(linpred[i])); 
    }

    (*fnew) = -(*ytlinpred) + sumlog;

  } else {

    (*ytlinpred)= 0;
    for (i=0; i< n; i++) { linpred[i]= 0; ypred[i]= 0.5; }

    (*f) = - n * log(.5);
    //neglogl0_logreg(fnew, th, sel, thlength, lm, funargs);

  }

}


//Obtain gradient and hessian wrt j
void negloglgradhess_logreg(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, n= *(lm->n), idxj;
  double *ypred, *ytX= lm->ytX, *x= lm->x;

  ypred= (*funargs)["ypred"];

  idxj= *(lm->n) * sel[j];

  (*grad)= -ytX[sel[j]];
  (*hess)= 0;
  for (i=0; i<n; i++) {

    (*grad) += ypred[i] * x[idxj + i];
    (*hess) += ypred[i] * (1-ypred[i]) * x[idxj + i] * x[idxj + i];

  }

}


void negloglgrad_logreg(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, n= *(lm->n), idxj;
  double *ypred, *ytX= lm->ytX, *x= lm->x;

  ypred= (*funargs)["ypred"];

  idxj= *(lm->n) * sel[j];

  (*grad)= -ytX[sel[j]];
  for (i=0; i<n; i++) { 
    (*grad) += ypred[i] * x[idxj + i]; 
  }

}


//Obtain hessian
void negloglhess_logreg(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, k, n= *(lm->n), nvars= *thlength, idxj, idxk;
  double *ypred, *x= lm->x;

  //linpred= (*funargs)["linpred"];
  ypred= (*funargs)["ypred"];
  //ytlinpred= (*funargs)["ytlinpred"];

  for (j=1; j<=nvars; j++) {
    idxj= *(lm->n) * sel[j-1];

    hess[j][j]= 0;
    for (i=0; i<n; i++) { hess[j][j] += ypred[i] * (1-ypred[i]) * x[idxj + i] * x[idxj + i]; }

    for (k=1; k<j; k++) {
      idxk= *(lm->n) * sel[k-1];

      hess[j][k]= 0;
      for (i=0; i<n; i++) { hess[j][k] += ypred[i] * (1-ypred[i]) * x[idxj + i] * x[idxk + i]; }
      hess[k][j]= hess[j][k];
    }

  }

}


/* EVALUATE LOG-LIKELIHOOD AND DERIVATIVES FOR BETA= MLE AT THE INTERCEPT-ONLY MODEL */

//Function to evaluate the negative loglikelihood
void neglogl0_logreg(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 

  double meany= (*(lm->sumy)) / (*(lm->n) + .0);

  (*f) = - (*(lm->n) + .0) * (meany * log(meany/(1-meany)) + log(1-meany));

}


//Obtain gradient and hessian wrt j
void negloglgradhess0_logreg(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  double *colsumsx= lm->colsumsx, *ytX= lm->ytX, meany= (*(lm->sumy)) / (*(lm->n) + .0);

  (*grad)= -ytX[sel[j]] + meany * colsumsx[sel[j]];

  (*hess)= meany * (1-meany) * (lm->XtX)->at(sel[j],sel[j]);

}

//Obtain hessian
void negloglhess0_logreg(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, nvars= *thlength; 
  double meany= (*(lm->sumy)) / (*(lm->n) + .0), vary= meany * (1-meany);

  for (i=1; i<=nvars; i++) {
    hess[i][i]= vary * (lm->XtX)->at(sel[i-1],sel[i-1]);
    for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= vary * (lm->XtX)->at(sel[i-1],sel[j-1]); }
  }

}


/* LOG-LIKELIHOOD AND DERIVATIVES AT BETA=0 */
void neglogl00_logreg(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 

  (*f) = (*(lm->n) + .0) * log(2.0);
}


void negloglgradhess00_logreg(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  double *colsumsx= lm->colsumsx, *ytX= lm->ytX;

  (*grad)= -ytX[sel[j]] + 0.5 * colsumsx[sel[j]];
  (*hess)= 0.25 * (lm->XtX)->at(sel[j],sel[j]);

}


void negloglhess00_logreg(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, nvars= *thlength;

  for (i=1; i<=nvars; i++) {
    hess[i][i]= 0.25 * (lm->XtX)->at(sel[i-1],sel[i-1]);
    for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= 0.25 * (lm->XtX)->at(sel[i-1],sel[j-1]); }
  }

}



//AUXILIARY FUNCTIONS TO EVALUATE LOG-JOINT AND DERIVATIVES. THESE CAN BE COPY/PASTED WITH MINIMAL ADAPTATION, IF ONE WISHES TO IMPLEMENT OTHER LIKELIHOODS (SEE POISSON EXAMPLE BELOW)

//Log joint under pMOM - gZellner prior
void fjoint_logreg_pmomgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_logreg, &pmomgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_logreg_pmomgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_logreg, &pmomgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_logreg_pmomgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_logreg, &pmomgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_logreg_pmomgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_logreg, &pmomgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_logreg_pmomgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_logreg, &pmomgzell_hess, hess, th, sel, thlength, lm, funargs);
}


//Log joint under peMOM - gZellner prior
void fjoint_logreg_pemomgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_logreg, &pemomgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_logreg_pemomgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_logreg, &pemomgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_logreg_pemomgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_logreg, &pemomgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_logreg_pemomgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_logreg, &pemomgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_logreg_pemomgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_logreg, &pemomgzell_hess, hess, th, sel, thlength, lm, funargs);
}


//Log joint under gZellner - gZellner prior
void fjoint_logreg_gzellgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_logreg, &gzellgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_logreg_gzellgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_logreg, &gzellgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_logreg_gzellgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_logreg, &gzellgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_logreg_gzellgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_logreg, &gzellgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_logreg_gzellgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_logreg, &gzellgzell_hess, hess, th, sel, thlength, lm, funargs);
}

void fjoint0_logreg_gzellgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl0_logreg, &gzellgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointgh0_logreg_gzellgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess0_logreg, &gzellgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh0_logreg_gzellgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess0_logreg, &gzellgzell_hess, hess, th, sel, thlength, lm, funargs);
}





//*************************************************************************************
// POISSON REGRESSION
//*************************************************************************************


//Function to evaluate the negative loglikelihood, and update funargs
void neglogl_poisson(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 
  int i, n= *(lm->n), nvars= *thlength;
  double *ypred, *linpred, *ytlinpred, *ytX= lm->ytX, sumypred=0, *sumlogyfact= lm->sumlogyfact;

  ypred= (*funargs)["ypred"];
  linpred= (*funargs)["linpred"];
  ytlinpred= (*funargs)["ytlinpred"];
  (*ytlinpred)= 0;

  if (*thlength >0) {

    for (i=0; i< nvars; i++) { (*ytlinpred) += ytX[sel[i]] * th[i]; }
    Aselvecx(lm->x, th, linpred, 0, n-1, sel, &nvars); //Returns linpred= x[,sel] %*% th
    for (i=0; i< n; i++) {
      ypred[i]= exp(linpred[i]);
      sumypred += ypred[i];
    }

    (*f) = -(*ytlinpred) + sumypred + *sumlogyfact;

  } else {

    (*ytlinpred)= 0;
    for (i=0; i< n; i++) { linpred[i]= 0; ypred[i]= 1; }
    neglogl0_poisson(f, th, sel, thlength, lm, funargs);

  }

}



//Update the negative loglikelihood due to changing th[j] into thjnew, and update funargs
void negloglupdate_poisson(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  int i, n= *(lm->n), idxj;
  double *linpred, *ypred, *ytlinpred, *ytX= lm->ytX, *x= lm->x, sumypred=0, thdif, *sumlogyfact= lm->sumlogyfact;

  linpred= (*funargs)["linpred"];
  ypred= (*funargs)["ypred"];
  ytlinpred= (*funargs)["ytlinpred"];

  if (*thlength >0) {

    //Update inner product of ytX and linear predictor
    idxj= n * sel[j];
    thdif= *thjnew - th[j];
    (*ytlinpred) += ytX[sel[j]] * thdif; 

    //Update linear predictor
    for (i=0; i< n; i++) { 
      linpred[i] += x[i + idxj] * thdif;
      ypred[i]= exp(linpred[i]);
      sumypred += ypred[i];
    }

    (*fnew) = -(*ytlinpred) + sumypred + *sumlogyfact;

  } else {

    (*ytlinpred)= 0;
    for (i=0; i< n; i++) { linpred[i]= 0; ypred[i]= 1; }
    neglogl0_poisson(fnew, th, sel, thlength, lm, funargs);

  }

}


//Obtain gradient and hessian wrt j
void negloglgradhess_poisson(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, n= *(lm->n), idxj;
  double *ypred, *ytX= lm->ytX, *x= lm->x;

  ypred= (*funargs)["ypred"];

  idxj= *(lm->n) * sel[j];

  (*grad)= -ytX[sel[j]];
  (*hess)= 0;
  for (i=0; i<n; i++) {

    (*grad) += ypred[i] * x[idxj + i];
    (*hess) += ypred[i] * x[idxj + i] * x[idxj + i];

  }

}


void negloglgrad_poisson(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, n= *(lm->n), idxj;
  double *ypred, *ytX= lm->ytX, *x= lm->x;

  ypred= (*funargs)["ypred"];

  idxj= *(lm->n) * sel[j];

  (*grad)= -ytX[sel[j]];
  for (i=0; i<n; i++) { 
    (*grad) += ypred[i] * x[idxj + i]; 
  }

}


//Obtain hessian
void negloglhess_poisson(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, k, n= *(lm->n), nvars= *thlength, idxj, idxk;
  double *ypred, *x= lm->x;

  //linpred= (*funargs)["linpred"];
  ypred= (*funargs)["ypred"];
  //ytlinpred= (*funargs)["ytlinpred"];

  for (j=1; j<=nvars; j++) {
    idxj= *(lm->n) * sel[j-1];

    hess[j][j]= 0;
    for (i=0; i<n; i++) { hess[j][j] += ypred[i] * x[idxj + i] * x[idxj + i]; }

    for (k=1; k<j; k++) {
      idxk= *(lm->n) * sel[k-1];

      hess[j][k]= 0;
      for (i=0; i<n; i++) { hess[j][k] += ypred[i] * x[idxj + i] * x[idxk + i]; }
      hess[k][j]= hess[j][k];
    }

  }

}



/* EVALUATE LOG-LIKELIHOOD AND DERIVATIVES FOR BETA= MLE AT THE INTERCEPT-ONLY MODEL */
//Function to evaluate the negative loglikelihood
void neglogl0_poisson(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 
  double meany= (*(lm->sumy)) / (*(lm->n) + .0);
  (*f)= - (*(lm->sumy)) * (log(meany) - 1) + *(lm->sumlogyfact);
}


//Obtain gradient and hessian wrt j
void negloglgradhess0_poisson(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  double *colsumsx= lm->colsumsx, *ytX= lm->ytX, meany= (*(lm->sumy)) / (*(lm->n) + .0);

  (*grad)= -ytX[sel[j]] + meany * colsumsx[sel[j]];
  (*hess)= meany * (lm->XtX)->at(sel[j],sel[j]);

}

//Obtain hessian
void negloglhess0_poisson(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, nvars= *thlength;
  double meany= (*(lm->sumy)) / (*(lm->n) + .0);

  for (i=1; i<=nvars; i++) {
    hess[i][i]= meany * (lm->XtX)->at(sel[i-1],sel[i-1]);
    for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= meany * (lm->XtX)->at(sel[i-1],sel[j-1]); }
  }

}


/* LOG-LIKELIHOOD GRADIENT AND HESSIANS AT th=0 (USED BY ALA) */
void neglogl00_poisson(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) { 
  (*f) = (*(lm->n) + .0) + *(lm->sumlogyfact);
}


void negloglgradhess00_poisson(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  double *colsumsx= lm->colsumsx, *ytX= lm->ytX;

  (*grad)= -ytX[sel[j]] + colsumsx[sel[j]];
  (*hess)= (lm->XtX)->at(sel[j],sel[j]);

}

void negloglhess00_poisson(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, nvars= *thlength;

  for (i=1; i<=nvars; i++) {
    hess[i][i]= (lm->XtX)->at(sel[i-1],sel[i-1]);
    for (j=1; j<i; j++) { hess[i][j]= hess[j][i]= (lm->XtX)->at(sel[i-1],sel[j-1]); }
  }

}




//AUXILIARY FUNCTIONS TO EVALUATE LOG-JOINT AND DERIVATIVES

//Log joint under pMOM - gZellner prior
void fjoint_poisson_pmomgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_poisson, &pmomgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_poisson_pmomgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_poisson, &pmomgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_poisson_pmomgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_poisson, &pmomgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_poisson_pmomgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_poisson, &pmomgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_poisson_pmomgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_poisson, &pmomgzell_hess, hess, th, sel, thlength, lm, funargs);
}


//Log joint under peMOM - gZellner prior
void fjoint_poisson_pemomgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_poisson, &pemomgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_poisson_pemomgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_poisson, &pemomgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_poisson_pemomgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_poisson, &pemomgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_poisson_pemomgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_poisson, &pemomgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_poisson_pemomgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_poisson, &pemomgzell_hess, hess, th, sel, thlength, lm, funargs);
}


//Log joint under gZellner - gZellner prior
void fjoint_poisson_gzellgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl_poisson, &gzellgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointu_poisson_gzellgzell(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint_update(&negloglupdate_poisson, &gzellgzell_log, fnew, thjnew, j, f, th, sel, thlength, lm, funargs);
}

void fjointg_poisson_gzellgzell(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_grad(&negloglgrad_poisson, &gzellgzell_grad, grad, j, th, sel, thlength, lm, funargs);
}

void fjointgh_poisson_gzellgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess_poisson, &gzellgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh_poisson_gzellgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess_poisson, &gzellgzell_hess, hess, th, sel, thlength, lm, funargs);
}

void fjoint0_poisson_gzellgzell(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  fjoint(&neglogl0_poisson, &gzellgzell_log, f, th, sel, thlength, lm, funargs);
}

void fjointgh0_poisson_gzellgzell(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_gradhess(&negloglgradhess0_poisson, &gzellgzell_gradhess, grad, hess, j, th, sel, thlength, lm, funargs);
}

void fjointh0_poisson_gzellgzell(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  fjoint_hess(&negloglhess0_poisson, &gzellgzell_hess, hess, th, sel, thlength, lm, funargs);
}


