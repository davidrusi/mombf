#include "survival.h"


//*************************************************************************************
// MARGINAL LIKELIHOOD FOR ACCELERATED FAILURE TIME MODELS
//*************************************************************************************

//Negative log-likelihood for AFT model with Normal errors
void negloglnormalAFT(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) {
  int i, nuncens, n= *((*lm).n), nvars= *thlength -1;
  double rho= th[*thlength -1], exprho= exp(rho), *ypred, *y= (*lm).y, sumres2, sumlogPhires, *res, *pnormres;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  (*f)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * rho);

  if (*thlength >1) {
    ypred= dvector(0,n);
    Aselvecx((*lm).x, th, ypred, 0, n-1, sel, &nvars); //Returns ypred= x[,sel] %*% th
    for (i=0, sumres2=0; i< nuncens; i++) { res[i]= exprho * y[i] - ypred[i]; sumres2 += res[i]*res[i]; } //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i] - ypred[i]; pnormres[i-nuncens]= pnormC(-res[i]); sumlogPhires += log(pnormres[i-nuncens]); } //Censored observations
    //for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i] - ypred[i]; sumlogPhires += log(pnormC(- res[i])); } //Censored observations
    free_dvector(ypred, 0,n);

  } else {

    for (i=0, sumres2=0; i< nuncens; i++) { res[i]= exprho * y[i]; sumres2 += res[i]*res[i]; }     //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i]; pnormres[i-nuncens]= pnormC(-res[i]); sumlogPhires += log(pnormres[i-nuncens]); } //Censored observations
    //for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i]; sumlogPhires += log(pnormC(- res[i])); } //Censored observations

  }

  (*f)= (*f) + 0.5 * sumres2 - sumlogPhires;

}

//Negative log-likelhood for AFT model with Normal errors
void negloglnormalAFTupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {

  int i, idxj, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], *y= (*lm).y, sumres2, sumlogPhires, *res, *pnormres, *x= (*lm).x, thdif, exprhodif;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  idxj= *((*lm).n) * sel[j];

  if (j < *thlength -1) { //updating a regression coefficient

    (*fnew)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * rho);
    thdif= th[j] - *thjnew;
    for (i=0, sumres2=0; i< nuncens; i++) { //Contribution from uncensored observations
      res[i] += x[i + idxj] * thdif; //Update res[i]= exprho * y[i] - ypred[i];
      sumres2 += res[i]*res[i];
    }
    for (i=nuncens, sumlogPhires=0; i< n; i++) { //Contribution from censored observations
      res[i] += x[i + idxj] * thdif;
      pnormres[i-nuncens]= pnormC(-res[i]);
      sumlogPhires += log(pnormres[i-nuncens]);
      //sumlogPhires += log(pnormC(- res[i]));
    }

  } else { //updating rho= log(residual precision)

    (*fnew)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * (*thjnew));
    exprhodif= exp(*thjnew) - exp(th[*thlength -1]);
    for (i=0, sumres2=0; i< nuncens; i++) { res[i] += y[i] * exprhodif; sumres2 += res[i]*res[i]; } //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i] += y[i] * exprhodif; pnormres[i-nuncens]= pnormC(-res[i]); sumlogPhires += log(pnormres[i-nuncens]); }  //Censored observations
    //for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i] += y[i] * exprhodif; sumlogPhires += log(pnormC(- res[i])); }  //Censored observations

  }

  (*fnew)= (*fnew) + 0.5 * sumres2 - sumlogPhires;

}

//Gradient and hessian wrt th[j] of negative log-likelihood for AFT model with Normal errors
void negloglnormalAFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, idxj, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], exprho, *y= (*lm).y, *x= (*lm).x, *res, *pnormres, ytres, *sumy2obs, sumy2D, r;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  sumy2obs= (*funargs)["sumy2obs"];
  idxj= *((*lm).n) * sel[j];
  (*grad)= (*hess)= 0;

  if (j < *thlength -1) { //updating a regression coefficient

    for (i=0; i< nuncens; i++) (*grad) -= res[i] * x[idxj +i]; //Uncensored observations
    (*hess)= ((*lm).XtXuncens)->at(sel[j],sel[j]);
    for (i=nuncens; i< n; i++) { //Censored observations
      r= dnormC(-res[i],0) / pnormres[i-nuncens]; //r= invmillsnorm(-res[i]);
      (*grad) -= r * x[idxj +i];
      (*hess) += x[i + idxj] * x[i + idxj] * r*(r-res[i]);
      //(*grad) -= invmillsnorm(-res[i]) * x[idxj +i]; (*hess) += x[i + idxj] * x[i + idxj] * infopropAFT(res[i]); //Old version, slower as it requires evaluating invmillsnorm twice
    }

  } else { //updating rho= log(residual precision)

    exprho= exp(rho); ytres= sumy2D= 0;
    for (i=0; i< nuncens; i++) ytres += res[i] * y[i]; //Uncensored observations
    for (i=nuncens; i< n; i++) {                       //Censored observations
      r= dnormC(-res[i],0) / pnormres[i-nuncens];       //r= invmillsnorm(-res[i]);
      ytres += r * y[i];
      sumy2D += y[i]*y[i] * r*(r-res[i]);
      //ytres += invmillsnorm(-res[i]) * y[i]; sumy2D += y[i]*y[i]*infopropAFT(res[i]); //Old version, slower as it requires evaluating invmillsnorm twice
    }
    (*grad)= -(*(*funargs)["nuncens"]) + exprho * ytres;
    (*hess)= exprho * ytres + exprho * exprho * (*sumy2obs + sumy2D);

  }

}


//Full hessian matrix H[1..nsel][1..nsel] of log-likelihood for AFT model with Normal errors (only upper-triangular elements are returned)

void negloglnormalAFThess(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, l, idxj, idxl, nuncens, n= *((*lm).n), nvars= *thlength -1;
  double rho= th[*thlength -1], exprho, *y= (*lm).y, *x= (*lm).x, *ytXuncens= (*lm).ytXuncens, *res, *pnormres, ytres, *sumy2obs, *D, sumD=0, sumy2D=0, xyD, r;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  sumy2obs= (*funargs)["sumy2obs"];
  D= dvector(0, n-nuncens);

  //hessian wrt log-residual precision
  exprho= exp(rho); ytres= 0;
  for (i=0; i< nuncens; i++) ytres += res[i] * y[i]; //Uncensored observations
  for (i=nuncens; i< n; i++) { //Censored observations
    r= dnormC(-res[i],0) / pnormres[i-nuncens]; //r= invmillsnorm(-res[i]);
    ytres += r * y[i];
    D[i-nuncens]= r * (r-res[i]);
    //ytres += invmillsnorm(-res[i]) * y[i]; D[i-nuncens]= infopropAFT(res[i]); //Old version, slower as it required evaluating invmillsnorm twice
    sumD += D[i-nuncens];
    sumy2D += y[i]*y[i]*D[i-nuncens];
  }
  hess[*thlength][*thlength]= exprho * ytres + exprho * exprho * (*sumy2obs + sumy2D);

  //hessian wrt regression coefficients
  for (j=0; j< nvars; j++) {
    idxj= n * sel[j];
    for (l=j; l < nvars; l++) {
      idxl= n * sel[l];
      hess[j+1][l+1]= ((*lm).XtXuncens)->at(sel[j],sel[l]);
      for (i=nuncens; i< n; i++)  hess[j+1][l+1] += x[i + idxj] * x[i + idxl] * D[i - nuncens];
    }
  }

  //hessian wrt (log-residual precision, regression coefficients)
  l= *thlength;
  for (j=0; j< l-1; j++) {
    hess[j+1][l]= -exprho * ytXuncens[sel[j]];
    for (i=nuncens, idxj= n*sel[j], xyD=0; i< n; i++) xyD += x[i + idxj] * y[i] * D[i-nuncens];
    hess[j+1][l] -= exprho * xyD;
  }

  free_dvector(D, 0, n-nuncens);
}




//Fast approximation to Negative log-likelihood for AFT model with Normal errors (uses apnorm, ainvmillsnorm in cstat.cpp)
void anegloglnormalAFT(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) {
  int i, nuncens, n= *((*lm).n), nvars= *thlength -1;
  double rho= th[*thlength -1], exprho= exp(rho), *ypred, *y= (*lm).y, sumres2, sumlogPhires, *res, *pnormres;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  (*f)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * rho);

  if (*thlength >1) {
    ypred= dvector(0,n);
    Aselvecx((*lm).x, th, ypred, 0, n-1, sel, &nvars); //Returns ypred= x[,sel] %*% th
    for (i=0, sumres2=0; i< nuncens; i++) { res[i]= exprho * y[i] - ypred[i]; sumres2 += res[i]*res[i]; } //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i] - ypred[i]; pnormres[i-nuncens]= apnorm(-res[i],false); sumlogPhires += log(pnormres[i-nuncens]); } //Censored observations
    free_dvector(ypred, 0,n);

  } else {

    for (i=0, sumres2=0; i< nuncens; i++) { res[i]= exprho * y[i]; sumres2 += res[i]*res[i]; }     //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i]= exprho * y[i]; pnormres[i-nuncens]= apnorm(-res[i],false); sumlogPhires += log(pnormres[i-nuncens]); } //Censored observations

  }

  (*f)= (*f) + 0.5 * sumres2 - sumlogPhires;

}

//Same as anegloglnormalAFT, but assumes that regression coefficients th[0,...,*thlength-2]= 0
//NOTE: error log-variance rho not assumed to be zero, th[*thlength -1] is taken
void anegloglnormalAFT0(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs) {
  int i, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], exprho= exp(rho), *y= (*lm).y, sumres2, sumlogPhires, *res, *pnormres;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  (*f)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * rho);

  //Uncensored observations
  for (i=0, sumres2=0; i< nuncens; i++) { res[i]= exprho * y[i]; sumres2 += res[i]*res[i]; }

  //Censored observations
  for (i=nuncens, sumlogPhires=0; i< n; i++) { 
     res[i]= exprho * y[i]; 
     pnormres[i-nuncens]= apnorm(-res[i],false); 
     sumlogPhires += log(pnormres[i-nuncens]); 
  }

  (*f)= (*f) + 0.5 * sumres2 - sumlogPhires;

}


//Fast approximation to Negative log-likelhood for AFT model with Normal errors (uses apnorm, ainvmillsnorm in cstat.cpp)
void anegloglnormalAFTupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {

  int i, idxj, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], *y= (*lm).y, sumres2, sumlogPhires, *res, *pnormres, *x= (*lm).x, thdif, exprhodif;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  idxj= *((*lm).n) * sel[j];

  if (j < *thlength -1) { //updating a regression coefficient

    (*fnew)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * rho);
    thdif= th[j] - *thjnew;
    for (i=0, sumres2=0; i< nuncens; i++) { //Contribution from uncensored observations
      res[i] += x[i + idxj] * thdif; //Update res[i]= exprho * y[i] - ypred[i];
      sumres2 += res[i]*res[i];
    }
    for (i=nuncens, sumlogPhires=0; i< n; i++) { //Contribution from censored observations
      res[i] += x[i + idxj] * thdif;
      pnormres[i-nuncens]= apnorm(-res[i],false);
      sumlogPhires += log(pnormres[i-nuncens]);
    }

  } else { //updating rho= log(residual precision)

    (*fnew)= 0.5 * (*(*funargs)["nuncens"]) * (LOG_M_2PI - 2.0 * (*thjnew));
    exprhodif= exp(*thjnew) - exp(th[*thlength -1]);
    for (i=0, sumres2=0; i< nuncens; i++) { res[i] += y[i] * exprhodif; sumres2 += res[i]*res[i]; } //Uncensored observations
    for (i=nuncens, sumlogPhires=0; i< n; i++) { res[i] += y[i] * exprhodif; pnormres[i-nuncens]= apnorm(-res[i],false); sumlogPhires += log(pnormres[i-nuncens]); }  //Censored observations

  }

  (*fnew)= (*fnew) + 0.5 * sumres2 - sumlogPhires;

}

//Fast approximation to Gradient and hessian wrt th[j] of negative log-likelihood for AFT model with Normal errors (uses apnorm, ainvmillsnorm from cstat.cpp)
void anegloglnormalAFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, idxj, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], exprho, *y= (*lm).y, *x= (*lm).x, *res, *pnormres, ytres, *sumy2obs, sumy2D, r;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  sumy2obs= (*funargs)["sumy2obs"];
  idxj= *((*lm).n) * sel[j];
  (*grad)= (*hess)= 0;

  if (j < *thlength -1) { //updating a regression coefficient

    for (i=0; i< nuncens; i++) (*grad) -= res[i] * x[idxj +i]; //Uncensored observations
    (*hess)= ((*lm).XtXuncens)->at(sel[j],sel[j]);
    for (i=nuncens; i< n; i++) { //Censored observations
      //Obtain r= ainvmillsnorm(-res[i]);
      if (res[i]> 1.756506) {
        r= (res[i] + 1.0/(res[i]+2.0/(res[i]+3.0/(res[i]+4.0/(res[i]+5.0/(res[i]+11.5/(res[i] + 4.890096)))))));
      } else {
        r= dnormC(-res[i],0) / pnormres[i-nuncens];
      }
      (*grad) -= r * x[idxj +i];
      (*hess) += x[i + idxj] * x[i + idxj] * r*(r-res[i]);
      //(*grad) -= invmillsnorm(-res[i]) * x[idxj +i]; (*hess) += x[i + idxj] * x[i + idxj] * infopropAFT(res[i]); //Old version, slower as it requires evaluating invmillsnorm twice
    }

  } else { //updating rho= log(residual precision)

    exprho= exp(rho); ytres= sumy2D= 0;
    for (i=0; i< nuncens; i++) ytres += res[i] * y[i]; //Uncensored observations
    for (i=nuncens; i< n; i++) {                       //Censored observations
      //Obtain r= ainvmillsnorm(-res[i]);
      if (res[i]> 1.756506) {
        r= (res[i] + 1.0/(res[i]+2.0/(res[i]+3.0/(res[i]+4.0/(res[i]+5.0/(res[i]+11.5/(res[i] + 4.890096)))))));
      } else {
        r= dnormC(-res[i],0) / pnormres[i-nuncens];
      }
      ytres += r * y[i];
      sumy2D += y[i]*y[i] * r*(r-res[i]);
      //ytres += invmillsnorm(-res[i]) * y[i]; sumy2D += y[i]*y[i]*infopropAFT(res[i]); //Old version, slower as it requires evaluating invmillsnorm twice
    }
    (*grad)= -(*(*funargs)["nuncens"]) + exprho * ytres;
    (*hess)= exprho * ytres + exprho * exprho * (*sumy2obs + sumy2D);

  }

}


void anegloglnormalAFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, idxj, nuncens, n= *((*lm).n);
  double rho= th[*thlength -1], exprho, *y= (*lm).y, *x= (*lm).x, *res, *pnormres, ytres, sumy2D, r;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  //sumy2obs= (*funargs)["sumy2obs"];
  idxj= *((*lm).n) * sel[j];
  (*grad)= 0;

  if (j < *thlength -1) { //updating a regression coefficient

    for (i=0; i< nuncens; i++) (*grad) -= res[i] * x[idxj +i]; //Uncensored observations
    for (i=nuncens; i< n; i++) { //Censored observations
      //Obtain r= ainvmillsnorm(-res[i]);
      if (res[i]> 1.756506) {
        r= (res[i] + 1.0/(res[i]+2.0/(res[i]+3.0/(res[i]+4.0/(res[i]+5.0/(res[i]+11.5/(res[i] + 4.890096)))))));
      } else {
        r= dnormC(-res[i],0) / pnormres[i-nuncens];
      }
      (*grad) -= r * x[idxj +i];
    }

  } else { //updating rho= log(residual precision)

    exprho= exp(rho); ytres= sumy2D= 0;
    for (i=0; i< nuncens; i++) ytres += res[i] * y[i]; //Uncensored observations
    for (i=nuncens; i< n; i++) {                       //Censored observations
      //Obtain r= ainvmillsnorm(-res[i]);
      if (res[i]> 1.756506) {
        r= (res[i] + 1.0/(res[i]+2.0/(res[i]+3.0/(res[i]+4.0/(res[i]+5.0/(res[i]+11.5/(res[i] + 4.890096)))))));
      } else {
        r= dnormC(-res[i],0) / pnormres[i-nuncens];
      }
      ytres += r * y[i];
      sumy2D += y[i]*y[i] * r*(r-res[i]);
    }
    (*grad)= -(*(*funargs)["nuncens"]) + exprho * ytres;

  }

}


//Fast approx to Full hessian matrix H[1..nsel][1..nsel] of log-likelihood for AFT model with Normal errors (only upper-triangular elements are returned)

void anegloglnormalAFThess(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {

  int i, j, l, idxj, idxl, nuncens, n= *((*lm).n), nvars= *thlength -1;
  double rho= th[*thlength -1], exprho, *y= (*lm).y, *x= (*lm).x, *ytXuncens= (*lm).ytXuncens, *res, *pnormres, ytres, *sumy2obs, *D, sumD=0, sumy2D=0, xyD, r;

  nuncens= (int) (*(*funargs)["nuncens"] +.1);
  res= (*funargs)["residuals"];
  pnormres= (*funargs)["pnormres"];
  sumy2obs= (*funargs)["sumy2obs"];
  D= dvector(0, n-nuncens);

  //hessian wrt log-residual precision
  exprho= exp(rho); ytres= 0;
  for (i=0; i< nuncens; i++) ytres += res[i] * y[i]; //Uncensored observations
  for (i=nuncens; i< n; i++) { //Censored observations
    //Obtain r= ainvmillsnorm(-res[i]);
    if (res[i]> 1.756506) {
      r= (res[i] + 1.0/(res[i]+2.0/(res[i]+3.0/(res[i]+4.0/(res[i]+5.0/(res[i]+11.5/(res[i] + 4.890096)))))));
    } else {
      r= dnormC(-res[i],0) / pnormres[i-nuncens];
    }
    ytres += r * y[i];
    D[i-nuncens]= r * (r-res[i]);
    //ytres += invmillsnorm(-res[i]) * y[i]; D[i-nuncens]= infopropAFT(res[i]); //Old version, slower as it required evaluating invmillsnorm twice
    sumD += D[i-nuncens];
    sumy2D += y[i]*y[i]*D[i-nuncens];
  }
  hess[*thlength][*thlength]= exprho * ytres + exprho * exprho * (*sumy2obs + sumy2D);

  //hessian wrt regression coefficients
  for (j=0; j< nvars; j++) {
    idxj= n * sel[j];
    for (l=j; l < nvars; l++) {
      idxl= n * sel[l];
      hess[j+1][l+1]= ((*lm).XtXuncens)->at(sel[j],sel[l]);
      for (i=nuncens; i< n; i++)  hess[j+1][l+1] += x[i + idxj] * x[i + idxl] * D[i - nuncens];
    }
  }

  //hessian wrt (log-residual precision, regression coefficients)
  l= *thlength;
  for (j=0; j< l-1; j++) {
    hess[j+1][l]= -exprho * ytXuncens[sel[j]];
    for (i=nuncens, idxj= n*sel[j], xyD=0; i< n; i++) xyD += x[i + idxj] * y[i] * D[i-nuncens];
    hess[j+1][l] -= exprho * xyD;
  }

  free_dvector(D, 0, n-nuncens);
}


//Proportion of information in the AFT model contained in an observation censored z standard deviations after the mean
/*
double infopropAFT(double z) {
  double ans;
  ans= dnormC(z,0) / (1.0 - pnormC(z));
  ans= ans * (ans - z);
  return ans;
}
*/

double pmomgmomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= pmomgmomSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double pmomgmomSurvMarg(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood under pMOM(tau) + group MOM(taugroup) prior

    prod_j pMOM(beta_j; tau)  prod_j N(delta_j; 0, (taugroup/[ncol(S_j)+2]) n (S_j)^{-1})

    where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
    corresponding to coefficients for individual variables and grouped variables (respectively)

    S_j is the submatrix of S indicated by sel[0], sel[1] etc.
  */

  return SurvMargALA(sel, nsel, lm, 10); //priorcode=10 is pMOM + group pMOM

}

double gmomgmomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= gmomgmomSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double gmomgmomSurvMarg(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood under pMOM(tau) + group MOM(taugroup) prior

    prod_j pMOM(beta_j; (tau/3) n (S_j)^{-1})  prod_j gMOM(delta_j; 0, (taugroup/[ncol(S_j)+2]) n (S_j)^{-1})

    S_j is the submatrix of S indicated by sel[0], sel[1] etc.
  */

  return SurvMargALA(sel, nsel, lm, 50); //priorcode=50 is group pMOM + group pMOM

}

double gmomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= gmomgzellSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double gmomgzellSurvMarg(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood under group pMOM(tau) + group Zellner(taugroup) prior

    prod_j pMOM(beta_j; (tau/3) n (S_j)^{-1})  prod_j N(delta_j; 0, (taugroup/ncol(S_j)) n (S_j)^{-1})

    S_j is the submatrix of S indicated by sel[0], sel[1] etc.
  */

    return SurvMargALA(sel, nsel, lm, 53); //priorcode=53 is group pMOM + group Zellner

}


double pmomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= pmomgzellSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double pmomgzellSurvMarg(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood under pMOM(tau) + block Zellner(taugroup) prior

    prod_j pMOM(beta_j; tau)  prod_j N(delta_j; 0, (taugroup/ncol(S_j)) (S_j)^{-1})

    where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
    corresponding to coefficients for individual variables and grouped variables (respectively)

    S_j is the submatrix of S indicated by sel[0], sel[1] etc.
  */

  if (*(*lm).method ==2) {
    return SurvMargALA(sel, nsel, lm, 13); //priorcode=13 is pMOM + block Zellner
  } else {
    return SurvMarg(sel, nsel, lm, 13); //priorcode=13 is pMOM + block Zellner
  }

}


// peMOM on individual coef, group eMOM on groups
double pemomgemomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  Rf_error("peMOM + group eMOM not currently implemented for the AFT Normal model");
}


// peMOM on individual coef, block Zellner on groups
double pemomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= pemomgzellSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double pemomgzellSurvMarg(int *sel, int *nsel, lmObject *lm) {
  /*Marginal likelihood under peMOM(tau) + block Zellner(taugroup) prior

    prod_j pMOM(beta_j; tau)  prod_j N(delta_j; 0, (taugroup/ncol(S_j)) (S_j)^{-1})

    where beta=(beta_1,...,beta_p) and delta=(delta_1,...,delta_q) are subsets of th
    corresponding to coefficients for individual variables and grouped variables (respectively)

    S_j is the submatrix of S indicated by sel[0], sel[1] etc.
  */

  return SurvMarg(sel, nsel, lm, 33); //priorcode=33 is pMOM + block Zellner
}

// Zellner on individual coef, block Zellner on groups
double gzellgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold=nullptr, arma::mat *m=nullptr, arma::mat *cholVinv=nullptr) {
  int *selvec, nsel= sel->n_nonzero;
  double ans;
  selvec= ivector(0, nsel);
  spmat_to_ivector(selvec, nullptr, sel);
  ans= gzellgzellSurvMarg(selvec, &nsel, lm);
  free_ivector(selvec, 0, nsel);
  return ans;
}

double gzellgzellSurvMarg (int *sel, int *nsel, lmObject *lm) {

  if (*(*lm).method ==2) {
    return SurvMargALA(sel, nsel, lm, 43); //priorcode=43 is block Zellner + block Zellner
  } else {
    return SurvMarg(sel, nsel, lm, 43);
  }

}



double SurvMargALA(int *sel, int *nsel, lmObject *lm, int priorcode) {
  /*Marginal likelihood for AFT survival model under pMOM/peMOM + block Zellner(taugroup) prior

    priorcode indicates what prior is used. Currently implemented options

    10: pMOM + group MOM
    13: pMOM + group Zellner
    43: group Zellner + group Zellner
    50: group pMOM + group pMOM
    53: group pMOM + group Zellner
   */

  std::map<string, double *> funargs;
  bool momsingle, momgroup, converged;
  int i, nselgroupsint, cholSsize, *uncens, thlength= *nsel +1;
  double ans, nuncens, sumy2obs=0, *residuals, nselgroups, *nvarinselgroups, *firstingroup, *selgroups, *ldetSinv, *cholSini, *cholSinv, *Sinv, *thini, *thopt, fini, *y, *pnormres, *g, **H, **Hinv, **cholH, logdispersion, *delta;
  modselFunction *msfun;

  g= dvector(1,thlength); H= dmatrix(1,thlength,1,thlength); Hinv= dmatrix(1,thlength,1,thlength); cholH= dmatrix(1,thlength,1,thlength);
  thopt= dvector(0, *nsel); thini= dvector(0, *nsel);     delta= dvector(1, thlength);
  y= (*lm).y;

  //Initialize static elements in funargs (not changed by msfun)
  uncens= (*lm).uncens;
  for (i=0; (i< *((*lm).n) && (uncens[i]==1)); i++) { sumy2obs+= y[i] * y[i]; }
  nuncens= (double) i;
  funargs["nuncens"]= &nuncens; //number of uncensored observations
  funargs["sumy2obs"]= &sumy2obs; //sum of squares for uncensored observations, i.e. sum_{i: uncens[i]==1} y[i]^2

  nvarinselgroups= dvector(0, min_xy(*nsel, *((*lm).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*lm).ngroups))); selgroups= dvector(0, *nsel -1);
  findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*lm).nvaringroup, (*lm).ngroups); //copy subset of nvaringroup into nvarinselgroups
  funargs["nvarinselgroups"]= nvarinselgroups;
  funargs["firstingroup"]= firstingroup;
  funargs["nselgroups"]= &nselgroups;
  funargs["selgroups"]= selgroups;
  nselgroupsint= (int) (nselgroups +.1);

  //Obtain Cholesky decomp and determinant of prior scale covariances for each group
  ldetSinv= dvector(0, nselgroupsint); cholSini= dvector(0, nselgroupsint);
  cholSini_indexes(cholSini, &cholSsize, nselgroupsint, nvarinselgroups);
  cholSinv= dvector(0, cholSsize); Sinv= dvector(0, cholSsize);

  funargs["cholSini"]= cholSini; //cholSini[j] is the index in cholSinv at which Sinv_j starts
  gzell_Sinv_byprior(Sinv, cholSinv, ldetSinv, &nselgroupsint, nvarinselgroups, sel, cholSini, (*lm).XtX, (*lm).n, (*lm).tau, (*lm).taugroup, priorcode);
  funargs["ldetSinv"]= ldetSinv; funargs["cholSinv"]= cholSinv; funargs["Sinv"]= Sinv;

  //Initialize dynamic elements in funargs (changed by msfun)
  residuals= dvector(0, *((*lm).n));
  pnormres= dvector(0, *((*lm).n) - nuncens);
  funargs["residuals"]= residuals;
  funargs["pnormres"]= pnormres;

  //Assign functions to evaluate log-posterior, update log-posterior, gradient and hessians
  msfun= new modselFunction(sel, thlength, lm, NULL);

  //ALA to integrated likelihood under the base Normal prior
  msfun->funupdate= &fgzellgzellSurvupdate;  //objective function
  msfun->gradhessUniv= &fgzellgzell_AFTgradhess; msfun->hess= &fgzellgzellhess_AFT; msfun->gradUniv= &fgzellgzell_AFTgrad; //derivatives
  msfun->ftol= 0.001; msfun->thtol= 0.001;

  //Initialize. If not previously computed, find the optimal value of log-error dispersion, given theta=0
  for (i=0; i< thlength; i++) thini[i]= 0;

  if (*((*lm).usethinit) == 2) {
    //take previously-computed error log-variance parameter
    thini[*nsel]= ((*lm).thinit)[*((*lm).p)];
    //evaluate fini at regression coef th=0 and initialize funargs
    msfun->fun= &fgzellgzellSurv0;
    msfun->evalfun(&fini, thini, &funargs);
    msfun->fun= &fgzellgzellSurv;
  } else { 
    logdispersion= 0;
    //evaluate fini at regression coef th=0 and initialize funargs
    msfun->fun= &fgzellgzellSurv0; 
    msfun->evalfun(&fini, thini, &funargs); 
    //optimize error log-variance parameter
    msfun->fun= &fgzellgzellSurv;
    msfun->Newtonuniv(&logdispersion, *nsel, &fini, &converged, thini, &funargs, 5); //fini returns f at optimal log dispersion
    thini[*nsel]= ((*lm).thinit)[*((*lm).p)]= logdispersion;
    (*((*lm).usethinit))= 2;
  } 

  ans= msfun->ALA(thini, &fini, g, H, cholH, Hinv, true, true, 1.0, &funargs); //aprox marginal likelihood and return g, H, cholH and Hinv

  //If needed, add term corresponding to the non-local prior penalty
  momsingle= ((priorcode==10) || (priorcode==13) || (priorcode==50) || (priorcode==53)); //pMOM or groupMOM on single coef was set
  momgroup= ((priorcode==10) || (priorcode)==50); //groupMOM on groups of coef was set

 
  if (momsingle || momgroup) {

    //Compute thopt= thini - Hinv g
    Ax(Hinv,g,delta,1,thlength,1,thlength); 
    for (i=0; i<= *nsel; i++) thopt[i]= thini[i] - delta[i+1];

    gmompenalty_approx(momsingle, momgroup, thopt, Hinv, Sinv, exp(thopt[*sel]), thlength, *nsel, nselgroupsint, nvarinselgroups, firstingroup, cholSini);

  }

  //Free memory
  free_dvector(g,1,thlength); free_dmatrix(H,1,thlength,1,thlength); free_dmatrix(Hinv,1,thlength,1,thlength); free_dmatrix(cholH, 1,thlength,1,thlength);
  free_dvector(thopt, 0, *nsel); free_dvector(thini, 0, *nsel); free_dvector(delta, 1, thlength);
  free_dvector(residuals, 0, *((*lm).n));
  free_dvector(pnormres, 0, *((*lm).n) - nuncens);
  free_dvector(nvarinselgroups, 0, min_xy(*nsel, *((*lm).ngroups)));
  free_dvector(firstingroup, 0, min_xy(*nsel, *((*lm).ngroups)));
  free_dvector(selgroups, 0, *nsel -1);
  free_dvector(ldetSinv, 0, nselgroupsint); free_dvector(cholSini, 0, nselgroupsint); free_dvector(cholSinv, 0, cholSsize); free_dvector(Sinv, 0, cholSsize);
  delete msfun;

  return ans;
}




double SurvMarg(int *sel, int *nsel, lmObject *lm, int priorcode) {
  /*Marginal likelihood for AFT survival model under pMOM/peMOM + block Zellner(taugroup) prior

    priorcode indicates what prior is used. Currently implemented options

    13: pMOM + group Zellner
    32: peMOM + group eMOM
    33: peMOM + group Zellner
    43: group Zellner + group Zellner
   */

  std::map<string, double *> funargs;
  bool posdef, orthoapprox=false, converged;
  int i, nselgroupsint, cholSsize, *uncens, thlength= *nsel +1;
  double ans, nuncens, sumy2obs=0, *residuals, nselgroups, *nvarinselgroups, *firstingroup, *selgroups, *ldetSinv, *cholSini, *cholSinv, *Sinv, *thini, *thopt, fini, fopt, *y, *pnormres, *g, **H, **Hinv, **cholH;
  modselFunction *msfun;

  g= dvector(1,thlength); H= dmatrix(1,thlength,1,thlength); Hinv= dmatrix(1,thlength,1,thlength); cholH= dmatrix(1,thlength,1,thlength);
  thopt= dvector(0, *nsel); thini= dvector(0, *nsel);
  y= (*lm).y;

  //For MOM priors, if ALA is specified then approximate the mean of products via the product of means
  if ((priorcode==10) | (priorcode == 13)) {
      if ((*(*lm).method ==2) | ((*(*lm).method == -1) & ((*nsel)>0)))  { orthoapprox= true; }
  }

  //Initialize static elements in funargs (not changed by msfun)
  uncens= (*lm).uncens;
  for (i=0; (i< *((*lm).n) && (uncens[i]==1)); i++) { sumy2obs+= y[i] * y[i]; }
  nuncens= (double) i;
  funargs["nuncens"]= &nuncens; //number of uncensored observations
  funargs["sumy2obs"]= &sumy2obs; //sum of squares for uncensored observations, i.e. sum_{i: uncens[i]==1} y[i]^2

  nvarinselgroups= dvector(0, min_xy(*nsel, *((*lm).ngroups))); firstingroup= dvector(0, min_xy(*nsel, *((*lm).ngroups))); selgroups= dvector(0, *nsel -1);
  findselgroups(nvarinselgroups, firstingroup, &nselgroups, selgroups, sel, nsel, (*lm).nvaringroup, (*lm).ngroups); //copy subset of nvaringroup into nvarinselgroups
  funargs["nvarinselgroups"]= nvarinselgroups;
  funargs["firstingroup"]= firstingroup;
  funargs["nselgroups"]= &nselgroups;
  funargs["selgroups"]= selgroups;
  nselgroupsint= (int) (nselgroups +.1);

  //Obtain Cholesky decomp and determinant of prior scale covariances for each group
  ldetSinv= dvector(0, nselgroupsint); cholSini= dvector(0, nselgroupsint);
  cholSini_indexes(cholSini, &cholSsize, nselgroupsint, nvarinselgroups);
  funargs["cholSini"]= cholSini; //cholSini[j] is the index in cholSinv at which Sinv_j starts

  cholSinv= dvector(0, cholSsize); Sinv= dvector(0, cholSsize);
  gzell_Sinv(Sinv, cholSinv, ldetSinv, &nselgroupsint, nvarinselgroups, sel, cholSini, (*lm).XtX, (*lm).tau, (*lm).taugroup, orthoapprox);
  funargs["ldetSinv"]= ldetSinv; funargs["cholSinv"]= cholSinv; funargs["Sinv"]= Sinv;

  //Initialize dynamic elements in funargs (changed by msfun)
  residuals= dvector(0, *((*lm).n));
  pnormres= dvector(0, *((*lm).n) - nuncens);
  funargs["residuals"]= residuals;
  funargs["pnormres"]= pnormres;

  //Assign functions to evaluate log-posterior, update log-posterior, gradient and hessians
  msfun= new modselFunction(sel, thlength, lm, NULL);

  //Initialize posterior mode
  msfun->fun= &fgzellgzellSurv; msfun->funupdate= &fgzellgzellSurvupdate; msfun->gradhessUniv= &fgzellgzell_AFTgradhess; msfun->hess= &fgzellgzellhess_AFT; //Zell
  msfun->gradUniv= &fgzellgzell_AFTgrad;
  msfun->ftol= 0.001; msfun->thtol= 0.001;
  if (*((*lm).optim_maxit) >= 0) msfun->maxiter= *((*lm).optim_maxit);
   
  for (i=0; i< thlength; i++) thini[i]= 0;
  msfun->evalfun(&fini, thini, &funargs); //call evalfun for its side effect of initializing funargs
  msfun->hess(H, thini, sel, &thlength, lm, &funargs);
  inv_posdef(H, thlength, Hinv, &posdef);
  for (i=0; i< thlength; i++) { msfun->gradUniv(g+1+i, i, thini, sel, &thlength, lm, &funargs); g[i+1]= -g[i+1]; }
  Ax(Hinv,g,thini-1,1,thlength,1,thlength);
   
  //Stored posterior mode under previously visited model
  if (*((*lm).usethinit) == 2) {
    for (i=0; i< *nsel; i++) { thopt[i]= ((*lm).thinit)[sel[i]]; }
    thopt[*nsel]= ((*lm).thinit)[*((*lm).p)];
    msfun->evalfun(&fini, thini, &funargs);
    msfun->evalfun(&fopt, thopt, &funargs);
    if (fopt < fini) { for (i=0; i< *nsel; i++) thini[i]= thopt[i]; }
  }
   
  //Optimize and approximate the integrated likelihood
  if (priorcode != 43) {
    if (priorcode == 13) {
      if (!orthoapprox) {
        msfun->fun= &fpmomgzellSurv;
        msfun->funupdate= &fpmomgzellSurvupdate;
        msfun->gradUniv= &fpmomgzell_AFTgrad;
        msfun->gradhessUniv= &fpmomgzell_AFTgradhess;
        msfun->hess= &fpmomgzellhess_AFT;
      }
    } else if (priorcode==33) {
      msfun->fun= &fpemomgzellSurv;
      msfun->funupdate= &fpemomgzellSurvupdate;
      msfun->gradUniv= &fpemomgzell_AFTgrad;
      msfun->gradhessUniv= &fpemomgzell_AFTgradhess;
      msfun->hess= &fpemomgzellhess_AFT;
    } else {
      Rf_error("priorcode in SurvMarg not recognized\n");
    }
  }
   
  if ((priorcode != 43) && (!((priorcode == 13) && orthoapprox))) {   //Avoid exact zeroes (0 prior density under non-local priors)
    for (i=0; i< *nsel; i++) {
      if (fabs(thini[i]) < 1.0e-5) {
        double fminus, fplus;
        thini[i]= -1.0e-5; msfun->evalfun(&fminus, thini, &funargs);
        thini[i]=  1.0e-5; msfun->evalfun(&fplus, thini, &funargs);;
        if (fminus<=fplus) { thini[i]= -1.0e-5; } else { thini[i]= 1.0e-5; }
      }
    }
  }
   
  if (*nsel >=15) {
    msfun->cdaNewton(thopt, &fopt, &converged, thini, &funargs, 5);
  } else {
    msfun->Newton(thopt, &fopt, &converged, thini, &funargs, 5);
    if (!converged) msfun->cdaNewton(thopt, &fopt, &converged, thini, &funargs, 5);
  }
   
  ans= msfun->laplaceapprox(thopt, &fopt, H, cholH, true, &funargs); //Laplace approx (also returns H and cholH)
  //ans= msfun->laplaceapprox(thopt, &fopt, &funargs); //Laplace approx
   
   
  if ((priorcode == 13) && orthoapprox) { //orthogonal approx to posterior expectation of MOM penalty
    double pen;
   
    inv_posdef(H, thlength, Hinv, &posdef, cholH); //compute Hinv

    pen= pmompenalty_approx(thopt, Hinv, (*lm).tau, nselgroupsint, nvarinselgroups, firstingroup);
    ans += pen;
  }
     
  //Store optimal value for use in subsequent calls
  if (*((*lm).usethinit) > 0) {
    int iall;
    for (iall=0; iall< sel[0]; iall++) ((*lm).thinit)[iall]= 0;
    for (i=0; i< *nsel; i++) {
      ((*lm).thinit)[sel[i]]= thopt[i];
      if (i< *nsel -1) { for (iall=sel[i]+1; iall< sel[i+1]; iall++) ((*lm).thinit)[iall]= 0; }
    }
    ((*lm).thinit)[*((*lm).p)]= thopt[*nsel];
    (*((*lm).usethinit))= 2; //next time SurvMarg is called it will initialize at (*lm).thinit
  }


  //Free memory
  free_dvector(g,1,thlength); free_dmatrix(H,1,thlength,1,thlength); free_dmatrix(Hinv,1,thlength,1,thlength); free_dmatrix(cholH, 1,thlength,1,thlength);
  free_dvector(thopt, 0, *nsel);
  free_dvector(thini, 0, *nsel);
  free_dvector(residuals, 0, *((*lm).n));
  free_dvector(pnormres, 0, *((*lm).n) - nuncens);
  free_dvector(nvarinselgroups, 0, min_xy(*nsel, *((*lm).ngroups)));
  free_dvector(firstingroup, 0, min_xy(*nsel, *((*lm).ngroups)));
  free_dvector(selgroups, 0, *nsel -1);
  free_dvector(ldetSinv, 0, nselgroupsint); free_dvector(cholSini, 0, nselgroupsint); free_dvector(cholSinv, 0, cholSsize); free_dvector(Sinv, 0, cholSsize);
  delete msfun;

  return ans;
}



//Evaluate negative log-likelihood + log-prior (pMOM + group MOM) and initialize funargs
void fpmomgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens=0;

  anegloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  //negloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  dmomgzell(&priordens, th, (*lm).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  (*f) -= priordens;
}


//Evaluate negative log-likelihood + log-prior (peMOM + group MOM) and initialize funargs
void fpemomgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens=0;

  anegloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  //negloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  demomgzell(&priordens, th, (*lm).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  (*f) -= priordens;
}

//Evaluate negative log-likelihood + log-prior (group Zellner + group Zellner) and initialize funargs
void fgzellgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens=0;

  anegloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  //negloglnormalAFT(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  dgzellgzell(&priordens, th, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  (*f) -= priordens;
}


//Same as fgzllgzellSurv assuming that regression coef th=0 (error log-dispersion not assumed to be 0), and initialize funargs
void fgzellgzellSurv0(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double priordens=0;

  anegloglnormalAFT0(f, th, sel, thlength, lm, funargs); //evaluate -log(likelihood), initialize funargs
  dgzellgzell(&priordens, th, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  (*f) -= priordens;
}


//Update log-likelihood and funargs due to changing th[j] into thjnew
void fpmomgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double thtmp, priordens=0;

  anegloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  //negloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  thtmp= th[j]; th[j]= *thjnew;
  dmomgzell(&priordens, th, (*lm).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  th[j]= thtmp;
  (*fnew) -= priordens;
}

void fpemomgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double thtmp, priordens=0;

  anegloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  //negloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  thtmp= th[j]; th[j]= *thjnew;
  demomgzell(&priordens, th, (*lm).tau, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  th[j]= thtmp;
  (*fnew) -= priordens;
}

void fgzellgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs) {
  double thtmp, priordens=0;

  anegloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  //negloglnormalAFTupdate(fnew,thjnew,j,f,th,sel,thlength,lm,funargs); //update -log(likelihood) and funargs["residuals"]
  thtmp= th[j]; th[j]= *thjnew;
  dgzellgzell(&priordens, th, (*funargs)["nvarinselgroups"], (*funargs)["nselgroups"], (*funargs)["ldetSinv"], (*funargs)["cholSinv"], (*funargs)["cholSini"], true);
  priordens += dinvgammaC(exp(-2.0*th[*thlength -1]), *((*lm).alpha)/2.0, *((*lm).lambda)/2.0, 1) + log(2.0) - 2.0*th[*thlength -1];
  th[j]= thtmp;
  (*fnew) -= priordens;
}


//Gradient and hessian
void fpmomgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood
  //negloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pmomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad; (*hess) -= priorhess;
}

void fpmomgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgrad(grad, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pmomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad;
}


void fpemomgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood
  //negloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pemomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad; (*hess) -= priorhess;
}

void fpemomgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgrad(grad, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pemomgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad;
}

void fgzellgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood
  //negloglnormalAFTgradhess(grad, hess, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  gzellgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad; (*hess) -= priorhess;
}

void fgzellgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  double priorgrad, priorhess;

  anegloglnormalAFTgrad(grad, j, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  gzellgzellig_gradhess(&priorgrad, &priorhess, j, th, sel, thlength, lm, funargs); //contribution from the log-prior

  (*grad) -= priorgrad;
}


void fpmomgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  //int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  //double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  anegloglnormalAFThess(hess, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pmomgzellig_hess(hess, th, sel, thlength, lm, funargs);

}

void fpemomgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  //int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  //double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  anegloglnormalAFThess(hess, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  pemomgzellig_hess(hess, th, sel, thlength, lm, funargs);

}


void fgzellgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs) {
  /* HESSIAN FOR AFT LOG-LIKELIHOOD Xobs + LOG-LIKELIHOOD Xcens evaluated at any th   */
  //int j, k, kk, l, idxini, ngroups, ningroup, firstingroup;
  //double priorgrad, priorhess, *Sinv= (*funargs)["Sinv"], *nvaringroup= (*funargs)["nvarinselgroups"], *cholSini= (*funargs)["cholSini"];

  anegloglnormalAFThess(hess, th, sel, thlength, lm, funargs); //contribution from the log-likelihood

  gzellgzellig_hess(hess, th, sel, thlength, lm, funargs);

}



