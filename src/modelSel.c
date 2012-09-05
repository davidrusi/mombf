#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "cstat.h"
#include "modelSel.h"
#include "do_mombf.h"

//Global variables defined for minimization/integration routines
struct marginalPars f2opt_pars, f2int_pars;


//*************************************************************************************
//SETTING PRIOR & MARGINALS
//*************************************************************************************

pt2margFun set_marginalFunction(int *prCoef, int *knownphi) {
  //Returns pointer to function to compute the marginal density of the data for a given model indicator
  // - prCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
  // - knownphi: 1 if residual variance phi is know, 0 otherwise. 
  // Note: if phi known, when actually calling the returned pt2margFun, phi must be set in the parameter of type struct marginalPars *
  pt2margFun ans=NULL;
  if (*prCoef==0) {
    if (*knownphi==1) { ans= pmomMarginalKC; } else { ans= pmomMarginalUC; }
  } else if (*prCoef==1) {
    if (*knownphi==1) { ans= pimomMarginalKC; } else { ans= pimomMarginalUC; }
  } else if (*prCoef==2) {
    if (*knownphi==1) { ans= pemomMarginalKC; } else { ans= pemomMarginalUC; }
  }
  return ans;
}

pt2margFun set_priorFunction(int *prDelta) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - prDelta: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2margFun ans=NULL;
  if (*prDelta==0) { ans= unifPrior; } else if ((*prDelta==1) || (*prDelta==2)) { ans= binomPrior; }
  return ans;
}

pt2modavgPrior set_priorFunction_modavg(int *priorModel) {
  //Returns pointer to function to compute the prior probability of a model indicator
  // - priorModel: 0 for uniform, 1 for binomial, 2 for beta-binomial
  pt2modavgPrior ans=NULL;
  if (*priorModel==0) { ans= unifPrior_modavg; } else if (*priorModel==1) { ans= binomPrior_modavg; } else if (*priorModel==2) { ans= betabinPrior_modavg; }
  return ans;
}


//********************************************************************************************
// GENERAL ALGEBRA
//********************************************************************************************

//multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]
//Note: A is indexed at 1. x and sel are indexed at 0. ans is indexed at 1.
void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<=fi;_i++) {
    for (_j=_i, ans[_i]=0; _j<=fi; _j++) { ans[_i]+= A[_i][_j] * x[sel[_j-1]]; }
    for (_j= 1; _j<_i; _j++) { ans[_i]+= A[_j][_i] * x[sel[_j-1]]; }
  } 
}

//multiply symmetric A[1..ncolA][1..ncolA] (formatted as vector) with x[1..nsel]. 
//Use only selected elems in A and all elems in x[1..nsel]
void Asel_x(double *A, int ncolA, double *x, int nsel, int *sel, double *ans) {
  int _i, _j;
  for (_i=1;_i<= nsel;_i++) {
    for (_j=1, ans[_i]=0; _j<= nsel; _j++) { ans[_i]+= A[sel[_j]*ncolA+sel[_i]] * x[_j]; }
  } 
}


//Add constant ct to diagonal elements in XtX[sel,sel]. XtX[0..p-1][0..p-1] is formatted as a vector indexed at 0, V[sel[0]..sel[nsel-1]][sel[0]..sel[nsel-1]] as a matrix indexed at 1, sel is indexed at 0
//Note: Only diagonal & upper-diagonal elements in V are set.
void addct2XtX(double *ct, double *XtX, int *sel, int *nsel, int *p, double **V) {
  int i,j;
  for (i=1;i<=(*nsel);i++) { V[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]] + (*ct); }
  for (i=1;i<=(*nsel);i++) {
    for (j=i+1;j<=(*nsel);j++) {
      V[i][j]= XtX[sel[j-1]*(*p) + sel[i-1]];
    }
  }
}


//*************************************************************************************
// MODEL AVERAGING ROUTINES
//*************************************************************************************

void set_modavgPars(struct modavgPars *pars, int *n, int *p1, int *p2, int *isbinary, int *ybinary, double *y, double *sumy2, double *x1, double *x2, double *XtX, double *ytX, double *cholS2, double *S2inv, double *cholS2inv, double *colsumx1sq, double *alpha, double *lambda, int *priorCoef, int *r, double *tau1, double *tau2, int *priorTau1, double *atau1, double *btau1, int *priorModel, double *prModelpar) {
  (*pars).n= n;
  (*pars).p1= p1;
  (*pars).p2= p2;
  (*pars).isbinary= isbinary;
  (*pars).ybinary= ybinary;
  (*pars).y= y;
  (*pars).sumy2= sumy2;
  (*pars).x1= x1;
  (*pars).x2= x2;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).cholS2= cholS2;
  (*pars).S2inv= S2inv;
  (*pars).cholS2inv= cholS2inv;
  (*pars).colsumx1sq= colsumx1sq;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).priorCoef= priorCoef;
  (*pars).r= r;
  (*pars).tau1= tau1;
  (*pars).tau2= tau2;
  (*pars).priorTau1= priorTau1;
  (*pars).atau1= atau1;
  (*pars).btau1= btau1;
  (*pars).priorModel= priorModel;
  (*pars).prModelpar= prModelpar;
}


//MH within Gibbs scheme to sample from the joint posterior of (theta,delta) under a product MOM prior and a linear regression model
//Input
// - pars: data, pre-computed quantities and prior parameters. See definition of struct modavgPars for details
// - niter: number of MCMC iterations
// - thinning: only 1 out of each thinning iterations is kept
// - burnin: burning
// - niniModel: number of variables initially in the model
// - iniModel: initial model (vector of length p1)
// - iniCoef1: initial coefficient values for variables being selected (length p1)
// - iniCoef2: initial coefficient values for variables always in the model (length p2)
// - iniPhi: initial residual variance value
// - iniOthers: initial values for other variables. Currently this indicates initial tau value, and is only use if (*pars).priorTau1 !=0.
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
//Output
// - postModel: MCMC saves for variable inclusion indicators
// - margpp: marginal posterior prob for inclusion of each covariate (uses MH acceptance prob, which is more precise than simply averaging inclusion indicators)
// - postCoef1: MCMC saves for regression coefficients of variables being selected
// - postCoef2: MCMC saves for regression coefficients of variables which are always in the model
// - postPhi: MCMC saves for residual variance
// - postOther: MCMC saves for other parameters. Currently saves tau values (pMOM prior precision) when (*pars).priorTau1 != 0.
SEXP pmomLM_I(SEXP postModel, SEXP margpp, SEXP postCoef1, SEXP postCoef2, SEXP postPhi, SEXP postOther, SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP XtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar) {
  struct modavgPars pars;
  SEXP ans;

  set_modavgPars(&pars,INTEGER(n),INTEGER(p1),INTEGER(p2),INTEGER(isbinary),INTEGER(ybinary),REAL(y),REAL(sumy2),REAL(x1),REAL(x2),REAL(XtX),REAL(ytX),REAL(cholS2),REAL(S2inv),REAL(cholS2inv),REAL(colsumx1sq),REAL(alpha),REAL(lambda),INTEGER(priorCoef),INTEGER(r),REAL(tau1),REAL(tau2),INTEGER(priorTau1),REAL(atau1),REAL(btau1),INTEGER(priorModel),REAL(prModelpar));
  pmomLM(INTEGER(postModel), REAL(margpp), REAL(postCoef1), REAL(postCoef2), REAL(postPhi), REAL(postOther), &pars, INTEGER(niter), INTEGER(thinning), INTEGER(burnin), INTEGER(niniModel), INTEGER(iniModel), REAL(iniCoef1), REAL(iniCoef2), REAL(iniPhi), REAL(iniOthers), INTEGER(verbose));

  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;
}


void pmomLM(int *postModel, double *margpp, double *postCoef1, double *postCoef2, double *postPhi, double *postOther, struct modavgPars *pars, int *niter, int *thinning, int *burnin, int *niniModel, int *iniModel, double *iniCoef1, double *iniCoef2, double *iniPhi, double *iniOthers, int *verbose) {
  int i, j, k, ilow, iupper, savecnt, niterthin, niter10, nsel= *niniModel, *curModel, newdelta, n=*(*pars).n, p1=*(*pars).p1, p2=*(*pars).p2, psn, resupdate, isbinary=*(*pars).isbinary;
  double *res, *partialres, sumres2, sumpartialres2, newcoef, *curCoef1, *curCoef2, curPhi, *linpred1, *linpred2, pinclude, *temp;
  if (*verbose) Rprintf("Running MCMC");
  niterthin= floor((*niter - *burnin +.0)/(*thinning +.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  if (*burnin >0) { ilow= - *burnin; savecnt=0; iupper= *niter - *burnin +1; } else { ilow=0; savecnt=1; iupper= *niter; }
  //Initialize
  curModel= ivector(0,p1); curCoef1= dvector(0,p1); curCoef2= dvector(0,p2); linpred1= dvector(0,n); linpred2= dvector(0,n);
  res= dvector(0, n); partialres= dvector(0,n);
  for (i=0; i<p1; i++) { margpp[i]= 0; curModel[i]= postModel[i*niterthin]= iniModel[i]; curCoef1[i]= postCoef1[i*niterthin]= iniCoef1[i]; }
  for (i=0; i<p2; i++) { curCoef2[i]= postCoef2[i*niterthin]= iniCoef2[i]; }
  if (isbinary) { curPhi= postPhi[0]= 1.0; } else { curPhi= postPhi[0]= *iniPhi; }
  postOther[0]= iniOthers[0];
  Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);
  Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
  if (isbinary) sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
  for (i=0, sumres2=0; i<n; i++) { res[i]= partialres[i]= (*pars).y[i] - linpred1[i] - linpred2[i]; sumres2+= res[i]*res[i]; }
  sumpartialres2= sumres2;
  //MCMC iterations
  for (i=ilow; i< iupper; i++) {
    //Sample (curCoef1,curModel)
    for (j=0; j< *(*pars).p1; j++) {
      if (curModel[j]) {
	for (k=0, sumpartialres2=0; k<n; k++) { partialres[k]= res[k] + curCoef1[j] * ((*pars).x1[n*j+k]); sumpartialres2+= partialres[k]*partialres[k]; }
      } 
      MHTheta1pmom(&newdelta, &newcoef, &pinclude, &resupdate, res, partialres, &sumres2, &sumpartialres2, j, &nsel, curModel, curCoef1, &curPhi, pars);
      if (newdelta > curModel[j]) { nsel++; } else if (newdelta < curModel[j]) { nsel--; }
      curModel[j]= newdelta; curCoef1[j]= newcoef; 
      if (i>=0) margpp[j]+= pinclude;
      if (resupdate) { temp= partialres; partialres= res; res=temp; }
    }
    //Sample curCoef2
    for (k=0; k<n; k++) res[k]+= linpred2[k];
    simTheta2(curCoef2, res, &curPhi, pars);
    Avecx((*pars).x2, curCoef2, linpred2, 0, n-1, 0, p2-1);
    for (k=0, sumres2=0; k<n; k++) { res[k]-= linpred2[k]; sumres2+= res[k]*res[k]; }
    //Sample phi
    if (isbinary==0) { curPhi= simPhipmom(&nsel, curModel, curCoef1, curCoef2, &sumres2, pars); }
    //Sample tau
    if (*(*pars).priorTau1 !=0) { *(*pars).tau1= simTaupmom(&nsel, curModel, curCoef1, &curPhi, pars); }
    //Sample latent variables (only for probit model)
    if (isbinary) {
      Avecx((*pars).x1, curCoef1, linpred1, 0, n-1, 0, p1-1);  //update linpred1 (linpred2 already updated)
      sample_latentProbit((*pars).y,res,&sumres2,(*pars).ybinary,linpred1,linpred2,pars);
    }
    //Save values
    if ((i>0) && (i % (*thinning))==0) {
      for (j=0; j<p1; j++) {
        psn= niterthin*j+savecnt;
        postModel[psn]= curModel[j];
        postCoef1[psn]= curCoef1[j];
      }
      for (j=0; j<p2; j++) postCoef2[niterthin*j+savecnt]= curCoef2[j];
      postPhi[savecnt]= curPhi;
      if (*(*pars).priorTau1 !=0) postOther[savecnt]= *(*pars).tau1;
      savecnt++;
    }
    if ((*verbose ==1) && (i%niter10)==0) Rprintf(".");
  } //end MCMC for
  if (iupper>ilow) { for (j=0; j< p1; j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*verbose ==1) Rprintf("Done.\n");
  free_ivector(curModel,0,p1); free_dvector(curCoef1,0,p1); free_dvector(curCoef2,0,p2); free_dvector(linpred1,0,n); free_dvector(linpred2,0,n);
  free_dvector(res,0,n); free_dvector(partialres,0,n);
}


//Sample from the posterior of latent variables in probit model given the regression coefficients (i.e. the linear predictor)
//Input:
// - ybinary: response variable (1: success; 0: failure)
// - linpred1: linear predictor for current regression coefficients associated to variables under selection 
// - linpred2: linear predictor associated to adjustment variables
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - y: sampled values for the latent variables
// - res: sampled residuals, i.e. y - linpred1 - linpred2
// - sumres2: sum(res^2)
// - (*pars).ytX: updated t(y) %*% x1
// - (*pars).sumy2: updated sum(y^2)
void sample_latentProbit(double *y, double *res, double *sumres2, int *ybinary, double *linpred1, double *linpred2, struct modavgPars *pars) {
  int i;
  double linpred, plinpred, u;
  for (i=0, *sumres2=0, *(*pars).sumy2=0; i< *(*pars).n; i++) {
    linpred= linpred1[i] + linpred2[i];
    plinpred= pnormC(-linpred,0,1);
    if (ybinary[i]) {
      u= plinpred + (1.0-plinpred) * runif();  //u ~ Unif(plinpred,1)
    } else {
      u= plinpred * runif(); //u ~ Unif(0,plinpred)
    }
    res[i]= qnormC(u,0,1);
    (*sumres2)+= res[i]*res[i];
    y[i]= linpred + res[i];
    (*(*pars).sumy2)+= y[i]*y[i];
  }
  Atvecx((*pars).x1,y,(*pars).ytX,0,*(*pars).p1 -1,0,*(*pars).n -1); //update ytX=Xty
}

//Univariate MH update of (coef,delta)
//Input:
// - res: vector with residuals. Only used if curModel[j]==0
// - partialres: vector with partial residuals removing the variable from the model. Only used if curModel[j]==1
// - sumres2: sum(res*res)
// - sumpartialres2: sum(partialres*partialres)
// - j: index of the variable for which move is to be proposed
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - newdelta: new value
// - newcoef: new coefficient
// - pinclude: probability of including the variable in the model
// - res: on input, vector with residuals given the current coefficient value. On output, residuals given the updated coefficient value.
void MHTheta1pmom(int *newdelta, double *newcoef, double *pinclude, int *resupdate, double *res, double *partialres, double *sumres2, double *sumpartialres2, int j, int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int n= *(*pars).n, logscale=1, nsel0, nsel1, deltaprop, nu, i;
  double m1, *xj, m0, logbf, logpratio, thetaprop, m, S, propPars[5], lhood, lprior, lprop, lambda, num, den, sqrtPhi=sqrt(*curPhi);
  pt2modavgPrior priorFunction= NULL;
  *resupdate= 0;
  xj= (*pars).x1+j*n; //pointer to variable j in x1
  priorFunction= set_priorFunction_modavg((*pars).priorModel);
  //Propose delta
  if (curModel[j]) {
    m1= pmomMargKuniv(partialres, xj, sumpartialres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(partialres, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel -1; nsel1= *nsel;
  } else {
    m1= pmomMargKuniv(res, xj, sumres2, (*pars).colsumx1sq+j, &n, curPhi, (*pars).tau1, (*pars).r, &logscale);
    m0= dnormC_jvec(res, *(*pars).n, 0, sqrtPhi, 1);
    nsel0= *nsel; nsel1= *nsel +1;
  }
  logbf= m0-m1;
  logpratio= priorFunction(curModel, &nsel0, pars) - priorFunction(curModel, &nsel1, pars); //we use curModel in both cases as priorFunction currently only depends on nb vars
  *pinclude= 1.0/(1.0+exp(logbf+logpratio));
  if (runif() < *pinclude) { deltaprop=1; } else { deltaprop=0; }
  //Propose coef
  nu= sqrt(n);
  if ((curModel[j]==0) && (deltaprop==0)) {  //proposal is to keep variable out of the model
    *newdelta=0; *newcoef=0;
  } else {
    S= (*pars).colsumx1sq[j] + 1.0/(*(*pars).tau1);
    if (curModel[j]) { 
      for (i=0, m=0; i<n; i++) m+= xj[i]*partialres[i]; 
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, partialres, xj, &m1, &nu);
    } else { 
      for (i=0, m=0; i<n; i++) m+= xj[i]*res[i]; 
      m= m/S;
      proposalpmom(propPars, &m, &S, curPhi, (*pars).r, (*pars).tau1, &n, res, xj, &m1, &nu);
    }
    if (curModel[j] && deltaprop) {  //proposal is to keep variable in the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, lhood=0; i<n; i++) {
        partialres[i]-= thetaprop*xj[i];
        lhood+= dnormC(partialres[i],0,sqrtPhi,1) - dnormC(res[i],0,sqrtPhi,1);
      }
      lprior= dmomNorm(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1) - dmomNorm(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lprop= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) - dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1);
      lambda= exp(lhood+lprior+lprop);
    } else if ((curModel[j]==0) && deltaprop) { //proposal is to add variable to the model
      thetaprop= rtmixC(propPars, propPars+2, propPars+4, nu, 2);
      for (i=0, num=0; i<n; i++) {
        partialres[i]= res[i] - thetaprop*xj[i];
        num+= dnormC(partialres[i],0,sqrtPhi,1);
      }
      num+= dmomNorm(thetaprop,0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      den= dtmixC(thetaprop,propPars,propPars+2,propPars+4,nu,2,1) + m1;
      lambda= exp(num-den);
    } else if (curModel[j] && (deltaprop==0)) { // proposal is to drop variable from the model
      thetaprop=0;
      num= dtmixC(curCoef1[j],propPars,propPars+2,propPars+4,nu,2,1) + m1;
      for (i=0, den=0; i<n; i++) { den+= dnormC(res[i],0,sqrtPhi,1); }
      den+= dmomNorm(curCoef1[j],0,*(*pars).tau1,*curPhi,*(*pars).r,1);
      lambda= exp(num-den);
    }
    if (runif()<lambda) {
      *newdelta=deltaprop; *newcoef= thetaprop; 
      *resupdate= 1;  //signal that res and partialres have to be interchanged after exiting the function
      for (i=0, *sumres2=0; i< n; i++) (*sumres2)+= partialres[i]*partialres[i];
    } else {
      *newdelta=curModel[j]; *newcoef= curCoef1[j];
    }
  }
}

//Find parameters for univariate pmom proposal distribution (2 component mixture of T distributions)
//   Posterior: N(e; xj*theta; phi*I) * pmom(theta; phi, r, tau) / m1
//   Mixture: w1 * T_nu(theta;mu1,sigma21) + (1-w1) * T_nu(theta;mu2,sigma22)  (sigma21, sigma22 denote variances)
//Input:
// - m, S: posterior location & scale parameters
// - phi: residual variance
// - r: product MOM power parameter
// - tau1: product MOM prior dispersion parameter
// - e: response variable
// - xj: predictor
// - m1: normalization constant
// - nu: desired degrees of freedom
//Output: means in propPars[0:1], SD in propPars[2:3], weights in propPars[4:5]
void proposalpmom(double *propPars, double *m, double *S, double *phi, int *r, double *tau1, int *n, double *e, double *xj, double *m1, int *nu) {
  int i;
  double eps, fmode, sqrtPhi=sqrt(*phi), temp, temp2, doubler=2*(*r), ct2;
  //Find modes
  eps= sqrt((*m)*(*m) + 8.0*(*r)*(*phi)/(*S));
  propPars[0]= .5*(*m - eps); propPars[1]= .5*(*m + eps);
  //Find density at the mode
  for (i=0, fmode=0; i< *n; i++) fmode+= dnormC(e[i],propPars[1]*xj[i],sqrtPhi,1);
  fmode+= dmomNorm(propPars[1],0,*tau1,*phi,*r,1) - *m1;
  fmode= exp(fmode);
  //Proposal variances
  temp= (*S)/(*phi);
  propPars[2]= sqrt(1.0/(temp + doubler/(propPars[0]*propPars[0]))); propPars[3]= sqrt(1.0/(temp + doubler/(propPars[1]*propPars[1])));
  temp2= .5*(*nu); temp= temp2+.5;
  //Weights
  ct2= exp(gamln(&temp) - .5*log(*nu) - gamln(&temp2) - .5*log(M_PI*propPars[3]*propPars[3]));
  propPars[4]= max_xy(0,(fmode-ct2)/(dnormC(propPars[1],propPars[0],propPars[2],0) - ct2));
  propPars[5]= 1-propPars[4];
}

//Univariate marginal density under a product MOM prior (known variance case)
// integral N(y; x*theta, phi*I) * (theta^2/(tau*phi))^r * N(theta; 0; tau*phi) / (2r-1)!! d theta
// - y: response variable (must be a vector)
// - x: design matrix (must be a vector)
// - sumy2: sum(y*y)
// - n: length of y
// - phi: residual variance
// - tau: prior variance parameter
// - logscale: if set to 1 the log of the integral is returned
double pmomMargKuniv(double *y, double *x, double *sumy2, double *sumxsq, int *n, double *phi, double *tau, int *r, int *logscale) {
  int i; double ans, m, s, I, doubler=2.0*(*r);
  s= *sumxsq + 1.0/(*tau);
  for (i=0, m=0; i< *n; i++) { m+= y[i]*x[i]; }
  m/= s;
  I= log(mnorm(doubler,m,sqrt(*phi/s)));
  ans= I -.5*(*sumy2 - s*m*m)/(*phi) - .5*(*n)*log(2*M_PI*(*phi)) - .5*(log(s)+log(*tau)) - ldoublefact(doubler-1) - (*r)*log((*tau)*(*phi));
  if (*logscale ==0) ans=exp(ans);
  return(ans);
}

//Sample from conditional posterior of coefficients for variables not undergoing selection
//Input:
// - partialres: partial residuals obtained by removing adjustment variables from the model
// - phi: current value of the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output:
// - theta2: sample from conditional posterior of theta2 given the data and all other parameters
void simTheta2(double *theta2, double *partialres, double *phi, struct modavgPars *pars) {
  int i, j; double *tmp, *m, **cholS, sqrtPhi=sqrt(*phi);
  //m= S2inv * t(x2) * partialres
  tmp= dvector(0,*(*pars).p2); m= dvector(0,*(*pars).p2); cholS= dmatrix(1,*(*pars).p2,1,*(*pars).p2);
  Atvecx((*pars).x2, partialres, tmp, 0, *(*pars).p2 -1, 0, *(*pars).n -1); 
  Avecx((*pars).S2inv, tmp, m, 0, *(*pars).p2, 0, *(*pars).p2);
  //S= S2inv * phi
  for (i=0; i< *(*pars).p2; i++) { for (j=0; j< *(*pars).p2; j++) { cholS[i+1][j+1]= sqrtPhi * (*pars).cholS2inv[i+j*(*(*pars).p2)]; } }
  //Generate theta2 ~ N(m,S)
  rmvnormC(theta2-1,*(*pars).p2,m-1,cholS);
  free_dvector(tmp,0,*(*pars).p2); free_dvector(m,0,*(*pars).p2); free_dmatrix(cholS,1,*(*pars).p2,1,*(*pars).p2);
}


//Sample from conditional posterior of residual variance phi under a product MOM prior
//Input:
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curCoef2: current values of the regression coefficients for variables not undergoing selection
// - ssr: residual sum of squares for current coefficient values
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of phi given the data and all other parameteres
double simPhipmom(int *nsel, int *curModel, double *curCoef1, double *curCoef2, double *ssr, struct modavgPars *pars) {
  int i; double a, b, sumth1, sumth2;
  a= *(*pars).alpha + *(*pars).n + (2*(*(*pars).r)+1)*(*nsel) + *(*pars).p2;
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  for (i=0, sumth2=0; i< *(*pars).p2; i++) { sumth2+= curCoef2[i]*curCoef2[i]; }
  b= *(*pars).lambda + sumth1/(*(*pars).tau1) + sumth2/(*(*pars).tau2) + *ssr;
  return(1.0/rgammaC(.5*a,.5*b));
}

//Sample from conditional posterior of tau given all other parameters
//Input:
// - nsel: number of variables undergoing selection currently in the model
// - curModel: vector of 0's and 1's indicating which variables are currently in the model
// - curCoef1: current values of the regression coefficients for variables undergoing selection
// - curPhi: current value for the residual variance
// - pars: data, pre-computed quantities and prior parameters. See struct modavgPars for details.
//Output: random draw from conditional posterior of tau given the data and all other parameters
double simTaupmom(int *nsel, int *curModel, double *curCoef1, double *curPhi, struct modavgPars *pars) {
  int i; double a, b, sumth1;
  a= *(*pars).atau1 + (2*(*(*pars).r)+1)*(*nsel);
  for (i=0, sumth1=0; i< *(*pars).p1; i++) { if (curModel[i]==1) sumth1+= curCoef1[i]*curCoef1[i]; }
  b= *(*pars).btau1 + sumth1/(*curPhi);
  return(1.0/rgammaC(.5*a,.5*b));
}


//********************************************************************************************
// GENERAL MARGINAL DENSITY CALCULATION ROUTINES
//********************************************************************************************

void set_marginalPars(struct marginalPars *pars, int *n,int *p,double *y,double *sumy2,double *x,double *XtX,double *ytX,int *method,int *B,double *alpha,double *lambda,double *phi,double *tau,int *r,double *prDeltap,double *parprDeltap, int *logscale) {
  (*pars).n= n;
  (*pars).p= p;
  (*pars).y= y;
  (*pars).sumy2= sumy2;
  (*pars).x= x;
  (*pars).XtX= XtX;
  (*pars).ytX= ytX;
  (*pars).method= method;
  (*pars).B= B;
  (*pars).alpha= alpha;
  (*pars).lambda= lambda;
  (*pars).phi= phi;
  (*pars).tau= tau;
  (*pars).r= r;
  (*pars).prDeltap= prDeltap;
  (*pars).parprDeltap= parprDeltap;
  (*pars).logscale= logscale;
}

void set_f2opt_pars(double *m, double **S, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel) {
  f2opt_pars.m= m;
  f2opt_pars.S= S;
  f2opt_pars.sumy2= sumy2;
  f2opt_pars.XtX= XtX;
  f2opt_pars.ytX= ytX;
  f2opt_pars.alpha= alpha;
  f2opt_pars.lambda= lambda;
  f2opt_pars.phi= phi;
  f2opt_pars.tau= tau;
  f2opt_pars.r= r;
  f2opt_pars.n= n;
  f2opt_pars.p= p;
  f2opt_pars.sel= sel;
  f2opt_pars.nsel= nsel;
}

void set_f2int_pars(double *XtX, double *ytX, double *tau, int *n, int *p, int *sel, int *nsel, double *y, double *sumy2, int *method, int *B, double *alpha, double *lambda, int *logscale) {
  f2int_pars.XtX= XtX;
  f2int_pars.ytX= ytX;
  f2int_pars.tau= tau;
  f2int_pars.n= n;
  f2int_pars.p= p;
  f2int_pars.sel= sel;
  f2int_pars.nsel= nsel;
  f2int_pars.y= y;
  f2int_pars.sumy2= sumy2;
  f2int_pars.method= method;
  f2int_pars.B= B;
  f2int_pars.alpha= alpha;
  f2int_pars.lambda= lambda;
  f2int_pars.logscale= logscale;
}



double LogFactorial(int);

/* Get the next n-tuple in lexicographical order with values 0, 1, ..., base-1.
   Return true if successful, false if unable to increment further. */
int GetNextTuple(int* tuple, int n, int base)
{
  int j = 0;
  while (j < n && tuple[j] == base-1)
    {
      tuple[j] = 0;
      j++;
    }
  if (j < n)
    tuple[j]++;
  return (j < n);
}

int BinomialCoefficient(int power, int nu) {
  int ans;
  /* power can only be 2 or 4 */
  /* 0 <= nu <= power */
  if (power == 2) {
    ans= 1 + nu%2;
  } else if (power == 4) {
      if (nu == 0 || nu == 4)
	ans= 1;
      if (nu == 1 || nu == 3)
	ans= 4;
      if (nu == 2)
	ans= 6;
  } else ans= 0; /* should never happen */
  return ans;
}

/* NB: The equation in the paper is off by one.
   The following code, taken from the paper's definition of kappa,
   acually computes kappa + 1. */
double one_plus_kappa
(
 double dof, /* degrees of freedom */
 int r/* summation index */
 )
{
  double product = 1.0;
  int i;

  if (r == 0)
    return 1.0;

  for (i = 1; i <= r; i++)
    product *= 0.5*dof - i;
  return pow(0.5*dof- 1.0, r) / product;
}


/*
    Expectation of prod (x_i)^(2*power), where x_i ~ multivariate T_dof(mu,sigma). Set dof= -1 for x_i ~ N(mu,sigma)
    mu    = mean, a vector of length n
    sigma = covariance, n by n matrix
    dof   = degrees of freedom, -1 signals normal case
    power = 2 or 4, exponent on components
    Note: This function was written by John D. Cook
*/
double mvtexpect(const double* mu, const double* sigma, int n, int power, double dof)
{
  double product = 1.0;
  double sum = 0.0;
  double covariance_term, mean_term, temp;
  int index_sum;

  int s = power*n;
  int half_s = s/2; 
  int j, k, r, half_power;

  int* nu_index; /* indices nu_0 through nu_{n-1} */
  nu_index = ivector(0,n);

  for (r = 0; r <= half_s; r++)
    {

      for (j = 0; j < n; j++)
	nu_index[j] = 0;

      do
	{
	  product = 1.0;

	  index_sum = 0;
	  for (j = 0; j < n; j++)
	    index_sum += nu_index[j];
	  if (index_sum % 2)
	    product *= -1;

	  for (j = 0; j < n; j++)
	    product *= BinomialCoefficient(power, nu_index[j]); 

	  if (dof > 0.0)
	    product *= one_plus_kappa(dof, r);
	  /* else normal case and multiplicative term is 1 */

	  covariance_term = 0.0;
	  half_power = power/2;
	  for (j = 0; j < n; j++)
	    {
	      temp = 0.0; /* double_sum accumulates the product sigma*h */
	      for (k = 0; k < n; k++)
		{
		  /* sigma is stored by columns, so sigma_{ij} = sigma[i + j*n] */
		  temp += sigma[j + k*n]*(half_power - nu_index[k]);
		}
	      covariance_term += (half_power - nu_index[j])*temp;
	    }
	  product *= pow(0.5*covariance_term, r);

	  mean_term = 0.0;
	  for (j = 0; j < n; j++)
	    mean_term += (half_power - nu_index[j])*mu[j];
	  product *= pow(mean_term, s - 2*r); 
	  product /= exp(LogFactorial(r) + LogFactorial(s - 2*r));

	  sum += product;

	} while( GetNextTuple(nu_index, n, power+1) );
    }
  free_ivector(nu_index,0,n);
  return sum;
}

double LogFactorial(int n) {
  double ans;
  if (n > 254) {
      double x = n + 1;
      ans= (x - 0.5)*log(x) - x + 0.5*log(2*M_PI) + 1.0/(12.0*x);
  } else {
        double lf[] = {
            0.000000000000000,
            0.000000000000000,
            0.693147180559945,
            1.791759469228055,
            3.178053830347946,
            4.787491742782046,
            6.579251212010101,
            8.525161361065415,
            10.604602902745251,
            12.801827480081469,
            15.104412573075516,
            17.502307845873887,
            19.987214495661885,
            22.552163853123421,
            25.191221182738683,
            27.899271383840894,
            30.671860106080675,
            33.505073450136891,
            36.395445208033053,
            39.339884187199495,
            42.335616460753485,
            45.380138898476908,
            48.471181351835227,
            51.606675567764377,
            54.784729398112319,
            58.003605222980518,
            61.261701761002001,
            64.557538627006323,
            67.889743137181526,
            71.257038967168000,
            74.658236348830158,
            78.092223553315307,
            81.557959456115029,
            85.054467017581516,
            88.580827542197682,
            92.136175603687079,
            95.719694542143202,
            99.330612454787428,
            102.968198614513810,
            106.631760260643450,
            110.320639714757390,
            114.034211781461690,
            117.771881399745060,
            121.533081515438640,
            125.317271149356880,
            129.123933639127240,
            132.952575035616290,
            136.802722637326350,
            140.673923648234250,
            144.565743946344900,
            148.477766951773020,
            152.409592584497350,
            156.360836303078800,
            160.331128216630930,
            164.320112263195170,
            168.327445448427650,
            172.352797139162820,
            176.395848406997370,
            180.456291417543780,
            184.533828861449510,
            188.628173423671600,
            192.739047287844900,
            196.866181672889980,
            201.009316399281570,
            205.168199482641200,
            209.342586752536820,
            213.532241494563270,
            217.736934113954250,
            221.956441819130360,
            226.190548323727570,
            230.439043565776930,
            234.701723442818260,
            238.978389561834350,
            243.268849002982730,
            247.572914096186910,
            251.890402209723190,
            256.221135550009480,
            260.564940971863220,
            264.921649798552780,
            269.291097651019810,
            273.673124285693690,
            278.067573440366120,
            282.474292687630400,
            286.893133295426990,
            291.323950094270290,
            295.766601350760600,
            300.220948647014100,
            304.686856765668720,
            309.164193580146900,
            313.652829949878990,
            318.152639620209300,
            322.663499126726210,
            327.185287703775200,
            331.717887196928470,
            336.261181979198450,
            340.815058870798960,
            345.379407062266860,
            349.954118040770250,
            354.539085519440790,
            359.134205369575340,
            363.739375555563470,
            368.354496072404690,
            372.979468885689020,
            377.614197873918670,
            382.258588773060010,
            386.912549123217560,
            391.575988217329610,
            396.248817051791490,
            400.930948278915760,
            405.622296161144900,
            410.322776526937280,
            415.032306728249580,
            419.750805599544780,
            424.478193418257090,
            429.214391866651570,
            433.959323995014870,
            438.712914186121170,
            443.475088120918940,
            448.245772745384610,
            453.024896238496130,
            457.812387981278110,
            462.608178526874890,
            467.412199571608080,
            472.224383926980520,
            477.044665492585580,
            481.872979229887900,
            486.709261136839360,
            491.553448223298010,
            496.405478487217580,
            501.265290891579240,
            506.132825342034830,
            511.008022665236070,
            515.890824587822520,
            520.781173716044240,
            525.679013515995050,
            530.584288294433580,
            535.496943180169520,
            540.416924105997740,
            545.344177791154950,
            550.278651724285620,
            555.220294146894960,
            560.169054037273100,
            565.124881094874350,
            570.087725725134190,
            575.057539024710200,
            580.034272767130800,
            585.017879388839220,
            590.008311975617860,
            595.005524249382010,
            600.009470555327430,
            605.020105849423770,
            610.037385686238740,
            615.061266207084940,
            620.091704128477430,
            625.128656730891070,
            630.172081847810200,
            635.221937855059760,
            640.278183660408100,
            645.340778693435030,
            650.409682895655240,
            655.484856710889060,
            660.566261075873510,
            665.653857411105950,
            670.747607611912710,
            675.847474039736880,
            680.953419513637530,
            686.065407301994010,
            691.183401114410800,
            696.307365093814040,
            701.437263808737160,
            706.573062245787470,
            711.714725802289990,
            716.862220279103440,
            722.015511873601330,
            727.174567172815840,
            732.339353146739310,
            737.509837141777440,
            742.685986874351220,
            747.867770424643370,
            753.055156230484160,
            758.248113081374300,
            763.446610112640200,
            768.650616799717000,
            773.860102952558460,
            779.075038710167410,
            784.295394535245690,
            789.521141208958970,
            794.752249825813460,
            799.988691788643450,
            805.230438803703120,
            810.477462875863580,
            815.729736303910160,
            820.987231675937890,
            826.249921864842800,
            831.517780023906310,
            836.790779582469900,
            842.068894241700490,
            847.352097970438420,
            852.640365001133090,
            857.933669825857460,
            863.231987192405430,
            868.535292100464630,
            873.843559797865740,
            879.156765776907600,
            884.474885770751830,
            889.797895749890240,
            895.125771918679900,
            900.458490711945270,
            905.796028791646340,
            911.138363043611210,
            916.485470574328820,
            921.837328707804890,
            927.193914982476710,
            932.555207148186240,
            937.921183163208070,
            943.291821191335660,
            948.667099599019820,
            954.046996952560450,
            959.431492015349480,
            964.820563745165940,
            970.214191291518320,
            975.612353993036210,
            981.015031374908400,
            986.422203146368590,
            991.833849198223450,
            997.249949600427840,
            1002.670484599700300,
            1008.095434617181700,
            1013.524780246136200,
            1018.958502249690200,
            1024.396581558613400,
            1029.838999269135500,
            1035.285736640801600,
            1040.736775094367400,
            1046.192096209724900,
            1051.651681723869200,
            1057.115513528895000,
            1062.583573670030100,
            1068.055844343701400,
            1073.532307895632800,
            1079.012946818975000,
            1084.497743752465600,
            1089.986681478622400,
            1095.479742921962700,
            1100.976911147256000,
            1106.478169357800900,
            1111.983500893733000,
            1117.492889230361000,
            1123.006317976526100,
            1128.523770872990800,
            1134.045231790853000,
            1139.570684729984800,
            1145.100113817496100,
            1150.633503306223700,
            1156.170837573242400,
	  };
        ans= lf[n];
  }
  return ans;
}



//********************************************************************************************
// MODEL SELECTION ROUTINES
//********************************************************************************************

//modelSelectionC: Gibbs sampler for model selection in linear regression for several choices of prior distribution
//Input parameters
// - knownphi: is residual variance phi known?
// - priorCoef: 0 for product MOM, 1 for product iMOM, 2 for product eMOM
// - priorDelta: 0 for uniform, 1 for binomial, 2 for binomial with beta hyper-prior for success prob
// - niter: number of Gibbs iterations
// - ndeltaini: length of deltaini
// - deltaini: vector with indexes of covariates initially in the model (both deltaini and its indexes must be indexed at 0)
// - verbose: set verbose==1 to print iteration progress every 10% of the iterations
// - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
//Output
// - postSample: matrix with niter rows and p columns with posterior samples for covariate inclusion/exclusion (formatted as a vector in column order)
// - postOther: matrix with niter rows and nOther columns with posterior samples for other model parameters
// - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
// - postMode: model with highest posterior probability amongst all those visited
// - postModeProb: unnormalized posterior prob of posterior mode (log scale)
// - postProb: unnormalized posterior prob of each visited model (log scale)

SEXP modelSelectionCI(SEXP SpostSample, SEXP SpostOther, SEXP Smargpp, SEXP SpostMode, SEXP SpostModeProb, SEXP SpostProb, SEXP Sknownphi, SEXP SpriorCoef, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale);
  modelSelectionC(INTEGER(SpostSample), REAL(SpostOther), REAL(Smargpp), INTEGER(SpostMode), REAL(SpostModeProb), REAL(SpostProb), INTEGER(Sknownphi), INTEGER(SpriorCoef), INTEGER(SpriorDelta), INTEGER(Sniter), INTEGER(Sthinning), INTEGER(Sburnin), INTEGER(Sndeltaini), INTEGER(Sdeltaini), INTEGER(Sverbose), &pars);

  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;
}

void modelSelectionC(int *postSample, double *postOther, double *margpp, int *postMode, double *postModeProb, double *postProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *thinning, int *burnin, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, k, *sel, *selnew, *selaux, nsel, nselnew, niter10, niterthin, savecnt, ilow, iupper;
  double currentJ, newM, newP, newJ, ppnew, u;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi);
  priorFunction= set_priorFunction(prDelta);
 
  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose ==1) Rprintf("Running Gibbs sampler");
  niterthin= floor((*niter - *burnin +.0)/(*thinning+.0));
  if (*niter >10) { niter10= *niter/10; } else { niter10= 1; }
  for (j=0; j< *(*pars).p; j++) { margpp[j]= 0; }
  nsel= *ndeltaini;
  for (j=0; j< nsel; j++) { sel[j]= deltaini[j]; postSample[deltaini[j]*niterthin]= postMode[deltaini[j]]= 1; }
  if ((*prDelta)==2) { postOther[0]= *(*pars).prDeltap; }
  currentJ= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);
  postProb[0]= *postModeProb= currentJ;
  if (*burnin >0) { ilow=-(*burnin); savecnt=0; iupper= *niter - *burnin +1; } else { ilow=1; savecnt=1; iupper= *niter; } //if no burnin, start at i==1 & save initial value

  //Iterate
  for (i=ilow; i< iupper; i++) {
    for (j=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      newM= marginalFunction(selnew,&nselnew,pars);
      newP= priorFunction(selnew,&nselnew,pars);
      newJ= newM + newP;
      if (newJ > *postModeProb) {   //update posterior mode
        *postModeProb= newJ;
        for (k=0; k< *(*pars).p; k++) { postMode[k]= 0; }
        for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; } 
      }
      ppnew= 1.0/(1.0+exp(currentJ-newJ));
      if (i>=0) { if (nselnew>nsel) { margpp[j]+= ppnew; } else { margpp[j]+= (1-ppnew); } }
      u= runif();
      if (u<ppnew) {  //update model indicator
        selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; currentJ= newJ;
      }
    }  //end j for
    if ((i>0) && ((i%(*thinning))==0)) {
      if ((*prDelta)==2) { 
        *(*pars).prDeltap= rbetaC(nsel + (*pars).parprDeltap[0], *(*pars).p - nsel + (*pars).parprDeltap[1]);
        postOther[savecnt]= *(*pars).prDeltap;
      }
      for (j=0; j<nsel; j++) { postSample[sel[j]*niterthin+savecnt]= 1; }
      postProb[savecnt]= currentJ;
      savecnt++;
    }
    if ((*verbose==1) && ((i%niter10)==0)) { Rprintf("."); }
  }
  if (iupper>ilow) { for (j=0; j< *(*pars).p; j++) { margpp[j] /= (iupper-imax_xy(0,ilow)+.0); } } //from sum to average
  if (*verbose==1) Rprintf(" Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}

//greedyVarSelC: greedy search for posterior mode in variable selection.
//               Similar to Gibbs sampling, except that deterministic updates are made iff there is an increase in post model prob
//               The scheme proceeds until no variable is included/excluded or niter iterations are reached
// Input arguments: same as in modelSelectionC.

SEXP greedyVarSelCI(SEXP SpostMode, SEXP SothersMode, SEXP SpostModeProb, SEXP Sknownphi, SEXP SpriorCoef, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose) {
  int logscale=1;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars, INTEGER(Sn), INTEGER(Sp), REAL(Sy), REAL(Ssumy2), REAL(Sx), REAL(SXtX), REAL(SytX), INTEGER(Smethod), INTEGER(SB), REAL(Salpha),REAL(Slambda), REAL(Sphi), REAL(Stau), INTEGER(Sr), REAL(SprDeltap), REAL(SparprDeltap), &logscale);
  greedyVarSelC(INTEGER(SpostMode),REAL(SothersMode),REAL(SpostModeProb),INTEGER(Sknownphi),INTEGER(SpriorCoef),INTEGER(SpriorDelta),INTEGER(Sniter),INTEGER(Sndeltaini),INTEGER(Sdeltaini),INTEGER(Sverbose),&pars);
  PROTECT(ans = allocVector(REALSXP, 1));
  *REAL(ans)= 1.0;
  UNPROTECT(1);
  return ans;

}

void greedyVarSelC(int *postMode, double *othersMode, double *postModeProb, int *knownphi, int *prCoef, int *prDelta, int *niter, int *ndeltaini, int *deltaini, int *verbose, struct marginalPars *pars) {
  int i, j, *sel, *selnew, *selaux, nsel, nselnew, nchanges;
  double newJ, a, b;
  pt2margFun marginalFunction=NULL, priorFunction=NULL; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  marginalFunction= set_marginalFunction(prCoef, knownphi);
  priorFunction= set_priorFunction(prDelta);
  sel= ivector(0,*(*pars).p); selnew= ivector(0,*(*pars).p);

  //Initialize
  if (*verbose==1) Rprintf("Greedy searching posterior mode... ");
  for (j=0, nsel=*ndeltaini; j< nsel; j++) { sel[j]= deltaini[j]; postMode[deltaini[j]]= 1; }
  *postModeProb= marginalFunction(sel,&nsel,pars) + priorFunction(sel,&nsel,pars);

  //Iterate
  for (i=0, nchanges=1; (i< *niter) && (nchanges>0); i++) {
    for (j=0, nchanges=0; j< *(*pars).p; j++) {
      sel2selnew(j,sel,&nsel,selnew,&nselnew); //copy sel into selnew, adding/removing jth element
      newJ= marginalFunction(selnew,&nselnew,pars) + priorFunction(selnew,&nselnew,pars);
      if (newJ > *postModeProb) {
        *postModeProb= newJ;  //update post mode prob
        if (postMode[j]==0) { postMode[j]= 1; } else { postMode[j]= 0; }  //update post mode
        //for (k=0; k< *(*pars).p; k++) { postMode[k]= 0; }
        //for (k=0; k< nselnew; k++) { postMode[selnew[k]]= 1; } 
        selaux= sel; sel=selnew; selnew=selaux; nsel=nselnew; //update model indicator
        nchanges++;
      }
    } //end j for
    if ((*prDelta)==2) { 
      a= nsel + (*pars).parprDeltap[0]; b= *(*pars).p - nsel + (*pars).parprDeltap[1];
      if ((a>1) && (b>1)) { *(*pars).prDeltap= (a-1)/(a+b-2); } else { *(*pars).prDeltap= a/(a+b); }
    }
  }
  othersMode[0]= *(*pars).prDeltap;
  if (*verbose==1) Rprintf("Done.\n");

  free_ivector(sel,0,*(*pars).p); free_ivector(selnew,0,*(*pars).p);
}



void sel2selnew(int newelem,int *sel,int *nsel,int *selnew,int *nselnew) {
//Copy sel into selnew. 
// - If j in sel, don't copy it in selnew and set nselnew=nsel-1. 
// - If j not in sel, add it to selnew and set nselnew=nsel+1.
  int i, found;
  for (i=0, found=0; (i< *nsel) && (found==0); i++) { selnew[i]= sel[i]; found= (newelem==sel[i]); }
  if (found==0) { //add newelem
    selnew[i]= newelem; (*nselnew)= (*nsel)+1;
  } else {  //remove new elem
    for (i=i; i< *nsel; i++) { selnew[i-1]= sel[i]; }
    (*nselnew)= (*nsel)-1;
  }
}


//********************************************************************************************
// PRIORS ON MODEL SPACE (always return on log scale)
//********************************************************************************************

double unifPrior(int *sel, int *nsel, struct marginalPars *pars) { return 0.0; }
double unifPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) { return 0.0; }

//nsel ~ Binom(p,prDeltap)
double binomPrior(int *sel, int *nsel, struct marginalPars *pars) {
  return dbinomial(*nsel,*(*pars).p,*(*pars).prDeltap,1);
} 
double binomPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return dbinomial(*nsel,*(*pars).p1,(*pars).prModelpar[0],1);
} 

//nsel ~ Beta-Binomial(prModelpar[0],prModelPar[1])
double betabinPrior_modavg(int *sel, int *nsel, struct modavgPars *pars) {
  return bbPrior(*nsel,*(*pars).p1,(*pars).prModelpar[0],(*pars).prModelpar[1],1);
}

//*************************************************************************************
// PRODUCT MOM ROUTINES
//*************************************************************************************

double f2opt_mom(double *th) {
  return fmomNegC_non0(th+1,f2opt_pars.m+1,f2opt_pars.S,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.r,f2opt_pars.n,f2opt_pars.nsel);
}

//Note: th and m are assumed to be indexed at 0; S indexed at 1
double fmomNegC_non0(double *th, double *m, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i;
  double ans, sumlogth, *z;
  z= dvector(0,*nsel);
  for (i=0, sumlogth=0; i<(*nsel); i++) { sumlogth+= log(th[i]*th[i]); z[i]= th[i]-m[i]; }
  ans= .5*quadratic_xtAx(z-1,S,1,*nsel)/(*phi) - (*r +.0)*sumlogth;
  //ans= .5*quadratic_xtAselx(z, XtXplusct, p, nsel, sel)/(*phi) - (*r +.0)*sumlogth;
  free_dvector(z,0,*nsel);
  return ans;
}

void fppmomNegC_non0(double **ans, double *th, double **S, double *phi, double *tau, int *r, int *n, int *nsel) {
  int i, j;
  for (i=1; i<=(*nsel); i++) { ans[i][i]= S[i][i]/(*phi) + 2.0*(*r)/(th[i]*th[i]); }
  for (i=1; i<=(*nsel); i++) { for (j=i+1; j<=(*nsel); j++) { ans[i][j]= ans[j][i]= S[i][j]/(*phi); } }
}

void momIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *n, int *nsel, double *m, double **S, double *detS, double *phi, double *tau, int *r, int *logscale) {
  int i, emptyint, iter, maxit=50;
  double emptydouble, ftol= 1.0e-2, **Vopt, detVopt, **dirth;

  Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  set_f2opt_pars(m,S,&emptydouble,&emptydouble,&emptydouble,&emptydouble,&emptydouble,phi,tau,r,n,nsel,&emptyint,nsel);

  //Minimization
  for (i=1; i<=(*nsel); i++) { thopt[i]= m[i]; }  //init
  ddiag(dirth,1,*nsel);
  minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_mom, maxit);

  //Laplace approx
  fppmomNegC_non0(Vopt,thopt,S,phi,tau,r,n,nsel);
  invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);

  (*ILaplace)= -(*fopt) + .5*(log(*detS)-log(detVopt)- (*nsel)*log(*phi)) ;

  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
  free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
}


//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ N(m,Sinv)
double MC_mom_normal(double *m,double **Sinv,int *r,int *nsel, int *B) {
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvnormC(thsim,*nsel,m,cholSinv);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac); 
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}

//Monter Carlo evaluation of E(prod(z^(2*r))), where z ~ T_nu(m,Sinv)
double MC_mom_T(double *m,double **Sinv,int *nu,int *r,int *nsel, int *B) {
  int i;
  double **cholSinv, *thsim, ans, normfac;

  thsim= dvector(1,*nsel);
  cholSinv= dmatrix(1,*nsel,1,*nsel);
  choldc(Sinv,*nsel,cholSinv); //compute cholesky decomposition
  normfac= rsumlogsq(m,r,nsel);
  for (i=0, ans=0; i<(*B); i++) {
    rmvtC(thsim,*nsel,m,cholSinv,*nu);
    ans+= exp(rsumlogsq(thsim,r,nsel) - normfac); 
  }
  ans= log(ans/(*B +.0)) + normfac;

  free_dvector(thsim,1,*nsel);
  free_dmatrix(cholSinv,1,*nsel,1,*nsel);
  return ans;
}


// PRODUCT MOMENT MARGINALS
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - r: MOM power parameter
// - method==0 for Laplace; method==1 for Monte Carlo; method==2 for plug-in
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale) {
  struct marginalPars pars;
  double *rans, emptydouble=0;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),&emptydouble,REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),&emptydouble,&emptydouble,REAL(Sphi),REAL(Stau),INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale));
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalKC(INTEGER(Ssel),INTEGER(Snsel),&pars);
  UNPROTECT(1);
  return ans;
}

// Function to compute r * sum(log(th^2))
double rsumlogsq(double *th, int *r, int *nsel) {
  int i; double ans;
  for (i=1, ans=0; i<=(*nsel); i++) { ans+= log(th[i]*th[i]); }
  ans*= (*r);
  return(ans);
}

double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int i,j;
  double *m, s, **S, **Sinv, detS, num, den, logtau= log(*(*pars).tau), tauinv= 1.0/(*(*pars).tau), logphi= log(*(*pars).phi), ans, *thopt, **Voptinv, fopt;

  if (*nsel ==0) {
    m= dvector(1,1);
    m[1]=0; s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m[1],s,1);
    free_dvector(m,1,1);
  } else {
    m= dvector(1,*nsel);
    S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);

    num= -.5*(*(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel))/(*(*pars).phi);
    den= .5*((*(*pars).n +.0)*(LOG_M_2PI+logphi) + log(detS) + (*nsel)*logtau) + (*nsel)*(*(*pars).r)*(logtau+logphi+ldoublefact(2*(*(*pars).r)-1));
    if (*(*pars).method ==0) { //Laplace
      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,(*pars).phi,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);
    } else if (*(*pars).method ==1) { //MC
      for (i=1; i<=(*nsel); i++) { Sinv[i][i]= (*(*pars).phi)*Sinv[i][i]; for (j=i+1; j<=(*nsel); j++) { Sinv[i][j]=Sinv[j][i]= (*(*pars).phi)*Sinv[i][j]; } }
      ans= MC_mom_normal(m,Sinv,(*pars).r,nsel,(*pars).B);
    } else if (*(*pars).method ==2) { //Plug-in
      ans= rsumlogsq(m,(*pars).r,nsel);
    }
    ans+= num - den;
    free_dvector(m,1,*nsel);
    free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;  
}


SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  double *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,INTEGER(Sn),INTEGER(Sp),REAL(Sy),REAL(Ssumy2),REAL(Sx),REAL(SXtX),REAL(SytX),INTEGER(Smethod),INTEGER(SB),REAL(Salpha),REAL(Slambda),&emptydouble,REAL(Stau),INTEGER(Sr),&emptydouble,&emptydouble,INTEGER(Slogscale));
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pmomMarginalUC(INTEGER(Ssel), INTEGER(Snsel), &pars);
  UNPROTECT(1);
  return ans;
}


double pmomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, nu;
  double num, den, ans, term1, *m, **S, **Sinv, detS, *thopt, **Voptinv, fopt, phiadj, tauinv= 1.0/(*(*pars).tau), nuhalf, alphahalf=.5*(*(*pars).alpha), lambdahalf=.5*(*(*pars).lambda);
  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*LOG_M_PI + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
  } else {
    m= dvector(1,*nsel); S= dmatrix(1,*nsel,1,*nsel); Sinv= dmatrix(1,*nsel,1,*nsel);
    addct2XtX(&tauinv,(*pars).XtX,sel,nsel,(*pars).p,S);
    invdet_posdef(S,*nsel,Sinv,&detS);
    Asym_xsel(Sinv,*nsel,(*pars).ytX,sel,m);
    nuhalf= (*(*pars).r)*(*nsel) + .5*(*(*pars).n + *(*pars).alpha);
    nu= 2*nuhalf;

    num= gamln(&nuhalf) + alphahalf*log(lambdahalf) + nuhalf*(log(2) - log(*(*pars).lambda + *(*pars).sumy2 - quadratic_xtAx(m,S,1,*nsel)));
    den= (*nsel)*ldoublefact(2*(*(*pars).r)-1.0) + .5*(*(*pars).n * LOG_M_2PI + log(detS)) + (*nsel)*(.5 + *(*pars).r)*log(*(*pars).tau) + gamln(&alphahalf);
    if (*(*pars).method ==0) { //Laplace
      thopt= dvector(1,*nsel); Voptinv= dmatrix(1,*nsel,1,*nsel);
      phiadj= (nu+.0)/(nu-2.0);
      momIntegralApproxC(&ans,thopt,Voptinv,&fopt,(*pars).n,nsel,m,S,&detS,&phiadj,(*pars).tau,(*pars).r,(*pars).logscale);
      free_dvector(thopt,1,*nsel); free_dmatrix(Voptinv,1,*nsel,1,*nsel);
    } else if (*(*pars).method ==1) {  //MC
      term1= (*(*pars).lambda + *(*pars).sumy2 - quadratic_xseltAxsel((*pars).ytX,Sinv,1,nsel,sel))/(nu+.0);
      for (i=1; i<= *nsel; i++) { for (j=i; j<= *nsel; j++) { Sinv[i][j]= Sinv[j][i]= Sinv[i][j]*term1; } } //Vinv matrix
      ans= MC_mom_T(m,Sinv,&nu,(*pars).r,nsel,(*pars).B);
    } else if (*(*pars).method ==2) {  //Plug-in
      ans= rsumlogsq(m,(*pars).r,nsel);
    }
    ans+= num - den;
    free_dvector(m,1,*nsel); free_dmatrix(S,1,*nsel,1,*nsel); free_dmatrix(Sinv,1,*nsel,1,*nsel);
  }
  if (*(*pars).logscale !=1) { ans= exp(ans); }
  return ans;
}


//********************************************************************************************
// PRODUCT IMOM ROUTINES
//********************************************************************************************

//fimomNeg: minus log integrand of the function needed to compute product iMOM marginal density of the data conditional under a given linear model
//
// fimomNegC
// Input
// - th: theta value at which to evaluate the function (includes coef for all covariates, even those excluded in the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - n: sample size (length of y)
// - p: number of columns in XtX
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// Output: scalar evaluating the function for a single value of theta
double fimomNegC(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[sel[i]];
    th2= th[sel[i]] * th[sel[i]];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xseltAselxsel(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


double f2opt_imom(double *th) {
  double ans;
  ans= fimomNegC_non0(th+1,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.phi,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomNegC_non0(double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
//same as fimomNegC but loops over all elem in th (i.e. th has length *nsel and contains non-zero elements only). th is indexed at 0.
  int i;
  double ans, ytXth, sumlogth, suminvth, th2;
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(quadratic_xtAselx(th, XtX, p, nsel, sel) - 2*ytXth)/(*phi) + (*tau)*(*phi)*suminvth + sumlogth;
  return ans;
}


//Hessian of fimomNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)][1:(*nsel)])
// - th: th[1:(*nsel)] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomNegC_non0(double **ans, double *th, double *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2;

  for (i=1; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]]/(*phi) + 6.0*(*tau)*(*phi)/(th2*th2) - 2.0/th2;
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= XtX[sel[i-1]*(*p)+sel[j-1]]/(*phi);
    }
  }
}


void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *logscale) {
  int iter, maxit=50, emptyint;
  double **V, **Vinv, ftol= 1.0e-2, **dirth, **Vopt, detVopt, emptydouble=0, **emptymatrix;

  V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel); Vopt= dmatrix(1,*nsel,1,*nsel); dirth= dmatrix(1,*nsel,1,*nsel);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to XtX diagonal, store in V
  inv_posdef_upper(V,*nsel,Vinv);
  Asym_xsel(Vinv,*nsel,ytX,sel,thopt);  //product Vinv * selected elements in ytX
  ddiag(dirth,1,*nsel);
  set_f2opt_pars(&emptydouble,emptymatrix,&emptydouble,XtX,ytX,&emptydouble,&emptydouble,phi,tau,&emptyint,n,p,sel,nsel);
  //Minimization
  minimize(thopt, dirth, *nsel, ftol, &iter, fopt, f2opt_imom, maxit);

  //Laplace approx
  fppimomNegC_non0(Vopt,thopt,XtX,ytX,phi,tau,n,p,sel,nsel);
  invdet_posdef(Vopt,*nsel,Voptinv,&detVopt);
  (*ILaplace)= -(*fopt) - 0.5*log(detVopt);

  free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel); free_dmatrix(Vopt,1,*nsel,1,*nsel); free_dmatrix(dirth,1,*nsel,1,*nsel);
  free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
}


//Product iMOM marginal density for known phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - phi: residual variance
// - tau: prior dispersion parameter
// - method: method to approximate the marginal. method==0 for Laplace approximation (may underestimate the true value), method==1 for Importance Sampling Monte Carlo 
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale

SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *XtX=REAL(SXtX), *ytX=REAL(SytX), *phi=REAL(Sphi), *tau=REAL(Stau), *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,&emptydouble,XtX,ytX,method,B,&emptydouble,&emptydouble,phi,tau,&r,&emptydouble,&emptydouble,logscale);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalKC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  int one=1;
  double k, ans, m, s, ILaplace, *thopt, **Voptinv, fopt;
  thopt= dvector(1,*nsel);
  Voptinv= dmatrix(1,*nsel,1,*nsel);
  if ((*nsel)==0) {
    m= 0;
    s= sqrt(*(*pars).phi);
    ans= dnormC_jvec((*pars).y,*(*pars).n,m,s,1);
  } else {
    imomIntegralApproxC(&ILaplace,thopt,Voptinv,&fopt,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,&one);
    k= .5*((*nsel)*log(*(*pars).tau) - (*(*pars).sumy2)/(*(*pars).phi) - (*(*pars).n +.0) *LOG_M_2PI - (*(*pars).n - *nsel)*log(*(*pars).phi) - (*nsel)*LOG_M_PI);
    if ((*(*pars).method)==0) {
      ans= k + ILaplace;
    } else {
      ans= k + IS_imom(thopt,Voptinv,sel,nsel,(*pars).n,(*pars).p,(*pars).XtX,(*pars).ytX,(*pars).phi,(*pars).tau,(*pars).B);
    }
  }
  if ((*(*pars).logscale)!=1) { ans= exp(ans); }
  free_dvector(thopt,1,*nsel);
  free_dmatrix(Voptinv,1,*nsel,1,*nsel);
  return(ans);
}


//Evaluation of iMOM integral via Importance Sampling (result is returned in log-scale)
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, double *XtX, double *ytX, double *phi, double *tau, int *B) {
  int i,j;
  double *sdprop, **Vprop, *sopt, **cholVprop, **cholVpropinv, detVpropinv, *mprop, *thsim, *logr, maxlogr, ans;

  sdprop= dvector(1,*nsel); sopt= dvector(1,*nsel);
  mprop= dvector(1,*nsel); thsim= dvector(1, *nsel);
  logr= dvector(0,999);
  Vprop= dmatrix(1,*nsel,1,*nsel); cholVprop= dmatrix(1,*nsel,1,*nsel); cholVpropinv= dmatrix(1,*nsel,1,*nsel);

  for (i=1; i<=(*nsel); i++) {
    mprop[i]= 0;
    sopt[i]= sqrt(Voptinv[i][i]);
    sdprop[i]= .5*fabs(thopt[i] + 2*dsign(thopt[i])*sopt[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i; j<=(*nsel); j++) {
      Vprop[i][j]= Vprop[j][i]= sdprop[i]*sdprop[j]*Voptinv[i][j]/(sopt[i]*sopt[j]);
    }
  }
  choldc(Vprop,*nsel,cholVprop);
  choldc_inv(Vprop,*nsel,cholVpropinv);
  detVpropinv= choldc_det(cholVpropinv, *nsel);
  rmvtC(thsim, *nsel, mprop, cholVprop, 1);
  maxlogr= logr[0]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
  for (i=1;i<1000;i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    logr[i]= -fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1);
    if (logr[i]>maxlogr) { maxlogr= logr[i]; }
  }
  for (i=0, ans=0; i<1000; i++) { ans+= exp(logr[i]-maxlogr+500); }
  for (i=1000;i<(*B);i++) {
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    ans+= exp(-fimomNegC_non0(thsim+1,XtX,ytX,phi,tau,n,p,sel,nsel) - dmvtC(thsim,*nsel,mprop,cholVpropinv,detVpropinv,1,1) -maxlogr+500);
  }
  ans= log(ans/(.0+ (*B))) + maxlogr-500;

  free_dvector(sdprop,1,*nsel); free_dvector(sopt,1,*nsel);
  free_dvector(mprop, 1,*nsel); free_dvector(thsim, 1, *nsel);
  free_dvector(logr,0,999);
  free_dmatrix(Vprop,1,*nsel,1,*nsel); free_dmatrix(cholVprop,1,*nsel,1,*nsel); free_dmatrix(cholVpropinv,1,*nsel,1,*nsel);
  return(ans);
}




double f2opt_imomU(double *th) {
  //last element in th corresponds to eta=log(tau), i.e. log residual variance
  double ans;
  ans= fimomUNegC_non0(th+1,f2opt_pars.sumy2,f2opt_pars.XtX,f2opt_pars.ytX,f2opt_pars.alpha,f2opt_pars.lambda,f2opt_pars.tau,f2opt_pars.n,f2opt_pars.p,f2opt_pars.sel,f2opt_pars.nsel);
  return(ans);
}

double fimomUNegC_non0(double *th, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
//loops over all elem in th (i.e. th has length *nsel+1 and contains non-zero elements only). th is indexed at 0.
//Note: last element in th corresponds to eta=log(tau), i.e. log residual variance
  int i;
  double ans, ytXth, sumlogth, suminvth, th2, eta, phi;
  eta= th[*nsel]; phi= exp(eta);
  for (i=0, ytXth=0, sumlogth=0, suminvth=0; i<(*nsel); i++) {
    ytXth+= ytX[sel[i]] * th[i];
    th2= th[i] * th[i];
    suminvth+= 1/th2;
    sumlogth+= log(th2);
  }
  ans= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth + sumlogth + .5*eta*(*n - *nsel + *alpha);
  return ans;
}


//Hessian of fimomUNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)+1][1:(*nsel)+1])
// - th: th[1:(*nsel)+1] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void fppimomUNegC_non0(double **ans, double *th, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *n, int *p, int *sel, int *nsel) {
  int i, j;
  double th2, eta, phi, suminvth, ytXth, *XtXth;

  XtXth= dvector(1,*nsel);
  eta= th[*nsel +1]; phi= exp(eta);
  Asel_x(XtX,*p,th,*nsel,sel-1,XtXth);
  for (i=1, ytXth=0, suminvth=0; i<=(*nsel); i++) {
    th2= th[i]*th[i];
    ans[i][i]= XtX[sel[i-1]*(*p)+sel[i-1]]/phi + 6.0*(*tau)*phi/(th2*th2) - 2.0/th2;
    ans[i][*nsel+1]= ans[*nsel+1][i]= -2.0*(*tau)*phi/(th2*th[i]) - (XtXth[i]-ytX[sel[i-1]])/phi;
    ytXth+= ytX[sel[i-1]] * th[i];
    suminvth+= 1/(th[i]*th[i]);
  }
  for (i=1; i<=(*nsel); i++) {
    for (j=i+1; j<=(*nsel); j++) {
      ans[i][j]= ans[j][i]= XtX[sel[i-1]*(*p)+sel[j-1]]/phi;
    }
  }
  ans[*nsel+1][*nsel+1]= .5*(*lambda + *sumy2 - 2*ytXth + quadratic_xtAselx(th+1, XtX, p, nsel, sel))/phi + (*tau)*phi*suminvth;
  free_dvector(XtXth,1,*nsel);
}


void imomUIntegralApproxC(double *ILaplace, double *thopt, int *sel, int *nsel, int *n, int *p, double *sumy2, double *XtX, double *ytX, double *alpha, double *lambda, double *tau, int *logscale) {
  int iter, maxit=50, emptyint;
  double ftol= 1.0e-4, **dirth, **Vopt, **Voptinv, detVopt, emptydouble=0, **emptymatrix, fopt;

  Vopt= dmatrix(1,*nsel +1,1,*nsel +1); Voptinv= dmatrix(1,*nsel +1,1,*nsel +1);
  dirth= dmatrix(1,*nsel +1,1,*nsel +1);
  emptymatrix= dmatrix(1,1,1,1);
  //Initialize
  ddiag(dirth,1,*nsel +1);
  set_f2opt_pars(&emptydouble,emptymatrix,sumy2,XtX,ytX,alpha,lambda,&emptydouble,tau,&emptyint,n,p,sel,nsel);

  //Minimization
  minimize(thopt, dirth, *nsel +1, ftol, &iter, &fopt, f2opt_imomU, maxit);

  //Laplace approx
  fppimomUNegC_non0(Vopt,thopt,sumy2,XtX,ytX,alpha,lambda,tau,n,p,sel,nsel);
  invdet_posdef(Vopt,*nsel +1,Voptinv,&detVopt);
  (*ILaplace)= -fopt - 0.5*log(detVopt) + .5*(*nsel)*log(2.0*(*tau));

  free_dmatrix(Vopt,1,*nsel +1,1,*nsel +1); free_dmatrix(Voptinv,1,*nsel +1,1,*nsel +1); 
  free_dmatrix(dirth,1,*nsel +1,1,*nsel+1); free_dmatrix(emptymatrix,1,1,1,1);
  if ((*logscale)!=1) { (*ILaplace)= exp(*ILaplace); }
}



//Product iMOM marginal density for unknown phi
// Input:
// - sel: model indicator. Vector of length p indicating the index of the variables in the model (starting the indexing at 0)
// - nsel: length of sel
// - n: sample size (length of y)
// - p: number of columns in XtX
// - y: observed response vector (length n)
// - sumy2: sum of y*y
// - x: design matrix (includes all covariates, even those excluded under the current model)
// - XtX: X'X where X is the design matrix (includes all covariates, even those excluded under the current model)
// - ytX: vector of length p containing y'X (where y is the length n response vector)
// - tau: prior dispersion parameter
// - method: method to approximate the integral for known phi. 0 for Laplace approx which may underestimate true value, 1 for exact evaluation which can be very computationally expensive. 2 ('Hybrid') integrates wrt phi numerically (Romberg) and wrt theta via Laplace approx (Laplace error is adjusted based on exact evaluation for single phi value)
// - B: number of Monte Carlo samples. Ignored unless method==1.
// - logscale: if set to 1 result is returned in log scale
// - alpha, lambda: prior for phi (residual variance) is Inverse Gamma (.5*alpha,.5*lambda)
double f2int_imom(double phi) {
  double ans, *inputphi= f2int_pars.phi;
  f2int_pars.phi= &phi;
  ans= pimomMarginalKC(f2int_pars.sel,f2int_pars.nsel,&f2int_pars) * dinvgammaC(phi,.5*(*f2int_pars.alpha),.5*(*f2int_pars.lambda));
  f2int_pars.phi= inputphi;
  return(ans);
}

SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda) {
  int *sel=INTEGER(Ssel), *nsel=INTEGER(Snsel), *n=INTEGER(Sn), *p=INTEGER(Sp), *method=INTEGER(Smethod), *B=INTEGER(SB), *logscale=INTEGER(Slogscale), r=1;
  double *y=REAL(Sy), *sumy2=REAL(Ssumy2), *x=REAL(Sx), *XtX=REAL(SXtX), *ytX=REAL(SytX), *tau=REAL(Stau), *alpha=REAL(Salpha), *lambda=REAL(Slambda), *rans, emptydouble=0;
  struct marginalPars pars;
  SEXP ans;

  set_marginalPars(&pars,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,&emptydouble,tau,&r,&emptydouble,&emptydouble,logscale);
  PROTECT(ans = allocVector(REALSXP, 1));
  rans = REAL(ans);
  *rans= pimomMarginalUC(sel, nsel, &pars);
  UNPROTECT(1);
  return ans;
}


double pimomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  int i, j, zero=0, one=1;
  double ans, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, adj, *inputphi, num, den, term1, alphahalf=.5*(*(*pars).alpha);

  if (*nsel ==0) {
    term1= .5*(*(*pars).n + *(*pars).alpha);
    num= .5*(*(*pars).alpha)*log(*(*pars).lambda) + gamln(&term1);
    den= .5*(*(*pars).n)*LOG_M_PI + gamln(&alphahalf);
    ans= num -den - term1*log(*(*pars).lambda + *(*pars).sumy2);
    if ((*(*pars).logscale)!=1) ans= exp(ans);
  } else {
    if ((*(*pars).method)==1) {  //MC + Univariate integration
      set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
      adj= (qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
      if ((*(*pars).logscale)==1) ans= log(ans);
    } else {                     //Laplace or Hybrid-Laplace MC
      V= dmatrix(1,*nsel,1,*nsel); 
      Vinv= dmatrix(1,*nsel,1,*nsel);
      thest= dvector(1,*nsel+1);
     
      addct2XtX((*pars).tau,(*pars).XtX,sel,nsel,(*pars).p,V); //add tau to diagonal elem of XtX
      inv_posdef_upper(V,*nsel,Vinv);
      Asym_xsel(Vinv,*nsel,(*pars).ytX,sel,thest);
      for (i=0, sumer2=0; i<(*(*pars).n); i++) {
        for (j=1, ypred=0; j<=(*nsel); j++) { ypred += (*pars).x[i + (*(*pars).n)*sel[j-1]] * thest[j]; }
        er= (*pars).y[i] - ypred;
        sumer2+= er*er;
      }
      phiest= (sumer2 + (*(*pars).lambda))/(*(*pars).alpha + *(*pars).n);
      if ((*(*pars).method)==0) {  //Laplace
        thest[*nsel +1]= log(phiest);
        imomUIntegralApproxC(&ans,thest,sel,nsel,(*pars).n,(*pars).p,(*pars).sumy2,(*pars).XtX,(*pars).ytX,(*pars).alpha,(*pars).lambda,(*pars).tau,&one);
	ans= ans + alphahalf*log(.5*(*(*pars).lambda)) - .5*(*(*pars).n)*LOG_M_2PI - gamln(&alphahalf);
        if ((*(*pars).logscale)!=1) ans= exp(ans);
      } else if ((*(*pars).method)==2) {  //Hybrid Laplace - MC - Univariate integration
        set_f2int_pars((*pars).XtX,(*pars).ytX,(*pars).tau,(*pars).n,(*pars).p,sel,nsel,(*pars).y,(*pars).sumy2,(*pars).method,(*pars).B,(*pars).alpha,(*pars).lambda,&zero);
        inputphi= (*pars).phi; (*pars).phi= &phiest; 
        (*(*pars).method)= 1; //IS evaluation of marginal for phi=phiest
        intmc= pimomMarginalKC(sel, nsel, pars); 
        (*(*pars).method)= 0; //Laplace approx for phi=phiest
        intlapl= pimomMarginalKC(sel, nsel, pars); 
        (*pars).phi= inputphi; (*(*pars).method)= 2;  //reset input values for phi, method
        if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
        adj= intmc/intlapl;
        f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx
        ans= adj * (qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
        if ((*(*pars).logscale)==1) ans= log(ans);
      }
      free_dmatrix(V,1,*nsel,1,*nsel); 
      free_dmatrix(Vinv,1,*nsel,1,*nsel);
      free_dvector(thest,1,*nsel+1);
    }
  }
  return(ans);
}

/*
double pimomMarginalUC(int *sel, int *nsel, int *n, int *p, double *y, double *sumy2, double *x, double *XtX, double *ytX, double *tau, int *method, int *B, int *logscale, double *alpha, double *lambda) {
  int i, j, zero=0, one=1;
  double ans, er, sumer2, **V, **Vinv, *thest, ypred, phiest, intmc, intlapl, adj;

  set_f2int_pars(XtX,ytX,tau,n,p,sel,nsel,y,sumy2,method,B,alpha,lambda);
  if ((*method)==2) {
    V= dmatrix(1,*nsel,1,*nsel); Vinv= dmatrix(1,*nsel,1,*nsel);
    thest= dvector(1,*nsel);

    addct2XtX(tau,XtX,sel,nsel,p,V); //add tau to diagonal elem of XtX
    inv_posdef_upper(V,*nsel,Vinv);
    Asym_xsel(Vinv,*nsel,ytX,sel,thest);
    for (i=0, sumer2=0; i<(*n); i++) {
      for (j=1, ypred=0; j<=(*nsel); j++) { ypred += x[i + (*n)*sel[j-1]] * thest[j]; }
      er= y[i] - ypred;
      sumer2+= er*er;
    }
    phiest= (sumer2 + (*lambda))/(*alpha + *n);
    intmc= pimomMarginalKC(sel,nsel,n,p,y,sumy2,XtX,ytX,&phiest,tau,&one,B,&zero); //IS evaluation of marginal for phi=phiest
    intlapl= pimomMarginalKC(sel,nsel,n,p,y,sumy2,XtX,ytX,&phiest,tau,&zero,B,&zero); //Laplace approx for phi=phiest
    if (intlapl==0) { intmc+= 1.0e-300; intlapl+= 1.0e-300; } //avoid numerical zero
    adj= intmc/intlapl;
    f2int_pars.method= &zero;  //set method to eval marginal for known phi to Laplace approx

    free_dmatrix(V,1,*nsel,1,*nsel); free_dmatrix(Vinv,1,*nsel,1,*nsel);
    free_dvector(thest,1,*nsel);
  } else {
    adj= 1.0;
  }

  ans= adj * (qromo(f2int_imom,0.0,100,midpnt) + qromo(f2int_imom,100,1.0e30,midinf));
  if ((*logscale)==1) ans= log(ans);
  return(ans);
}
*/



//*************************************************************************************
// Product eMOM routines
//*************************************************************************************

double pemomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}

double pemomMarginalUC(int *sel, int *nsel, struct marginalPars *pars) {
  return 0.0;
}

