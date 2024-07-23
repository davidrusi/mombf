// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//Include other headers
#include "cstat.h"
#include "crossprodmat.h"
#include "modelSel_regression.h"
#include "modselIntegrals.h"
#include "modselFunction.h"
#include "Polynomial.h"
#include "nlpMarg.h"

#include <map>
#include <string>

//Syntax calling from R, see nlpMarginal.R
//nlpMarginalCI(sel,       nsel,    familyint,           prior,          priorgr,       n,       p,       y,       uncens,       sumy2,       x,       colsumsx,       XtX,       ytX,       method,       hesstype,       optimMethod, optim_maxit,       B,       alpha,       lambda,       tau,       taugroup,       taualpha,       fixatanhalpha,       r,       groups,       ngroups,       nvaringroup,       constraints,      invconstraints, Dmat, logscale)

// [[Rcpp::export]]
SEXP nlpMarginalCI(SEXP Sknownphi, SEXP Ssel, SEXP Snsel, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Ssumy, SEXP Ssumlogyfact, SEXP Sx, SEXP Scolsumsx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Sadjoverdisp, SEXP Shesstype, SEXP SoptimMethod, SEXP Soptim_maxit, SEXP Sthinit, SEXP Susethinit, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP Sa, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP SDmat, SEXP Slogscale) {
  int i, j, idxj, nuncens, *isgroup, *nconstraints, *ninvconstraints, ngroupsconstr=0, p= INTEGER(Sp)[0], usethinit= INTEGER(Susethinit)[0], maxvars_empty= -1;
  double *rans, *ytXuncens=NULL, emptydouble=0, *thinit;
  intptrvec constraints, invconstraints;
  crossprodmat *XtX, *XtXuncens=NULL, *Pmat;
  lmObject *lm;
  //struct marginalPars pars;
  SEXP ans;

  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);

  isgroup= ivector(0, p);
  nconstraints= ivector(0,INTEGER(Sngroups)[0]); ninvconstraints= ivector(0,INTEGER(Sngroups)[0]);
  countConstraints(nconstraints, &constraints, ninvconstraints, &invconstraints, &ngroupsconstr, isgroup, INTEGER(Sngroups), INTEGER(Snvaringroup), Sconstraints, Sinvconstraints);

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],p,true);
  Pmat= new crossprodmat(REAL(SDmat), p, p, false);

  if (LENGTH(Suncens)>0) { //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
    int n=INTEGER(Sn)[0], *uncens= INTEGER(Suncens);
    double *pty= REAL(Sy), *ptx= REAL(Sx);
    for (nuncens=0; (nuncens<n) && (uncens[nuncens]==1); nuncens++) { }
    XtXuncens= new crossprodmat(REAL(Sx), INTEGER(Sn)[0], p, false, nuncens, 0);
    ytXuncens= dvector(0,p);
    for (j=0; j< p; j++) { for (i=0, ytXuncens[j]=0, idxj=j*n; i< nuncens; i++) { ytXuncens[j] += pty[i] * ptx[i + idxj]; } }
  } else { nuncens= INTEGER(Sn)[0]; }

  thinit= dvector(0, p);
  if (usethinit != 3) {
    for (j=0; j<= p; j++) { thinit[j]= 0; }
  } else {
    for (j=0; j<= p; j++) { thinit[j]= REAL(Sthinit)[j]; }
  }

  lm= new lmObject(INTEGER(SpriorCoef), INTEGER(SpriorGroup), INTEGER(Sfamily), INTEGER(Sn), &nuncens, INTEGER(Sp), REAL(Sy), INTEGER(Suncens), REAL(Ssumy2), REAL(Ssumy), REAL(Ssumlogyfact), REAL(Sx), REAL(Scolsumsx), XtX, REAL(SytX), INTEGER(Smethod), INTEGER(Sadjoverdisp), INTEGER(Shesstype), INTEGER(SoptimMethod), INTEGER(Soptim_maxit), &usethinit, thinit, INTEGER(SB), REAL(Salpha),REAL(Slambda), INTEGER(Sknownphi), &emptydouble, REAL(Stau), REAL(Staugroup), REAL(Staualpha), REAL(Sfixatanhalpha), INTEGER(Sr), REAL(Sa), REAL(SDmat), Pmat, &emptydouble, &emptydouble, &emptydouble, &emptydouble, &maxvars_empty, INTEGER(Slogscale), &emptydouble, INTEGER(Sgroups), isgroup, INTEGER(Sngroups), 0, INTEGER(Snvaringroup), 0, 0, XtXuncens,ytXuncens);

  arma::SpMat<short> *model;
  model= new arma::SpMat<short>(INTEGER(Sp)[0],1);
  ivector_to_spmat(INTEGER(Ssel), INTEGER(Snsel), model); //Convert integer vector storing indexes of non-zero rows into arma::SpMat with 1's in those indexes
  (*rans)= nlpMarginal(model, lm);

  delete model;
  delete lm;
  //delete_marginalPars(&pars);
  delete XtX;
  delete Pmat;
  free_dvector(thinit, 0, p);
  UNPROTECT(1);
  return ans;
}

double nlpMarginal(arma::SpMat<short> *model, lmObject *lm) {
  double ans;
  pt2margFun marginalFunction; //same as double (*marginalFunction)(int *, int *, lmObject *);

  marginalFunction = set_marginalFunction(lm);
  ans = marginalFunction(model, lm, nullptr, nullptr, nullptr, nullptr);
  return(ans);
}
