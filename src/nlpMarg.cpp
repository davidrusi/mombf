// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//Include other headers
#include "cstat.h"
#include "crossprodmat.h"
#include "modelSel.h"
#include "modselIntegrals.h"
#include "modselFunction.h"
#include "Polynomial.h"
#include "nlpMarg.h"

#include <map>
#include <string>

SEXP nlpMarginalCI(SEXP Ssel, SEXP Snsel, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Slogscale) {
  int i, j, idxj, nuncens, emptyint=0, *isgroup, *nconstraints, *ninvconstraints, ngroupsconstr=0;
  double *rans, *ytXuncens=NULL, emptydouble=0;
  intptrvec constraints, invconstraints;
  crossprodmat *XtX, *XtXuncens=NULL;
  struct marginalPars pars;
  SEXP ans;

  PROTECT(ans = Rf_allocVector(REALSXP, 1));
  rans = REAL(ans);

  isgroup= ivector(0, INTEGER(Sp)[0]);
  nconstraints= ivector(0,INTEGER(Sngroups)[0]); ninvconstraints= ivector(0,INTEGER(Sngroups)[0]);
  countConstraints(nconstraints, &constraints, ninvconstraints, &invconstraints, &ngroupsconstr, isgroup, INTEGER(Sngroups), INTEGER(Snvaringroup), Sconstraints, Sinvconstraints);

  XtX= new crossprodmat(REAL(SXtX),INTEGER(Sn)[0],INTEGER(Sp)[0],true);
  if (LENGTH(Suncens)>0) { //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
    int n=INTEGER(Sn)[0], *uncens= INTEGER(Suncens);
    double *pty= REAL(Sy), *ptx= REAL(Sx);
    for (nuncens=0; (nuncens<n) && (uncens[nuncens]==1); nuncens++) { }
    XtXuncens= new crossprodmat(REAL(Sx), INTEGER(Sn)[0], INTEGER(Sp)[0], false, nuncens, 0);
    ytXuncens= dvector(0,INTEGER(Sp)[0]);
    for (j=0; j< INTEGER(Sp)[0]; j++) { for (i=0, ytXuncens[j]=0, idxj=j*n; i< nuncens; i++) { ytXuncens[j] += pty[i] * ptx[i + idxj]; } }
  } else { nuncens= INTEGER(Sn)[0]; }

  set_marginalPars(&pars, INTEGER(Sn), &nuncens, INTEGER(Sp), REAL(Sy), INTEGER(Suncens), REAL(Ssumy2), REAL(Sx), XtX, REAL(SytX), INTEGER(Smethod), INTEGER(Shesstype), INTEGER(SoptimMethod), &emptyint, &emptydouble, INTEGER(SB), REAL(Salpha),REAL(Slambda), &emptydouble, REAL(Stau), REAL(Staugroup), REAL(Staualpha), REAL(Sfixatanhalpha), INTEGER(Sr), &emptydouble, &emptydouble, &emptydouble, &emptydouble, INTEGER(Slogscale), &emptydouble, INTEGER(Sgroups), 0, INTEGER(Sngroups), 0, INTEGER(Snvaringroup), 0, 0, XtXuncens,ytXuncens);
  *rans= nlpMarginal(INTEGER(Ssel), INTEGER(Snsel), INTEGER(Sfamily), INTEGER(SpriorCoef), INTEGER(SpriorGroup), &pars);
  delete XtX;
  UNPROTECT(1);
  return ans;
}

double nlpMarginal(int *sel, int *nsel, int *family, int *prCoef, int *prGroup, struct marginalPars *pars) {
  int priorcode, knownphi=0;
  double ans;
  pt2margFun marginalFunction; //same as double (*marginalFunction)(int *, int *, struct marginalPars *);

  priorcode = mspriorCode(prCoef, prGroup, pars);
  marginalFunction = set_marginalFunction(&priorcode, &knownphi, family, pars);
  ans = marginalFunction(sel, nsel, pars);
  return(ans);
}
