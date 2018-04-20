/*
 * DO_MOMBF.H - Public include for .Call interface to model selection routines
 */

#ifndef DO_MOMBF_H
#define DO_MOMBF_H      1

#include <Rdefines.h>
#include "R_ext/Rdynload.h"



/*
 * Function Prototypes
 */
extern "C" {

  SEXP bsplineCI(SEXP x, SEXP degree, SEXP knots);

  SEXP eprod_I(SEXP m, SEXP S, SEXP n, SEXP power, SEXP dof);

  SEXP pmomLM_I(SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP XtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar);

  //Model selection
  SEXP modelSelectionEnumCI(SEXP Snmodels, SEXP Smodels, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose);

  SEXP modelSelectionGibbsCI(SEXP SpostModeini, SEXP SpostModeiniProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose);

  SEXP greedyVarSelCI(SEXP Sknownphi, SEXP SpriorCoef, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sincludevars, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose);

  //Non-local prior marginal likelihoods
  SEXP pmomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale);
  SEXP pmomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Sr, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda);

  SEXP pimomMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale);
  SEXP pimomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda);

  SEXP pemomMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Smethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda);

  //Integrated likelihood for linear models under Zellner's prior
  SEXP zellnerMarginalKI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP SXtX, SEXP SytX, SEXP Sphi, SEXP Stau, SEXP Slogscale);
  SEXP zellnerMarginalUI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Slogscale, SEXP Salpha, SEXP Slambda);

  //Integrated likelihoods under skew normal, laplace and asymmetric laplace residuals
  SEXP nlpMarginalSkewNormI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP Smethod, SEXP SoptimMethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP SprCoef);
  SEXP nlpMarginalAlaplI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP Ssymmetric, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP SprCoef);
}

static R_CallMethodDef callMethods[]  = {
  {"bsplineCI", (DL_FUNC) &bsplineCI, 3},
  {"eprod_I", (DL_FUNC) &eprod_I, 5},
  {"pmomLM_I", (DL_FUNC) &pmomLM_I, 36},
  {"modelSelectionEnumCI", (DL_FUNC) &modelSelectionEnumCI, 27},
  {"modelSelectionGibbsCI", (DL_FUNC) &modelSelectionGibbsCI, 32},
  {"greedyVarSelCI", (DL_FUNC) &greedyVarSelCI, 27},
  {"pmomMarginalKI", (DL_FUNC) &pmomMarginalKI, 14},
  {"pmomMarginalUI", (DL_FUNC) &pmomMarginalUI, 16},
  {"pimomMarginalKI", (DL_FUNC) &pimomMarginalKI, 13},
  {"pimomMarginalUI", (DL_FUNC) &pimomMarginalUI, 15},
  {"pemomMarginalUI", (DL_FUNC) &pemomMarginalUI, 15},
  {"zellnerMarginalKI", (DL_FUNC) &zellnerMarginalKI, 11},
  {"zellnerMarginalUI", (DL_FUNC) &zellnerMarginalUI, 13},
  {"nlpMarginalSkewNormI", (DL_FUNC) &nlpMarginalSkewNormI, 20},
  {"nlpMarginalAlaplI", (DL_FUNC) &nlpMarginalAlaplI, 22},
  {NULL, NULL, 0}
};

void R_init_mombf(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
}


#endif /* DO_MOMBF_H */
