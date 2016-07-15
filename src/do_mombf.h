/*
 * DO_MOMBF.H - Public include for .Call interface to model selection routines
 */

#ifndef DO_MOMBF_H
#define DO_MOMBF_H      1

#include <Rdefines.h>

/*
 * Function Prototypes
 */
extern "C" {

  SEXP eprod_I(SEXP m, SEXP S, SEXP n, SEXP power, SEXP dof);

  SEXP pmomLM_I(SEXP postModel, SEXP margpp, SEXP postCoef1, SEXP postCoef2, SEXP postPhi, SEXP postOther, SEXP niter, SEXP thinning, SEXP burnin, SEXP niniModel, SEXP iniModel, SEXP iniCoef1, SEXP iniCoef2, SEXP iniPhi, SEXP iniOthers, SEXP verbose, SEXP n, SEXP p1, SEXP p2, SEXP isbinary, SEXP ybinary, SEXP y, SEXP sumy2, SEXP x1, SEXP x2, SEXP XtX, SEXP ytX, SEXP cholS2, SEXP S2inv, SEXP cholS2inv, SEXP colsumx1sq, SEXP alpha, SEXP lambda, SEXP priorCoef, SEXP r, SEXP tau1, SEXP tau2, SEXP priorTau1, SEXP atau1, SEXP btau1, SEXP priorModel, SEXP prModelpar);

  //Model selectio
  SEXP modelSelectionCI(SEXP SpostSample, SEXP SpostOther, SEXP Smargpp, SEXP SpostMode, SEXP SpostModeProb, SEXP SpostProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sniter, SEXP Sthinning, SEXP Sburnin, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shess, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose);

  SEXP greedyVarSelCI(SEXP SpostMode, SEXP SpostModeProb, SEXP Sknownphi, SEXP Sfamily, SEXP SpriorCoef, SEXP Sniter, SEXP Sndeltaini, SEXP Sdeltaini, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shess, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Sphi, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP SpriorDelta, SEXP SprDeltap, SEXP SparprDeltap, SEXP Sverbose);

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
  SEXP nlpMarginalSkewNormI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP Smethod, SEXP SoptimMethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP SprCoef);
  SEXP nlpMarginalAlaplI(SEXP Ssel, SEXP Snsel, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Stau, SEXP Staualpha, SEXP Sr, SEXP Ssymmetric, SEXP Smethod, SEXP Shess, SEXP SoptimMethod, SEXP SB, SEXP Slogscale, SEXP Salpha, SEXP Slambda, SEXP SprCoef);
}

#endif /* DO_MOMBF_H */
