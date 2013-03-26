/*
 * DO_MOMBF.H - Public include for .Call interface to model selection routines
 */

#ifndef DO_MOMBF_H
#define DO_MOMBF_H      1

#ifdef __cplusplus
extern "C" {
#endif

#include <Rdefines.h>


/*
 * Function Prototypes
 */
extern SEXP pmomLM_I(
  SEXP vnt_postModel,
  SEXP vnt_margpp,
  SEXP vnt_postCoef1,
  SEXP vnt_postCoef2,
  SEXP vnt_postPhi,
  SEXP vnt_postOther,
  SEXP vnt_niter,
  SEXP vnt_thinning,
  SEXP vnt_burnin,
  SEXP vnt_niniModel,
  SEXP vnt_iniModel,
  SEXP vnt_iniCoef1,
  SEXP vnt_iniCoef2,
  SEXP vnt_iniPhi,
  SEXP vnt_iniOthers,
  SEXP vnt_verbose,        /* :TODO: Move to end */
  SEXP vnt_n,
  SEXP vnt_p1,
  SEXP vnt_p2,
  SEXP vnt_isbinary,
  SEXP vnt_ybinary,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_x1,
  SEXP vnt_x2,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_cholS2,
  SEXP vnt_S2inv,
  SEXP vnt_cholS2inv,
  SEXP vnt_colsumx1sq,
  SEXP vnt_alpha,
  SEXP vnt_lambda,
  SEXP vnt_priorCoef,
  SEXP vnt_r,
  SEXP vnt_tau1,
  SEXP vnt_tau2,
  SEXP vnt_priorTau1,
  SEXP vnt_atau1,
  SEXP vnt_btau1,
  SEXP vnt_priorModel,
  SEXP vnt_prModelpar);

extern SEXP modelSelectionCI(
  SEXP vnt_postSample,
  SEXP vnt_postOther,
  SEXP vnt_margpp,
  SEXP vnt_postMode,
  SEXP vnt_postModeProb,
  SEXP vnt_postProb,
  SEXP vnt_knownphi,
  SEXP vnt_priorCoef,
  SEXP vnt_niter,
  SEXP vnt_thinning,
  SEXP vnt_burnin,
  SEXP vnt_ndeltaini,
  SEXP vnt_deltaini,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_x,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_alpha,
  SEXP vnt_lambda,
  SEXP vnt_phi,
  SEXP vnt_tau,
  SEXP vnt_r,
  SEXP vnt_priorDelta,
  SEXP vnt_prDeltap,
  SEXP vnt_parprDeltap,
  SEXP vnt_verbose);

extern SEXP greedyVarSelCI(
  SEXP vnt_postMode,
  SEXP vnt_postModeProb,
  SEXP vnt_knownphi,
  SEXP vnt_priorCoef,
  SEXP vnt_niter,
  SEXP vnt_ndeltaini,
  SEXP vnt_deltaini,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_x,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_alpha,
  SEXP vnt_lambda,
  SEXP vnt_phi,
  SEXP vnt_tau,
  SEXP vnt_r,
  SEXP vnt_priorDelta,
  SEXP vnt_prDeltap,
  SEXP vnt_parprDeltap,
  SEXP vnt_verbose);

extern SEXP pmomMarginalKI(
  SEXP vnt_sel,
  SEXP vnt_nsel,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_phi,
  SEXP vnt_tau,
  SEXP vnt_r,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_logscale);

extern SEXP pmomMarginalUI(
  SEXP vnt_sel,
  SEXP vnt_nsel,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_x,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_tau,
  SEXP vnt_r,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_logscale,
  SEXP vnt_alpha,
  SEXP vnt_lambda);

extern SEXP pimomMarginalKI(
  SEXP vnt_sel,
  SEXP vnt_nsel,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_phi,
  SEXP vnt_tau,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_logscale);

extern SEXP pimomMarginalUI(
  SEXP vnt_sel,
  SEXP vnt_nsel,
  SEXP vnt_n,
  SEXP vnt_p,
  SEXP vnt_y,
  SEXP vnt_sumy2,
  SEXP vnt_x,
  SEXP vnt_XtX,
  SEXP vnt_ytX,
  SEXP vnt_tau,
  SEXP vnt_method,
  SEXP vnt_B,
  SEXP vnt_logscale,
  SEXP vnt_alpha,
  SEXP vnt_lambda);

#ifdef __cplusplus
}
#endif

#endif /* DO_MOMBF_H */

