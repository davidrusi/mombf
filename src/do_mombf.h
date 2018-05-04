/*
 * DO_MOMBF.H - Public include for .Call interface to model selection routines
 */

#ifndef DO_MOMBF_H
#define DO_MOMBF_H      1

#include <Rdefines.h>
#include "R_ext/Rdynload.h"




static R_CallMethodDef callMethods[]  = {
  {"bsplineCI", (DL_FUNC) &bsplineCI, 3},
  {"mnormCI", (DL_FUNC) &mnormCI, 3},
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
  {"normalmixGibbsCI", (DL_FUNC) &normalmixGibbsCI, 12},
  {NULL, NULL, 0}
};

void R_init_mombf(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
   R_useDynamicSymbols(info, FALSE);
}


#endif /* DO_MOMBF_H */
