/*
 * DO_MOMBF.H - Public include for .Call interface to model selection routines
 */

#ifndef DO_MOMBF_H
#define DO_MOMBF_H      1

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include <Rdefines.h>
#include "modelSel.h"
#include "mixtures.h"
#include "cstat.h"

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
/*extern SEXP bsplineCI(SEXP, SEXP, SEXP);
extern SEXP eprod_I(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP greedyVarSelCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mnormCI(SEXP, SEXP, SEXP);
extern SEXP modelSelectionEnumCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP modelSelectionGibbsCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlpMarginalAlaplI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP nlpMarginalSkewNormI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP normalmixGibbsCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pemomMarginalUI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pimomMarginalKI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pimomMarginalUI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pmomLM_I(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pmomMarginalKI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pmomMarginalUI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnlpCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnlpPostCI_lm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rnorm_truncMultCI(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rtmvnormCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rtmvnormProdCI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP zellnerMarginalKI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP zellnerMarginalUI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
*/

/*
static const R_CallMethodDef CallEntries[] = {
    {"bsplineCI",             (DL_FUNC) &bsplineCI,              3},
    {"eprod_I",               (DL_FUNC) &eprod_I,                5},
    {"greedyVarSelCI",        (DL_FUNC) &greedyVarSelCI,        28},
    {"mnormCI",               (DL_FUNC) &mnormCI,                3},
    {"modelSelectionEnumCI",  (DL_FUNC) &modelSelectionEnumCI,  27},
    {"modelSelectionGibbsCI", (DL_FUNC) &modelSelectionGibbsCI, 33},
    {"nlpMarginalAlaplI",     (DL_FUNC) &nlpMarginalAlaplI,     24},
    {"nlpMarginalSkewNormI",  (DL_FUNC) &nlpMarginalSkewNormI,  22},
    {"normalmixGibbsCI",      (DL_FUNC) &normalmixGibbsCI,      13},
    {"pemomMarginalUI",       (DL_FUNC) &pemomMarginalUI,       15},
    {"pimomMarginalKI",       (DL_FUNC) &pimomMarginalKI,       13},
    {"pimomMarginalUI",       (DL_FUNC) &pimomMarginalUI,       15},
    {"pmomLM_I",              (DL_FUNC) &pmomLM_I,              36},
    {"pmomMarginalKI",        (DL_FUNC) &pmomMarginalKI,        14},
    {"pmomMarginalUI",        (DL_FUNC) &pmomMarginalUI,        16},
    {"rnlpCI",                (DL_FUNC) &rnlpCI,                 9},
    {"rnlpPostCI_lm",         (DL_FUNC) &rnlpPostCI_lm,         11},
    {"rnorm_truncMultCI",     (DL_FUNC) &rnorm_truncMultCI,      5},
    {"rtmvnormCI",            (DL_FUNC) &rtmvnormCI,             7},
    {"rtmvnormProdCI",        (DL_FUNC) &rtmvnormProdCI,         9},
    {"zellnerMarginalKI",     (DL_FUNC) &zellnerMarginalKI,     13},
    {"zellnerMarginalUI",     (DL_FUNC) &zellnerMarginalUI,     15},
    {NULL, NULL, 0}
};

void R_init_mombf(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
*/

#endif /* DO_MOMBF_H */
