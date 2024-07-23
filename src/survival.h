#ifndef SURVIVAL_H
#define SURVIVAL_H 1


#include <RcppArmadillo.h>
#include <list>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "crossprodmat.h"
#include "covariancemat.h"
#include "cstat.h"
#include "modelSel_regression.h"
#include "modselFunction.h"
using namespace std;



//*************************************************************************************
// MARGINAL LIKELIHOOD FOR ACCELERATED FAILURE TIME MODELS
//*************************************************************************************

//log-likelihood of Normal AFT model and its derivatives
void negloglnormalAFT(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs);
void negloglnormalAFTupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void negloglnormalAFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void loglnormalAFThess(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);

//Approx log-likelihood of Normal AFT model and its derivatives (based on apnorm, ainvmillsnorm)
void anegloglnormalAFT(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs);
void anegloglnormalAFT0(double *f, double *th, int *sel, int *thlength, lmObject *lm,  std::map<string, double *> *funargs);
void anegloglnormalAFTupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void anegloglnormalAFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void anegloglnormalAFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void aloglnormalAFThess(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);



// Computation of marginal likelihoods
double SurvMargALA(int *sel, int *nsel, lmObject *lm, int priorcode);  //same as SurvMarg, using ALA
double SurvMarg(int *sel, int *nsel, lmObject *lm, int priorcode);  //wrapper function calling the function corresponding to the specified prior

double pmomgmomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double pmomgmomSurvMarg(int *sel, int *nsel, lmObject *lm); // pMOM on individual coef, group MOM on groups
double pemomgemomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv); // peMOM on individual coef, group eMOM on groups

double gmomgmomSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double gmomgmomSurvMarg(int *sel, int *nsel, lmObject *lm);
double gmomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double gmomgzellSurvMarg(int *sel, int *nsel, lmObject *lm);
double pmomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double pmomgzellSurvMarg(int *sel, int *nsel, lmObject *lm); // pMOM on individual coef, block Zellner on groups
double pemomgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double pemomgzellSurvMarg(int *sel, int *nsel, lmObject *lm); // peMOM on individual coef, block Zellner on groups
double gzellgzellSurvMarg(arma::SpMat<short> *sel, lmObject *lm, arma::mat *cholV_old, arma::SpMat<short> *modelold, arma::mat *m, arma::mat *cholVinv);
double gzellgzellSurvMarg (int *sel, int *nsel, lmObject *lm); // Zellner on individual coef, block Zellner on groups


//Evaluate log-posterior
void fpmomgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void fpemomgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void fgzellgzellSurv(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void fgzellgzellSurv0(double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);

void fpmomgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void fpemomgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);
void fgzellgzellSurvupdate(double *fnew, double *thjnew, int j, double *f, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double *> *funargs);


//Evaluate log-posterior gradient & hessian
void fpmomgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fpmomgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fpemomgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fpemomgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fgzellgzell_AFTgradhess(double *grad, double *hess, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fgzellgzell_AFTgrad(double *grad, int j, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);

void fpmomgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fpemomgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);
void fgzellgzellhess_AFT(double **hess, double *th, int *sel, int *thlength, lmObject *lm, std::map<string, double*> *funargs);




#endif /* SURVIVAL_H */
