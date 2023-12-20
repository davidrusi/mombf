#ifndef MODSELINTEGRALS
#define MODSELINTEGRALS 1

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <map>
#include <string>
#include "modelSel_regression.h"
#include "ggm.h"
#include "cstat.h"
using namespace std;


/***********************************************************************************/
/* Integrated likelihoods for regression models                                    */
/***********************************************************************************/

class modselIntegrals {

public:

  modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars);  //initialize logjoint to fun, maxVars to nvars
  
  ~modselIntegrals();

  double getJoint(int *sel, int *nsel, struct marginalPars *pars); //Return logjoint(). Uses logjointSaved if available, else adds result to logjointSaved

  double maxIntegral; //Stores value of largest integral

  string maxModel; //Stores model with largest integral, e.g. "10001" 

private:

  int maxVars; //Maximum number of covariates
  char *zerochar;  //Store model id (vars in the model) in character format, e.g. "00000"
  pt2margFun marginalFunction;  //Function computing log(marginal likelihood)
  pt2margFun priorFunction;     //Function computing log(model prior)
  std::map<string, double> logjointSaved; //Saves previously computed logjoint
  long unsigned int maxsave;

};


//*************************************************************************************
// typedefs
//*************************************************************************************

//pointer to function to compute log joint (log marginal likelihood + log model prior prob) for a given row
//GGMrow_marg is an example of such a function
typedef void(*pt2GGM_rowmarg)(double *, arma::mat *, arma::mat *, arma::SpMat<short> *, unsigned int, ggmObject *, arma::mat *);  


/************************************************************************************/
/* Integrated likelihoods for Gaussian graphical models with precision matrix Omega */
/*                                                                                  */
/* Models are defined by non-zero entries in column cold, and are conditional on    */
/* a given value of Omegainv, the inverse of Omega[-colid,-colid]                   */
/************************************************************************************/


class modselIntegrals_GGM {

public:

  modselIntegrals_GGM(pt2GGM_rowmarg jointFunction, ggmObject *ggm, unsigned int colid, arma::mat *Omegainv); 
  
  ~modselIntegrals_GGM();

  void getJoint(double *logjoint, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, bool postSample); //Return logjoint() and posterior sample for off-diagonal and diagonal elements 

  double maxIntegral; //Stores value of largest integral

  string maxModel; //Stores model with largest integral, e.g. "10001" 

  int nvars; //number of variables (ncol(Omega) - 1)

private:

  pt2GGM_rowmarg jointFunction; //Function computing log(marginal likelihood) + log(model prior)
  ggmObject *ggm;  //Object storing info about the Gaussian graphical model and the algorithms to be used
  unsigned int colid; //column of Omega for which log-joints are being calculated
  arma::mat *Omegainv; //inverse of Omega[-colid,-colid]

  //int maxVars; //Maximum number of covariates
  char *zerochar;  //Store model id (vars in the model) in character format, e.g. "00000"
  std::map<string, double> logjointSaved; //saves log-joint for each previously computed model
  std::map<string, arma::mat *> meanSaved; //saves posterior mean for each previously computed model
  std::map<string, arma::mat *> cholVSaved; //save Cholesky decomp of the posterior covariance for each previously computed model

  long unsigned int maxsave; //if size of logjointSaved, meanSaved, cholVsaved  >= maxsave, save only models with non-negligible post prob vs maxModel

  void get_Omegainv_model(arma::mat *Omegainv_model, arma::SpMat<short> *model);

};


#endif

