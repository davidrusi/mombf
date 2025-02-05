#ifndef MODSELINTEGRALS
#define MODSELINTEGRALS 1

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <map>
#include <string>
#include "modelSel_regression.h"
#include "cstat.h"
using namespace std;



//*************************************************************************************
// CLASS ggmObject
//*************************************************************************************

class ggmObject {

public:

  //Constructor and destructor

  ggmObject(arma::mat *y, List prCoef, List prModel, List samplerPars, bool use_tempering, bool computeS);
  ggmObject(ggmObject *ggm);
  ~ggmObject();

  //PUBLIC METHODS PROVIDED BY THE CLASS

  int n;    //sample size nrow(y)
  int ncol; //number of variables ncol(y)

  int burnin;  //number of MCMC burnin iterations
  double pbirth;  //probability of birth move, only used when sampler is "birthdeath" or "LIT"
  double pdeath;  //probability of death move, only used when sampler is "birthdeath"
  double pswap;   //probability of a swap move= 1 - prob birth - prob death
  double log_pbirth, log_pdeath, log_pswap; 
  double lbound_death; //In the LIT sampler, a lower bound on the log-proposal probability of a death move
  double ubound_death; //In the LIT sampler, an upper bound on the log-proposal probability of a death move
  double lbound_birth; //In the LIT sampler, a lower bound on the log-proposal probability of a birth move
  double ubound_birth; //In the LIT sampler, an upper bound on the log-proposal probability of a birth move

  int updates_per_iter; //an iteration consists of choosing updates_per_iter columns at random, and proposing updates_per_column updates for each column
  int updates_per_column; //see updates_per_iter
  int niter; //number of MCMC iterations
  double tempering; //tempering parameter in global proposal
  double truncratio; //truncation ratio in global proposal. If prob(model) < prob(top model) / truncratio, then prob(model) = prob(top model) / truncratio

  arma::mat S; //t(y) * y

  double prCoef_lambda; //Prior on diagonal entries Omega_{jj} ~ Exp(lambda)
  double prCoef_tau; //Prior on off-diagonal Omega_{jk} | Omega_{jk} != 0 ~ N(0, tau)
  //List prCoef;  //prior on parameters

  std::string priorlabel; //Label for model space prior. Currently only "binomial" is possible, P(Omega_{jk} != 0) = priorPars_p
  double priorPars_p;
  //List prModel; //prior on model

  std::string sampler; //MCMC sampler type, e.g. Gibbs, birth-death
  //List samplerPars; //posterior sampler parameters

  double prob_global; //proposal probability of almost-global update
  bool global_regression; //use almost-global regression based proposal?
  bool global_insample;  //use almost-global in-sample based proposal?
  bool use_tempering;
  bool verbose;

};



/***********************************************************************************/
/* Integrated likelihoods for regression models                                    */
/***********************************************************************************/

struct LM_logjoint {
  double logjoint; //saves log-joint posterior probability for a previously computed model (inclusion/exclusion of covariates in the regression model)
  arma::mat *mean; //saves posterior mean for a previously computed model
  arma::mat *cholVinv; //save Cholesky decomposition of the posterior precision matrix (e.g. given by X^T X + prior precision in the linear regression case)
};

class modselIntegrals {

public:

  modselIntegrals(pt2margFun marfun, pt2modelpriorFun priorfun, int nvars, lmObject *lm);  //initialize logjoint to fun, maxVars to nvars
  
  ~modselIntegrals();

  double getJoint(arma::SpMat<short> *model, lmObject *lm, arma::SpMat<short> *modelold); //Return logjoint(). Uses logjointSaved if available, else adds result to logjointSaved

  std::string getModelid(arma::SpMat<short> *model); //Return string with model id, e.g. "100010"

  double maxIntegral; //Stores value of largest integral

  string maxModel; //Stores model with largest integral, e.g. "10001" 

private:

  lmObject *lm; //lmObject storing data, prior parameters etc.
  int maxVars; //Maximum number of covariates
  char *zerochar;  //Store model id (vars in the model) in character format, e.g. "00000"
  pt2margFun marginalFunction;  //Function computing log(marginal likelihood)
  pt2modelpriorFun priorFunction;     //Function computing log(model prior)
  std::map<string, LM_logjoint> logjointSaved; //saves log-joint, posterior mean and Cholesky decomp for each previously considered model
  //std::map<string, double> logjointSaved; //Saves previously computed logjoint
  long unsigned int maxsave;

};


//*************************************************************************************
// typedefs
//*************************************************************************************

//pointer to function to compute log joint (log marginal likelihood + log model prior prob) for a given row
//GGMrow_marg is an example of such a function
typedef void(*pt2GGM_rowmarg)(double *, arma::mat *, arma::mat *, arma::mat *, arma::SpMat<short> *, unsigned int, ggmObject *, arma::mat *, arma::mat *, arma::SpMat<short> *);  


/************************************************************************************/
/* Integrated likelihoods for Gaussian graphical models with precision matrix Omega */
/*                                                                                  */
/* Models are defined by non-zero entries in column cold, and are conditional on    */
/* a given value of Omegainv, the inverse of Omega[-colid,-colid]                   */
/************************************************************************************/


struct GGM_logjoint {
  double logjoint; //saves log-joint posterior probability for a previously computed model (presence/absence of edges in the GGM)
  arma::mat *mean; //saves posterior mean for a previously computed model
  arma::mat *cholV; //save Cholesky decomp of the posterior covariance for a previously computed model
  arma::mat *cholVinv; //save Cholesky decomposition of the inverse posterior covariance (given by X^T X in the regression proposal)
};


class modselIntegrals_GGM {

public:

  modselIntegrals_GGM(pt2GGM_rowmarg jointFunction, ggmObject *ggm, unsigned int colid, arma::mat *Omegainv); 
  
  ~modselIntegrals_GGM();

  void getJoint(double *logjoint, arma::mat *mean_offdiag, double *mean_diag, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, arma::SpMat<short> *modelold, bool postSample); //Return logjoint() and posterior sample for off-diag & diagonal elements 

  void getMode(double *logjoint, arma::mat *mode_offdiag, double *mode_diag, arma::SpMat<short> *model, arma::SpMat<short> *modelold); //Return logjoint() and posterior mode for off-diag & diagonal elements

  std::string getModelid(arma::SpMat<short> *model); //Convert model to a string, e.g. "10001"

  double maxIntegral; //Stores value of largest integral

  string maxModel; //Stores model with largest integral, e.g. "10001" 

  int nvars; //number of variables (ncol(Omega) - 1)

private:

  pt2GGM_rowmarg jointFunction; //Function computing log(marginal likelihood) + log(model prior)
  ggmObject *ggm;  //Object storing info about the Gaussian graphical model and the algorithms to be used
  unsigned int colid; //column of Omega for which log-joints are being calculated
  arma::mat *Omegainv; //inverse of Omega[-colid,-colid]

  char *zerochar;  //Store model id (vars in the model) in character format, e.g. "00000"

  std::map<string, GGM_logjoint> logjointSaved; //saves log-joint, posterior mean and Cholesky decomp for each previously considered model

  long unsigned int maxsave; //if size of logjointSaved, meanSaved, cholVsaved  >= maxsave, save only models with non-negligible post prob vs maxModel

  void get_Omegainv_model(arma::mat *Omegainv_model, arma::SpMat<short> *model);

};


#endif

