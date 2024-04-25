#ifndef GGM_H
#define GGM_H 1

// Important: this definition ensures Armadillo enables SuperLU
// Commented out as the configuration of SuperLU is system-dependent
//#define ARMA_USE_SUPERLU 1

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


#include "crossprodmat.h"
#include "cstat.h"


using namespace Rcpp;
using namespace std;



//*************************************************************************************
// CLASS ggmObject
//*************************************************************************************

class ggmObject {

public:

  //Constructor and destructor

  ggmObject(arma::mat *y, List prCoef, List prModel, List samplerPars, bool use_tempering, bool computeS);
  ~ggmObject();

  //PUBLIC METHODS PROVIDED BY THE CLASS

  int n;    //sample size nrow(y)
  int ncol; //number of variables ncol(y)

  int burnin;  //number of MCMC burnin iterations
  double pbirth;  //probability of birth move, ignored unless sampler is "birthdeath"
  int nbirth; //number of birth/death updates to perform when updating each column of the precision matrix
  int niter; //number of MCMC iterations
  double tempering; //tempering parameter in almost-parallel proposal

  arma::mat S; //t(y) * y

  double prCoef_lambda; //Prior on diagonal entries Omega_{jj} ~ Exp(lambda)
  double prCoef_tau; //Prior on off-diagonal Omega_{jk} | Omega_{jk} != 0 ~ N(0, tau)
  //List prCoef;  //prior on parameters

  std::string priorlabel; //Label for model space prior. Currently only "binomial" is possible, P(Omega_{jk} != 0) = priorPars_p
  double priorPars_p;
  //List prModel; //prior on model

  std::string sampler; //MCMC sampler type, e.g. Gibbs, birth-death
  //List samplerPars; //posterior sampler parameters

  bool use_tempering;
  bool verbose;

//private:

//  arma::mat *y;

};



//*************************************************************************************
// FUNCTIONS
//*************************************************************************************


List modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_Gibbs(arma::sp_mat *samples, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini);

//void GGM_parallel_propdensity(arma::mat *propdens, double *dpropini, std::vector<arma::sp_mat> *samples, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_Gibbs_parallel(std::vector<arma::SpMat<short>> *models, ggmObject *ggm, arma::sp_mat *Omegaini, arma::mat *model_logprop);

void GGM_parallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_samples, arma::mat *propdens, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini);

void update_Omegaini(arma::sp_mat *Omegaini, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag);

arma::mat get_invOmega_j(arma::sp_mat *Omega, int j);

void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob);

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob);

void save_ggmsample_col(arma::sp_mat *ans, arma::SpMat<short> *model, double *sample_diag, arma::mat *sample_offdiag, int col2save, unsigned int colid);

void save_ggmmodel_col(arma::SpMat<short> *ans, arma::SpMat<short> *model, int col2save, unsigned int colid);

void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model);

double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm);

//Matrix manipulation
void spmatsym_save2flat(arma::sp_mat *ans, arma::sp_mat *A, int col2store); //copy symmetric sp_mat in flat format to A(,col2store)

void spmat_rowcol2zero(arma::sp_mat *A, int colid); //Set row and colum colid of A to 0

void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int j); //drop row & column j from A

void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model); //copy A[model,model] into Aout, excluding column excludecol

void symmat2vec(arma::vec *Aflat, arma::mat *A); //flatten symmetric matrix A, in column-wise order


#endif /* MODELSEL_H */
