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
  ggmObject(ggmObject *ggm);
  ~ggmObject();

  //PUBLIC METHODS PROVIDED BY THE CLASS

  int n;    //sample size nrow(y)
  int ncol; //number of variables ncol(y)

  int burnin;  //number of MCMC burnin iterations
  bool fullscan; //if true, an MCMC iteration consists of updating all columns. If false, an MCMC iteration consists of updating a randomly chosen column
  double pbirth;  //probability of birth move, ignored unless sampler is "birthdeath"
  double pdeath;  //probability of death move, ignored unless sampler is "birthdeath"
  int updates_per_column; //number of birth/death updates to perform when updating each column of the precision matrix
  int niter; //number of MCMC iterations
  double tempering; //tempering parameter in parallel proposal
  double truncratio; //truncation ratio in parallel proposal. If prob(model) < prob(top model) / truncratio, then prob(model) = prob(top model) / truncratio

  arma::mat S; //t(y) * y

  double prCoef_lambda; //Prior on diagonal entries Omega_{jj} ~ Exp(lambda)
  double prCoef_tau; //Prior on off-diagonal Omega_{jk} | Omega_{jk} != 0 ~ N(0, tau)
  //List prCoef;  //prior on parameters

  std::string priorlabel; //Label for model space prior. Currently only "binomial" is possible, P(Omega_{jk} != 0) = priorPars_p
  double priorPars_p;
  //List prModel; //prior on model

  std::string sampler; //MCMC sampler type, e.g. Gibbs, birth-death
  //List samplerPars; //posterior sampler parameters

  double prob_parallel; //proposal probability of almost-parallel update
  bool parallel_regression; //use almost-parallel regression based proposal?
  bool parallel_insample;  //use almost-parallel in-sample based proposal?
  bool use_tempering;
  bool verbose;

};



//*************************************************************************************
// FUNCTIONS
//*************************************************************************************


List modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_Gibbs(arma::sp_mat *samples, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_parallel_proposal(std::vector<arma::SpMat<short>> *models, std::vector<std::vector<double>> *model_logprop, std::vector<std::map<string, double>> *map_logprob, double *logprop_modelini, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_parallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, std::vector<std::map<string, double>> *map_logprob, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_onlyparallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_CDA(arma::sp_mat *Omega, ggmObject *ggm);

void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob, double *modelini_logprob);

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int *number_accept, int *number_proposed, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob, double *modelini_logprob);

void GGM_birthdeath_proposal(arma::SpMat<short> *modelnew, int *idx_update, bool *birth, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, bool setmodelnew);

void GGM_birthdeathswap_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, double *pdeath, bool setmodelnew);

void niter_GGM_proposal(int *niter_prop, int *burnin_prop, int *niter, int *burnin, int *p);

void unique_model_logprob(arma::SpMat<short> *uniquemodels, std::vector<double> *unique_logprob, std::map<string, double> *map_logprob, arma::SpMat<short> *models, arma::mat *models_logprob, double *maxratio, double *logprobini);

std::string getModelid(arma::SpMat<short> *model, char *zerochar);

void update_Omega(arma::sp_mat *Omega, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag);

arma::mat get_invOmega_j(arma::sp_mat *Omega, int j);
void update_invOmega_submat(arma::mat *Omega_submat_inv, arma::sp_mat *Omega, int *oldcol, int *newcol);
void mapindexes_submat(int *mapforw, int *coldif, int *col1, int *col2, int *p);
void mapindexes_submat(int *mapforw, int *mapback, int *coldif, int *col1, int *col2, int *p);


void save_ggmsample_col(arma::sp_mat *ans, arma::SpMat<short> *model, double *sample_diag, arma::mat *sample_offdiag, int col2save, unsigned int colid);

void save_ggmmodel_col(arma::SpMat<short> *ans, arma::SpMat<short> *model, int col2save, unsigned int colid);

void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model);

void GGMrow_marg_regression(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model);

double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm);

//Matrix manipulation
void spmatsym_save2flat(arma::sp_mat *ans, arma::sp_mat *A, int col2store); //copy symmetric sp_mat in flat format to A(,col2store)

bool checkNonZeroDiff(const arma::SpMat<short>* A, const arma::SpMat<short>* B, int maxdif);

void spmat_rowcol2zero(arma::sp_mat *A, int colid); //Set row and colum colid of A to 0

void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int j); //drop row & column j from A

void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model); //copy A[model,model] into Aout

void symmat2vec(arma::vec *Aflat, arma::mat *A); //flatten symmetric matrix A, in column-wise order


#endif /* MODELSEL_H */
