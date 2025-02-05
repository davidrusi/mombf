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

#include "modselIntegrals.h"
#include "crossprodmat.h"
#include "cstat.h"


using namespace Rcpp;
using namespace std;


//*************************************************************************************
// FUNCTIONS
//*************************************************************************************

//Functions implementing MCMC

List modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_MHwithinGibbs(arma::sp_mat *samples, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_global_proposal(std::vector<arma::SpMat<short>> *models, std::vector<std::vector<double>> *model_logprop, std::vector<std::map<string, double>> *map_logprob, double *logprop_modelini, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_MHwithinGibbs_global(arma::sp_mat *postSample, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, std::vector<std::map<string, double>> *map_logprob, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_MHwithinGibbs_onlyglobal(arma::sp_mat *postSample, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_CDA(arma::sp_mat *Omega, ggmObject *ggm);

void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob, double *modelini_logprob);

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::vec *margpp, arma::Col<int> *margppcount, int *number_accept, int *number_proposed, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, bool *use_LIT, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob, double *modelini_logprob);

void update_margpp_raoblack(arma::vec *margpp, double ppnew, arma::SpMat<short> *model, arma::SpMat<short> *modelnew);

void GGM_birthdeath_proposal(arma::SpMat<short> *modelnew, int *idx_update, bool *birth, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, bool setmodelnew);

void GGM_birthdeathswap_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, double *pdeath, double *pswap, bool setmodelnew);

void GGM_LIT_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, ggmObject *ggm, modselIntegrals_GGM *ms, bool setmodelnew);
void dprop_LIT_birth_GGM(std::vector<double> *proposal_kernel, std::vector<int> *indexes_birth, arma::SpMat<short> *model, ggmObject *ggm, modselIntegrals_GGM *ms);
void dprop_LIT_death_GGM(std::vector<double> *proposal_kernel, std::vector<int> *indexes_death, int *colid, arma::SpMat<short> *model, ggmObject *ggm, modselIntegrals_GGM *ms);

void niter_GGM_proposal(int *niter_prop, int *burnin_prop, int *niter, int *burnin, int *p);

void unique_model_logprob(arma::SpMat<short> *uniquemodels, std::vector<double> *unique_logprob, std::map<string, double> *map_logprob, arma::SpMat<short> *models, arma::mat *models_logprob, double *maxratio, double *logprobini);

std::string getModelid(arma::SpMat<short> *model, char *zerochar);

//Updating inverses
void update_Omega(arma::sp_mat *Omega, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag);
void addto_Omega(arma::mat *Omega, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag);


arma::mat get_invOmega_j(arma::sp_mat *Omega, int j);
void update_invOmega_submat(arma::mat *Omega_submat_inv, arma::sp_mat *Omega, int *oldcol, int *newcol);
void mapindexes_submat(int *mapforw, int *coldif, int *col1, int *col2, int *p);
void mapindexes_submat(int *mapforw, int *mapback, int *coldif, int *col1, int *col2, int *p);

//Saving output
void save_ggmsample_col(arma::sp_mat *ans, arma::SpMat<short> *model, double *sample_diag, arma::mat *sample_offdiag, int col2save, unsigned int colid);
void save_ggmmodel_col(arma::SpMat<short> *ans, arma::SpMat<short> *model, int col2save, unsigned int colid);
void save_postmean(arma::mat *postmean, arma::Mat<int> *postmeancount, arma::SpMat<short> *model, arma::mat *mean_offdiag, double *mean_diag, unsigned int colid);


//Obtaining marginal likelihoods
void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::mat *cholU, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv, arma::mat *cholU_old, arma::SpMat<short> *modelold);
void get_Omegainv_model(arma::mat *Omegainv_model, arma::mat *Omegainv, arma::SpMat<short> *model, unsigned int colid);

void GGMrow_marg_regression(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::mat *cholXtX, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model, arma::mat *cholXtX_old, arma::SpMat<short> *modelold);

double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm);

//Matrix manipulation
void spmatsym_save2flat(arma::sp_mat *ans, arma::sp_mat *A, int col2store); //copy symmetric sp_mat in flat format to A(,col2store)

bool checkNonZeroDiff(const arma::SpMat<short>* A, const arma::SpMat<short>* B, int maxdif);

void spmat_rowcol2zero(arma::sp_mat *A, int colid); //Set row and colum colid of A to 0

void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int j); //drop row & column j from A

void symmat2vec(arma::vec *Aflat, arma::mat *A); //flatten symmetric matrix A, in column-wise order


#endif /* MODELSEL_H */
