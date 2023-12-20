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

  ggmObject(arma::mat *y, List prCoef, List prModel, List samplerPars, bool computeS);
  ~ggmObject();

  //PUBLIC METHODS PROVIDED BY THE CLASS

  int n();    //sample size nrow(y)
  int ncol(); //number of variables ncol(y)

  CharacterVector sampler(); //sampler type
  int niter();
  int burnin();
  double pbirth();  //probability of birth move, ignored unless sampler is "birthdeath"
  int nbirth(); //number of birth/death updates to perform when updating each column of the precision matrix

  arma::mat S; //t(y) * y

  List prCoef;  //prior on parameters
  List prModel; //prior on model
  List samplerPars; //posterior sampler parameters

  bool verbose;

private:

  arma::mat *y;

};



//*************************************************************************************
// FUNCTIONS
//*************************************************************************************


List modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_Gibbs(arma::sp_mat *ans, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_Gibbs_parallel(std::list<arma::sp_mat> *ans, ggmObject *ggm, arma::sp_mat *Omegaini);

arma::mat get_invOmega_j(arma::sp_mat *Omega, int j);

void GGM_Gibbs_singlecol(arma::sp_mat *ans, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest);

void GGM_birthdeath_singlecol(arma::sp_mat *ans, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest);

void save_ggmsample_col(arma::sp_mat *ans, arma::SpMat<short> *model, double *sample_diag, arma::mat *sample_offdiag, int col2save, unsigned int colid);

void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model);

double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm);

//Matrix manipulation
void spmatsym_save2flat(arma::sp_mat *ans, arma::sp_mat *A, int col2store); //copy symmetric sp_mat in flat format to A(,col2store)

void spmat_rowcol2zero(arma::sp_mat *A, int colid); //Set row and colum colid of A to 0

void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int j); //drop row & column j from A

void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model); //copy A[model,model] into Aout, excluding column excludecol

void symmat2vec(arma::vec *Aflat, arma::mat *A); //flatten symmetric matrix A, in column-wise order


#endif /* MODELSEL_H */
