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

  arma::mat S; //t(y) * y

  List prCoef;  //prior on parameters
  List prModel; //prior on model
  List samplerPars; //posterior sampler parameters

private:

  arma::mat *y;

};


//*************************************************************************************
// typedefs
//*************************************************************************************

//pointer to function to compute log joint (log marginal likelihood + log model prior prob) for a given row
//GGMrow_marg is an example of such a function
typedef void(*pt2GGM_rowmarg)(double *, arma::mat *, arma::mat *, arma::SpMat<short> *, unsigned int, ggmObject *);  



//*************************************************************************************
// FUNCTIONS
//*************************************************************************************


//void print_mat( mat_type A );
//void print_mat(arma::mat *A);
//void print_spmat(arma::sp_mat *A);

arma::sp_mat modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_Gibbs(arma::sp_mat *ans, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_Gibbs_singlecol(arma::sp_mat *ans, int iter, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest);

void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm);

double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm);

//Matrix manipulation
void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int *j); //drop row & column j from A
void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model); //copy A[model,model] into Aout, excluding column excludecol



#endif /* MODELSEL_H */
