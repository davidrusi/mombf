#ifndef GGM_H
#define GGM_H 1

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

  ggmObject(NumericMatrix y, List prCoef, List prModel, List samplerPars);
  ~ggmObject();

  //PUBLIC METHODS PROVIDED BY THE CLASS

  int nrow(); //number of variables
  int ncol(); //number of individuals

  CharacterVector sampler(); //sampler type
  int niter();
  int burnin();

private:

  NumericMatrix y;

  List prCoef;  //prior on parameters
  List prModel; //prior on model
  List samplerPars; //posterior sampler parameters

};



//*************************************************************************************
// FUNCTIONS
//*************************************************************************************

arma::sp_mat modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini);

void GGM_Gibbs(arma::sp_mat *ans, ggmObject *ggm, arma::sp_mat *Omegaini);

void GGM_Gibbs_singlerow(arma::sp_mat *ans, int iter, int rowid, ggmObject *ggm, arma::sp_mat *Omegaini);



#endif /* MODELSEL_H */
