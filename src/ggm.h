#ifndef GGM_H
#define GGM_H 1

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace std;

arma::sp_mat rowSumsC(NumericMatrix x, List prCoef, List prModel, List samplerPars);

arma::sp_mat modelSelectionGGMC(Rcpp::NumericMatrix y, List prCoef, List prModel, List samplerPars);



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


#endif /* MODELSEL_H */
