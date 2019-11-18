// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef COVARIANCEMAT
#define COVARIANCEMAT 1

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;


// THIS CLASS HANDLES STORAGE OF COVARIANCE MATRICES
// IT ALLOWS INCREMENTALLY POPULATING THE MATRIX AS ELEMENTS ARE NEEDED, SO ONE DOESN'T NEED TO PRE-COMPUTE AND ALLOCATE MEMORY FOR ALL ENTRIES.
// WHEN AN ENTRY IS REQUIRED IT IS COMPUTED ON-THE-FLY (OUTSIDE COVARIANCEMAT CLASS) AND STORED INTERNALLY INTO A SPARSE MATRIX USING CLASS SpMat FROM ARMADILLO
//
// NOTE: all indexes in this class start at 0, e.g. A(0,1) returns element in row 0, column 1; A(0) returns element in row 0 column 0
//

class covariancemat {

public:

  covariancemat(int ncolx);

  ~covariancemat();

  double at(int i, int j);  //Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
  bool computed_at(int i, int j);
  void set(int i, int j, double value);

private:

  int ncolx;
  arma::sp_mat XtXs;  //equivalent to SpMat<double> XtXs
  arma::SpMat<short> XtXcomputed; //bool entries indicating if XtX has been computed

};

#endif

