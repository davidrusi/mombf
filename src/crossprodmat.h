// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef CROSSPRODMAT
#define CROSSPRODMAT 1

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rinternals.h>
#include <map>
#include <string>
#include "cstat.h"
using namespace std;


// THE PURPOSE OF THIS CLASS IS TO STORE CROSS-PRODUCT MATRICES
// MAIN FEATURE: AVOID RESERVING MEMORY FOR & PRE-COMPUTING ALL ENTRIES.
// WHEN AN ENTRY IS REQUIRED IT IS COMPUTED ON-THE-FLY AND STORED INTERNALLY INTO A SPARSE MATRIX
//
// SECONDARY FEATURE: ACCESS DENSE MATRICES STORED AS VECTORS USING MATRIX ACCESSORS A[i,j]

// Let x be an (nrowx,ncolx) matrix, stored as a vector (column1,column2,...)
// The cross-product matrix XtX has (i,j) element equal to the inner product between columns i and j in x
//
// Example 1. Create object from x
// A= crossprodmat(x, nrow, ncol);  //no inner-products are computed
// double a= A(0,1);  //inner product between columns (0,1) is computed, stored internally and saved into a
// double b= A(0,1) + 1;  //inner product is not re-computed, as it was stored earlier
//
// Example 2. Create object from XtX
// A= crossprodmat(XtX, nrow, ncol);  //store pointer to pre-computed inner-products in XtX
// double a= A(0,1);  //access element (0,1) in XtX
// double b= A(0,1) + 1;  //inner product is not re-computed, as it was stored earlier

//NOTE: all indexes in this class start at 0, e.g. A(0,1) returns element in row 0, column 1; A(0) returns element in row 0 column 0

class crossprodmat {

public:

  crossprodmat(double *mymat, int *nrowx, int *ncolx, bool dense);
  crossprodmat(double *mymat, int nrowx, int ncolx, bool dense);

  ~crossprodmat();

  double operator() (const int i, const int j);  //Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1

  double operator() (const int k);  //Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow

private:

  double *x;
  int *nrowx;
  int *ncolx;
  bool dense; //if true then matrix is stored in XtXd, else in XtXs
  double *XtXd;
  arma::sp_mat XtXs;  //equivalent to SpMat<double> XtXs
  arma::SpMat<short> XtXcomputed; //bool entries indicating if XtX has been computed

};

#endif

