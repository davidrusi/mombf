// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef CROSSPRODMAT
#define CROSSPRODMAT 1

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;


// THIS CLASS HANDLES CROSS-PRODUCTS MATRICES XtX (INNER-PRODUCTS BETWEEN COLUMNS of x)
// IT ALLOWS INCREMENTALLY POPULATING THE MATRIX AS ELEMENTS ARE NEEDED, SO ONE DOESN'T NEED TO PRE-COMPUTE AND ALLOCATE MEMORY FOR ALL ENTRIES.
// WHEN AN ENTRY IS REQUIRED IT IS COMPUTED ON-THE-FLY AND STORED INTERNALLY INTO A SPARSE MATRIX USING CLASS SpMat FROM ARMADILLO
//
// Let x be an (nrowx,ncolx) matrix stored as a vector (in column-order, i.e. column 1, column 2, etc)
// The cross-product matrix XtX has (i,j) element equal to the inner product between columns i and j in x
//
// NOTE: all indexes in this class start at 0, e.g. A(0,1) returns element in row 0, column 1; A(0) returns element in row 0 column 0
//
/* EXAMPLE: matrix x has 4 rows, 3 columns
      bool posdef;
      int i; double *x, *cholXtX, detXtX;
      x= dvector(0,12);
      for (i=0; i<12; i++) x[i]= (double) i;

      // OPTION 1. PRE-COMPUTE XtX
      double *XtX;
      XtX= dvector(0,9);
      XtX[0]=  14; XtX[3]=  38; XtX[6]=  62;
      XtX[1]=  38; XtX[4]= 126; XtX[7]= 214;
      XtX[2]=  62; XtX[5]= 214; XtX[8]= 366;

      crossprodmat *Ad;
      Ad= new crossprodmat(XtX,4,3,true);  //true indicates that 1st argument is the pre-computed XtX
      double tmp= Ad->at(0,0);  //returns the pre-computed XtX(0,0)
      cholXtX= dvector(0,2);    //store Cholesky decomp of XtX[1:2,1:2] requires 3 elements
      XtX->choldc(1,2,*cholXtX,&detXtX,&posdef);
      delete Ad; free_dvector(x,0,9); free_dvector(XtX,0,9); free_dvector(cholXtX, 0,2);

      // OPTION 2. COMPUTE XtX ON-THE-FLY
      crossprodmat *As;
      As= new crossprodmat(x,4,3,false); //false indicates that 1st argument is the matrix x
      double tmp= As->at(0,0); //computes and stores internally the inner-product, then returns it
      tmp= As->at(0,0); //returns the stored value
      delete As; free_dvector(x,0,9); free_dvector(XtX,0,9);
*/



//CLASS crossprodmat using goodies from Rcpp
class crossprodmatRcpp {

public:

  crossprodmatRcpp(NumericMatrix mymat, bool dense); //if dense==true, mymat is pointer to pre-computed XtX; if dense==false, mymat is pointer to x

  ~crossprodmatRcpp();

  double at(int i, int j);  //Access XtX(i,j). Careful: row & column indexes start at 0

private:

  NumericMatrix x;
  bool dense; //if true then matrix is stored in XtXd, else in XtXs
  arma::mat XtXd;
  arma::sp_mat XtXs;  //equivalent to SpMat<double> XtXs
  arma::SpMat<short> XtXcomputed; //bool entries indicating if XtX has been computed

};


//CLASS crossprodmat using basic C++
class crossprodmat {

public:

  crossprodmat(double *mymat, int nrowx, int ncolx, bool dense, int nuserows, int *userows); //if dense==true, mymat is pointer to pre-computed XtX; if dense==false, mymat is pointer to x
  crossprodmat(double *mymat, int nrowx, int ncolx, bool dense, int nuserows, int userowsini);
  crossprodmat(double *mymat, int nrowx, int ncolx, bool dense);

  ~crossprodmat();

  double at(int i, int j);  //Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
  double at(int k);  //Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
  arma::mat * submat(std::vector<int> *row_indexes, std::vector<int> *col_indexes); //Return dense submatrix with rows / cols in row_indexes / col_indexes

  void choldc(int idxini, int idxfi, double *cholXtX, double *detXtX, bool *posdef); //Cholesky decomp and determinant

private:

  double *x;           //X stored as a vector
  int nrowx;
  int ncolx;
  int *userows; //optional slot indicating the indexes of the rows in x to be used when computing XtX. That is XtX= t(x[userows,]) %*% x[userows,]
  int nuserows; //number of rows in x to be used when computing XtX
  int userowsini; //if userows not provided, use x[userowsini : userowsini+nuserows-1,] to compute XtX (default userowsini=0)
  bool dense; //if true then matrix is stored in XtXd, else in XtXs
  double *XtXd;
  arma::sp_mat XtXs;  //equivalent to SpMat<double> XtXs
  arma::SpMat<short> XtXcomputed; //bool entries indicating if XtX has been computed

};

#endif

