#include "crossprodmat.h"
#include "cstat.h"
using namespace std;
using namespace arma;


//Constructor
// If dense=true, mymat is assumed to point to a pre-computed dense matrix containing XtXd
// If dense=false, mymat is assumed to be the x matrix. No XtX entries are pre-computed at the time of creation
crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense) {

  this->nrowx= nrowx;
  this->ncolx= ncolx;

  if (dense) {
    this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    arma::SpMat<short> (this->XtXcomputed)= arma::SpMat<short>(nrowx, ncolx);
  }

}


//Class destructor
crossprodmat::~crossprodmat() { }



//Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
double crossprodmat::at(const int i, const int j) {
//double crossprodmat::operator()(const int i, const int j) {
  
  if (dense) {

    return XtXd[i + j * nrowx];

  } else {

    if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed

      int k; double ans= 0;
      int iini= i * nrowx;
      int jini= j * nrowx;

      for (k=0; k< nrowx; k++) ans += x[k + iini] * x[k + jini];

      XtXcomputed.at(i,j)= 1;
      XtXs.at(i,j)= ans;
    }

    return XtXs.at(i,j);

  }
}


//Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
double crossprodmat::at(const int k) {
//double crossprodmat::operator() (const int k) {

  return (*this).at(k % nrowx, k/nrowx);
  //return (*this).operator()(k % nrowx, k/nrowx);

}



