#include "crossprodmat.h"
#include "cstat.h"
using namespace std;
using namespace arma;


//Constructor
// If dense=true, mymat is assumed to point to a pre-computed dense matrix containing XtXd
// If dense=false, mymat is assumed to be the x matrix. No XtX entries are pre-computed at the time of creation
crossprodmat::crossprodmat(double *mymat, int *nrowx, int *ncolx, bool dense) {

  this->nrowx= nrowx;
  this->ncolx= ncolx;

  if (dense) {
    this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    arma::SpMat<short> (this->XtXcomputed)= arma::SpMat<short>(*nrowx, *ncolx);
  }

}

crossprodmat::crossprodmat(double *XtX, int nrowx, int ncolx, bool dense) { crossprodmat(XtX, &nrowx, &ncolx, dense); }


//Class destructor
crossprodmat::~crossprodmat() { }



//Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
double crossprodmat::operator() (const int i, const int j) {

  if (dense) {

    return XtXd[i + (j-1) * (*nrowx)];

  } else {

    if (XtXcomputed.at(i,j) == 1) {  //if this entry has been already computed

      return XtXs.at(i,j);

    } else {  //else compute it and store into XtXs

      int k; double ans= 0;
      int iini= (i-1) * (*nrowx);
      int jini= (j-1) * (*nrowx);

      for (k=0; k< *nrowx; k++) ans += x[k + iini] * x[k + jini];

      XtXcomputed.at(i,j)= 1;
      XtXs.at(i,j)= ans;
      return ans;
    }

  }
}


//Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
double crossprodmat::operator() (const int k) {

  return (*this).operator()(k % (*nrowx), k/(*nrowx));

}



