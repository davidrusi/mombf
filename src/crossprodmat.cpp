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
    (this->XtXs)= arma::sp_mat(nrowx, ncolx);
    (this->XtXcomputed)= arma::SpMat<short>(nrowx, ncolx);
  }

}


//Class destructor
crossprodmat::~crossprodmat() { }



//Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
double crossprodmat::at(int i, int j) {

  if (dense) {

    return XtXd[i + j * ncolx];

  } else {

    if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed

      int iini, jini, k; double ans= 0;
      for (k=0, iini=i*nrowx, jini=j*nrowx; k< nrowx; k++) ans += x[k + iini] * x[k + jini];

      XtXcomputed(i,j)= 1;
      XtXs(i,j)= ans;
    }

    return XtXs.at(i,j);

  }
}


//Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
double crossprodmat::at(int k) {

  if (dense) {

    return XtXd[k];

  } else {

    int i= k % ncolx;
    int j= k / ncolx;

    if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed

      int iini, jini, k; double ans= 0;
      for (k=0, iini=i*nrowx, jini=j*nrowx; k< nrowx; k++) ans += x[k + iini] * x[k + jini];

      XtXcomputed(i,j)= 1;
      XtXs(i,j)= ans;
    }

    return XtXs.at(i,j);

  }


}



void crossprodmat::choldc(int idxini, int idxfi, double *cholx, double *detx, bool *posdef) {
  /*Cholesky decomposition of XtX[idxini..idxfi][idxini..idxfi], where input a is crossprodmat and output cholx a vector
    chol(a) is stored into a vector cholx in column order (1st column, 2nd column, etc).
    This means that the element (i,j) of chol(a) is stored into cholx[i + (j-1) (n-j/2)]

    Input: idxini, idxfi: first and last row/column indexes
    Ouput: cholx contains the Cholesky decomposition, detx the determinant of XtX[idxini..idxfi][idxini..idxfi]
  */
  int i,j,k, n=idxfi-idxini, kidx;
  double sum, *p, max_a;

  *posdef= true;
  *detx= 1.0;
  p= dvector(1,n);
  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      sum= this->at(idxini+i-1,idxini+j-1);
      for (k=i-1; k>=idxini; k--) { kidx= (k-1)*(n-k/2); sum -= cholx[i + kidx] * cholx[j + kidx]; }
      //for (sum=a[i][j],k=i-1; k>=idxini; k--) sum -= cholx[i][k]*cholx[j][k];
      if (i == j) {
	if (sum <= 0.0) *posdef= false;
	cholx[i + (i-1) * (n-i/2)]=sqrt(sum);  //cholx[i][i]=sqrt(sum);
	(*detx) *= sum;
      } else {
	max_a=max_xy(fabs(cholx[i + (i-1) * (n-i/2)]), 1e-10);  //max_a=max_xy(fabs(cholx[i][i]), 1e-10);
	cholx[j + (i-1)*(n-i/2)]= sum/max_a; //cholx[j][i]=sum/max_a;
      }
    }
  }
  free_dvector(p, 1,n);

}


