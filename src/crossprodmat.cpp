#include "crossprodmat.h"
#include "cstat.h"
using namespace std;
using namespace arma;




/***************************************************************************************************/
/* CLASS crossprodmat Rcpp goodies                                                                 */
/***************************************************************************************************/

//if dense==true, mymat is pointer to pre-computed XtX; if dense==false, mymat is pointer to x
crossprodmatRcpp::crossprodmatRcpp(NumericMatrix mymat, bool dense) {

  if (dense) {
    this->XtXd = arma::mat(mymat.begin(), mymat.nrow(), mymat.ncol(), false);
    //this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    (this->XtXs)= arma::sp_mat(mymat.ncol(), mymat.ncol());
    (this->XtXcomputed)= arma::SpMat<short>(mymat.ncol(), mymat.ncol());
  }

}

//Class destructor
crossprodmatRcpp::~crossprodmatRcpp() { }


//Access XtX(i,j)
double crossprodmatRcpp::at(int i, int j) { 

//adapt code below
  if (dense) {

    return XtXd.at(i,j);

  } else {

    if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed

      NumericVector coli= this->x(_,i), colj= this->x(_,j);
      XtXs(i,j)= sum(coli * colj);
      //XtXs(i,j)= arma::as_scalar(coli.t() * colj);
      XtXcomputed(i,j)= 1; 

    }

    return XtXs.at(i,j);

  }

}


/***************************************************************************************************/
/* CLASS crossprodmat USING BASIC C++                                                              */
/***************************************************************************************************/

/*Constructor

 INPUT

 - mymat: If dense==true, mymat is assumed to point to a pre-computed dense matrix containing XtXd (an ncolx * ncolx matrix)
          If dense==false, mymat is assumed to be the x matrix (an nrowx * ncolx matrix). No XtX entries are pre-computed at the time of creation

 - nrowx: number of rows in x
 - ncolx: number of columns in x (equivalently, number of rows & columns in XtX)
 - dense: indicator if a pre-computed XtX was provided

*/
crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense, int nuserows, int *userows) {

  this->nrowx= nrowx;
  this->ncolx= ncolx;
  this->userowsini= userows[0];
  this->nuserows= nuserows;
  this->userows= userows;

  if (dense) {
    this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    (this->XtXs)= arma::sp_mat(ncolx, ncolx);
    (this->XtXcomputed)= arma::SpMat<short>(ncolx, ncolx);
  }

}

crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense, int nuserows, int userowsini) {

  this->nrowx= nrowx;
  this->ncolx= ncolx;
  this->userowsini= userowsini;
  this->nuserows= nuserows;
  this->userows= NULL;

  if (dense) {
    this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    (this->XtXs)= arma::sp_mat(ncolx, ncolx);
    (this->XtXcomputed)= arma::SpMat<short>(ncolx, ncolx);
  }

}

crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense) {

  this->nrowx= nrowx;
  this->ncolx= ncolx;
  this->userowsini= 0;
  this->nuserows= nrowx;
  this->userows= NULL;

  if (dense) {
    this->XtXd= mymat;
    this->dense= true;
  } else {
    this->x= mymat;
    this->dense= false;
    (this->XtXs)= arma::sp_mat(ncolx, ncolx);
    (this->XtXcomputed)= arma::SpMat<short>(ncolx, ncolx);
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
      if (this->userows ==NULL) {
        for (k= this->userowsini, iini=i*nrowx, jini=j*nrowx; k< this->nuserows+this->userowsini; k++) ans += x[k + iini] * x[k + jini];
      } else {
        for (k= 0, iini=i*nrowx, jini=j*nrowx; k< this->nuserows; k++) ans += x[(this->userows[k]) + iini] * x[(this->userows[k]) + jini];
      }

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

    int i= k % ncolx, j= k / ncolx;

    if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed

      int iini, jini, k; double ans= 0;
      if (this->userows ==NULL) {

        for (k= this->userowsini, iini=i*nrowx, jini=j*nrowx; k< this->nuserows+this->userowsini; k++) ans += x[k + iini] * x[k + jini];

      } else {

        for (k= 0, iini=i*nrowx, jini=j*nrowx; k< this->nuserows; k++) ans += x[(this->userows[k]) + iini] * x[(this->userows[k]) + jini];

      }

      XtXcomputed(i,j)= 1;
      XtXs(i,j)= ans;
    }

    return XtXs.at(i,j);

  }


}

//Return dense submatrix with rows / cols in row_indexes / col_indexes
arma::mat * crossprodmat::submat(std::vector<int> *row_indexes, std::vector<int> *col_indexes) {
  int i, j, rowi, nrows= row_indexes->size(), ncols= col_indexes->size();
  arma::mat *ans;
  ans= new arma::mat(nrows, ncols);

  for (i=0; i<nrows; i++) {
    rowi= row_indexes->at(i);
    ans->at(i,i)= this->at(rowi,rowi);
    for (j=0; j<i; j++) ans->at(i,j)= ans->at(j,i)= this->at(rowi,col_indexes->at(j));
  }

  return ans;
}



void crossprodmat::choldc(int idxini, int idxfi, double *cholXtX, double *detXtX, bool *posdef) {
  /*Cholesky decomposition of XtX[idxini..idxfi][idxini..idxfi], where input a is crossprodmat and output cholXtX a vector
    chol(a) is stored into a vector cholXtX in column order (1st column, 2nd column, etc).

    - Element (i,i) of chol(XtX) is stored into cholXtX[ii], where ii= (i-1)*n - (i-1)*(i-2)/2;
    - For i>j, element (i,j) of chol(XtX) is stored into cholXtX[jj + i - j], where jj= (j-1)*n - (j-1)*(j-2)/2;

    Input: idxini, idxfi: first and last row/column indexes
    Ouput: cholXtX[0,...,idxfi-idxini] contains the Cholesky decomposition, detXtX the determinant of XtX[idxini..idxfi][idxini..idxfi]
  */
  int i,ii,j,k,kk, n=idxfi-idxini+1;
  double sum, *p, max_a;

  *posdef= true;
  *detXtX= 1.0;
  p= dvector(1,n);
  for (i=1;i<=n;i++) {
    ii= (i-1)*n - (i-1)*(i-2)/2;
    for (j=i;j<=n;j++) {
      sum= this->at(idxini+i-1,idxini+j-1);
      for (k=i-1; k>=1; k--) { kk= (k-1)*n - (k-1)*(k-2)/2; sum -= cholXtX[kk + i-k] * cholXtX[kk + j-k]; } //sum -= cholXtX[i][k]*cholXtX[j][k];
      if (i == j) {
        if (sum <= 0.0) *posdef= false;
        cholXtX[ii]=sqrt(sum);  //cholXtX[i][i]=sqrt(sum);
        (*detXtX) *= sum;
      } else {
        max_a=max_xy(fabs(cholXtX[ii]), 1e-10);  //max_a=max_xy(fabs(cholXtX[i][i]), 1e-10);
        cholXtX[ii + j-i]= sum/max_a; //cholXtX[j][i]=sum/max_a;
      }
    }
  }
  free_dvector(p, 1,n);

}


