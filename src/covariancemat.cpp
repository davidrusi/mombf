#include <stdexcept>
#include "covariancemat.h"


/*Constructor

 INPUT

 - ncolx: number of columns in x (equivalently, number of rows & columns in XtX)

TODO: store block determinants. A vector of length ngroups?
*/
covariancemat::covariancemat(int ncolx) {

  this->ncolx= ncolx;

  (this->XtXs)= arma::sp_mat(ncolx, ncolx);
  (this->XtXcomputed)= arma::SpMat<short>(ncolx, ncolx);
}

//Class destructor
covariancemat::~covariancemat() { }



//Access element with matrix-type index, e.g. A(0,1) is element in row 0, column 1
double covariancemat::at(int i, int j) {

  if (XtXcomputed.at(i,j) == 1) {
    return XtXs.at(i,j);
  } else {
    throw std::runtime_error("covariancemat value not yet computed");
  }
}

// Check if element has already been computed with matrix-type index,
// e.g. A(0,1) is element in row 0, column 1
bool covariancemat::computed_at(int i, int j) {

  if (XtXcomputed.at(i,j) == 0) {  //if this entry has not been already computed
    return false;
  } else {
    return true;
  }
}

// Set element with matrix-type index and mark it as computed,
// e.g. A(0,1) is element in row 0, column 1
void covariancemat::set(int i, int j, double value) {

  XtXcomputed(i,j)= 1;
  XtXs(i,j)= value;
}

