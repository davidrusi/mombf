#include "ggm.h"


//*************************************************************************************
// CLASS ggmObject
//*************************************************************************************

//Constructor
ggmObject::ggmObject(NumericMatrix y, List prCoef, List prModel, List samplerPars) {

  this->y = y;
  this->prCoef = prCoef;
  this->prModel = prModel;
  this->samplerPars = samplerPars;

}

//Destructor
ggmObject::~ggmObject() {

}

int ggmObject::nrow() {
  return (this->y).nrow();
}

int ggmObject::ncol() {
  return (this->y).ncol();
}


//Posterior sampler type
CharacterVector ggmObject::sampler() { 
  return this->samplerPars["sampler"]; 
}

int ggmObject::niter() {
  return this->samplerPars["niter"];
}

int ggmObject::burnin() {
  return this->samplerPars["burnin"];
}


// Example of returning an sp_mat object
//arma::sp_mat armaEx(S4 mat, bool show) {
//  IntegerVector dims = mat.slot("Dim");
//   arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
//   arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
//   arma::vec x = Rcpp::as<arma::vec>(mat.slot("x"));
//   int nrow = dims[0], ncol = dims[1];
//   arma::sp_mat ans(i, p, x, nrow, ncol);
//   Rprintf("Before %d; After %d\n", prova, ans2(2,1));
//   //if (show) Rcpp::Rcout << res << std::endl;
//   return ans;
//}





// [[Rcpp::export]]
arma::sp_mat modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars) {
  ggmObject *ggm;
  ggm= new ggmObject(y, prCoef, prModel, samplerPars);

  int niter= ggm->niter(), p= ggm->ncol();
  int npars= p*(p+1)/2;

  crossprodmatRcpp S(y, false);
  double s11= S.at(0,0);
  double s12= S.at(0,1);
  Rprintf("s11=%f, s12=%f\n", s11, s12);

  //Test output
  arma::sp_mat ans(niter, npars);
  ans(1,2)= 1;

  delete ggm;

  return ans;
}



