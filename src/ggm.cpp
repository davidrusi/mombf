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
arma::sp_mat modelSelectionGGMC(NumericMatrix y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
  ggmObject *ggm;
  ggm= new ggmObject(y, prCoef, prModel, samplerPars);

  int niter= ggm->niter(), p= ggm->ncol();
  int npars= p*(p+1)/2;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), zigzag("zigzag");

  //Create output matrix
  arma::sp_mat ans(niter, npars);

  if (sampler == Gibbs) {

    GGM_Gibbs(&ans, ggm, &Omegaini);

  } else if (sampler == zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //Test crossprodmat
  //crossprodmatRcpp S(y, false);
  //double s11= S.at(0,0);
  //double s12= S.at(0,1);
  //Rprintf("s11=%f, s12=%f\n", s11, s12);

  delete ggm;

  return ans;
}


void GGM_Gibbs(arma::sp_mat *ans, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, *rowidx;
  arma::SpMat<short> model(ggm->ncol(), ggm->ncol());

  //Initialize order in which rows of Omega will be sampled
  rowidx= ivector(0, ggm->ncol() -1);
  for (i=0; i< ggm->ncol(); i++) rowidx[i]= i; 

  //Initialize model
  arma::sp_mat::const_iterator it;
  for (it= Omegaini->begin(); it != Omegaini->end(); ++it) model[it.row(), it.col()]= 1;
    
  for (i=0; i <= ggm->niter(); i++) {

    samplei_wr(rowidx, ggm->ncol(), ggm->ncol()); //permute row indexes

    for (j=0; j <= ggm->ncol(); j++) {
      GGM_Gibbs_singlerow(ans, i, rowidx[j], ggm, Omegaini); //update row given by rowidx[j]
    }

  }

  free_ivector(rowidx, 0, ggm->ncol() -1);

}

void GGM_Gibbs_singlerow(arma::sp_mat *ans, int iter, int rowid, ggmObject *ggm, arma::sp_mat *Omegaini) {


  //double *prova= Omegaini->at(rowid);
  //arma::uword ptr= (Omegaini->col_ptrs)[rowid];
  //double *Omegaini_rowid= Omegaini->at(ptr);
  
  //arma::colvec Omegaini_rowid= Omegaini->col(rowid); //careful, this may copy the whole row

  arma::SpSubview_col<double> Omegaini_rowid= Omegaini->col(rowid); //careful, this compiles but may copy the whole row
  
  Rprintf("Omegaini row %d is %f, %f, %f", rowid, Omegaini_rowid[0], Omegaini_rowid[1], Omegaini_rowid[2]);
  
}

