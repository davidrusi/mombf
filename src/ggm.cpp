#include "ggm.h"


//*************************************************************************************
// sp_mat usage examples (from Rcpp gallery)
//*************************************************************************************

/*
//Select one column from a sparse matrix (.row selects a row)
arma::sp_mat col_slice(const arma::sp_mat& x, const int n) {
    return x.col(n - 1);
}

//Create sp_mat by selecting a column of an sp_mat, then using iterator to compute a sum
double sum_by_col(const arma::sp_mat& x) {
    double result = 0;
    for (size_t i = 0; i < x.n_cols; i++) {
        arma::sp_mat col(x.col(i));
        for (arma::sp_mat::iterator j = col.begin(); j != col.end(); ++j) {
            result += *j * (i + 1);
        }
    }
    return result;
}

//Iterate over non-zero elements of a sparse matrix
double sum_by_iterator(const arma::sp_mat& x) {
    double result = 0;
    for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
        result += *i * (i.row() + 1);
    }
    return result;
}

//Identify non-zero elements in a column
IntegerVector findAdjacentStates(sp_mat adjacency, int col) {
    IntegerVector out;
    sp_mat::const_col_iterator start = adjacency.begin_col(col);
    sp_mat::const_col_iterator end = adjacency.end_col(col);
    for ( sp_mat::const_col_iterator i = start; i != end; ++i )
    {
        out.push_back(i.row());
    }
    return out;
}

*/

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

int ggmObject::n() {
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


template <typename mat_type>
void print_mat( mat_type A ) {
  int k, kk, nrow= A->n_rows, ncol= A->n_cols;
  double tmp;
  for (k=0; k < nrow; k++) {
    Rprintf("Row %d: ",k);
    for (kk=0; kk < ncol; kk++) {
      tmp= A->at(k,kk);
      Rprintf("%f, ", tmp);
    }
    Rprintf("\n");
  }
}

//Print matrix on screen
//void print_mat(arma::mat *A) {
//  int k, kk, nrow= A->n_rows, ncol= A->n_cols;
//  for (k=0; k < nrow; k++) {
//    Rprintf("Row %d: ",k);
//    for (kk=0; kk < ncol; kk++) Rprintf("%f, ", A->at(k,kk));
//    Rprintf("\n");
//  }
//}
// 
// 
//void print_spmat(arma::sp_mat *A) {
//  int k, kk, nrow= A->n_rows, ncol= A->n_cols;
//  double tmp;
//  for (k=0; k < nrow; k++) {
//    Rprintf("Row %d: ",k);
//    for (kk=0; kk < ncol; kk++) {
//      tmp= A->at(k,kk);
//      Rprintf("%f, ", tmp);
//    }
//    Rprintf("\n");
//  }
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

  delete ggm;

  return ans;
}


//Drop row and column j from sparse matrix A, return the new matrix in A_minusj
void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int *j) {
  int itcol, itrow;
  arma::sp_mat::const_iterator it;
  for (it= A->begin(); it != A->end(); ++it) {
    itcol= it.col(); itrow= it.row();
    if (itcol == *j) continue; else if (itcol > *j) itcol--;
    if (itrow == *j) continue; else if (itrow > *j) itrow--;
    (*A_minusj)(itrow, itcol)= (*A_minusj)(itcol,itrow)= A->at(it.row(),it.col());
  }
}

/*Gibbs sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)
  - Omegaini: initial value of Omega

OUTPUT
  - ans: sparse matrix where each row corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then ans has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

*/

void GGM_Gibbs(arma::sp_mat *ans, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, *colidx, p= ggm->ncol();

  //Initialize order in which rows of Omega will be sampled
  colidx= ivector(0, p-1);
  for (i=0; i< p; i++) colidx[i]= i; 

  for (i=0; i <= ggm->niter(); i++) {

    samplei_wr(colidx, p, p); //permute row indexes

    for (j=0; j < p; j++) {  //for each column

      //Create Omegaini[-j,-j]
      arma::sp_mat Omega_j(p-1, p-1);
      spmat_droprowcol(&Omega_j, Omegaini, &j);

      //Invert Omegaini[-j,-j]
      arma::mat I= arma::eye(p-1, p-1);
      arma::mat invOmega_j= arma::spsolve(Omega_j, I, "lapack"); //lapack solver, slower but more portable
      //arma::mat invOmega_j= arma::spsolve(Omega_j, I, "superlu"); //superLU solver, faster but requires -lsuperlu compiler flag
      //print_mat(&invOmega_j);  //check

      arma::sp_mat Omegacol= Omegaini->col(colidx[j]);
      GGM_Gibbs_singlecol(ans, i, colidx[j], ggm, &Omegacol,  &invOmega_j); //update row given by colidx[j]
    }

  }

  free_ivector(colidx, 0, p-1);

}

/*Gibbs sampling for column colid of Omega in a Gaussian graphical model

INPUT
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - ans: Omegacol after applying a Gibbs update to its entries
*/

void GGM_Gibbs_singlecol(arma::sp_mat *ans, int iter, int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest) {

  int p= ggm->ncol();
  arma::sp_mat::const_iterator it;

  //Identify non-zero entries in colid
  //IntegerVector sel;
  //for (it= Omegacol->begin(); it != Omegacol->end(); ++it) sel.push_back(it.row());

  //Initialize model and modelnew
  //print_mat <arma::sp_mat *> (Omegacol); //debug
  arma::SpMat<short> model(p, 1), modelnew(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model(it.row(), it.col())= modelnew(it.row(), it.col())= 1;

  //print_mat <arma::SpMat<short> *> (&model); //debug

  //Obtain log-marginal + log-prior for current model

  //For each entry j in colid, obtain log-marginal + log-prior for birth/death of j
  
}


//Compute log-marginal likelihood for model specifying non-zero entries in column colid of Omega, given inverse of Omega[-colid,-colid]
double GGMrow_marg(arma::sp_mat *model, ggmObject *ggm, int *colid, arma::sp_mat *Omega, arma::mat *Omega_colid_inv) {

  return 0.0;

}

//double GGMrow_margupdate(int *changevar, ggmObject *ggm, int *colid, arma::sp_mat *Omega, arma::mat *Omega_colid_inv) {
// 
//}
