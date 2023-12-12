
#include "modselIntegrals.h"
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
ggmObject::ggmObject(arma::mat *y, List prCoef, List prModel, List samplerPars, bool computeS=true) {

  this->y = y;
  this->prCoef = prCoef;
  this->prModel = prModel;
  this->samplerPars = samplerPars;

  if (computeS) {
    this->S= (*y).t() * (*y);
  }

}

//Destructor
ggmObject::~ggmObject() {

}

int ggmObject::n() {
  return this->y->n_rows;
  //return (this->y).nrow();
}

int ggmObject::ncol() {
  return this->y->n_cols;
  //return (this->y).ncol();
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


/*Print matrix on screen. The function works for any matrix that can be printed as a double

Usage examples:

arma::sp_mat A(5,2);
print_mat <arma::sp_mat *> (&A);

arma::SpMat<short> As(5, 1);
print_mat <arma::SpMat<short> *> (&As);

*/
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





// [[Rcpp::export]]
arma::sp_mat modelSelectionGGMC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, true);

  int niter= ggm->niter(), p= ggm->ncol(), burnin= ggm->burnin();
  int npars= p*(p+1)/2;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), zigzag("zigzag");

  //Create output matrix
  arma::sp_mat ans(npars, niter - burnin);

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

//Copy A[model,model] into Aout
void copy_submatrix(arma::mat *Aout, arma::mat *A, arma::SpMat<short> *model) {
  int i, j;
  arma::SpMat<short>::iterator it, it2;

  for (it= model->begin(), i=0; it != model->end(); ++it, i++) {
    for (it2= model->begin(), j=0; it2 != model->end(); ++it2, j++) {
        Aout->at(i,j)= A->at(it.row(), it2.row());
    }
  }
}


/*Gibbs sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - ans: sparse matrix where each row corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then ans has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_Gibbs(arma::sp_mat *ans, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, *colidx, p= ggm->ncol(), iter= -1, burnin= ggm->burnin();
  arma::sp_mat ans_row(p, 1);
  arma::sp_mat::iterator it;

  //Initialize order in which rows of Omega will be sampled
  colidx= ivector(0, p-1);
  for (i=0; i< p; i++) colidx[i]= i; 

  for (i=0; i < ggm->niter(); i++) {

    //samplei_wr(colidx, p, p); //permute row indexes

    for (j=0; j < p; j++) {  //for each column

      //Create Omegaini(-colidx[j],-colidx[j])
      arma::sp_mat Omega_j(p-1, p-1);
      spmat_droprowcol(&Omega_j, Omegaini, colidx + j);
      //Consider replacing spmat_droprowcol by shed_row, shed_col, and then insert_row, insert_col at end of for j loop

      //Invert Omegaini(-colidx[j],-colidx[j])
      arma::mat I= arma::eye(p-1, p-1);
      arma::mat invOmega_j= arma::spsolve(Omega_j, I, "lapack"); //lapack solver, slower but more portable
      //arma::mat invOmega_j= arma::spsolve(Omega_j, I, "superlu"); //superLU solver, faster but requires -lsuperlu compiler flag

      arma::sp_mat Omegacol= Omegaini->col(colidx[j]);
      if (i >= burnin) iter++;
      GGM_Gibbs_singlecol(&ans_row, iter, iter, (unsigned int) colidx[j], ggm, &Omegacol,  &invOmega_j); //update row given by colidx[j]

      //Update Omegaini
      for (it= ans_row.begin(); it != ans_row.end(); ++it) {
        k= it.row(); //ans_row->at(k,0) stores entry (colid, k) of ans
        Omegaini->at(colidx[j],k)= Omegaini->at(k,colidx[j])= ans_row.at(k,0);
      }     

      //Copy from Omegaini into ans(,iter)
      if (i >= burnin) save_spmat2_flatspmat(ans, Omegaini, iter);

    }

  }

  free_ivector(colidx, 0, p-1);

}


/* Copy sparse matrix into flattened sparse symmetric matrix

  INPUT
  - ans: a sparse matrix A.
  - col2store: column of ans into which A has to stored in a flattened format

  OUTPUT: A is copied in a flattened way into ans(,col2store).
          If A is p x p then ans must have p * (p+1) / 2 rows.
          Entries in ans(,col2store) are ordered A(1,1), A(2,1), A(2,2), A(3,1)...

*/
void save_spmat2_flatspmat(arma::sp_mat *ans, arma::sp_mat *A, int col2store) {
  int k, l, pos;
  arma::sp_mat::iterator it;

  for (it= A->begin(); it != A->end(); ++it) {
    k= it.row(); l= it.col();
    if (k <= l) {
      pos= (k+1)*(k+2)/2 - 1 + l - k;
      ans->at(pos, col2store)= A->at(k,l);
    }
  }

}


/* Copy sparse row (1-column sparse matrix) into flattened sparse symmetric matrix

  INPUT
    - ans_row: stores column colid of matrix A as a 1-column sparse matrix
    - colid: column of A stored in ans_row
    - iter: column of ans in which to copy A[,colid]

  OUTPUT: ans[,iter] stores the matrix A in flat format, A(1,1), A(2,1), A(2,2), A(3,1), ...


*/
//void save_spvec2_flatspmat(arma::sp_mat *ans, arma::sp_mat *ans_row, int colid, int iter) {
// 
//  int pos, posj= colidx[j] * (colid+1) / 2;
//  arma::sp_mat::iterator it;
// 
//  for (it= ans_row->begin(); it != ans_row->end(); ++it) {
//    k= it.row(); //ans_row->at(k,0) stores entry (colid, k) of ans
//    //Stored as flat vector (colid,k) goes to column start index + abs(colid - k)
//    //Column start index is k * (k+1)/2 if k < colid, else it's posj obtained above
//    if (colid > k) pos= k * (k+1) / 2 + k - colid; else pos= posj + colid - k;
//    ans->at(pos,iter)= ans_row->at(k,0);
// 
//  }     
// 
//}

/*Gibbs sampling for column colid of Omega in a Gaussian graphical model

INPUT
  - iterini, iterfi: sampled values of Omega are stored in ans[iterini:iterfi,]. Negative entries in iterini:iterfi are not stored (burnin)
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - ans: each column of ans contains Omegacol after applying a Gibbs update to its entries
*/

void GGM_Gibbs_singlecol(arma::sp_mat *ans, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest) {

  int i, j, p= ggm->ncol();
  unsigned int k;
  double mcurrent, mnew, ppnew, sample_diag, samplenew_diag;
  arma::mat *sample_offdiag, *samplenew_offdiag;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  arma::SpMat<short>::iterator it_short, it2_short;
  //arma::SpMat<short> model(p, 1), modelnew(p, 1), *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model);

  for (i= iterini; i <= iterfi; i++) {
    //For each entry j in colid, obtain log-marginal + log-prior for adding/removing j
    for (j= 0; j < p; j++) {

      if (j == (int) colid) continue;  //diagonal entry is always in

      if (model->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;
      samplenew_offdiag= new arma::mat(modelnew->n_nonzero - 1, 1);

      ms->getJoint(&mnew, samplenew_offdiag, &samplenew_diag, modelnew);

      ppnew = exp(mnew - mcurrent);
      ppnew /= (1.0 + ppnew);

      if (runif() < ppnew) { //if new model is accepted

        model_tmp_ptr= model;
        model= modelnew;
        modelnew= model_tmp_ptr;

        delete sample_offdiag;
        sample_offdiag= samplenew_offdiag;
        sample_diag= samplenew_diag;

      } else {

        delete samplenew_offdiag;

      }

      modelnew->at(j,0)= model->at(j,0);

    } //end j for (iteration over columns)

    //Copy current model into ans[,i] as flat vector
    if (i >= 0) {

      for (it_short= model->begin(), j=0; it_short != model->end(); ++it_short) {
        k= it_short.row();
        if (k != colid) { ans->at(k,0)= sample_offdiag->at(j,0); j++; } else ans->at(k,0)= sample_diag;
      }

    }

  } //end i for (Gibbs iterations)

  delete sample_offdiag;
  delete model;
  delete modelnew;
  
}





/* Compute log-joint (log-marginal likelihood + log prior) for model specifying non-zero entries in column colid of Omega, given the entries selected by model of the inverse of Omega[-colid,-colid]

  INPUT
  - model: entries that are non-zero
  - colid: column id
  - ggm: object storing info about the Gaussian graphical model
  - Omegainv_model: entries selected by model from Omega[-colid,-colid]

  OUTPUT
  - logjoint: log-marginal + log-prior for model
  - m: u1 ~ multivariate Normal with mean m and covariance Uinv, where u1= -Omega[-colid,colid] are the (negative) off-diagonal entries
  - cholUinv: Cholesky decomposition of Uinv, i.e. Uinv= t(cholUinv) * cholUinv

*/
void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model) {

  unsigned int npar= model->n_nonzero -1;
  arma::vec tau = as<arma::vec>(ggm->prCoef["tau"]); //Prior is Omega_{jk} | Omega_{jk} != 0 ~ N(0, tau)
  arma::SpMat<short> model_offdiag(ggm->ncol() -1, 1);

  if (npar > 0) { //if some off-diagonal entry is non-zero

    unsigned int i, j;
    arma::vec lambda= as<arma::vec>(ggm->prCoef["lambda"]); //Prior is Omega_{jj} ~ Exp(lambda)
    double tauinv= 1.0/tau[0];
    arma::mat U(npar,npar), s(npar, 1);
    arma::SpMat<short>::iterator it;

    //model_offdiag indicates non-zero off-diagonal elements
    //copy sample covariance(model_offdiag, colid) into s
    for (it= model->begin(), i=0; it != model->end(); ++it) {
      if (it.row() == colid) continue;
      model_offdiag.at(i, 0)= model->at(it.row(), 0);
      s.at(i,0)= ggm->S.at(it.row(), colid);
      i++;
    }
     
    //Create U matrix
    double ct= lambda[0] + ggm->S.at(colid, colid);
    for (i=0; i<npar; i++) {
      U.at(i,i)= ct * Omegainv_model->at(i,i) + tauinv;
      for (j=i+1; j<npar; j++) U.at(i,j)= U.at(j,i)= ct * Omegainv_model->at(i,j);
    }
     
    double logdetUinv= 0;
    arma::mat Uinv(npar,npar);
    choldcinv_det(&Uinv, cholUinv, &logdetUinv, &U);
    //arma::mat Uinv= inv_sympd(U);
     
    (*m)= Uinv * s;
     
    (*logjoint) = 0.5 * (arma::as_scalar(m->t() * U * (*m)) - ((double) npar) * log(tau[0]) + logdetUinv);

  } else { //if all off-diagonal entries are zero

    (*logjoint) = - ((double) npar) * log(tau[0]);

  }

  (*logjoint) += logprior_GGM(&model_offdiag, ggm);

}


double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm) {
  double ans;
  string priorlabel= as<string> (ggm->prModel["priorlabel"]);

  if (priorlabel == "binomial") {
//  if (ggm->prModel["priorlabel"] == "binomial") {
    double npar= (double) model->n_nonzero; 
    double p= as<double>(ggm->prModel["priorPars.p"]); 
    //arma::vec p= as<arma::vec>(ggm->prModel["priorPars.p"]); 

    ans= npar * log(p) + ((double) (model->n_rows) - npar) * log(1-p);

  } else Rf_error("This model prior is not implemented\n");
  
  return ans;
}


//double GGMrow_margupdate(int *changevar, ggmObject *ggm, int *colid, arma::sp_mat *Omega, arma::mat *Omega_colid_inv) {
// 
//}
