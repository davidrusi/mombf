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

  arma::vec v = as<arma::vec>(samplerPars["verbose"]);
  if (v[0] == 1) this->verbose= true; else this->verbose= false;


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


double ggmObject::pbirth() {
  return this->samplerPars["pbirth"];
}


int ggmObject::nbirth() {
  return this->samplerPars["nbirth"];
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
List modelSelectionGGMC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, true);

  int niter= ggm->niter(), p= ggm->ncol(), burnin= ggm->burnin();
  int npars= p*(p+1)/2;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (sampler == Gibbs);
  bool use_birthdeath= (sampler == birthdeath);
  bool use_zigzag= (sampler == zigzag);

  //Create output matrix
  arma::sp_mat postSample(npars, niter - burnin);

  arma::mat margpp = arma::zeros(p,p); 
  arma::Mat<int> margppcount;
  margppcount.zeros(p, p);
  arma::vec margppflat(p * (p+1)/2);

  //Obtain posterior samples
  if (use_gibbs || use_birthdeath) {

    GGM_Gibbs(&postSample, &margpp, &margppcount, ggm, &Omegaini);

    margpp = margpp / margppcount;
    for (int i=0; i < p; i++) margpp.at(i,i)= 1;
    symmat2vec(&margppflat, &margpp);

  } else if (use_zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //Free memory and return output
  delete ggm;

  List ret;
  ret["postSample"] = postSample;
  ret["margpp"] = margppflat;
  return ret;

}


//arma::sp_mat modelSelectionGGMC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
//  ggmObject *ggm;
//  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, true);
// 
//  int niter= ggm->niter(), p= ggm->ncol(), burnin= ggm->burnin();
//  int npars= p*(p+1)/2;
// 
//  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
//  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
//  bool use_gibbs= (sampler == Gibbs);
//  bool use_birthdeath= (sampler == birthdeath);
//  bool use_zigzag= (sampler == zigzag);
// 
//  //Create output matrix
//  arma::sp_mat ans(npars, niter - burnin);
// 
//  if (use_gibbs || use_birthdeath) {
// 
//    GGM_Gibbs(&ans, ggm, &Omegaini);
// 
//  } else if (use_zigzag) {
//    
//    Rprintf("zigzag will be implemented soon\n");
// 
//  } else Rf_error("This sampler type is not currently implemented\n");
// 
//  delete ggm;
// 
//  return ans;
//}




/*Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - ans: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then ans has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_Gibbs(arma::sp_mat *ans, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, p= ggm->ncol(), iter= 0, burnin= ggm->burnin(), niter= ggm->niter(), niter10;
  arma::sp_mat::iterator it;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (sampler == Gibbs);
  bool use_birthdeath= (sampler == birthdeath);

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs requires the sampler to be Gibbs or birthdeath");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  for (i=0; i < niter; i++) {

    for (j=0; j < p; j++) {  //for each column

      //Rprintf("i=%d, j=%d\n", i, j); //debug
      //Omegaini->print("Omegaini"); //debug

      arma::mat invOmega_j= get_invOmega_j(Omegaini, j);
      arma::sp_mat Omegacol= Omegaini->col(j);
      arma::sp_mat ans_row(p, 1);
      arma::vec margpp_row= arma::zeros(p);
      arma::Col<int> margppcount_row;
      margppcount_row.zeros(p);

      if (use_gibbs) {
        GGM_Gibbs_singlecol(&ans_row, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&ans_row, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j); //update row given by j
      }

      //Update Omegaini
      spmat_rowcol2zero(Omegaini, j); //set row & column j in Omegaini to zero

      for (it= ans_row.begin(); it != ans_row.end(); ++it) {
        k= it.row(); //ans_row->at(k,0) stores entry (j, k) of ans
        Omegaini->at(j,k)= Omegaini->at(k,j)= ans_row.at(k,0);
      }

      //Copy posterior marginal inclusion probabilities
      for (k=0; k<p; k++) { 
        margpp->at(k,j) += margpp_row.at(k); 
        margpp->at(j,k)= margpp->at(k,j);
        margppcount->at(k,j) += margppcount_row.at(k);
        margppcount->at(j,k)= margppcount->at(k,j);
      }

    } //end for j

    //Copy from Omegaini into ans(,iter)
    if (i >= burnin) {
      if (i >= burnin) spmatsym_save2flat(ans, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &niter10);

  } //end for i

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}


void GGM_Gibbs_new(arma::sp_mat *ans, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, p= ggm->ncol(), iter= 0, burnin= ggm->burnin(), niter= ggm->niter(), niter10, colupdated;
  arma::sp_mat::iterator it;
  arma::sp_mat ans_row(p, 1), newvals, oldvals, difvals;
  arma::mat invOmega_j, invOmega_j_safe;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (sampler == Gibbs);
  bool use_birthdeath= (sampler == birthdeath);

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs requires the sampler to be Gibbs or birthdeath");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  invOmega_j= get_invOmega_j(Omegaini, 0);

  for (i=0; i < niter; i++) {

    for (j=0; j < p; j++) {  //for each column

      //Rprintf("i=%d, j=%d\n", i, j); //debug
      //Omegaini->print("Omegaini"); //debug
      //Update invOmega_j
      if ((i>0) | (j > 0)) {

        if (p >= 50) {

          if (j>0) { colupdated= j-1;  } else { colupdated= p-2; } //index of the column that was last updated
          oldvals= Omegaini->col(j); //previous values in Omega_j
          oldvals.shed_row(colupdated);
          newvals= Omegaini->col(colupdated); //new values in Omega[,j]
          newvals.shed_row(j);
          if (j==0) { //put last entry in newvals as first
            arma::sp_mat first= newvals.row(newvals.n_rows - 1);
            newvals.shed_row(newvals.n_rows - 1);
            newvals= join_cols(first, newvals);
          }
          arma::sp_mat difvals= newvals - oldvals;
          //oldvals.print("oldvals"); //debug
          //newvals.print("newvals"); //debug
          //difvals.print("difvals"); //debug
	   
          //Less accurate matrix inversion, quadratic cost
          if (j>0) {
            symmat_inv_colupdate(&invOmega_j, &difvals, j-1); 
          } else {
            symmat_inv_colupdate(&invOmega_j, &difvals, 0);
            //to do: permute last column into first
          }
          invOmega_j_safe= get_invOmega_j(Omegaini, j); //more accurate, cubic cost
	   
          //double approxerror= norm(invOmega_j - invOmega_j_safe, "fro"); //approx error in rank 1  
          //Rprintf("i=%d, j=%d, approx error= %f\n", i, j, approxerror; //debug
          //invOmega_j.print("invOmega_j (rank 1 update)"); //debug
          //invOmega_j_safe.print("invOmega_j (full update)"); //debug

        } else {

          invOmega_j= get_invOmega_j(Omegaini, j); //more accurate, cubic cost

        }

      }

      arma::sp_mat Omegacol= Omegaini->col(j);
      arma::vec margpp_row= arma::zeros(p);
      arma::Col<int> margppcount_row;
      margppcount_row.zeros(p);

      if (use_gibbs) {
        GGM_Gibbs_singlecol(&ans_row, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&ans_row, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j); //update row given by j
      }

      //Update Omegaini
      spmat_rowcol2zero(Omegaini, j); //set row & column j in Omegaini to zero

      for (it= ans_row.begin(); it != ans_row.end(); ++it) {
        k= it.row(); //ans_row->at(k,0) stores entry (j, k) of ans
        Omegaini->at(j,k)= Omegaini->at(k,j)= ans_row.at(k,0);
      }

      //Copy posterior marginal inclusion probabilities
      for (k=0; k<p; k++) { 
        margpp->at(k,j) += margpp_row.at(k); 
        margpp->at(j,k)= margpp->at(k,j);
        margppcount->at(k,j) += margppcount_row.at(k);
        margppcount->at(j,k)= margppcount->at(k,j);
      }

    } //end for j

    //Copy from Omegaini into ans(,iter)
    if (i >= burnin) {
      if (i >= burnin) spmatsym_save2flat(ans, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &niter10);

  } //end for i

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}




// [[Rcpp::export]]
List GGM_Gibbs_parallelC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
/* Interface for GGM_Gibbs_parallel to be called from R */

  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, true);

  int j, niter= ggm->niter(), p= ggm->ncol(), burnin= ggm->burnin();

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (sampler == Gibbs);
  bool use_birthdeath= (sampler == birthdeath);
  bool use_zigzag= (sampler == zigzag);

  //Create output matrix
  std::list<arma::sp_mat> ans;
  //List ans;
  for (j=0; j<p; j++) {
    ans.push_back(arma::sp_mat(p, niter-burnin));
  }

  //Obtain posterior samples
  if (use_gibbs || use_birthdeath) {

    GGM_Gibbs_parallel(&ans, ggm, &Omegaini);

  } else if (use_zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //Free memory and return output
  delete ggm;

  std::list<arma::sp_mat>::iterator it;
  List ret(p);
  for (it= ans.begin(), j=0; it != ans.end(); ++it, j++) { 
    ret[j]= (*it);
  }

  //return Rcpp::List::create(ans);
  return ret;


}

/*Parallel Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n. Each column is sampled independently from its conditional posterior, given the current value of Omega, i.e. the algorithm targets a pseudo-posterior rather than the actual posterior

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - ans: list where entry j stores posterior samples for Omega[,j], in a sparse matrix format (p rows and number of columns equal to the number of iterations)

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_Gibbs_parallel(std::list<arma::sp_mat> *ans, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int j, p= ggm->ncol(), burnin= ggm->burnin(), niter= ggm->niter(), p10;
  std::list<arma::sp_mat>::iterator it;

  std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (sampler == Gibbs);
  bool use_birthdeath= (sampler == birthdeath);

  arma::vec margpp_row= arma::zeros(p);
  arma::Col<int> margppcount_row;
  margppcount_row.zeros(p);

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs_parallel requires the sampler to be Gibbs or birthdeath");

  if (p >10) { p10= p / 10; } else { p10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  for (it= ans->begin(), j=0; it != ans->end(); ++it, j++) { //for each column  

    arma::mat invOmega_j= get_invOmega_j(Omegaini, j);
    arma::sp_mat Omegacol= Omegaini->col(j);
    arma::sp_mat *ans_row= &(*it);

    if (use_gibbs) {
      GGM_Gibbs_singlecol(ans_row, &margpp_row, &margppcount_row, -burnin, niter-burnin-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j);
    } else {
      GGM_birthdeath_singlecol(ans_row, &margpp_row, &margppcount_row, -burnin, niter-burnin-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j);
    }

  } //end for j

  if (ggm->verbose) print_iterprogress(&j, &p, &p10);

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}


/* Compute the inverse of Omega[-j,-j] */
arma::mat get_invOmega_j(arma::sp_mat *Omega, int j) {
  int p= Omega->n_rows;
  //Create Omega(-colidx[j],-colidx[j])
  arma::sp_mat Omega_j(p-1, p-1);
  spmat_droprowcol(&Omega_j, Omega, j);
  //To do: consider replacing spmat_droprowcol by shed_row, shed_col, and then insert_row, insert_col at end of for j loop

  //Invert Omega(-colidx[j],-colidx[j])
  arma::mat I= arma::eye(p-1, p-1);
  arma::mat invOmega_j= arma::spsolve(Omega_j, I, "lapack"); //lapack solver, slower but more portable
  //arma::mat invOmega_j= arma::spsolve(Omega_j, I, "superlu"); //superLU solver, faster but requires -lsuperlu compiler flag
  return invOmega_j;
}



/* Copy sparse matrix into flattened sparse symmetric matrix

  INPUT
  - ans: a sparse matrix A.
  - col2store: column of ans into which A has to stored in a flattened format

  OUTPUT: A is copied in a flattened way into ans(,col2store).
          If A is p x p then ans must have p * (p+1) / 2 rows.
          Entries in ans(,col2store) are ordered A(1,1), A(2,1), A(2,2), A(3,1)...

*/
void spmatsym_save2flat(arma::sp_mat *ans, arma::sp_mat *A, int col2store) {
  int k, l, pos;
  arma::sp_mat::iterator it;

  for (it= A->begin(); it != A->end(); ++it) {
    k= it.row(); l= it.col();
    if (k < l) {
      pos= l*(l+1)/2 + k;
    } else if (k == l) {
      pos= (l+1)*(l+2)/2 - 1;
    } else {
      continue;
    }
    ans->at(pos, col2store)= A->at(k,l);
  }

}


/*Gibbs sampling for column colid of Omega in a Gaussian graphical model

INPUT
  - iterini, iterfi: sampled values of Omega are stored in ans[iterini:iterfi,]. Negative entries in iterini:iterfi are not stored (burnin)
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - ans: each column of ans contains Omegacol after applying a Gibbs update to its entries

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

IMPORTANT: at input ans should only contain zeroes

*/

void GGM_Gibbs_singlecol(arma::sp_mat *ans, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest) {

  int i, j, p= ggm->ncol();
  double mcurrent, mnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  //sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false);

  for (i= iterini; i <= iterfi; i++) {
    //For each entry j in colid, obtain log-marginal + log-prior for adding/removing j
    for (j= 0; j < p; j++) {

      if (j == (int) colid) continue;  //diagonal entry is always in

      if (model->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;

      ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);

      ppnew = exp(mnew - mcurrent);
      ppnew /= (1.0 + ppnew);

      if (model->at(j,0) == 0) margpp->at(j) += ppnew; else margpp->at(j) += 1 - ppnew;
      margppcount->at(j) ++;

      if (runif() < ppnew) { //if new model is accepted

        model_tmp_ptr= model;
        model= modelnew;
        modelnew= model_tmp_ptr;

      }

      modelnew->at(j,0)= model->at(j,0);

    } //end j for (iteration over columns)

    sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
    ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, true); //update parameter values

    //Copy current model into ans[,i] as flat vector
    if (i >= 0) save_ggmsample_col(ans, model, &sample_diag, sample_offdiag, i, colid);  

  } //end i for (Gibbs iterations)

  delete sample_offdiag;
  delete model;
  delete modelnew;
  delete ms;
  
}



/*Birth-death sampling for column colid of Omega in a Gaussian graphical model

INPUT
  - iterini, iterfi: sampled values of Omega are stored in ans[iterini:iterfi,]. Negative entries in iterini:iterfi are not stored (burnin)
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - ans: each column of ans contains Omegacol after applying a total of (iterfi - iterini + 1) birth-death updates

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

IMPORTANT: at input ans should only contain zeroes

*/

void GGM_birthdeath_singlecol(arma::sp_mat *ans, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest) {

  bool birth;
  int i, j, p= ggm->ncol(), nbirth= ggm->nbirth(), idx_update;
  double pbirth= ggm->pbirth(), mcurrent, mnew, dpropcurrent, dpropnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false);

  for (i= iterini; i <= iterfi; i++) {

      for (j=1; j<=nbirth; j++) {  //Make nbirth birth/death updates

        //Make birth/death proposal
        arma::SpMat<short> model_colid= (*model);
        model_colid.shed_row(colid);
        rbirthdeath(&idx_update, &birth, &model_colid, pbirth);
         
        //Store birth/death proposal into modelnew_colid and modelnew
        arma::SpMat<short> modelnew_colid= model_colid;
        if (birth) modelnew_colid.at(idx_update,0)= 1; else modelnew_colid.at(idx_update,0)= 0;
         
        if (idx_update >= (int) colid) idx_update++;
        if (birth) modelnew->at(idx_update,0)= 1; else modelnew->at(idx_update,0)= 0;
         
        dpropnew= dbirthdeath(&modelnew_colid, &model_colid, pbirth, true);
        dpropcurrent= dbirthdeath(&model_colid, &modelnew_colid, pbirth, true);
         
        //model->print("model"); //debug
        //Rprintf("idx_update= %d; birth= %B\n", idx_update, birth); //debug
        //modelnew->print("modelnew"); //debug
        //Rprintf("dpropnew= %f; dpropcurrent= %f\n", dpropnew, dpropcurrent); //debug
         
        //Obtain posterior sample for modelnew
        ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);
         
        ppnew = exp(mnew - mcurrent + dpropcurrent - dpropnew);
        ppnew /= (1.0 + ppnew);

        //Update posterior marginal inclusion probabilities
        if (birth) margpp->at(idx_update) += ppnew; else margpp->at(idx_update) += 1 - ppnew;
        margppcount->at(idx_update) ++;
         
        if (runif() < ppnew) { //if new model is accepted
         
          model_tmp_ptr= model;
          model= modelnew;
          modelnew= model_tmp_ptr;
         
        }
         
        modelnew->at(idx_update,0)= model->at(idx_update,0);

    } //end for j

    sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
    ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, true); //update parameter values

    //Copy current model into ans[,i] as flat vector
    if (i >= 0) save_ggmsample_col(ans, model, &sample_diag, sample_offdiag, i, colid);  

  } //end i for (Gibbs iterations)

  delete sample_offdiag;
  delete model;
  delete modelnew;
  delete ms;
  
}


/* Save sample_diag and sample_diag into ans[,col2save] */
void save_ggmsample_col(arma::sp_mat *ans, arma::SpMat<short> *model, double *sample_diag, arma::mat *sample_offdiag, int col2save, unsigned int colid) {

  int j;
  unsigned int k;
  arma::SpMat<short>::iterator it_short;

  for (it_short= model->begin(), j=0; it_short != model->end(); ++it_short) {
    k= it_short.row();
    if (k != colid) { ans->at(k,col2save)= sample_offdiag->at(j,0); j++; } else ans->at(k,col2save)= *sample_diag;
  }

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


/************************************************************************************************

 MATRIX MANIPULATION

************************************************************************************************/


//Set row and colum colid of A to 0
void spmat_rowcol2zero(arma::sp_mat *A, int colid) {
  int j;
  long unsigned int k;
  std::vector<int> idx2delete;
  arma::sp_mat::iterator it;
  //Find indexes of rows to be deleted
  for (it= A->begin_col(colid); it != A->end_col(colid); ++it) idx2delete.push_back(it.row());
  //Set entries to zero
  for (k=0; k < idx2delete.size(); k++) {
    j= idx2delete.at(k);
    A->at(j,colid)= A->at(colid,j)= 0;
  }
}

//Drop row and column j from sparse matrix A, return the new matrix in A_minusj
void spmat_droprowcol(arma::sp_mat *A_minusj, arma::sp_mat *A, int j) {
  int itcol, itrow;
  arma::sp_mat::const_iterator it;
  for (it= A->begin(); it != A->end(); ++it) {
    itcol= it.col(); itrow= it.row();
    if (itcol == j) continue; else if (itcol > j) itcol--;
    if (itrow == j) continue; else if (itrow > j) itrow--;
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


//Flatten symmetric matrix A, in column-wise order
void symmat2vec(arma::vec *Aflat, arma::mat *A) {
  unsigned int i, j, k=0, p= A->n_rows;
  if (Aflat->size() != p * (p+1)/2) Rf_error("Error in symmat2vec: matrix dimensions don't match");

  for (i=0; i<p; i++) {
    for (j=0; j<=i; j++) {
      Aflat->at(k)= A->at(i,j);
      k++;
    }
  }

}
