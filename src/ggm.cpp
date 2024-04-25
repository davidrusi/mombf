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
ggmObject::ggmObject(arma::mat *y, List prCoef, List prModel, List samplerPars, bool use_tempering, bool computeS=true) {

  //this->y = y;
  //this->prCoef = prCoef;
  //this->prModel = prModel;
  //this->samplerPars = samplerPars;

  //Data summaries
  this->n= y->n_rows;
  this->ncol= y->n_cols;

  if (computeS) {
    this->S= (*y).t() * (*y);
  } else {
    Rf_error("Error in ggmObject. Currently computeS must be set to true");
  }


  //Parameters of prior distribution  
  this->prCoef_lambda= as<double>(prCoef["lambda"]);
  this->prCoef_tau = as<double>(prCoef["tau"]);

  this->priorlabel= as<string> (prModel["priorlabel"]);
  this->priorPars_p= as<double>(prModel["priorPars.p"]);

  //MCMC sampler parameters
  CharacterVector samplerR= samplerPars["sampler"];
  std::string samplerC = Rcpp::as<std::string>(samplerR);
  this->sampler= samplerC;

  this->burnin= as<int>(samplerPars["burnin"]);
  this->nbirth= as<int>(samplerPars["nbirth"]);
  this->niter= as<int>(samplerPars["niter"]);
  this->tempering= as<double>(samplerPars["tempering"]);
  this->pbirth= as<double>(samplerPars["pbirth"]);
  this->use_tempering= use_tempering;

  //Set print progress iteration to true/false
  arma::vec v = as<arma::vec>(samplerPars["verbose"]);
  if (v[0] == 1) this->verbose= true; else this->verbose= false;

}

//Destructor
ggmObject::~ggmObject() {

}

//int ggmObject::n() { return this->y->n_rows; }

//int ggmObject::ncol() { return this->y->n_cols; }

//Posterior sampler type
//CharacterVector ggmObject::sampler() { return this->samplerPars["sampler"]; }

//int ggmObject::niter() { return this->samplerPars["niter"]; }

//int ggmObject::burnin() {  return this->samplerPars["burnin"]; }

//double ggmObject::pbirth() return this->samplerPars["pbirth"]; }

// int ggmObject::nbirth() { return this->samplerPars["nbirth"]; }

// double ggmObject::tempering() { return this->samplerPars["tempering"]; }





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
  bool use_tempering= false;
  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, use_tempering, true);

  int p= ggm->ncol, npars= p*(p+1)/2;

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_zigzag= (ggm->sampler == zigzag);

  //Create output matrix
  arma::sp_mat postSample(npars, ggm->niter - ggm->burnin);

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




/*Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - samples: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then samples has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_Gibbs(arma::sp_mat *samples, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, p= ggm->ncol, iter= 0, burnin= ggm->burnin, niter= ggm->niter, niter10;
  arma::sp_mat::iterator it;

  //std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  arma::SpMat<short> *models= NULL;
  arma::mat *model_logprob= NULL;

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs requires the sampler to be Gibbs or birthdeath");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  for (i=0; i < niter; i++) {

    for (j=0; j < p; j++) {  //for each column

      arma::mat invOmega_j= get_invOmega_j(Omegaini, j);
      arma::sp_mat Omegacol= Omegaini->col(j);
      arma::sp_mat samples_row(p, 1);
      arma::vec margpp_row= arma::zeros(p);
      arma::Col<int> margppcount_row;
      margppcount_row.zeros(p);

      if (use_gibbs) {
        GGM_Gibbs_singlecol(&samples_row, models, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprob); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&samples_row, models, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprob); //update row given by j
      }

      //Update Omegaini
      spmat_rowcol2zero(Omegaini, j); //set row & column j in Omegaini to zero

      for (it= samples_row.begin(); it != samples_row.end(); ++it) {
        k= it.row(); //samples_row->at(k,0) stores entry (j, k) of samples
        Omegaini->at(j,k)= Omegaini->at(k,j)= samples_row.at(k,0);
      }

      //Copy posterior marginal inclusion probabilities
      for (k=0; k<p; k++) { 
        margpp->at(k,j) += margpp_row.at(k); 
        margpp->at(j,k)= margpp->at(k,j);
        margppcount->at(k,j) += margppcount_row.at(k);
        margppcount->at(j,k)= margppcount->at(k,j);
      }

    } //end for j

    //Copy from Omegaini into samples(,iter)
    if (i >= burnin) {
      if (i >= burnin) spmatsym_save2flat(samples, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &niter10);

  } //end for i

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}



//Not working yet. Attempt at improving GGM_Gibbs using rank 1 updates to obtain invOmega_j
void GGM_Gibbs_new(arma::sp_mat *samples, arma::mat *margpp, arma::Mat<int> *margppcount, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, p= ggm->ncol, iter= 0, burnin= ggm->burnin, niter= ggm->niter, niter10, colupdated;
  arma::sp_mat::iterator it;
  arma::sp_mat samples_row(p, 1), newvals, oldvals, difvals;
  arma::mat invOmega_j, invOmega_j_safe;

  //std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  arma::SpMat<short> *models= NULL;
  arma::mat *model_logprob= NULL;

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs requires the sampler to be Gibbs or birthdeath");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  invOmega_j= get_invOmega_j(Omegaini, 0);

  for (i=0; i < niter; i++) {

    for (j=0; j < p; j++) {  //for each column

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
        GGM_Gibbs_singlecol(&samples_row, models, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprob); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&samples_row, models, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprob); //update row given by j
      }

      //Update Omegaini
      spmat_rowcol2zero(Omegaini, j); //set row & column j in Omegaini to zero

      for (it= samples_row.begin(); it != samples_row.end(); ++it) {
        k= it.row(); //samples_row->at(k,0) stores entry (j, k) of samples
        Omegaini->at(j,k)= Omegaini->at(k,j)= samples_row.at(k,0);
      }

      //Copy posterior marginal inclusion probabilities
      for (k=0; k<p; k++) { 
        margpp->at(k,j) += margpp_row.at(k); 
        margpp->at(j,k)= margpp->at(k,j);
        margppcount->at(k,j) += margppcount_row.at(k);
        margppcount->at(j,k)= margppcount->at(k,j);
      }

    } //end for j

    //Copy from Omegaini into samples(,iter)
    if (i >= burnin) {
      if (i >= burnin) spmatsym_save2flat(samples, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &niter10);

  } //end for i

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}




// [[Rcpp::export]]
List GGM_Gibbs_parallelC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
/* Interface for GGM_Gibbs_parallel to be called from R 

   Output: list with p+2 elements. 
     - Elements j=0,...,p-1 (proposal_samples) are arma::sp_mat storing proposed samples for Omega[,j] given Omegaini[-j,-j]
     - Element p (propdens) is a matrix where column j stores the proposal density for the proposed samples of Omega[,j]
     - Element p+1 (postSample) are posterior samples for Omega, obtained from the proposals in proposal_samples
     - Element p+2 (prop_accept) is the proportion of accepted proposals
*/

  bool use_tempering= true;
  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, use_tempering, true);

  int j, niter= ggm->niter, p= ggm->ncol, burnin= ggm->burnin, npars= p*(p+1)/2;
  double *dpropini, prop_accept;

  //std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_zigzag= (ggm->sampler == zigzag);

  arma::mat propdens(p, niter-burnin);

  //Allocate memory
  dpropini= dvector(0, p-1);

  //Allocate vector where to store proposed models
  std::vector<arma::SpMat<short>> proposal_samples;
  for (j=0; j<p; j++) {
    proposal_samples.push_back(arma::SpMat<short>(p, niter-burnin));
  }

  //Propose models, store their proposal probability into propdens
  if (use_gibbs || use_birthdeath) {

    GGM_Gibbs_parallel(&proposal_samples, ggm, &Omegaini, &propdens);

  } else if (use_zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //MCMC using independent proposal MH to combine the chains
  ggm->use_tempering= false;
  arma::sp_mat postSample(npars, niter - burnin);
  GGM_parallel_MH_indep(&postSample, &prop_accept, &proposal_samples, &propdens, dpropini, ggm, &Omegaini);

  //Return output
  List ret(p+3);
  //List ret(2);
  ret[0]= postSample;
  ret[1]= prop_accept;

  //Store proposals into returned list (ret)
  std::vector<arma::SpMat<short>>::iterator it;
  for (it= proposal_samples.begin(), j=2; it != proposal_samples.end(); ++it, j++) { 
    ret[j]= (*it);
  }

  ret[p+2]= propdens;
 
   //Free memory and return output
  delete ggm;
  free_dvector(dpropini, 0, p-1);

  return ret;


}



/* Evaluate proposal density for GGM samples proposed in parallel

  Input
  - proposal_samples: vector where j=0,...,p-1 is an sp_mat with p rows and niter columns, storing the proposed values for Omega[,j] given Omegaini[-j,-j]
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)
  - Omegaini: initial guess at Omega used to propose new values for Omega

  Output
  - propdens: a matrix with j=0,...,p-1 rows and niter columns. propdens[j,i] is the proposal density for ith sample in proposal_samples[j]
  - dpropini: log-proposal density for Omegaini[,j] given Omegaini[-j,-j]

*/
/*void GGM_parallel_propdensity(arma::mat *propdens, double *dpropini, std::vector<arma::sp_mat> *proposal_samples, ggmObject *ggm, arma::sp_mat *Omegaini) {
  int i, j, nsamples;
  arma::vec lambda= as<arma::vec>(ggm->prCoef["lambda"]); //Prior is Omega_{jj} ~ Exp(lambda)
  double lambdahalf= 0.5 * lambda[0];
  std::vector<arma::sp_mat>::iterator it;

  for (it= proposal_samples->begin(), j=0; it != proposal_samples->end(); ++it, j++) { //for each column  

    double b= -0.5 * (ggm->S).at(j,j) + lambdahalf;
    arma::sp_mat Omega_j= (*Omegaini);

    arma::sp_mat colj= Omega_j.col(j);
    Omega_j.shed_row(j);
    Omega_j.shed_col(j);
    arma::mat invOmega_j= get_invOmega_j(Omegaini, j); //to do: use low-rank update

    colj.shed_row(j);
    dpropini[j]= b * arma::as_scalar(colj.t() * invOmega_j * colj);

    nsamples= (*it).n_cols;

    for (i= 0; i < nsamples; i++) {
        arma::sp_mat samplei= (*it).col(i);
        samplei.shed_row(j); 
        propdens->at(j,i)= b * arma::as_scalar(samplei.t() * invOmega_j * samplei);
    }

  }

}
*/


/* Obtain posterior samples for Omega via Metropolis-Hastings, using the parallel proposals

INPUT

  - proposal_samples: vector with parallel proposals in its entries 0,...,p-1. The proposals are model indicators (precision matrix entries being zero/non-zero)
  - propdens: matrix with log-proposal densities. Element (i,j) is the density for proposed sample i in column j
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT

  - postSample: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then postSample has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.
  - prop_accept: proportion of accepted proposals

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned
  - dpropini: vector where entry j is the log-proposal density for Omegaini[,j] given Omegaini[-j,-j]. Updated at each iteration, the value at the last iteration is returned

*/
void GGM_parallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_samples, arma::mat *propdens, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, iter, newcol, number_accept= 0, proposal_idx, burnin= ggm->burnin, niter= burnin + postSample->n_cols, p= ggm->ncol;
  double dpostnew, dpostold, dpropnew, dpropold, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::SpMat<short>> modelold(p);
  arma::SpMat<short> modelnew;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  //Initialize MCMC state
  /*
  for (j=0; j < p; j++) {

    //Init modelold
    arma::sp_mat Omegaini_colj= Omegaini->col(j);
    arma::SpMat<short> model(p,1);
    arma::sp_mat::iterator it;
    for (it= Omegaini_colj.begin(); it != Omegaini_colj.end(); ++it) {
      model.at(it.row(),0)= 1;
    }
    modelold[j]= model;

    //Init dpropini
    arma::mat invOmega_j= get_invOmega_j(Omegaini, j); //to do: use low-rank update

    modselIntegrals_GGM *ms;
    ms= new modselIntegrals_GGM(marfun, ggm, j, &invOmega_j);
    ms->getJoint(dpropini + j, sample_offdiag, &sample_diag, &(modelold[j]), false); //log-proposal for modelold[i]
    delete ms;

  }
  */
  for (newcol=0; newcol < p; newcol++) {
    proposal_idx= runifdisc(0, niter - burnin -1); //index of proposal
    modelold[newcol]= ((*proposal_samples)[newcol]).col(proposal_idx); //proposed model
    dpropini[newcol]= propdens->at(newcol, proposal_idx);  //log-proposal for proposed model

    //Init Omegaini[,newcol]
    arma::mat invOmega_newcol= get_invOmega_j(Omegaini, newcol); //to do: use low-rank update
    modselIntegrals_GGM *ms;
    ms= new modselIntegrals_GGM(marfun, ggm, newcol, &invOmega_newcol);
    sample_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
    ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &(modelold[newcol]), true); //sample Omega[,newcol]
    update_Omegaini(Omegaini, &newcol, &sample_diag, &modelnew, sample_offdiag); //copy sampled values into Omegaini[,newcol]

    delete ms;
    delete sample_offdiag;
  }


  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {
    
    newcol= runifdisc(0, p-1); //column to update
    proposal_idx= runifdisc(0, niter - burnin -1); //index of proposal
    
    modelnew= ((*proposal_samples)[newcol]).col(proposal_idx); //proposed value

    arma::mat invOmega_newcol= get_invOmega_j(Omegaini, newcol); //to do: use low-rank update

    modselIntegrals_GGM *ms;
    ms= new modselIntegrals_GGM(marfun, ggm, newcol, &invOmega_newcol);
    ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &modelnew, false); //log-posterior for proposed model
    ms->getJoint(&dpostold, sample_offdiag, &sample_diag, &(modelold[newcol]), false); //log-posterior for current model

    dpropnew= propdens->at(newcol, proposal_idx);  //log-proposal for proposed model
    dpropold= dpropini[newcol]; //log-proposal for current model

    ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);
    //Rprintf("newcol=%d", newcol); //debug
    //modelnew.print("modelnew"); //debug
    //(modelold[newcol]).print("modelold[newcol]"); //debug
    //Rprintf("dpostnew=%f, dpostold=%f, dpropnew=%f, dpropold=%f", dpostnew, dpostold, dpropnew, dpropold); //debug

    if ((ppnew > 1) | (runifC() < ppnew)) { //if update is accepted

      if (i >= burnin) number_accept++;

      dpropini[newcol]= dpropnew;
      modelold[newcol]= modelnew;

      //Sample Omega[,newcol]
      sample_offdiag= new arma::mat(modelnew.n_nonzero - 1, 1);
      ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &modelnew, true); //sample Omega[,newcol]

      update_Omegaini(Omegaini, &newcol, &sample_diag, &modelnew, sample_offdiag); //copy sampled values into Omegaini[,newcol]

      delete sample_offdiag;

    }

    //Copy from Omegaini into postSample[,iter]
    if (i >= burnin) {
      spmatsym_save2flat(postSample, Omegaini, iter);
      iter++;
    }

    delete ms;

  }

  (*prop_accept)= (number_accept + 0.0) / (iter + 0.0);

}

/*Update symmetric Omegaini[,newcol] by storing Omegaini[newcol,newcol]= sample_diag and storing sample_offdiag into the off-diagonal entries indicated by modelnew
  Input
  - newcol: column to be updated
  - sample_diag: new diagonal entry to be stored
  - modelnew: sparse column vector (dimension nrow(Omegaini) x 1) indicating what off-diagonal entries in Omegaini[,newcol] will be non-zero
  - sample_offdiag: dense vector of length equal to modelnew->n_nonzero storing the non-zero entries
  Output
  - Omegaini: row and column newcol are updated with the values in sample_diag and sample_offdiag. The rest of Omegaini is not altered

*/
void update_Omegaini(arma::sp_mat *Omegaini, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag) {
  int j, k;

  //Update Omegaini
  spmat_rowcol2zero(Omegaini, *newcol); //Set row and colum newcol to 0

  //Copy (sample_offdiag, sample_diag) to Omegaini
  Omegaini->at(*newcol, *newcol)= *sample_diag; //save diagonal element
  arma::SpMat<short>::iterator it;
  for (it= modelnew->begin(), j=0; it != modelnew->end(); ++it) {
    k= it.row();
    if (k != *newcol) { //skip diagonal element (already saved earlier)
      Omegaini->at(k, *newcol)= Omegaini->at(*newcol, k)= sample_offdiag->at(j,0); 
      j++;
    } 
  }

}


/* Old version where the proposal provided samples for Omega, rather than proposing only model indicators 
void GGM_parallel_MH_indep(arma::sp_mat *postSample, std::vector<double> *prop_accept, std::vector<arma::sp_mat> *proposal_samples, arma::mat *propdens, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, k, iter, newcol, proposal_idx, niter= postSample->n_cols, burnin= ggm->burnin, p= ggm->ncol;
  double dpostnew, dpostold, dpropnew, dpropold, diagnew, ppnew;
  arma::vec lambda= as<arma::vec>(ggm->prCoef["lambda"]); //Prior is Omega_{jj} ~ Exp(lambda)
  double lambdahalf= 0.5 * lambda[0];
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::sp_mat> u1old(p);
  arma::sp_mat u1new;
  double a= 0.5 * (double) ggm->n() + 1.0;
  std::vector<int> number_accept(p, 0);
  
  //Pre-compute constants & initialize MCMC state
  for (i=0; i < p; i++) {
    bnew[i]= -0.5 * (ggm->S).at(i,i) + lambdahalf;

    u1old[i]= Omegaini->col(i);
    diagcur[i]= (u1old[i]).at(i,0);
    (u1old[i]).shed_row(i);
  }

  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {

    newcol= runifdisc(0, p-1); //column to update
    proposal_idx= runifdisc(0, niter - burnin -1); //index of proposal
    
    u1new= ((*proposal_samples)[newcol]).col(proposal_idx); //proposed value
    u1new.shed_row(newcol); //remove diagonal entry

    arma::mat invOmega_newcol= get_invOmega_j(Omegaini, newcol); //to do: use low-rank update

    double ssnew= arma::as_scalar(u1new.t() * invOmega_newcol * u1new);
    dpostnew= bnew[newcol] * ssnew; //log-posterior for proposed value
    dpropnew= propdens->at(newcol, proposal_idx);  //log-proposal for proposed value

    dpropold= dpropini[newcol]; //log-proposal for current value
    dpostold= bnew[newcol] * arma::as_scalar((u1old[newcol]).t() * invOmega_newcol * u1old[newcol]); //log-posterior for current value

    ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);
    if ((ppnew > 1) | (runif() < ppnew)) { //if update is accepted

      if (i >= burnin) number_accept[newcol]= number_accept[newcol] + 1;
      u1old[newcol]= u1new;
      dpropini[newcol]= dpropnew;
      diagnew= rgammaC(a, -bnew[newcol]) + ssnew; //sample diagonal value from full conditional posterior

      //Update Omegaini
      spmat_rowcol2zero(Omegaini, newcol); //Set row and colum newcol to 0

      //Copy (u1new, diagnew) to Omegaini
      Omegaini->at(newcol, newcol)= diagnew;
      arma::sp_mat::iterator it;
      for (it= u1new.begin(); it != u1new.end(); ++it) {
        k= it.row();
        if (k >= newcol) k++;
        Omegaini->at(k, newcol)= Omegaini->at(newcol, k)= (*it);
      }


    }

    //Copy from Omegaini into postSample[iter,]
    if (i >= burnin) {
      spmatsym_save2flat(postSample, Omegaini, iter);
      iter++;
    }


  }

  for (i=0; i<p; i++) prop_accept->at(i)= (number_accept[i] + 0.0) / (iter + 0.0);

}
*/


/*Parallel Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n. 

Each column is sampled independently from its conditional posterior, given the current value of Omega, i.e. the algorithm targets a pseudo-posterior rather than the actual posterior

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)
  - Omegaini: initial value of Omega
  - model_logprop: log-posterior probability of the proposed model

OUTPUT
  - models: list where entry j stores indicators of Omega[,j] being zero. Formatted as a sparse matrix format (p rows and number of columns equal to the number of iterations)

*/

void GGM_Gibbs_parallel(std::vector<arma::SpMat<short>> *models, ggmObject *ggm, arma::sp_mat *Omegaini, arma::mat *model_logprop) {

  int j, p= ggm->ncol, burnin= ggm->burnin, niter= ggm->niter, p10;
//  std::vector<arma::SpMat<short>>::iterator it;

  //std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs_parallel requires the sampler to be Gibbs or birthdeath");

  if (p >10) { p10= p / 10; } else { p10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  #pragma omp parallel for default(none) shared(p, p10, use_gibbs, burnin, niter, ggm, Omegaini, models, model_logprop)
    for (j=0; j < p; j++) { //for each column

    arma::SpMat<short> *models_row= &(models->at(j));
    arma::mat invOmega_j= get_invOmega_j(Omegaini, j);
    arma::sp_mat Omegacol= Omegaini->col(j);
    arma::sp_mat *samples_row= NULL;
    arma::vec *margpp_row= NULL;
    arma::Col<int> *margppcount_row= NULL;

    if (use_gibbs) {
      GGM_Gibbs_singlecol(samples_row, models_row, margpp_row, margppcount_row, -burnin, niter-burnin-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprop);
    } else {
      GGM_birthdeath_singlecol(samples_row, models_row, margpp_row, margppcount_row, -burnin, niter-burnin-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, model_logprop);
    }

    //Rprintf("Column j=%d \n\n Proposed models \n",j); //debug
    //models_row->print("models_row"); //debug
    if (ggm->verbose) print_iterprogress(&j, &p, &p10);

  } //end for each colum

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
  - iterini, iterfi: sampled values of Omega are stored in samples[iterini:iterfi,]. Negative entries in iterini:iterfi are not stored (burnin)
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - samples: if not NULL, each column of samples will store Omegacol after applying a Gibbs update to its entries

  - models: if not NULL, stores model indicator corresponding to samples

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - model_logprob: model_logprop[,j] stores the log posterior probability of model sampled at each iteration, other columns are unchanged

IMPORTANT: if not NULL, at input samples and models should only contain zeroes

*/
void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = NULL) {

  int i, j, p= ggm->ncol;
  double mcurrent, mnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  std::default_random_engine generator((colid + 1));  //multi-thread version. Set seed according to colid
//  std::default_random_engine generator((omp_get_thread_num() + 1));  //multi-thread version. Set seed according to thread number
//  std::default_random_engine generator; //single thread version. Do not set seed
  std::uniform_real_distribution<double> unifdist(0.0,1.0);

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false);

  for (i= iterini; i <= iterfi; i++) {
    //For each entry j in colid, obtain log-marginal + log-prior for adding/removing j
    for (j= 0; j < p; j++) {

      if (j == (int) colid) continue;  //diagonal entry is always in

      if (model->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;

      ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);

      ppnew = exp(mnew - mcurrent);
      ppnew /= (1.0 + ppnew);

      if (margpp != NULL) {
        if (model->at(j,0) == 0) margpp->at(j) += ppnew; else margpp->at(j) += 1 - ppnew;
      }
      if (margppcount != NULL) margppcount->at(j) ++;

      double u= unifdist(generator);

      if (u < ppnew) { //if new model is accepted

        model_tmp_ptr= model;
        model= modelnew;
        modelnew= model_tmp_ptr;

      }

      modelnew->at(j,0)= model->at(j,0);

    } //end j for (iteration over columns)

    if (samples != NULL) {

      sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, true); //update parameter values. DO NOT USE IN A MULTI-THREAD ENVIRONMENT (WRITES TO SHARED RANDOM SEED)
      if (i >= 0) save_ggmsample_col(samples, model, &sample_diag, sample_offdiag, i, colid); //copy current model into samples[,i] as flat vector 
      delete sample_offdiag;

    } else {

      ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false); //obtain marginal likelihood, no posterior sample

    }

    if ((models != NULL) && (i >= 0)) save_ggmmodel_col(models, model, i, colid); //copy model indicator into models[,i] as flat vector

    if ((i>=0) && (model_logprob != NULL)) model_logprob->at(colid, i)= mcurrent;

  } //end i for (Gibbs iterations)

  delete model;
  delete modelnew;
  delete ms;
  
}



/*Birth-death sampling for column colid of Omega in a Gaussian graphical model

INPUT
  - iterini, iterfi: sampled values of Omega are stored in samples[iterini:iterfi,]. Negative entries in iterini:iterfi are not stored (burnin)
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters
  - Omegacol: current value of Omega[,colid]
  - invOmega_rest: inverse of Omega[-colid,-colid]

OUTPUT
  - samples: if not NULL, each column of samples will store Omegacol after applying a Gibbs update to its entries

  - models: if not NULL, stores model indicator corresponding to samples

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - model_logprob: model_logprop[,j] stores the log posterior probability of model sampled at each iteration, other columns are unchanged

IMPORTANT: if not NULL, at input samples and models should only contain zeroes

*/

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = NULL) {

  bool birth;
  int i, j, p= ggm->ncol, nbirth= ggm->nbirth, idx_update;
  double pbirth= ggm->pbirth, mcurrent, mnew, dpropcurrent, dpropnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  std::default_random_engine generator((colid + 1));  //multi-thread version. Set seed according to colid
//  std::default_random_engine generator((omp_get_thread_num() + 1));  //multi-thread version. Set seed according to thread number
//  std::default_random_engine generator; //single thread version. Do not set seed
  std::uniform_real_distribution<double> unifdist(0.0,1.0);

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false);

  for (i= iterini; i <= iterfi; i++) {

      for (j=1; j<=nbirth; j++) {  //Make nbirth birth/death updates

        //Make birth/death proposal
        arma::SpMat<short> model_colid= (*model);
        model_colid.shed_row(colid);
        rbirthdeath(&idx_update, &birth, &model_colid, pbirth);
         
        //If proposed birth/death move was possible (not all entries in model were already alive/dead)
        if (idx_update != -1) {  

          //Store birth/death proposal into modelnew_colid and modelnew
          arma::SpMat<short> modelnew_colid= model_colid;

          if (birth) modelnew_colid.at(idx_update,0)= 1; else modelnew_colid.at(idx_update,0)= 0;

          if (idx_update >= (int) colid) idx_update++;
          if (birth) modelnew->at(idx_update,0)= 1; else modelnew->at(idx_update,0)= 0;
           
          dpropnew= dbirthdeath(&modelnew_colid, &model_colid, pbirth, true);
          dpropcurrent= dbirthdeath(&model_colid, &modelnew_colid, pbirth, true);
           
          //Obtain marginal likelihood for modelnew
          ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);
           
          ppnew = exp(mnew - mcurrent + dpropcurrent - dpropnew);
          ppnew /= (1.0 + ppnew);
	   
          //Update posterior marginal inclusion probabilities
          if (margpp != NULL) {
            if (birth) margpp->at(idx_update) += ppnew; else margpp->at(idx_update) += 1 - ppnew;
          }
          if (margppcount != NULL) margppcount->at(idx_update) ++;
           
          if (unifdist(generator) < ppnew) { //if new model is accepted
           
            model_tmp_ptr= model;
            model= modelnew;
            modelnew= model_tmp_ptr;
           
          }
           
          modelnew->at(idx_update,0)= model->at(idx_update,0);

        }

    } //end for j

    if (samples != NULL) {

      sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, true); //obtain marginal likelihood & posterior sample
      if (i >= 0) save_ggmsample_col(samples, model, &sample_diag, sample_offdiag, i, colid); //copy current model into samples[,i] as flat vector
      delete sample_offdiag;

    } else {

      ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false); //obtain marginal likelihood, no posterior sample

    }

    if ((models != NULL) && (i >= 0)) save_ggmmodel_col(models, model, i, colid); //copy model indicator into models[,i] as flat vector

    if ((i>=0) && (model_logprob != NULL)) model_logprob->at(colid, i)= mcurrent;

  } //end i for (Gibbs iterations)

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


/* Save model indicator into ans[,col2save] */
void save_ggmmodel_col(arma::SpMat<short> *ans, arma::SpMat<short> *model, int col2save, unsigned int colid) {

  int j;
  unsigned int k;
  arma::SpMat<short>::iterator it_short;

  for (it_short= model->begin(), j=0; it_short != model->end(); ++it_short) {
    k= it_short.row();
    if (k != colid) { ans->at(k,col2save)= 1; j++; } else ans->at(k,col2save)= 1;
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
  //arma::vec tau = as<arma::vec>(ggm->prCoef["tau"]); //Prior is Omega_{jk} | Omega_{jk} != 0 ~ N(0, tau)
  arma::SpMat<short> model_offdiag(ggm->ncol -1, 1);

  if (npar > 0) { //if some off-diagonal entry is non-zero

    unsigned int i, j;
    //arma::vec lambda= as<arma::vec>(ggm->prCoef["lambda"]); //Prior is Omega_{jj} ~ Exp(lambda)
    double tauinv= 1.0/ ggm->prCoef_tau;
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
    double ct= ggm->prCoef_lambda + ggm->S.at(colid, colid);
    for (i=0; i<npar; i++) {
      U.at(i,i)= ct * Omegainv_model->at(i,i) + tauinv;
      for (j=i+1; j<npar; j++) U.at(i,j)= U.at(j,i)= ct * Omegainv_model->at(i,j);
    }
     
    double logdetUinv= 0;
    arma::mat Uinv(npar,npar);

    choldcinv_det(&Uinv, cholUinv, &logdetUinv, &U);
    //arma::mat Uinv= inv_sympd(U);
     
    (*m)= Uinv * s;
     
    (*logjoint) = 0.5 * (arma::as_scalar(m->t() * U * (*m)) - ((double) npar) * log( ggm->prCoef_tau ) + logdetUinv);

  } else { //if all off-diagonal entries are zero

    (*logjoint) = - ((double) npar) * log( ggm->prCoef_tau );

  }

  (*logjoint) += logprior_GGM(&model_offdiag, ggm);

}


double logprior_GGM(arma::SpMat<short> *model, ggmObject *ggm) {
  double ans;
  //string priorlabel= as<string> (ggm->prModel["priorlabel"]);

  if (ggm->priorlabel == "binomial") {
//  if (ggm->prModel["priorlabel"] == "binomial") {
    double npar= (double) model->n_nonzero; 
    double p= ggm->priorPars_p;
    //double p= as<double>(ggm->prModel["priorPars.p"]); 

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
