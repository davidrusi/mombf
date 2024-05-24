#include "modselIntegrals.h"
#include "ggm.h"
#include <cmath>

//*************************************************************************************
// sp_mat usage examples (from Rcpp gallery)
//*************************************************************************************

// Create dense matrix from sparse matrix
//  arma::mat models_dense = arma::conv_to<arma::mat>::from(*models);

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
  this->fullscan= as<bool>(samplerPars["fullscan"]);
  this->tempering= as<double>(samplerPars["tempering"]);
  this->truncratio= as<double>(samplerPars["truncratio"]);
  this->pbirth= as<double>(samplerPars["pbirth"]);
  this->pdeath= as<double>(samplerPars["pdeath"]);
  this->use_tempering= use_tempering;

  this->prob_parallel= as<double>(samplerPars["prob_parallel"]);
  CharacterVector almost_parallelR= samplerPars["almost_parallel"];
  std::string almost_parallelC = Rcpp::as<std::string>(almost_parallelR);
  std::string regression("regression"), insample("in-sample");
  this->parallel_regression= (almost_parallelC == regression);
  this->parallel_insample= (almost_parallelC == insample);

  //Set print progress iteration to true/false
  arma::vec v = as<arma::vec>(samplerPars["verbose"]);
  if (v[0] == 1) this->verbose= true; else this->verbose= false;

}


//Constructor
ggmObject::ggmObject(ggmObject *ggm) {

  //Data summaries
  this->n= ggm->n;
  this->ncol= ggm->ncol;
  this->S= ggm->S;

  //Parameters of prior distribution  
  this->prCoef_lambda= ggm->prCoef_lambda;
  this->prCoef_tau = ggm->prCoef_tau;

  this->priorlabel= ggm->priorlabel;
  this->priorPars_p= ggm->priorPars_p;

  //MCMC sampler parameters
  this->sampler= ggm->sampler;

  this->burnin= ggm->burnin;
  this->nbirth= ggm->nbirth;
  this->niter= ggm->niter;
  this->fullscan= ggm->fullscan;
  this->tempering= ggm->tempering;
  this->truncratio= ggm->truncratio;
  this->pbirth= ggm->pbirth;
  this->pdeath= ggm->pdeath;
  this->use_tempering= ggm->use_tempering;

  this->prob_parallel= ggm->prob_parallel;
  this->parallel_regression= ggm->parallel_regression;
  this->parallel_insample= ggm->parallel_insample;

  //Set print progress iteration to true/false
  this->verbose= ggm->verbose;

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
/*
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
*/




// [[Rcpp::export]]
List modelSelectionGGMC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
  bool use_tempering= false;
  double prop_accept;
  ggmObject *ggm;
  ggm= new ggmObject(&y, prCoef, prModel, samplerPars, use_tempering, true);

  int p= ggm->ncol, npars= p*(p+1)/2;

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_zigzag= (ggm->sampler == zigzag);

  GGM_CDA(&Omegaini, ggm); //update Omegaini to local posterior mode

  //Create output matrix
  arma::sp_mat postSample(npars, ggm->niter - ggm->burnin);

  arma::mat margpp = arma::zeros(p,p); 
  arma::Mat<int> margppcount;
  margppcount.zeros(p, p);
  arma::vec margppflat(p * (p+1)/2);

  //Obtain posterior samples
  if (use_gibbs || use_birthdeath) {

    GGM_Gibbs(&postSample, &margpp, &margppcount, &prop_accept, ggm, &Omegaini);

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
  ret["prop_accept"]= prop_accept;

  return ret;

}




/*Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - samples: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then samples has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - prop_accept: proportion of accepted proposals

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_Gibbs(arma::sp_mat *samples, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, jold, k, *sequence, p= ggm->ncol, iter= 0, burnin= ggm->burnin, niter= ggm->niter, niter10, number_accept, number_proposed, total_accept=0, total_proposed=0;
  double *modelini_logprob= NULL;
  arma::sp_mat::iterator it;

  //std::string sampler = Rcpp::as<std::string>(ggm->sampler());
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  arma::SpMat<short> *models= NULL;
  arma::mat *model_logprob= NULL, *invOmega_j;
  invOmega_j = new arma::mat(p-1, p-1);

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_Gibbs requires the sampler to be Gibbs or birthdeath");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  if (ggm->fullscan) { sequence= ivector(0, p-1); } else { sequence= ivector(0, 0); }

  for (i=0; i < niter; i++) {

    if (ggm->fullscan) { //update all columns in each iteration (in random order)
      for (j=0; j<p; j++) sequence[j]= j;
      samplei_wr(sequence, p, p);
    } else {             //update one column in each iteration
      sequence[0]= runifdisc(0, p-1);
    }

    for (int jj=0; jj < p; jj++) {  //for each column

      if ((!ggm->fullscan) && (jj != 0)) break;
      j= sequence[jj];

      //Obtain inverse of Omegaini[-j,-j]
      if ((i==0) && (jj==0)) {
        (*invOmega_j)= get_invOmega_j(Omegaini, j); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_j, Omegaini, &jold, &j);
      }
      jold= j;

      arma::sp_mat Omegacol= Omegaini->col(j);
      arma::sp_mat samples_row(p, 1);
      arma::vec margpp_row= arma::zeros(p);
      arma::Col<int> margppcount_row;
      margppcount_row.zeros(p);

      if (use_gibbs) {
        GGM_Gibbs_singlecol(&samples_row, models, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  invOmega_j, model_logprob, modelini_logprob); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&samples_row, models, &margpp_row, &margppcount_row, &number_accept, &number_proposed, 0, 0, (unsigned int) j, ggm, &Omegacol,  invOmega_j, model_logprob, modelini_logprob); //update row given by j
        total_accept += number_accept;
        total_proposed += number_proposed;
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

  if (use_gibbs) {
    (*prop_accept)= 1; 
  } else {
    (*prop_accept) = (total_accept + .0) / (total_proposed + .0);
  }

  if (ggm->fullscan) { free_ivector(sequence, 0, p-1); } else { free_ivector(sequence, 0, 0); }
  if (ggm->verbose) { Rcout << "\r Done\n"; }

  delete invOmega_j;

}




// [[Rcpp::export]]
List modelSelectionGGM_parallelC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
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

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_zigzag= (ggm->sampler == zigzag);

  //Allocate memory
  dpropini= dvector(0, p-1);

  //Update Omegaini to local posterior mode
  bool parallel_regression= ggm->parallel_regression;
  ggm->parallel_regression= false;
  GGM_CDA(&Omegaini, ggm);
  ggm->parallel_regression= parallel_regression;

  //Vectors where to store proposal's models & their probabilities
  std::vector<arma::SpMat<short>> proposal_models(p);
  std::vector<std::vector<double>> proposal_logprob(p);
  std::vector<std::map<string, double>> map_logprob(p);

  //Propose models, store their proposal probability into proposal_logprob
  if (use_gibbs || use_birthdeath) {

    GGM_parallel_proposal(&proposal_models, &proposal_logprob, &map_logprob, dpropini, ggm, &Omegaini);

  } else if (use_zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //MCMC using independent proposal MH to combine the chains
  ggm->use_tempering= false;
  ggm->parallel_regression= ggm->parallel_insample= false;
  arma::sp_mat postSample(npars, niter - burnin);
  if (ggm->prob_parallel < 0.9999) {
    GGM_parallel_MH_indep(&postSample, &prop_accept, &proposal_models, &proposal_logprob, dpropini, &map_logprob, ggm, &Omegaini);
  } else {
    GGM_onlyparallel_MH_indep(&postSample, &prop_accept, &proposal_models, &proposal_logprob, dpropini, ggm, &Omegaini);
  }

  //Return output
  List ret(p+3);
  //List ret(2);
  ret[0]= postSample;
  ret[1]= prop_accept;

  //Store proposals into returned list (ret)
  std::vector<arma::SpMat<short>>::iterator it;
  for (it= proposal_models.begin(), j=2; it != proposal_models.end(); ++it, j++) { 
    ret[j]= (*it);
  }

  ret[p+2]= proposal_logprob;
 
   //Free memory and return output
  delete ggm;
  free_dvector(dpropini, 0, p-1);

  return ret;


}



/*Parallel Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n. 

Each column is sampled independently from its conditional posterior, given the current value of Omega, i.e. the algorithm targets a pseudo-posterior rather than the actual posterior

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)
  - Omegaini: initial value of Omega

OUTPUT
  - models: list where entry j stores indicators of Omega[,j] being zero. Formatted as a sparse matrix format (p rows and number of columns equal to the number of iterations)
  - model_logprop: log-posterior probability of the proposed model
  - map_logprob: its jth entry contains the log-proposal probabilities, stored as a map. Useful to query the proposal probability of any given model 

*/

void GGM_parallel_proposal(std::vector<arma::SpMat<short>> *models, std::vector<std::vector<double>> *model_logprop, std::vector<std::map<string, double>> *map_logprob, double *logprop_modelini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int j, p= ggm->ncol, p10, niter_prop, burnin_prop;
  double truncratio= ggm->truncratio;
  std::string Gibbs("Gibbs"), birthdeath("birthdeath");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);

  niter_GGM_proposal(&niter_prop, &burnin_prop, &(ggm->niter), &(ggm->burnin), &p); //if not using a full scan, less iterations are needed

  if (!use_gibbs && !use_birthdeath) Rf_error("GGM_parallel_proposal requires the sampler to be Gibbs or birthdeath");

  if (p >10) { p10= p / 10; } else { p10= 1; }

  if (ggm->verbose) Rprintf(" Running parallel proposals\n");

  //#pragma omp parallel for default(none) shared(p, p10, use_gibbs, burnin, niter, ggm, Omegaini, models, model_logprop)
  for (j=0; j < p; j++) { //for each column

    arma::SpMat<short> models_row(p, niter_prop - burnin_prop);
    arma::mat invOmega_j;
    if (!ggm->parallel_regression) invOmega_j= get_invOmega_j(Omegaini, j);
    arma::sp_mat Omegacol= Omegaini->col(j);
    arma::sp_mat *samples_row= NULL;
    arma::vec *margpp_row= NULL;
    arma::Col<int> *margppcount_row= NULL;
    arma::mat modelj_logprop(1, niter_prop - burnin_prop);
    int number_accept, number_proposed;

    if (use_gibbs) {
      GGM_Gibbs_singlecol(samples_row, &models_row, margpp_row, margppcount_row, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, &modelj_logprop, logprop_modelini+j);
    } else {
      GGM_birthdeath_singlecol(samples_row, &models_row, margpp_row, margppcount_row, &number_accept, &number_proposed, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, &modelj_logprop, logprop_modelini+j);
    }

    //Obtain unique visited models and their probabilities, truncated so probability of top model / any other model <= truncratio
    unique_model_logprob(&(models->at(j)), &(model_logprop->at(j)), &(map_logprob->at(j)), &models_row, &modelj_logprop, &truncratio, logprop_modelini+j);
    if (truncratio > 0) {
      logprop_modelini[j]= max_xy(logprop_modelini[j], (model_logprop->at(j))[0] - log(truncratio)); //truncate probabilities according to maxratio
    }

    if (ggm->verbose) print_iterprogress(&j, &p, &p10);

  } //end for each colum

  //Alternative version where memory is not shared across threads (but requires making local copies of ggm, models and model_logprop)
  /*
  #pragma omp parallel for default(none) shared(p, p10, use_gibbs, burnin, niter, ggm, Omegaini, models, model_logprop)
  for (j=0; j < p; j++) { //for each column

    ggmObject *ggm_local;
    ggm_local= new ggmObject(ggm);
    arma::SpMat<short> models_row(p, niter_prop - burnin_prop);
    arma::mat invOmega_j;
    if (!ggm->parallel_regression) invOmega_j= get_invOmega_j(Omegaini, j);
    arma::sp_mat Omegacol= Omegaini->col(j);
    arma::sp_mat *samples_row= NULL;
    arma::vec *margpp_row= NULL;
    arma::Col<int> *margppcount_row= NULL;
    arma::mat modelj_logprop(1, niter_prop - burnin_prop);
    double logpropini;
    int number_accept, number_proposed;

    if (use_gibbs) {
      GGM_Gibbs_singlecol(samples_row, &models_row, margpp_row, margppcount_row, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm_local, &Omegacol,  &invOmega_j, &modelj_logprop, &logpropini);
    } else {
      GGM_birthdeath_singlecol(samples_row, &models_row, margpp_row, margppcount_row, &number_accept, &number_proposed, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm_local, &Omegacol,  &invOmega_j, &modelj_logprop, &logpropini);
    }

    delete ggm_local;

    //Obtain unique visited models and their probabilities, truncated so probability of top model / any other model <= truncratio
    arma::SpMat<short> models_row_sorted;
    arma::mat modelj_logprop_sorted;
    unique_model_logprob(&(models->at(j)), &(model_logprop->at(j)), &(map_logprob->at(j)), &models_row, &modelj_logprop, &truncratio, logprop_modelini+j);
    if (truncratio > 0) {
      logprop_modelini[j]= max_xy(logprop_modelini[j], (model_logprop->at(j))[0] - log(truncratio)); //truncate probabilities according to maxratio
    }

    if (ggm->verbose) print_iterprogress(&j, &p, &p10);

    #pragma omp critical 
    {
      (*models)[j]= models_row;
      model_logprop->row(j)= modelj_logprop;
      logprop_modelini[j]= logpropini;
    }

    //int nthreads= omp_get_num_threads(); //debug

  } //end for each colum
*/

  if (ggm->verbose) { Rcout << "\r Done\n"; }

}



/*Coordinate ascent algorithm for the precision matrix Omega in a Gaussian graphical model

  The best model is found by a coordinate descent algorithm on each column Omega[,colid], given Omega[-colid,-colid], and iterating over columns.
  The non-zero entries in Omega are set to their posterior mode under the chosen model

INPUT
  - colid: index of the column to be updated (starts at 0)
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters

OUTPUT

  - Omegacol: at input, the current value of Omega. At output, a local posterior mode of Omega

*/
void GGM_CDA(arma::sp_mat *Omega, ggmObject *ggm) {

  const int maxit=10;
  int i, j, k, colid, oldcolid, p= ggm->ncol, changes, changes_colid;
  double mcurrent, mnew, sample_diag;
  arma::mat *sample_offdiag= NULL, *invOmega_rest;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *modelopt, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun= &GGMrow_marg;

  invOmega_rest = new arma::mat(p-1, p-1);

  k= changes= 1;
  
  while ((changes >=1) && (k <= maxit)) {

    changes= 0;
    for (colid= 0; colid < p; colid++) {

      if ((k==1) && (colid==0)) {
        (*invOmega_rest)= get_invOmega_j(Omega, colid); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_rest, Omega, &oldcolid, &colid);
      }
      oldcolid= colid;
      //arma::mat invOmega_rest= get_invOmega_j(Omega, colid);
      arma::sp_mat Omegacol= Omega->col(colid);

      //Initialize model and modelnew
      modelopt= new arma::SpMat<short>(p, 1);
      modelnew= new arma::SpMat<short>(p, 1);
      for (it= Omegacol.begin(); it != Omegacol.end(); ++it) modelopt->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

      //Obtain log-marginal + log-prior for current model
      ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

      ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, modelopt, false);

      i= changes_colid= 1;
      while ((changes_colid > 0) && (i <= maxit)) {  // (update model for column colid)

        changes_colid= 0;

        //For each entry j in colid, set entries to zero/non-zero if it improves the posterior model probability
        for (j= 0; j < p; j++) {

          if (j == (int) colid) continue;  //diagonal entry is always in

          if (modelopt->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;

          ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);

          if (mnew > mcurrent) { //if new model has higher posterior probability

            model_tmp_ptr= modelopt;
            modelopt= modelnew;
            modelnew= model_tmp_ptr;
            mcurrent = mnew;
            changes++;
            changes_colid++;

          }

          modelnew->at(j,0)= modelopt->at(j,0);

        } //end j for (iteration over columns)

        i++;

      } //end i for (update model for column colid)

      sample_offdiag= new arma::mat(modelopt->n_nonzero - 1, 1);

      ms->getMode(&mcurrent, sample_offdiag, &sample_diag, modelopt); //obtain posterior mode for parameters given model
      update_Omega(Omega, &colid, &sample_diag, modelopt, sample_offdiag); //store posterior mode into Omega

      delete sample_offdiag;
      delete modelopt;
      delete modelnew;
      delete ms;
    }

    k++;

  }

  delete invOmega_rest;

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

  - model_logprob: If model_logprop has >1 rows, model_logprop[colid,j] stores the log posterior probability of model sampled at iteration j, other rows are unchanged. 
                   If model_logprop only has 1 row, then log posterior prob are stored in logprop[0,j]

  - modelini_logprob: log posterior probability of initial model

IMPORTANT: if not NULL, at input samples and models should only contain zeroes

*/
void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = NULL, double *modelini_logprob = NULL) {

  int i, j, p= ggm->ncol, col2save;
  double mcurrent, mnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun;

  if (ggm->parallel_regression) {
    marfun= &GGMrow_marg_regression;
  } else {
    marfun= &GGMrow_marg;
  }

  if ((model_logprob != NULL) && (model_logprob->n_rows > 1)) col2save= colid; else col2save= 0;

  //Set random number generators
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

  if (modelini_logprob != NULL) (*modelini_logprob)= mcurrent;

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

    if ((i>=0) && (model_logprob != NULL)) model_logprob->at(col2save, i)= mcurrent;

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

  - model_logprob: If model_logprop has >1 rows, model_logprop[colid,j] stores the log posterior probability of model sampled at iteration j, other rows are unchanged. 
                   If model_logprop only has 1 row, then log posterior prob are stored in logprop[0,j]

  - number_accept: number of accepted proposals

  - number_proposals: number of proposals

  - modelini_logprob: log posterior probability of initial model

IMPORTANT: if not NULL, at input samples and models should only contain zeroes

*/

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::vec *margpp, arma::Col<int> *margppcount, int *number_accept, int *number_proposed, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = NULL, double *modelini_logprob = NULL) {

  int i, j, p= ggm->ncol, nbirth= ggm->nbirth, index_birth, index_death, movetype, col2save, colid_int= (int) colid;
  double pbirth= ggm->pbirth, pdeath= ggm->pdeath, mcurrent, mnew, dpropcurrent, dpropnew, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun;

  if (ggm->parallel_regression) {
    marfun= &GGMrow_marg_regression;
  } else {
    marfun= &GGMrow_marg;
  }

  (*number_accept)= (*number_proposed)= 0;
  if ((model_logprob != NULL) && (model_logprob->n_rows > 1)) col2save= colid; else col2save= 0;

  //Set random number generators
  std::default_random_engine generator((colid + 1));  //multi-thread version. Set seed according to colid
//  std::default_random_engine generator; //single thread version. Do not set seed
  std::uniform_real_distribution<double> unifdist(0.0,1.0);

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  ms->getJoint(&mcurrent, sample_offdiag, &sample_diag, model, false);

  if (modelini_logprob != NULL) (*modelini_logprob)= mcurrent;

  for (i= iterini; i <= iterfi; i++) {

      for (j=1; j<=nbirth; j++) {  //Make nbirth birth/death updates

      //GGM_birthdeath_proposal(modelnew, &idx_update, &birth, &dpropnew, &dpropcurrent, model, &colid_int, &pbirth, false);
      GGM_birthdeathswap_proposal(modelnew, &index_birth, &index_death, &movetype, &dpropnew, &dpropcurrent, model, &colid_int, &pbirth, &pdeath, false);
      (*number_proposed)++;

      if ((index_birth != -1) || (index_death != -1)) {          
          //Obtain marginal likelihood for modelnew
          ms->getJoint(&mnew, sample_offdiag, &sample_diag, modelnew, false);
           
          ppnew = exp(mnew - mcurrent + dpropcurrent - dpropnew);
          ppnew /= (1.0 + ppnew);
	            
          if (unifdist(generator) < ppnew) { //if new model is accepted
           
            model_tmp_ptr= model;
            model= modelnew;
            modelnew= model_tmp_ptr;
            (*number_accept)++;
           
          }
           
          if (index_birth != -1) modelnew->at(index_birth,0)= model->at(index_birth,0);
          if (index_death != -1) modelnew->at(index_death,0)= model->at(index_death,0);

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

    if ((i>=0) && (model_logprob != NULL)) model_logprob->at(col2save, i)= mcurrent;

  } //end i for (Gibbs iterations)

  delete model;
  delete modelnew;
  delete ms;
  
}

/* Propose GGM birth-death model update. Return proposal density for the new state and for the current state

  INPUT

  - model: current model, a p x 1 matrix indicating the non-zero parameters
  - colid: column of the GGM precision matrix Omega currently being updated
  - pbirth: probability of a birth move
  - setmodelnew: if false, modelnew is assumed to be equal to model at input, so that only the entry being borth/dead is updated. Else, all entries in modelnew are set

  OUTPUT

  - modelnew: model proposed by the birth-death move
  - idx_update: index of parameter being born/dead. That is, modelnew[idx_update,0] is different from model[idx_update,0], all other rows are equal. 
                If proposed move was impossible (e.g. a birth was attempted but all parameters were already alive) then idx_update == -1 is returned
  - birth: true if a birth move was proposed (i.e. modelnew[idx_update,0]==1), false otherwise
  - dpropnew: log of the proposal probability of modelnew
  - dpropcurrent: log of the proposal probability of model

*/
void GGM_birthdeath_proposal(arma::SpMat<short> *modelnew, int *idx_update, bool *birth, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, bool setmodelnew) {

  if (setmodelnew) (*modelnew)= *model;

  //Make birth/death proposal
  arma::SpMat<short> model_colid= (*model);
  model_colid.shed_row(*colid);
  rbirthdeath(idx_update, birth, &model_colid, *pbirth);
         
  //If proposed birth/death move was possible (not all entries in model were already alive/dead)
  if (*idx_update != -1) {  

    //Store birth/death proposal into modelnew_colid and modelnew
    arma::SpMat<short> modelnew_colid= model_colid;

    if (*birth) modelnew_colid.at(*idx_update,0)= 1; else modelnew_colid.at(*idx_update,0)= 0;

    if (*idx_update >= (int) *colid) (*idx_update)++;
    if (*birth) modelnew->at(*idx_update,0)= 1; else modelnew->at(*idx_update,0)= 0;
           
    (*dpropnew)= dbirthdeath(&modelnew_colid, &model_colid, *pbirth, true);
    (*dpropcurrent)= dbirthdeath(&model_colid, &modelnew_colid, *pbirth, true);

  } else {

    (*dpropnew)= (*dpropcurrent) -INFINITY;

  }

}


/* Propose GGM birth-death-swap model update. Return proposal density for the new state and for the current state

  INPUT

  - model: current model, a p x 1 matrix indicating the non-zero parameters
  - colid: column of the GGM precision matrix Omega currently being updated
  - pbirth: probability of a birth move
  - pdeath: probability of a death move (swap move has probability 1 - pbirth - pdeath)
  - setmodelnew: if false, modelnew is assumed to be equal to model at input, so that only the entries being borth/dead are updated. Else, all entries in modelnew are set

  OUTPUT

  - modelnew: model proposed by the birth-death move
  - index_birth: index of parameter being born. If no entries were born, then index_birth= -1
  - index_death: index of parameter being killed. If no entries were killed, then index_death= -1
  - movetype: 1 if a birth move was chosen, 2 if a death move was chosen, 3 if a swap move was chosen
  - dpropnew: log of the proposal probability of modelnew
  - dpropcurrent: log of the proposal probability of model

*/
void GGM_birthdeathswap_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, double *pdeath, bool setmodelnew) {

  if (setmodelnew) (*modelnew)= *model;

  //Make birth/death proposal
  arma::SpMat<short> model_colid= (*model);
  model_colid.shed_row(*colid);
  rbirthdeathswap(index_birth, index_death, movetype, &model_colid, *pbirth, *pdeath);
         
  //If proposed  move was possible (not all entries in model were already alive/dead)
  if ((*index_birth != -1) || (*index_death != -1)) {  

    //Store birth/death proposal into modelnew_colid and modelnew
    arma::SpMat<short> modelnew_colid= model_colid;

    if (*index_birth != -1) {
      modelnew_colid.at(*index_birth,0)= 1; 
      if (*index_birth >= (int) *colid) (*index_birth)++;
      modelnew->at(*index_birth,0)= 1;
    }

    if (*index_death != -1) {
      modelnew_colid.at(*index_death,0)= 0;
      if (*index_death >= (int) *colid) (*index_death)++;
      modelnew->at(*index_death,0)= 0;
    }

    (*dpropnew)= dbirthdeathswap(&modelnew_colid, &model_colid, *pbirth, *pdeath, true);
    (*dpropcurrent)= dbirthdeathswap(&model_colid, &modelnew_colid, *pbirth, *pdeath, true);

  } else {

    (*dpropnew)= (*dpropcurrent) -INFINITY;

  }

}


/* Determine number of iterations to use in GGM parallel proposals 

  The number of expected draws from each parallel proposal is m= niter / p, and its SD is s= sqrt(niter (1/p) (1 - 1/p))
  This function returns m + 5 s, in order to propose 10 times the expected number of samples.
  However, in cases where 10 * niter / p is too small (<1000) we return 1000, and if it's too large (>niter) we return niter.
  Also, if niter was small to start with (niter <= 1000), it is returned unaltered.

  Input: niter, burnin, p
  Output: niter_prop, burnin_prop. The ratio of burnin_prop/niter_prop is the same as burnin/niter

*/
void niter_GGM_proposal(int *niter_prop, int *burnin_prop, int *niter, int *burnin, int *p) {

  if (*niter <= 1000) {

    (*niter_prop)= *niter;

  } else {

    double m= (*niter + .0) / (*p + .0), s= sqrt(m * (1 - 1/(*p + 0)));
    int niter_reduced= ceil(m + 5 * s);
    if (niter_reduced < 1000) {
      (*niter_prop)= *niter;
    } else if (niter_reduced > *niter) {
      (*niter_prop)= *niter;
    } else {
      (*niter_prop)= niter_reduced;
    }
    
  }

  (*burnin_prop)= (int) (*burnin + .0) / (*niter + .0) * (*niter_prop + .0);

}


/* Given possibly repeated models and their unnormalized probabilities, return the unique models and their normalized probabilities, sorted in decreasing probability
   and truncated so that top model / any model's probability <= maxratio 

INPUT
  - models: matrix specifying a model in each column, where rows correspond to model parameters being included/excluded
  - models_logprob: unnormalized log-probabilities, formatted as a matrix with 1 row and as many columns as models
  - maxratio: see above. If maxratio < 0, it is ignored 

OUTPUT
  - uniquemodels: same as models, but repeated models have been removed and models are sorted in decreasing probability
  - unique_logprob: normalized probabilities corresponding to uniquemodels
  - map_logprob: stores both uniquemodels and uniquemodels_logprob in a map (unsorted)

INPUT/OUTPUT

  - logprob_ini: if not null, at input it is a single unnormalized log-probability and at output the normalized log-probability is returned

*/
void unique_model_logprob(arma::SpMat<short> *uniquemodels, std::vector<double> *unique_logprob, std::map<string, double> *map_logprob, arma::SpMat<short> *models, arma::mat *models_logprob, double *maxratio, double *logprobini) {

  int i, idx, nmodels= models->n_cols, nunique, nvars= models->n_rows;
  double maxprob= -INFINITY, logmaxratio;
  std::vector<int> modelids_column(0);
  std::vector<double> modelids_logprob(0);
  std::vector<std::string> svec(0);
  char *zerochar;

  zerochar = (char *) calloc(nvars+1, sizeof(char));
  for (i=0; i<nvars; i++) zerochar[i]= '0';

  map_logprob->clear(); //empty the map

  //Store identifier of unique models in modelids, their column in modelids_column, their probabilities in modelids_logprob
  for (i=0; i<nmodels; i++) {
    //Store modelid into string s
    arma::SpMat<short> mymodel= models->col(i);
    std::string s= getModelid(&mymodel, zerochar);

    //If modelid does not exist, add to indexes
    if (map_logprob->count(s) == 0) {
      (*map_logprob)[s]= models_logprob->at(0,i);
      svec.push_back(s);
      modelids_column.push_back(i);
      modelids_logprob.push_back(models_logprob->at(0,i));
      maxprob= max_xy(maxprob, models_logprob->at(0,i));
    }

  }

  //Copy posterior probabilities of unique models onto modelids_logprob, and sort them decreasingly
  nunique= modelids_logprob.size();
  std::vector<int> sortedIndexes = sorted_indexes(modelids_logprob, true);

  //Copy unique models onto uniquemodels, unique_logprob
  arma::SpMat<short> mymodels(nvars, nunique);
  (*uniquemodels)= mymodels;
  unique_logprob->resize(nunique);
  if (*maxratio > 0) logmaxratio= log(*maxratio);
  double sumprob= 0, logsumprob;
  for (i=0; i<nunique; i++) {
    idx= sortedIndexes[i];
    uniquemodels->col(i)= models->col(modelids_column[idx]);
    if (*maxratio <= 0) {
      (*unique_logprob)[i]= modelids_logprob[idx] - maxprob;
    } else {
      (*unique_logprob)[i]= max_xy(modelids_logprob[idx] - maxprob, - logmaxratio); //truncate probabilities according to maxratio
    }
    sumprob += exp((*unique_logprob)[i]);
  }
  logsumprob= log(sumprob);

  //Normalize probabilities
  for (i=0; i<nunique; i++) {
    (*unique_logprob)[i] -= logsumprob;
    (*map_logprob)[svec[sortedIndexes[i]]]= (*unique_logprob)[i];
  }

  if (logprobini != NULL) {
    if (*maxratio <= 0) {
      (*logprobini)= *logprobini - maxprob - logsumprob;
    } else {
      (*logprobini)= max_xy(*logprobini - maxprob, - logmaxratio) - logsumprob; //truncate probabilities according to maxratio
    }
  }

  free((char  *) zerochar);

}

/* Return model identifier as a string with 0's and 1's

  INPUT

  - model: sparse matrix with 1 column, where non-zero entries are converted into a 1 in the output string
  - zerochar: char of length equal to model->n_rows, storing all zeroes, e.g. "00000". That is, before calling getModelid one has run something like

      zerochar = (char *) calloc(model.n_rows + 1, sizeof(char));
      for (i=0; i<model.n_rows; i++) zerochar[i]= '0';

  OUTPUT: a string indicating the non-zero rows in model, e.g. "10100"

*/
std::string getModelid(arma::SpMat<short> *model, char *zerochar) {
  if (model->n_cols !=1) Rf_error("In getModelid, argument model must have 1 column");
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) zerochar[it.row()]= '1'; //Set zerochar to current model
  std::string s (zerochar);
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) zerochar[it.row()]= '0'; //Return zerochar to its original empty model status
  return s;
}


/* Obtain posterior samples for Omega via Metropolis-Hastings, using the parallel proposals

INPUT

  - proposal_models: vector where jth entry contains the proposal models for Omega[,j] (precision matrix entries being zero/non-zero), sorted in decreasing probability
  - proposal_logprob: vector where jth entry contains log-proposal densities for proposal_models[j]
  - proposal_modelids: string identifiers of the models in proposal_models (unsorted). Useful to quickly check if a model is in proposal_models
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT

  - postSample: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then postSample has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.
  - prop_accept: proportion of accepted proposals

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/
void GGM_parallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, std::vector<std::map<string, double>> *map_logprob, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, i10, iter, newcol, oldcol, *sequence, index_birth, index_death, movetype, number_accept= 0, proposal_idx, burnin= ggm->burnin, niter= burnin + postSample->n_cols, p= ggm->ncol;
  double dpostnew, dpostold, dpropnew, dpropold, ppnew, sample_diag;
  double prob_parallel= ggm->prob_parallel;
  bool parallel_draw;
  char *zerochar;
  zerochar = (char *) calloc(p+1, sizeof(char));
  for (i=0; i<p; i++) zerochar[i]= '0';
  std::map<string, double>::iterator itmap;

  arma::mat *sample_offdiag= NULL, *invOmega_newcol;
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::SpMat<short>> modelold(p);
  arma::SpMat<short> modelnew;
  arma::SpMat<short>::iterator it;
  pt2GGM_rowmarg marfun= &GGMrow_marg;
  std::vector<double*> cdf_proposal(p);
  std::vector<int> nproposal(p);
  std::string s;

  invOmega_newcol = new arma::mat(p-1, p-1);

  if (ggm->verbose) Rprintf(" Running MH-within-Gibbs\n");
  if (niter >10) { i10= niter / 10; } else { i10= 1; }

  //Initialize MCMC state to Omegaini
  for (newcol=0; newcol < p; newcol++) {

    //Store current model into modelold
    arma::SpMat<short> mymodel(p,1);
    arma::sp_mat mycolumn= Omegaini->col(newcol);
    for (arma::sp_mat::iterator it= mycolumn.begin(); it != mycolumn.end(); ++it) { mymodel.at(it.row(), 0)= 1; }
    modelold[newcol]= mymodel;

    //Pre-compute parallel proposal's cdf, so we can use rdisc_pcum later
    nproposal[newcol]= (proposal_models->at(newcol)).n_cols;
    cdf_proposal[newcol]= dvector(0, nproposal[newcol]);
    double pcum=0, maxlogprob= (*proposal_logprob)[newcol][0];
    for (i=0; i<nproposal[newcol]; i++) {
      cdf_proposal[newcol][i]= exp((*proposal_logprob)[newcol][i] - maxlogprob);
      pcum += cdf_proposal[newcol][i];
    }
    for (i=0; i<nproposal[newcol]; i++) cdf_proposal[newcol][i] /= pcum;
    cumsum(cdf_proposal[newcol], cdf_proposal[newcol], &(nproposal[newcol])); 

  }

  if (ggm->fullscan) { sequence= ivector(0, p-1); } else { sequence= ivector(0, 0); }

  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {

    if (ggm->fullscan) { //update all columns in each iteration (in random order)
      sequence= ivector(0, p-1);
      for (j=0; j<p; j++) sequence[j]= j;
      samplei_wr(sequence, p, p);
    } else {             //update one column in each iteration
      sequence[0]= runifdisc(0, p-1);
    }
    for (j=0; j<p; j++) { //iterate over columns in random order

      newcol= sequence[j];
      if ((!ggm->fullscan) && (j != 0)) break;

      parallel_draw= (runifC() < prob_parallel); //choose proposal type

      //If current model has zero probability under the parallel proposal, use the birth-death proposal
      //if (parallel_draw) {
      //  s= getModelid(&(modelold[newcol]), zerochar);
      //  itmap= (map_logprob->at(newcol)).find(s);
      //  parallel_draw= (itmap != (map_logprob->at(newcol)).end());
      //}
      
      if (parallel_draw && (dpropini[newcol] != -INFINITY)) { //use parallel proposal

        dpropold= dpropini[newcol]; //log-proposal for current model
        proposal_idx= rdisc_pcum(cdf_proposal[newcol], nproposal[newcol]);  //index of proposal
        dpropnew= (proposal_logprob->at(newcol))[proposal_idx]; //log-proposal for new value
        modelnew= ((*proposal_models)[newcol]).col(proposal_idx); //proposed value

      } else { //use birth-death proposal

        //GGM_birthdeath_proposal(&modelnew, &idx_update, &birth, &dpropnew, &dpropold, &(modelold[newcol]), &newcol, &(ggm->pbirth), true);
        GGM_birthdeathswap_proposal(&modelnew, &index_birth, &index_death, &movetype, &dpropnew, &dpropold, &(modelold[newcol]), &newcol, &(ggm->pbirth), &(ggm->pdeath), true);

      }

      //Obtain inverse of Omegaini[-newcol,-newcol]
      if ((i==0) && (j==0)) {
        (*invOmega_newcol)= get_invOmega_j(Omegaini, newcol); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_newcol, Omegaini, &oldcol, &newcol);
        //Uncomment the lines below to keep track of numerical error in rank 1 calculation (none found in the tried examples)
        //arma::mat invOmega_costly= get_invOmega_j(Omegaini, newcol); 
        //double approxerror= norm(*invOmega_newcol - invOmega_costly, "fro"); 
        //Rprintf("i=%d, j=%d, approx error= %f\n", i, j, approxerror); 
      }
      oldcol= newcol;

      modselIntegrals_GGM *ms;
      ms= new modselIntegrals_GGM(marfun, ggm, newcol, invOmega_newcol);
      ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &modelnew, false); //log-posterior for proposed model
      ms->getJoint(&dpostold, sample_offdiag, &sample_diag, &(modelold[newcol]), false); //log-posterior for current model

      ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);

      //if ((newcol==29)) { //debug
      //  bool proptype= parallel_draw && (dpropini[newcol] != -INFINITY);
      //  Rprintf("i=%d, newcol=%d, dpostold=%f, dpostnew=%f, dpropold=%f, dpropnew=%f, parallel_update=%d \n", i, newcol, dpostold, dpostnew, dpropold, dpropnew, proptype);
      //  (modelold[newcol]).print("modelold");
      //  modelnew.print("modelnew");
      //}

      if ((ppnew >= 1) | (runifC() < ppnew)) { //if update is accepted

        if (i >= burnin) number_accept++;
        modelold[newcol]= modelnew;
        //Store log parallel proposal density of updated model
        s= getModelid(&modelnew, zerochar);
        itmap= (map_logprob->at(newcol)).find(s);
        if (itmap != (map_logprob->at(newcol)).end()) {
          dpropini[newcol]= itmap->second;
        } else {
          dpropini[newcol]= -INFINITY;
        }

      }

      //Sample Omega[,newcol] given model
      sample_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &(modelold[newcol]), true); //sample Omega[,newcol]
      update_Omega(Omegaini, &newcol, &sample_diag, &(modelold[newcol]), sample_offdiag); //copy sampled values into Omegaini[,newcol]
      delete sample_offdiag;
      delete ms;

    } //end j for (iterate over columns)

    //Copy from Omegaini into postSample[,iter]
    if (i >= burnin) {
      spmatsym_save2flat(postSample, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &i10);

  } //end i for

  if (ggm->fullscan) { free_ivector(sequence, 0, p-1); } else { free_ivector(sequence, 0, 0); }

  (*prop_accept)= (number_accept + 0.0) / (iter * p + 0.0);

  //Free memory
  for (newcol=0; newcol < p; newcol++) free_dvector(cdf_proposal[newcol], 0, nproposal[newcol]);
  free((char  *) zerochar);
  delete invOmega_newcol;

}


/* CODE FOR MIXTURE PROPOSAL, EVALUATING DENSITY OF WHOLE MIXTURE

      if (runifC() < prob_parallel) { //use parallel proposal

        proposal_idx= rdisc_pcum(cdf_proposal[newcol], nproposal[newcol]);  //index of proposal
        dpropnew_parallel= (proposal_logprob->at(newcol))[proposal_idx]; //log-proposal under parallel proposal
        modelnew= ((*proposal_models)[newcol]).col(proposal_idx); //proposed value
        if (!onlyparallel) {
          //log-proposal under birth-death mode
          is_birthdeath= checkNonZeroDiff(&modelnew, &(modelold[newcol]), 1) && (modelnew.n_nonzero != (modelold[newcol]).n_nonzero);
          if (is_birthdeath) {
            dpropnew_birthdeath= dbirthdeath(&modelnew, &(modelold[newcol]), ggm->pbirth, true); 
            dpropold_birthdeath= dbirthdeath(&(modelold[newcol]), &modelnew, ggm->pbirth, true);
          } else {
            dpropnew_birthdeath= dpropold_birthdeath= -INFINITY;
          }
        }
        //        dpropold= dpropold_parallel;
        //        dpropnew= dpropnew_parallel;

      } else { //use birth-death proposal

        GGM_birthdeath_proposal(&modelnew, &idx_update, &birth, &dpropnew_birthdeath, &dpropold_birthdeath, &(modelold[newcol]), &newcol, &(ggm->pbirth), true);

        //Evaluate parallel proposal density for modelnew
        s= getModelid(&modelnew, zerochar);
        itmap= (map_logprob->at(newcol)).find(s);
        if (itmap != (map_logprob->at(newcol)).end()) { dpropnew_parallel= itmap->second; } else { dpropnew_parallel= -INFINITY; }

      }

      dpropold_parallel= dpropini_parallel[newcol]; //log-proposal under parallel proposal

      //Evaluate proposal mixture density
      max_dprop= max_xy(dpropold_parallel, dpropold_birthdeath);
      dpropold= max_dprop + log(prob_parallel * exp(dpropold_parallel - max_dprop) + prob_birthdeath * exp(dpropold_birthdeath - max_dprop));

      max_dprop= max_xy(dpropnew_parallel, dpropnew_birthdeath);
      dpropnew= max_dprop + log(prob_parallel * exp(dpropnew_parallel - max_dprop) + prob_birthdeath * exp(dpropnew_birthdeath - max_dprop));

*/

void GGM_onlyparallel_MH_indep(arma::sp_mat *postSample, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, i10, iter, newcol, oldcol, *sequence, number_accept= 0, proposal_idx, burnin= ggm->burnin, niter= burnin + postSample->n_cols, p= ggm->ncol;
  double dpostnew, dpostold, dpropnew, dpropold, ppnew, sample_diag;
  arma::mat *sample_offdiag= NULL, *invOmega_newcol;
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::SpMat<short>> modelold(p);
  arma::SpMat<short> modelnew;
  pt2GGM_rowmarg marfun= &GGMrow_marg;
  std::vector<double*> cdf_proposal(p);
  std::vector<int> nproposal(p);

  invOmega_newcol = new arma::mat(p-1, p-1);

  if (ggm->verbose) Rprintf(" Running MH-within-Gibbs\n");
  if (niter >10) { i10= niter / 10; } else { i10= 1; }

  //Initialize MCMC state to Omegaini
  for (newcol=0; newcol < p; newcol++) {

    arma::SpMat<short> mymodel(p,1);
    arma::sp_mat mycolumn= Omegaini->col(newcol);
    for (arma::sp_mat::iterator it= mycolumn.begin(); it != mycolumn.end(); ++it) { mymodel.at(it.row(), 0)= 1; }
    modelold[newcol]= mymodel;

    //Pre-compute proposal cdf, so we can use rdisc_pcum later
    nproposal[newcol]= (proposal_models->at(newcol)).n_cols;
    cdf_proposal[newcol]= dvector(0, nproposal[newcol]);
    double pcum=0, maxlogprob= (*proposal_logprob)[newcol][0];
    for (i=0; i<nproposal[newcol]; i++) {
      cdf_proposal[newcol][i]= exp((*proposal_logprob)[newcol][i] - maxlogprob);
      pcum += cdf_proposal[newcol][i];
    }
    for (i=0; i<nproposal[newcol]; i++) cdf_proposal[newcol][i] /= pcum;
    cumsum(cdf_proposal[newcol], cdf_proposal[newcol], &(nproposal[newcol])); 

  }

  if (ggm->fullscan) { sequence= ivector(0, p-1); } else { sequence= ivector(0, 0); }

  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {

    if (ggm->fullscan) { //update all columns in each iteration (in random order)
      sequence= ivector(0, p-1);
      for (j=0; j<p; j++) sequence[j]= j;
      samplei_wr(sequence, p, p);
    } else {             //update one column in each iteration
      sequence[0]= runifdisc(0, p-1);
    }
    for (j=0; j<p; j++) { //iterate over columns in random order

      newcol= sequence[j];
      if ((!ggm->fullscan) && (j != 0)) break;

      proposal_idx= rdisc_pcum(cdf_proposal[newcol], nproposal[newcol]);  //index of proposal
      dpropnew= (proposal_logprob->at(newcol))[proposal_idx]; //log-proposal for proposed model
      modelnew= ((*proposal_models)[newcol]).col(proposal_idx); //proposed value

      //Obtain inverse of Omegaini[-newcol,-newcol]
      if ((i==0) && (j==0)) {
        (*invOmega_newcol)= get_invOmega_j(Omegaini, newcol); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_newcol, Omegaini, &oldcol, &newcol);
      }
      oldcol= newcol;

      modselIntegrals_GGM *ms;
      ms= new modselIntegrals_GGM(marfun, ggm, newcol, invOmega_newcol);
      ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &modelnew, false); //log-posterior for proposed model
      ms->getJoint(&dpostold, sample_offdiag, &sample_diag, &(modelold[newcol]), false); //log-posterior for current model

      dpropold= dpropini[newcol]; //log-proposal for current model

      ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);

      if ((ppnew >= 1) | (runifC() < ppnew)) { //if update is accepted

        if (i >= burnin) number_accept++;
        dpropini[newcol]= dpropnew;
        modelold[newcol]= modelnew;

      }

      //Sample Omega[,newcol] given model
      sample_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      ms->getJoint(&dpostnew, sample_offdiag, &sample_diag, &(modelold[newcol]), true); //sample Omega[,newcol]
      update_Omega(Omegaini, &newcol, &sample_diag, &(modelold[newcol]), sample_offdiag); //copy sampled values into Omegaini[,newcol]
      delete sample_offdiag;
      delete ms;

    } //end j for (iterate over columns)

    //Copy from Omegaini into postSample[,iter]
    if (i >= burnin) {
      spmatsym_save2flat(postSample, Omegaini, iter);
      iter++;
    }

    if (ggm->verbose) print_iterprogress(&i, &niter, &i10);

  } //end i for

  if (ggm->fullscan) { free_ivector(sequence, 0, p-1); } else { free_ivector(sequence, 0, 0); }

  (*prop_accept)= (number_accept + 0.0) / (iter * p + 0.0);

  for (newcol=0; newcol < p; newcol++) free_dvector(cdf_proposal[newcol], 0, nproposal[newcol]);

  delete invOmega_newcol;

}



/*Update symmetric Omegaini[,newcol] by storing Omegaini[newcol,newcol]= sample_diag and storing sample_offdiag into the off-diagonal entries indicated by modelnew
  Input
  - newcol: column to be updated
  - sample_diag: new diagonal entry to be stored
  - modelnew: sparse column vector (dimension nrow(Omegaini) x 1) indicating what off-diagonal entries in Omegaini[,newcol] will be non-zero
  - sample_offdiag: dense vector of length equal to modelnew->n_nonzero storing the non-zero entries
  Output
  - Omega: row and column newcol are updated with the values in sample_diag and sample_offdiag. The rest of Omegaini is not altered

*/
void update_Omega(arma::sp_mat *Omega, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag) {
  int j, k;

  //Update Omegaini
  spmat_rowcol2zero(Omega, *newcol); //Set row and colum newcol to 0

  //Copy (sample_offdiag, sample_diag) to Omegaini
  Omega->at(*newcol, *newcol)= *sample_diag; //save diagonal element
  arma::SpMat<short>::iterator it;
  for (it= modelnew->begin(), j=0; it != modelnew->end(); ++it) {
    k= it.row();
    if (k != *newcol) { //skip diagonal element (already saved earlier)
      Omega->at(k, *newcol)= Omega->at(*newcol, k)= sample_offdiag->at(j,0); 
      j++;
    } 
  }

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


/* Obtain inverse of Omega[-newcol,-newcol], given inverse of Omega[-oldcol,-oldcol]

  INPUT/OUTPUT
  - Omega_submat_inv: at input this contains Omega[-oldcol,-oldcol]. At output, it contains Omega[-newcol,-newcol]

  INPUT
  - Omega: symmetric matrix Omega
  - oldcol: index of the column defining the submatrix Omega[-oldcol,-oldcol]
  - newcol: index of the column defining the submatrix Omega[-newcol,-newcol]

*/
void update_invOmega_submat(arma::mat *Omega_submat_inv, arma::sp_mat *Omega, int *oldcol, int *newcol) {

  int i, coldif, sel, *mapforw, p= Omega_submat_inv->n_cols;

  if (*oldcol != *newcol) {
    mapforw= ivector(0, p-1); 

    //Obtain map such that Omega[-oldcol,-oldcol] and Omega[-newcol,-newcol] differ only in column coldif
    //(after permuting the rows/columns of the latter according to the map)
    mapindexes_submat(mapforw, &coldif, oldcol, newcol, &p);

    //Obtain new values for updated entries in Omega[-oldcol,-oldcol]
    arma::sp_mat A_newcol(p,1);
    if (mapforw[coldif] < *newcol) sel= mapforw[coldif]; else sel= mapforw[coldif] + 1; //Omega[,sel] stores the new values
    for (i=0; i<p; i++) {
      if (mapforw[i] < *newcol) A_newcol.at(i,0)= Omega->at(mapforw[i],sel); else A_newcol.at(i,0)= Omega->at(mapforw[i]+1,sel); //permute values
    }

    //Update inverse
    update_inverse_permutedindex(Omega_submat_inv, mapforw, &A_newcol, &coldif); //store Omega_submat_inv[mapforw[i],mapforw[j]]

    free_ivector(mapforw, 0, p-1);
  }

}

/* Map indexes of two submatrices A1 and A2, each of which is obtained by dropping a column from the same square symmetric matrix A

  The mapping is such that A2[mapforw,mapforw] is equal to A1 except for row/column coldif (the row in A1 that is not in A2)

  INPUT

  - col1: A1 is obtained by removing row/column col1
  - col2: A2 is obtained by removing row/column col2
  - p: dimension of A1 and A2

  OUTPUT
  - mapforw: a vector of length p giving the mapping
  - mapback (only in variation defined below): a vector of length p giving the inverse mapping
             That is, A2= A2[mapforwd,mapforw][mapback,mapback] (also, A1[mapback,mapback] is equal to A2 except for one row)
  - coldif: A1 is equal to A2[mapforw, mapforw] except for row/column coldif

*/
void mapindexes_submat(int *mapforw, int *coldif, int *col1, int *col2, int *p) {
  int i, l, u;

  if (*col1 < *col2) {

    l= *col1;
    u= *col2;
    //Obtain forward mapping
    for (i=0; i<*p; i++) {
      if (i+1 == u) {
        mapforw[i]= l;
      } else if ((i>=l) && (i<u)) {
        mapforw[i]= i+1;
      } else {
        mapforw[i]= i;
      }
    }
    (*coldif)= *col2 - 1;

  } else if (*col1 > *col2) {

    int *mapback= ivector(0, *p -1);
    l= *col2;
    u= *col1;
    //Obtain forward mapping
    for (i=0; i<*p; i++) {
      mapforw[i]= i;
      if (i+1 == u) {
        mapback[i]= l;
      } else if ((i>=l) && (i<u)) {
        mapback[i]= i+1;
      } else {
        mapback[i]= i;
      }
    }
    iindexsort(mapback, mapforw, 0, *p - 1, 1); //Obtain inverse mapping
    (*coldif)= *col2;
    free_ivector(mapback, 0, *p -1);

  } else { //A1 == A2
    for (i=0; i<*p; i++) mapforw[i]= i;
    return;
  }

}


//Same, but also returns the reverse map in mapback
void mapindexes_submat(int *mapforw, int *mapback, int *coldif, int *col1, int *col2, int *p) {
  int i, l, u;

  if (*col1 < *col2) {

    l= *col1;
    u= *col2;
    //Obtain forward mapping
    for (i=0; i<*p; i++) {
      mapback[i]= i;
      if (i+1 == u) {
        mapforw[i]= l;
      } else if ((i>=l) && (i<u)) {
        mapforw[i]= i+1;
      } else {
        mapforw[i]= i;
      }
    }
    iindexsort(mapforw, mapback, 0, *p - 1, 1); //Obtain inverse mapping
    (*coldif)= *col2 - 1;

  } else if (*col1 > *col2) {

    l= *col2;
    u= *col1;
    //Obtain forward mapping
    for (i=0; i<*p; i++) {
      mapforw[i]= i;
      if (i+1 == u) {
        mapback[i]= l;
      } else if ((i>=l) && (i<u)) {
        mapback[i]= i+1;
      } else {
        mapback[i]= i;
      }
    }
    iindexsort(mapback, mapforw, 0, *p - 1, 1); //Obtain inverse mapping
    (*coldif)= *col2;

  } else { //A1 == A2
    for (i=0; i<*p; i++) mapforw[i]= mapback[i]= i;
    return;
  }

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
  arma::SpMat<short> model_offdiag(ggm->ncol -1, 1);

  if (npar > 0) { //if some off-diagonal entry is non-zero

    unsigned int i, j;
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



/* Compute log-joint (log-marginal likelihood + log prior) of regression proposal for model specifying non-zero entries in column colid of Omega

  The likelihood is            y[,colid] ~ N( y[,-colid] beta, phi I)
  The prior on parameters is   beta_gamma ~ N(0, tau phi I), where beta_gamma are the non-zero entries in beta
                               phi ~ InvGamma(1, lambda/2)
  The prior on the model       gamma ~ prod Bern(gamma_j; priorPars_p), where gamma_j = I(beta_j != 0)

  INPUT
  - model: entries that are non-zero
  - colid: column id
  - ggm: object storing info about the Gaussian graphical model
  - Omegainv_model: ignored, kept for compatibility with GGMrow_marg

  OUTPUT
  - logjoint: log-marginal + log-prior for model
  - m: ignored, kept for compatibility with GGMrow_marg
  - cholUinv: ignored, kept for compatibility with GGMrow_marg

*/
void GGMrow_marg_regression(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model) {

  unsigned int npar= model->n_nonzero -1;
  arma::SpMat<short> model_offdiag(ggm->ncol -1, 1);
  double alphahalf= 0.5; //prior on diagonal entries is Gamma(alpha=1, lambda/2)
  double num, den, sumy2= ggm->S.at(colid, colid);

  if (npar > 0) { //if some off-diagonal entry is non-zero

    unsigned int i;
    double tauinv= 1.0/ ggm->prCoef_tau, logdetXtXinv, nuhalf, ss;
    arma::uvec covariate_indexes(npar);
    arma::mat XtX(npar,npar), Xty(npar, 1), XtXinv(npar,npar), cholXtXinv(npar,npar);
    arma::SpMat<short>::iterator it;

    //Create XtX and Xty sufficient statistic matrices for model specified by model_offdiag
    //Note: model_offdiag indicates non-zero off-diagonal elements
    for (it= model->begin(), i=0; it != model->end(); ++it) {
      if (it.row() == colid) continue;
      model_offdiag.at(i, 0)= model->at(it.row(), 0);
      Xty.at(i,0)= ggm->S.at(it.row(), colid);
      covariate_indexes[i]= it.row();
      i++;
    }
     
    XtX = ggm->S.submat(covariate_indexes, covariate_indexes);
    for (i=0; i<npar; i++) XtX.at(i,i) += tauinv;

    //Inverse and determinant of XtX
    choldcinv_det(&XtXinv, &cholXtXinv, &logdetXtXinv, &XtX);   
    //invdet_posdef(&XtX, &XtXinv, &detXtX);

    arma::mat m= XtXinv * Xty;

    ss= ggm->prCoef_lambda + sumy2 - arma::as_scalar(m.t() * XtX * m);
    nuhalf= .5*ggm->n + alphahalf;

    num= gamln(&nuhalf) + alphahalf * log(0.5 * ggm->prCoef_lambda) + nuhalf * (log(2.0) - log(ss));
    den= .5*(ggm->n * LOG_M_2PI - logdetXtXinv) - .5 * (npar * log(tauinv)) + gamln(&alphahalf);
    (*logjoint) = num - den;

  } else { //if all off-diagonal entries are zero

    double term1= .5*ggm->n + alphahalf;
    num= alphahalf * log(ggm->prCoef_lambda) + gamln(&term1);
    den= .5 * (ggm->n) * (LOG_M_PI) + gamln(&alphahalf);
    (*logjoint) = num -den - term1 * log(ggm->prCoef_lambda + sumy2);

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

/* Check if sparse matrices A and B have <= maxdif different non-zero entries

INPUT
- A: Pointer to the first sparse matrix. It can be an arma::sp_mat or arma::SpMat<short>
- B: Pointer to the second sparse matrix. It can be an arma::sp_mat or arma::SpMat<short>
- maxdif: An integer representing the maximum allowable difference in the number of non-zero elements between matrices A and B.

OUTPUT: true if the number of non-zero elements in A-B is less than or equal to maxdif; otherwise, returns false.

NOTE: A and B are assumed to have the same dimensions, the function does not check said dimensions

*/

bool checkNonZeroDiff(const arma::SpMat<short>* A, const arma::SpMat<short>* B, int maxdif) {

    int nonZeroCount = 0;

    // Iterate over non-zero elements in matrix A
    for (arma::SpMat<short>::const_iterator itA = A->begin(); itA != A->end(); ++itA) {
        // Find the corresponding element in matrix B
        double valB = (*B)(itA.row(), itA.col());

        // Compare the elements
        if (*itA != valB) {
            nonZeroCount++;
            if (nonZeroCount > maxdif) return false;             // Check if the number of differences exceeds maxdif
        }
    }

    // Check for any additional non-zero elements in matrix B
    for (arma::SpMat<short>::const_iterator itB = B->begin(); itB != B->end(); ++itB) {
        // Find the corresponding element in matrix A
        double valA = (*A)(itB.row(), itB.col());

        // Increment the count for any additional non-zero elements in B
        if (valA != *itB) {
            nonZeroCount++;
            if (nonZeroCount > maxdif) return false;             // Check if the number of differences exceeds maxdif
        }
    }

    return true;
}

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
