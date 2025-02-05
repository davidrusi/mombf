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
  Within gdb, you should be able to print matrix A by typing: call print_matrix<arma::Mat<double> >(A)

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

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), LIT("LIT"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_LIT= (ggm->sampler == LIT);
  bool use_zigzag= (ggm->sampler == zigzag);

  //GGM_CDA(&Omegaini, ggm); //update Omegaini to local posterior mode

  //Create output matrix
  arma::sp_mat postSample(npars, ggm->niter - ggm->burnin);

  arma::mat postmean = arma::zeros(p,p);
  arma::Mat<int> postmeancount;
  postmeancount.zeros(p, p);
  arma::vec postmeanflat(p * (p+1)/2);
  arma::mat margpp = arma::zeros(p,p); 
  arma::Mat<int> margppcount;
  margppcount.zeros(p, p);
  arma::vec margppflat(p * (p+1)/2);

  //Obtain posterior samples
  if (use_gibbs || use_birthdeath || use_LIT) {

    GGM_MHwithinGibbs(&postSample, &postmean, &postmeancount, &margpp, &margppcount, &prop_accept, ggm, &Omegaini);

    //Rao-blackwellized estimators of posterior mean & inclusion probabilities
    postmean = postmean / postmeancount; 
    symmat2vec(&postmeanflat, &postmean);

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
  ret["postmean"] = postmeanflat;
  ret["margpp"] = margppflat;
  ret["prop_accept"]= prop_accept;

  return ret;

}




/* Metropolis-Hastings within Gibbs for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n

   The algorithm samples each column of Omega conditional on the rest (Gibbs). 
   To sample within each column of Omega, either Gibbs, birth-death-swap, or the LIT proposal of Zhang et al can be used 

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT
  - samples: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then samples has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.

  - postmean: sum of posterior means (Rao-Blackwellized). It should be initialized to zeroes at input

  - postmeancount: number of terms added in postmean. Posterior means are given by postmean/postmeancount. It should be initialized to zeroes at input

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - prop_accept: proportion of accepted proposals

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/

void GGM_MHwithinGibbs(arma::sp_mat *samples, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, jold, k, *sequence, iter= 0, niter10, number_accept, number_proposed, total_accept=0, total_proposed=0;
  int burnin= ggm->burnin, niter= ggm->niter, updates_per_iter= ggm->updates_per_iter, p= ggm->ncol;
  double *modelini_logprob= nullptr;
  arma::sp_mat::iterator it;

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), LIT("LIT");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_LIT= (ggm->sampler == LIT);
  arma::SpMat<short> *models= nullptr;
  arma::mat *model_logprob= nullptr, *invOmega_j;
  invOmega_j = new arma::mat(p-1, p-1);

  if (!use_gibbs && !use_birthdeath && !use_LIT) Rf_error("GGM_MHwithinGibbs requires the sampler to be 'Gibbs', 'birthdeath' or 'LIT'");

  if (niter >10) { niter10= niter / 10; } else { niter10= 1; }

  if (ggm->verbose) Rprintf(" Obtaining posterior samples\n");

  sequence= ivector(0, p-1);
  for (j=0; j<p; j++) sequence[j]= j;

  for (i=0; i < niter; i++) {

    samplei_wr(sequence, p, updates_per_iter); //consider updates_per_iter columns, in random order

    for (int jj=0; jj < updates_per_iter; jj++) {  //for each column

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
        GGM_Gibbs_singlecol(&samples_row, models, postmean, postmeancount, &margpp_row, &margppcount_row, 0, 0, (unsigned int) j, ggm, &Omegacol,  invOmega_j, model_logprob, modelini_logprob); //update row given by j
      } else {
        GGM_birthdeath_singlecol(&samples_row, models, postmean, postmeancount, &margpp_row, &margppcount_row, &number_accept, &number_proposed, 0, 0, (unsigned int) j, ggm, &use_LIT, &Omegacol,  invOmega_j, model_logprob, modelini_logprob); //update row given by j
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

  free_ivector(sequence, 0, p-1);
  if (ggm->verbose) { Rcout << "\r Done\n"; }

  delete invOmega_j;

}




// [[Rcpp::export]]
List modelSelectionGGM_globalC(arma::mat y, List prCoef, List prModel, List samplerPars, arma::sp_mat Omegaini) {
/* Interface for GGM_MHwithinGibbs_global and GGM_MHwithinGibbs_onlyglobal to be called from R 

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

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), LIT("LIT"), zigzag("zigzag");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_LIT= (ggm->sampler == LIT);
  bool use_zigzag= (ggm->sampler == zigzag);

  //Allocate memory
  dpropini= dvector(0, p-1);

  //Update Omegaini to local posterior mode
  //bool global_regression= ggm->global_regression;
  //ggm->global_regression= false;
  //GGM_CDA(&Omegaini, ggm);
  //ggm->global_regression= global_regression;

  //Vectors where to store proposal's models & their probabilities
  std::vector<arma::SpMat<short>> proposal_models(p);
  std::vector<std::vector<double>> proposal_logprob(p);
  std::vector<std::map<string, double>> map_logprob(p);

  //Propose models, store their proposal probability into proposal_logprob
  if (use_gibbs || use_birthdeath || use_LIT) {

    int updates_per_column= ggm->updates_per_column;
    ggm->updates_per_column= 1;
    GGM_global_proposal(&proposal_models, &proposal_logprob, &map_logprob, dpropini, ggm, &Omegaini);
    ggm->updates_per_column= updates_per_column;

  } else if (use_zigzag) {
    
    Rprintf("zigzag will be implemented soon\n");

  } else Rf_error("This sampler type is not currently implemented\n");

  //MCMC using independent proposal MH to combine the chains
  ggm->use_tempering= false;
  ggm->global_regression= ggm->global_insample= false;
  arma::sp_mat postSample(npars, niter - burnin);

  arma::mat postmean = arma::zeros(p,p);
  arma::Mat<int> postmeancount;
  postmeancount.zeros(p, p);
  arma::vec postmeanflat(p * (p+1)/2);
  arma::mat margpp = arma::zeros(p,p); 
  arma::Mat<int> margppcount;
  margppcount.zeros(p, p);
  arma::vec margppflat(p * (p+1)/2);
  
  if (ggm->prob_global < 0.9999) {
    GGM_MHwithinGibbs_global(&postSample, &postmean, &postmeancount, &margpp, &margppcount, &prop_accept, &proposal_models, &proposal_logprob, dpropini, &map_logprob, ggm, &Omegaini);
  } else {
    GGM_MHwithinGibbs_onlyglobal(&postSample, &postmean, &postmeancount, &margpp, &margppcount, &prop_accept, &proposal_models, &proposal_logprob, dpropini, ggm, &Omegaini);
  }

  //Rao-blackwellized estimators of posterior mean & inclusion probabilities
  postmean = postmean / postmeancount; 
  symmat2vec(&postmeanflat, &postmean);

  margpp = margpp / margppcount;
  for (int i=0; i < p; i++) margpp.at(i,i)= 1;
  symmat2vec(&margppflat, &margpp);

  //Return output
  List ret(p+5);
  ret[0]= postSample;
  ret[1]= postmeanflat;
  ret[2]= margppflat;
  ret[3]= prop_accept;

  //Store proposals into returned list (ret)
  std::vector<arma::SpMat<short>>::iterator it;
  for (it= proposal_models.begin(), j=4; it != proposal_models.end(); ++it, j++) { 
    ret[j]= (*it);
  }

  ret[p+4]= proposal_logprob;
 
   //Free memory and return output
  delete ggm;
  free_dvector(dpropini, 0, p-1);

  return ret;


}



/*global Gibbs and birth-death sampling for the precision matrix of a Gaussian graphical models, i.e. Omega where y_i ~ N(0,Omega^{-1}) for i=1,...,n. 

Each column is sampled independently from its conditional posterior, given the current value of Omega, i.e. the algorithm targets a pseudo-posterior rather than the actual posterior

INPUT
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)
  - Omegaini: initial value of Omega

OUTPUT
  - models: list where entry j stores indicators of Omega[,j] being zero. Formatted as a sparse matrix format (p rows and number of columns equal to the number of iterations)
  - model_logprop: log-posterior probability of the proposed model
  - map_logprob: its jth entry contains the log-proposal probabilities, stored as a map. Useful to query the proposal probability of any given model 

*/

void GGM_global_proposal(std::vector<arma::SpMat<short>> *models, std::vector<std::vector<double>> *model_logprop, std::vector<std::map<string, double>> *map_logprob, double *logprop_modelini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int j, p= ggm->ncol, p10, niter_prop, burnin_prop;
  double truncratio= ggm->truncratio;

  std::string Gibbs("Gibbs"), birthdeath("birthdeath"), LIT("LIT");
  bool use_gibbs= (ggm->sampler == Gibbs);
  bool use_birthdeath= (ggm->sampler == birthdeath);
  bool use_LIT= (ggm->sampler == LIT);

  niter_GGM_proposal(&niter_prop, &burnin_prop, &(ggm->niter), &(ggm->burnin), &p); //if not using a full scan, less iterations are needed

  if (!use_gibbs && !use_birthdeath && !use_LIT) Rf_error("GGM_global_proposal requires the sampler to be Gibbs, birthdeath or LIT");

  if (p >10) { p10= p / 10; } else { p10= 1; }

  if (ggm->verbose) Rprintf(" Running global proposals\n");

  //#pragma omp global for default(none) shared(p, p10, use_gibbs, burnin, niter, ggm, Omegaini, models, model_logprop)
  for (j=0; j < p; j++) { //for each column

    arma::SpMat<short> models_row(p, niter_prop - burnin_prop);
    arma::mat invOmega_j;
    if (!ggm->global_regression) invOmega_j= get_invOmega_j(Omegaini, j);
    arma::sp_mat Omegacol= Omegaini->col(j);
    arma::sp_mat *samples_row= nullptr;

    arma::vec *margpp_row= nullptr;
    arma::Col<int> *margppcount_row= nullptr;
    arma::mat modelj_logprop(1, niter_prop - burnin_prop);
    int number_accept, number_proposed;

    if (use_gibbs) {
      GGM_Gibbs_singlecol(samples_row, &models_row, nullptr, nullptr, margpp_row, margppcount_row, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm, &Omegacol,  &invOmega_j, &modelj_logprop, logprop_modelini+j);
    } else {
      GGM_birthdeath_singlecol(samples_row, &models_row, nullptr, nullptr, margpp_row, margppcount_row, &number_accept, &number_proposed, -burnin_prop, niter_prop-burnin_prop-1, (unsigned int) j, ggm, &use_LIT, &Omegacol,  &invOmega_j, &modelj_logprop, logprop_modelini+j);
    }

    //Obtain unique visited models and their probabilities, truncated so probability of top model / any other model <= truncratio
    unique_model_logprob(&(models->at(j)), &(model_logprop->at(j)), &(map_logprob->at(j)), &models_row, &modelj_logprop, &truncratio, logprop_modelini+j);
    if (truncratio > 0) {
      logprop_modelini[j]= max_xy(logprop_modelini[j], (model_logprop->at(j))[0] - log(truncratio)); //truncate probabilities according to maxratio
    }

    if (ggm->verbose) print_iterprogress(&j, &p, &p10);

  } //end for each colum

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
  arma::mat *sample_offdiag= nullptr, *invOmega_rest;
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

      ms->getJoint(&mcurrent, nullptr, nullptr, sample_offdiag, &sample_diag, modelopt, nullptr, false);

      i= changes_colid= 1;
      while ((changes_colid > 0) && (i <= maxit)) {  // (update model for column colid)

        changes_colid= 0;

        //For each entry j in colid, set entries to zero/non-zero if it improves the posterior model probability
        for (j= 0; j < p; j++) {

          if (j == (int) colid) continue;  //diagonal entry is always in

          if (modelopt->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;

          ms->getJoint(&mnew, nullptr, nullptr, sample_offdiag, &sample_diag, modelnew, modelopt, false);

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

      ms->getMode(&mcurrent, sample_offdiag, &sample_diag, modelopt, nullptr); //obtain posterior mode for parameters given model
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
  - samples: if not nullptr, each column of samples will store Omegacol after applying a Gibbs update to its entries

  - models: if not nullptr, stores model indicator corresponding to samples

  - postmean: sum of posterior mean of non-diagonal Omega entries. It should be initialized to zeroes at input

  - postmeancount: number of terms added in postmean. Posterior mean are given by postmean/postmeancount. It should be initialized to zero at input

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - model_logprob: If model_logprop has >1 rows, model_logprop[colid,j] stores the log posterior probability of model sampled at iteration j, other rows are unchanged. 
                   If model_logprop only has 1 row, then log posterior prob are stored in logprop[0,j]

  - modelini_logprob: log posterior probability of initial model

IMPORTANT: if not nullptr, at input samples and models should only contain zeroes

*/
void GGM_Gibbs_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::vec *margpp, arma::Col<int> *margppcount, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = nullptr, double *modelini_logprob = nullptr) {

  int i, j, p= ggm->ncol, col2save;
  double mcurrent, mnew, ppnew, sample_diag, mean_diag;
  arma::mat *sample_offdiag= nullptr, *mean_offdiag= nullptr;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun;

  if (ggm->global_regression) {
    marfun= &GGMrow_marg_regression;
  } else {
    marfun= &GGMrow_marg;
  }

  if ((model_logprob != nullptr) && (model_logprob->n_rows > 1)) col2save= colid; else col2save= 0;

  //Set random number generators
  unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count(); //single-thread version
  //unsigned seed = colid + 1; //multi-thread version. Set seed according to colid
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> unifdist(0.0,1.0);

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  ms->getJoint(&mcurrent, nullptr, nullptr, sample_offdiag, &sample_diag, model, nullptr, false);

  if (modelini_logprob != nullptr) (*modelini_logprob)= mcurrent;

  for (i= iterini; i <= iterfi; i++) {
    //For each entry j in colid, obtain log-marginal + log-prior for adding/removing j
    for (j= 0; j < p; j++) {

      if (j == (int) colid) continue;  //diagonal entry is always in

      if (model->at(j,0) != 0) modelnew->at(j,0) = 0; else modelnew->at(j,0) = 1;

      ms->getJoint(&mnew, nullptr, nullptr, sample_offdiag, &sample_diag, modelnew, model, false);

      ppnew = 1 / (1 + exp(mcurrent - mnew));

      if (margpp != nullptr) {
        if (model->at(j,0) == 0) margpp->at(j) += ppnew; else margpp->at(j) += 1 - ppnew;
      }
      if (margppcount != nullptr) margppcount->at(j) ++;

      double u= unifdist(generator);

      if (u < ppnew) { //if new model is accepted

        mcurrent= mnew;
        model_tmp_ptr= model;
        model= modelnew;
        modelnew= model_tmp_ptr;

      }

      modelnew->at(j,0)= model->at(j,0);

    } //end j for (iteration over columns)

    if (samples != nullptr) {

      mean_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      ms->getJoint(&mcurrent, mean_offdiag, &mean_diag, sample_offdiag, &sample_diag, model, nullptr, true); //update parameter values. DO NOT USE IN A MULTI-THREAD ENVIRONMENT (WRITES TO SHARED RANDOM SEED)
      if (i >= 0) {
        save_ggmsample_col(samples, model, &sample_diag, sample_offdiag, i, colid); //copy current model into samples[,i] as flat vector 

        if (postmean != nullptr) {
          save_postmean(postmean, postmeancount, model, mean_offdiag, &mean_diag, colid); //save posterior mean of non-diagonal entries
        }
      }
      delete mean_offdiag;
      delete sample_offdiag;

    }

    if ((models != nullptr) && (i >= 0)) save_ggmmodel_col(models, model, i, colid); //copy model indicator into models[,i] as flat vector

    if ((i>=0) && (model_logprob != nullptr)) model_logprob->at(col2save, i)= mcurrent;

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
  - samples: if not nullptr, each column of samples will store Omegacol after applying a Gibbs update to its entries

  - models: if not nullptr, stores model indicator corresponding to samples

  - postmean: sum of posterior mean of non-diagonal Omega entries. It should be initialized to zeroes at input

  - postmeancount: its element i indicates the number of times that column i in postmean was updated. It should be initialized to zeroes at input

  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input

  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input

  - model_logprob: If model_logprop has >1 rows, model_logprop[colid,j] stores the log posterior probability of model sampled at iteration j, other rows are unchanged. 
                   If model_logprop only has 1 row, then log posterior prob are stored in logprop[0,j]

  - number_accept: number of accepted proposals

  - number_proposals: number of proposals

  - modelini_logprob: log posterior probability of initial model

IMPORTANT: if not nullptr, at input samples and models should only contain zeroes

*/

void GGM_birthdeath_singlecol(arma::sp_mat *samples, arma::SpMat<short> *models, arma::mat *postmean, arma::Mat<int> *postmeancount,  arma::vec *margpp, arma::Col<int> *margppcount, int *number_accept, int *number_proposed, int iterini, int iterfi, unsigned int colid, ggmObject *ggm, bool *use_LIT, arma::sp_mat *Omegacol, arma::mat *invOmega_rest, arma::mat *model_logprob = nullptr, double *modelini_logprob = nullptr) {

  int i, j, p= ggm->ncol, updates_per_column= ggm->updates_per_column, index_birth, index_death, movetype, col2save, colid_int= (int) colid, ppcount=0;
  double pbirth= ggm->pbirth, pdeath= ggm->pdeath, pswap= ggm->pswap, mcurrent, mnew, dpropcurrent, dpropnew, ppnew, sample_diag, mean_diag;
  arma::mat *mean_offdiag= nullptr, *sample_offdiag= nullptr;
  arma::sp_mat::const_iterator it;
  arma::SpMat<short> *model, *modelnew, *model_tmp_ptr;
  modselIntegrals_GGM *ms;
  pt2GGM_rowmarg marfun;

  if (ggm->global_regression) {
    marfun= &GGMrow_marg_regression;
  } else {
    marfun= &GGMrow_marg;
  }

  (*number_accept)= (*number_proposed)= 0;
  if ((model_logprob != nullptr) && (model_logprob->n_rows > 1)) col2save= colid; else col2save= 0;

  //Set random number generators
  unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count(); //single-thread version
  //unsigned seed = colid + 1; //multi-thread version. Set seed according to colid
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> unifdist(0.0,1.0);

  //Initialize model and modelnew
  model= new arma::SpMat<short>(p, 1);
  modelnew= new arma::SpMat<short>(p, 1);
  for (it= Omegacol->begin(); it != Omegacol->end(); ++it) model->at(it.row(), it.col())= modelnew->at(it.row(), it.col())= 1;

  //Obtain log-marginal + log-prior for current model
  ms= new modselIntegrals_GGM(marfun, ggm, colid, invOmega_rest);

  ms->getJoint(&mcurrent, nullptr, nullptr, sample_offdiag, &sample_diag, model, nullptr, false);

  if (modelini_logprob != nullptr) (*modelini_logprob)= mcurrent;

  for (i= iterini; i <= iterfi; i++) {

    for (j=1; j<=updates_per_column; j++) {  //Make updates_per_column updates

      if (*use_LIT) { //Locally-informed proposal

        GGM_LIT_proposal(modelnew, &index_birth, &index_death, &movetype, &dpropnew, &dpropcurrent, model, &colid_int, ggm, ms, true);

      } else {  //Non-informed birth-death-swap proposal

        GGM_birthdeathswap_proposal(modelnew, &index_birth, &index_death, &movetype, &dpropnew, &dpropcurrent, model, &colid_int, &pbirth, &pdeath, &pswap, true);

      } 
      
      (*number_proposed)++;

      if ((index_birth != -1) || (index_death != -1)) {     //if proposed modelnew != model

          ms->getJoint(&mnew, nullptr, nullptr, sample_offdiag, &sample_diag, modelnew, model, false);  //obtain marginal likelihood for modelnew

          ppnew = exp(mnew - mcurrent + dpropcurrent - dpropnew); //probability of accepting modelnew

          if (margpp != nullptr) {
            update_margpp_raoblack(margpp, min_xy(1,ppnew), model, modelnew);
            ppcount++;
          }
         
          if ((ppnew > 1) || (unifdist(generator) < ppnew)) { //if new model is accepted
           
            mcurrent= mnew;
            model_tmp_ptr= model;
            model= modelnew;
            modelnew= model_tmp_ptr;
            (*number_accept)++;
           
          }
           
          if (index_birth != -1) modelnew->at(index_birth,0)= model->at(index_birth,0);
          if (index_death != -1) modelnew->at(index_death,0)= model->at(index_death,0);

      } else {   //if proposed modelnew == model

          ppnew= 1;
          if (margpp != nullptr) {
            update_margpp_raoblack(margpp, min_xy(1,ppnew), model, modelnew);
            ppcount++;
          }


      }

    } //end for j

    if (samples != nullptr) {

      mean_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      sample_offdiag= new arma::mat(model->n_nonzero - 1, 1);
      ms->getJoint(&mcurrent, mean_offdiag, &mean_diag, sample_offdiag, &sample_diag, model, nullptr, true); //obtain marginal likelihood & posterior sample
      if (i >= 0) {
        save_ggmsample_col(samples, model, &sample_diag, sample_offdiag, i, colid); //copy current model into samples[,i] as flat vector
        if (postmean != nullptr) {
          save_postmean(postmean, postmeancount, model, mean_offdiag, &mean_diag, colid); //save posterior mean of non-diagonal entries
        }
      }
      delete mean_offdiag;
      delete sample_offdiag;

    }

    if ((models != nullptr) && (i >= 0)) save_ggmmodel_col(models, model, i, colid); //copy model indicator into models[,i] as flat vector

    if ((i>=0) && (model_logprob != nullptr)) model_logprob->at(col2save, i)= mcurrent;

  } //end i for (Gibbs iterations)

  if (margpp != nullptr) for (j = 0; j < p; j++) margppcount->at(j) += ppcount;

  delete model;
  delete modelnew;
  delete ms;
  
}


/* Update Rao-Blackwellized estimate of edge inclusion probability 

  INPUT
  - ppnew: probability of accepting modelnew
  - model: model at the past iteration
  - modelnew: model being proposed for the current iteration

  INTPUT/OUTPUT
  - margpp: sum of marginal inclusion probabilities for edges 0,...,p-1. At output, the probability that each edge is present at the new iteration is added to margpp

*/

void update_margpp_raoblack(arma::vec *margpp, double ppnew, arma::SpMat<short> *model, arma::SpMat<short> *modelnew) {

  arma::SpMat<short> modeldif= *model - *modelnew;
  arma::SpMat<short>::iterator it;
  for (it = modeldif.begin(); it != modeldif.end(); ++it) {
    if (*it < 0) {  //Case model==0, modelnew==1
      (margpp->at(it.row())) += min_xy(ppnew, 1);
    } else {        //Case model==1, modelnew==0
      (margpp->at(it.row())) += 1 - min_xy(ppnew, 1);
    }
  }
  for (it = model->begin(); it != model->end(); ++it) {
    if (modelnew->at(it.row(),0) == 1) (margpp->at(it.row())) += 1;  //Case model==1, modelnew==1
  }

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
  - dpropnew: log of the proposal probability of modelnew. If the proposed move was not possible, dpropnew= -INFINITY is returned
  - dpropcurrent: log of the proposal probability of model. If the propsed move was not possible, dpropcurrent is not updated

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

    (*dpropnew)= -INFINITY;
    //(*dpropnew)= (*dpropcurrent) -INFINITY;

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
  - dpropnew: log of the proposal probability of modelnew. If the proposed move was not possible, dpropnew= -INFINITY is returned
  - dpropcurrent: log of the proposal probability of model. If the propsed move was not possible, dpropcurrent is not updated

*/
void GGM_birthdeathswap_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, double *pbirth, double *pdeath, double *pswap, bool setmodelnew) {

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

    (*dpropnew)= dbirthdeathswap(&modelnew_colid, &model_colid, *pbirth, *pdeath, *pswap, true);
    (*dpropcurrent)= dbirthdeathswap(&model_colid, &modelnew_colid, *pbirth, *pdeath, *pswap, true);

  } else {

    (*dpropnew)= -INFINITY;
    //(*dpropnew)= (*dpropcurrent) -INFINITY;

  }

}


/*  Propose GGM model update using the LIT algorithm of Zhou, Yang, Vats, Roberts & Rosenthal. 
    Also, return the proposal density for the new model and for the current model (reverse move)

  INPUT

  - model: current model, a p x 1 matrix indicating the non-zero parameters
  - colid: column of the GGM precision matrix Omega currently being updated
  - ggm: ggmObject with info about the graphical model
  - ms: modselIntegrals_GGM object providing methods to compute & store marginal likelihoods for each possible model
  - setmodelnew: if false, modelnew is assumed to be equal to model at input, so that only the entries being born/dead are updated. Else, all entries in modelnew are set

  OUTPUT

  - modelnew: model proposed by the birth-death move
  - index_birth: index of parameter being born. If no entries were born, then index_birth= -1
  - index_death: index of parameter being killed. If no entries were killed, then index_death= -1
  - movetype: 1 if a birth move was chosen, 2 if a death move was chosen. Swap moves are not implemented
  - dpropnew: log of the proposal probability of modelnew. If the proposed move was not possible, dpropnew= -INFINITY is returned
  - dpropcurrent: log of the proposal probability of model. If the propsed move was not possible, dpropcurrent is not updated
*/
void GGM_LIT_proposal(arma::SpMat<short> *modelnew, int *index_birth, int *index_death, int *movetype, double *dpropnew, double *dpropcurrent, arma::SpMat<short> *model, int *colid, ggmObject *ggm, modselIntegrals_GGM *ms, bool setmodelnew) {

  const int dnonzero= model->n_nonzero - 1, d= model->n_rows - 1, dzero= d - dnonzero;
  double u;

  if (setmodelnew) (*modelnew)= *model;

  //Make birth/death proposal
  //arma::SpMat<short> model_colid= (*model);
  //model_colid.shed_row(*colid);

  (*index_birth)= (*index_death)= -1; //return value if proposed update cannot be done
  u = runifC();
  if (u < ggm->pbirth) { //birth move

    (*movetype)= 1;
    if (dzero > 0) {  //if dzero==0, all entries are already alive

      std::vector<int> indexes_birth(dzero);
      std::vector<double> proposal_birth(dzero);
      dprop_LIT_birth_GGM(&proposal_birth, &indexes_birth, model, ggm, ms); //obtain proposal kernel for birth moves

      //Propose new model
      int idx_birth= rmultinomial(&proposal_birth);
      (*dpropnew) = log(proposal_birth.at(idx_birth));
      (*dpropnew) += ggm->log_pbirth; 
      (*index_birth)= indexes_birth.at(idx_birth);
      modelnew->at(*index_birth, 0)= 1;

      //Proposal density for reverse move
      int dnonzero_new= dnonzero + 1;
      std::vector<int> indexes_death(dnonzero_new);
      std::vector<double> proposal_death(dnonzero_new);

      dprop_LIT_death_GGM(&proposal_death, &indexes_death, colid, modelnew, ggm, ms); //obtain proposal density for reverse death moves

      int idx_death= 0;
      while (indexes_death[idx_death] != (*index_birth)) { idx_death++; }
      (*dpropcurrent) = log(proposal_death.at(idx_death));
      (*dpropcurrent) += ggm->log_pdeath;

    } else {

      (*dpropnew)= -INFINITY;

    }

  } else { //death move

    (*movetype)= 2;
    if (dnonzero>0) {  //if all entries are already zero, none can be killed

      std::vector<int> indexes_death(dnonzero);
      std::vector<double> proposal_death(dnonzero);
      dprop_LIT_death_GGM(&proposal_death, &indexes_death, colid, modelnew, ggm, ms); //obtain proposal kernel for death moves

      //Propose new model
      int idx_death= rmultinomial(&proposal_death);
      (*dpropnew) = log(proposal_death.at(idx_death));
      (*dpropnew) += ggm->log_pdeath;
      (*index_death)= indexes_death.at(idx_death);
      modelnew->at(*index_death,0)= 0;

      //Proposal density for reverse move
      int dzero_new= dzero + 1;
      std::vector<int> indexes_birth(dzero_new);
      std::vector<double> proposal_birth(dzero_new);

      dprop_LIT_birth_GGM(&proposal_birth, &indexes_birth, modelnew, ggm, ms); //obtain proposal density for reverse birth move

      int idx_birth= 0;
      while (indexes_birth[idx_birth] != (*index_death)) { idx_birth++; }
      (*dpropcurrent) = log(proposal_birth.at(idx_birth));
      (*dpropcurrent) += ggm->log_pbirth;

    } else {

      (*dpropnew)= -INFINITY;

    }

  }

}


/* Obtain proposal probabilities for birth proposal in LIT algorithm of Zhou, Yang, Vats, Roberts & Rosenthal

  INPUT
  - ggm: ggmObject with info about the graphical model
  - ms: modselIntegrals_GGM object providing methods to compute & store marginal likelihoods for each possible model

  OUTPUT
  - proposal_kernel: proposal pmf for each possible birth move
  - indexes_birth: indexes of the edge being born in each possible birth move (ranging from 0 to p-1, p being the total number of variables)

*/
void dprop_LIT_birth_GGM(std::vector<double> *proposal_kernel, std::vector<int> *indexes_birth, arma::SpMat<short> *model, ggmObject *ggm, modselIntegrals_GGM *ms) {

  arma::SpMat<short> modelnew= (*model);

  int j, modelid= 0, dzero= proposal_kernel->size();
  double dpostcur, dpostnew, proposal_sum= 0, proposal_max= -INFINITY;

  ms->getJoint(&dpostcur, nullptr, nullptr, nullptr, nullptr, model, nullptr, false); //log-posterior for current model

  for (j=0; modelid < dzero; j++) {

    if (model->at(j,0) == 0) {

      modelnew.at(j,0)= 1;
      ms->getJoint(&dpostnew, nullptr, nullptr, nullptr, nullptr, &modelnew, model, false); //log-posterior for proposed model
      modelnew.at(j,0)= 0;
      proposal_kernel->at(modelid)= max_xy(min_xy(dpostnew - dpostcur, ggm->ubound_birth), ggm->lbound_birth);
      proposal_max= max_xy(proposal_max, proposal_kernel->at(modelid));
      indexes_birth->at(modelid)= j;
      modelid++;

    }

  }

  //Normalize probabilities
  for (j=0; j < dzero; j++) proposal_sum += exp(proposal_kernel->at(j) - proposal_max);
  for (j=0; j < dzero; j++) proposal_kernel->at(j) = exp(proposal_kernel->at(j) - proposal_max) / proposal_sum;

}


/* Obtain proposal probabilities for death proposal in LIT algorithm of Zhou, Yang, Vats, Roberts & Rosenthal

  INPUT
  - ggm: ggmObject with info about the graphical model
  - ms: modselIntegrals_GGM object providing methods to compute & store marginal likelihoods for each possible model
  - colid: index of the column in the GGM being updated (starts at 0)

  OUTPUT
  - proposal_kernel: proposal pmf for each possible death move
  - indexes_death: indexes of the edge being born in each possible death move (ranging from 0 to p-1, p being the total number of variables)

*/
void dprop_LIT_death_GGM(std::vector<double> *proposal_kernel, std::vector<int> *indexes_death, int *colid, arma::SpMat<short> *model, ggmObject *ggm, modselIntegrals_GGM *ms) {

  arma::SpMat<short> modelnew= (*model);
  arma::SpMat<short>::iterator it;

  int j, modelid= 0;
  double dpostcur, dpostnew, proposal_sum= 0, proposal_max= -INFINITY;

  ms->getJoint(&dpostcur, nullptr, nullptr, nullptr, nullptr, model, nullptr, false); //log-posterior for current model

  for (it= model->begin(); it != model->end(); ++it) {

    j= it.row();
    if (j == *colid) continue; //diagonal term is always in

    modelnew.at(j,0)= 0;
    ms->getJoint(&dpostnew, nullptr, nullptr, nullptr, nullptr, &modelnew, model, false); //log-posterior for proposed model
    modelnew.at(j,0)= 1;

    proposal_kernel->at(modelid)= max_xy(min_xy(dpostnew - dpostcur, ggm->ubound_death), ggm->lbound_death);
    proposal_max= max_xy(proposal_max, proposal_kernel->at(modelid));
    indexes_death->at(modelid)= j;
    modelid++;

  }

  //Normalize probabilities
  int nalive= proposal_kernel->size();
  for (j=0; j < nalive; j++) proposal_sum += exp(proposal_kernel->at(j) - proposal_max);
  for (j=0; j < nalive; j++) proposal_kernel->at(j) = exp(proposal_kernel->at(j) - proposal_max) / proposal_sum;

}


/* Determine number of iterations to use in GGM global proposals 

  The number of expected draws from each global proposal is m= niter / p, and its SD is s= sqrt(niter (1/p) (1 - 1/p))
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

  if (logprobini != nullptr) {
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


/* Obtain posterior samples for Omega via Metropolis-within-Gibbs

INPUT

  - proposal_models: vector where jth entry contains the proposal models for Omega[,j] (precision matrix entries being zero/non-zero), sorted in decreasing probability
  - proposal_logprob: vector where jth entry contains log-proposal densities for proposal_models[j]
  - proposal_modelids: string identifiers of the models in proposal_models (unsorted). Useful to quickly check if a model is in proposal_models
  - ggm: object of class ggmObject storing y, the prior distribution and Gibbs sampling parameters (burnin, number of iterations...)

OUTPUT

  - postSample: sparse matrix where each column corresponds to a posterior sample. If the Gaussian precision matrix Omega is p x p, then postSample has p*(p+1)/2 columns storing Omega in column order. That is, the first p entries correspond to Omega[,1] (first column in Omega), the next p-1 entries to Omega[2:p,2], the next p-2 entries to Omega[3:p,3], and so on.
  - postmean: sum of posterior mean of Omega across iterations.
  - postmeancount: number of additions performed for each entry in postmean
  - margpp: sum of marginal posterior inclusion probabilities (Rao-Blackwellized). It should be initialized to zeroes at input
  - margppcount: number of terms added in margpp. Posterior inclusion prob are given by margpp/margppcount. It should be initialized to zeroes at input
  - prop_accept: proportion of accepted proposals

INPUT/OUTPUT

  - Omegaini: initial value of Omega, this is updated at each iteration and the value at the last iteration is returned

*/
void GGM_MHwithinGibbs_global(arma::sp_mat *postSample, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, std::vector<std::map<string, double>> *map_logprob, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, i10, iter, newcol, oldcol, *sequence, index_birth, index_death, movetype, number_accept= 0, proposal_idx;
  int burnin= ggm->burnin, niter= burnin + postSample->n_cols, p= ggm->ncol, updates_per_column= ggm->updates_per_column, updates_per_iter= ggm->updates_per_iter;
  double dpostnew, dpostold, dpropnew, dpropold=0, ppnew, sample_diag, mean_diag;
  double prob_global= ggm->prob_global;
  bool global_draw;
  char *zerochar;
  zerochar = (char *) calloc(p+1, sizeof(char));
  for (i=0; i<p; i++) zerochar[i]= '0';
  std::map<string, double>::iterator itmap;

  arma::mat *mean_offdiag= nullptr, *sample_offdiag= nullptr, *invOmega_newcol;
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::SpMat<short>> modelold(p);
  arma::SpMat<short> modelnew;
  arma::SpMat<short>::iterator it;
  pt2GGM_rowmarg marfun= &GGMrow_marg;
  std::vector<double*> cdf_proposal(p);
  std::vector<int> nproposal(p);
  std::vector<arma::vec> margppcol(p);
  std::vector<int> ppcount(p);
  std::string s;

  for (i = 0; i < p; i++) {
    margppcol[i] = arma::zeros<arma::vec>(p);
    ppcount[i]= 0;
  }
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

    //Pre-compute global proposal's cdf, so we can use rdisc_pcum later
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

  sequence= ivector(0, p-1);
  for (j=0; j<p; j++) sequence[j]= j;

  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {

    samplei_wr(sequence, p, updates_per_iter); //do updates_per_iter column updates (in random order)

    for (j=0; j<updates_per_iter; j++) { //iterate over columns in random order

      newcol= sequence[j];

      //Obtain inverse of Omegaini[-newcol,-newcol]
      if ((i==0) && (j==0)) {
        (*invOmega_newcol)= get_invOmega_j(Omegaini, newcol); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_newcol, Omegaini, &oldcol, &newcol);
      }
      oldcol= newcol;

      modselIntegrals_GGM *ms;
      ms= new modselIntegrals_GGM(marfun, ggm, newcol, invOmega_newcol);

      for (k=0; k < updates_per_column; k++) {  //perform multiple updates per column

        global_draw= (runifC() < prob_global); //choose proposal type
    
        if (global_draw && (dpropini[newcol] != -INFINITY)) { //use global proposal

          dpropold= dpropini[newcol]; //log-proposal for current model
          proposal_idx= rdisc_pcum(cdf_proposal[newcol], nproposal[newcol]);  //index of proposal
          dpropnew= (proposal_logprob->at(newcol))[proposal_idx]; //log-proposal for new value
          modelnew= ((*proposal_models)[newcol]).col(proposal_idx); //proposed value

        } else { //use birth-death proposal

          GGM_birthdeathswap_proposal(&modelnew, &index_birth, &index_death, &movetype, &dpropnew, &dpropold, &(modelold[newcol]), &newcol, &(ggm->pbirth), &(ggm->pdeath), &(ggm->pswap), true);

        }

        if (k==0) ms->getJoint(&dpostold, nullptr, nullptr, sample_offdiag, &sample_diag, &(modelold[newcol]), nullptr, false); //log-posterior for current model
        ms->getJoint(&dpostnew, nullptr, nullptr, sample_offdiag, &sample_diag, &modelnew, &(modelold[newcol]), false); //log-posterior for proposed model

        ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);

        if (margpp != nullptr) {
          update_margpp_raoblack(&(margppcol[newcol]), min_xy(1,ppnew), &(modelold[newcol]), &modelnew); 
          (ppcount[newcol])++;
        }

        if ((ppnew >= 1) || (runifC() < ppnew)) { //if update is accepted

          number_accept++;
          modelold[newcol]= modelnew;
          dpostold= dpostnew;
          //Store log global proposal density of updated model
          s= getModelid(&modelnew, zerochar);
          itmap= (map_logprob->at(newcol)).find(s);
          if (itmap != (map_logprob->at(newcol)).end()) {
            dpropini[newcol]= itmap->second;
          } else {
            dpropini[newcol]= -INFINITY;
          }

        }

      } //end k for (perform multiple updates per column)

      //Sample Omega[,newcol] given model
      mean_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      sample_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      ms->getJoint(&dpostnew, mean_offdiag, &mean_diag, sample_offdiag, &sample_diag, &(modelold[newcol]), nullptr, true); //sample Omega[,newcol]
      if (postmean != nullptr) {
        save_postmean(postmean, postmeancount, &(modelold[newcol]), mean_offdiag, &mean_diag, newcol); //save posterior mean of non-diagonal entries
      }
      update_Omega(Omegaini, &newcol, &sample_diag, &(modelold[newcol]), sample_offdiag); //copy sampled values into Omegaini[,newcol]
      delete mean_offdiag;
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

  (*prop_accept)= (number_accept + 0.0) / (i * updates_per_iter * updates_per_column + 0.0);

  //Rao-Blackwellized edge probabilities
  for (i = 0; i < p; i++) {
    for (j = 0; j <=i; j++) { margpp->at(i,j) += margppcol[i][j]; margppcount->at(i,j) += ppcount[i]; }
    for (j = i; j < p; j++) { margpp->at(j,i) += margppcol[i][j]; margppcount->at(j,i) += ppcount[i]; }
  }

  //Free memory
  free_ivector(sequence, 0, p-1);
  for (newcol=0; newcol < p; newcol++) free_dvector(cdf_proposal[newcol], 0, nproposal[newcol]);
  free((char  *) zerochar);
  delete invOmega_newcol;

}


void GGM_MHwithinGibbs_onlyglobal(arma::sp_mat *postSample, arma::mat *postmean, arma::Mat<int> *postmeancount, arma::mat *margpp, arma::Mat<int> *margppcount, double *prop_accept, std::vector<arma::SpMat<short>> *proposal_models, std::vector<std::vector<double>> *proposal_logprob, double *dpropini, ggmObject *ggm, arma::sp_mat *Omegaini) {

  int i, j, k, i10, iter, newcol, oldcol, *sequence, number_accept= 0, proposal_idx;
  int burnin= ggm->burnin, niter= burnin + postSample->n_cols, p= ggm->ncol, updates_per_column= ggm->updates_per_column, updates_per_iter= ggm->updates_per_iter;
  double dpostnew, dpostold, dpropnew, dpropold=0, ppnew, sample_diag, mean_diag;
  arma::mat *mean_offdiag= nullptr, *sample_offdiag= nullptr, *invOmega_newcol;
  std::vector<double> diagcur(p), bnew(p);
  std::vector<arma::SpMat<short>> modelold(p);
  arma::SpMat<short> modelnew;
  pt2GGM_rowmarg marfun= &GGMrow_marg;
  std::vector<double*> cdf_proposal(p);
  std::vector<int> nproposal(p);
  std::vector<arma::vec> margppcol(p);
  std::vector<int> ppcount(p);

  for (i = 0; i < p; i++) {
    margppcol[i] = arma::zeros<arma::vec>(p);
    ppcount[i]= 0;
  }
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

  sequence= ivector(0, p-1);
  for (j=0; j<p; j++) sequence[j]= j;

  //MCMC iterations
  for (i = 0, iter = 0; i < niter; i++) {

    samplei_wr(sequence, p, updates_per_iter); //do updates_per_iter column updates (in random order)

    for (j=0; j<updates_per_iter; j++) { //iterate over columns in random order

      newcol= sequence[j];

      //Obtain inverse of Omegaini[-newcol,-newcol]
      if ((i==0) && (j==0)) {
        (*invOmega_newcol)= get_invOmega_j(Omegaini, newcol); 
      } else { //use rank 1 update (quadratic cost in p)
        update_invOmega_submat(invOmega_newcol, Omegaini, &oldcol, &newcol);
      }
      oldcol= newcol;

      modselIntegrals_GGM *ms;
      ms= new modselIntegrals_GGM(marfun, ggm, newcol, invOmega_newcol);

      for (k=0; k < updates_per_column; k++) {  //perform multiple updates per column

        proposal_idx= rdisc_pcum(cdf_proposal[newcol], nproposal[newcol]);  //index of proposal
        dpropnew= (proposal_logprob->at(newcol))[proposal_idx]; //log-proposal for proposed model
        modelnew= ((*proposal_models)[newcol]).col(proposal_idx); //proposed value

        if (k==0) ms->getJoint(&dpostold, nullptr, nullptr, sample_offdiag, &sample_diag, &(modelold[newcol]), nullptr, false); //log-posterior for current model
        ms->getJoint(&dpostnew, nullptr, nullptr, sample_offdiag, &sample_diag, &modelnew, &(modelold[newcol]), false); //log-posterior for proposed model

        dpropold= dpropini[newcol]; //log-proposal for current model

        ppnew = exp(dpostnew - dpostold + dpropold - dpropnew);

        if (margpp != nullptr) {
          update_margpp_raoblack(&(margppcol[newcol]), min_xy(1,ppnew), &(modelold[newcol]), &modelnew); 
          (ppcount[newcol])++;
        }

        if ((ppnew >= 1) || (runifC() < ppnew)) { //if update is accepted

          if (i >= burnin) number_accept++;
          dpropini[newcol]= dpropnew;
          modelold[newcol]= modelnew;
          dpostold= dpostnew;

        }

      } //end k for (perform multiple updates per column)

      //Sample Omega[,newcol] given model
      mean_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      sample_offdiag= new arma::mat((modelold[newcol]).n_nonzero - 1, 1);
      ms->getJoint(&dpostnew, mean_offdiag, &mean_diag, sample_offdiag, &sample_diag, &(modelold[newcol]), nullptr, true); //sample Omega[,newcol]
      if (postmean != nullptr) {
        save_postmean(postmean, postmeancount, &(modelold[newcol]), mean_offdiag, &mean_diag, newcol); //save posterior mean of non-diagonal entries
      }
      update_Omega(Omegaini, &newcol, &sample_diag, &(modelold[newcol]), sample_offdiag); //copy sampled values into Omegaini[,newcol]
      delete mean_offdiag;
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

  (*prop_accept)= (number_accept + 0.0) / (i * updates_per_iter * updates_per_column + 0.0);

  //Rao-Blackwellized edge probabilities
  for (i = 0; i < p; i++) {
    for (j = 0; j <=i; j++) { margpp->at(i,j) += margppcol[i][j]; margppcount->at(i,j) += ppcount[i]; }
    for (j = i; j < p; j++) { margpp->at(j,i) += margppcol[i][j]; margppcount->at(j,i) += ppcount[i]; }
  }

  //Free memory
  for (newcol=0; newcol < p; newcol++) free_dvector(cdf_proposal[newcol], 0, nproposal[newcol]);

  free_ivector(sequence, 0, p-1);
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

/* Same as update_Omega, but (sample_diag, sample_offdiag) are added to Omegaini[,newcol] (rather than replacing its entries) */
void addto_Omega(arma::mat *Omega, int *newcol, double *sample_diag, arma::SpMat<short> *modelnew, arma::mat *sample_offdiag) {
  int j, k;

  //Copy (sample_offdiag, sample_diag) to Omegaini
  Omega->at(*newcol, *newcol)+= *sample_diag; //save diagonal element
  
  arma::SpMat<short>::iterator it;
  for (it= modelnew->begin(), j=0; it != modelnew->end(); ++it) {
    k= it.row();
    if (k != *newcol) { //skip diagonal element (already saved earlier)
      Omega->at(k, *newcol)+= Omega->at(*newcol, k)= sample_offdiag->at(j,0); 
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


/* Add the posterior mean of Omega[,colid] to postmean

  INPUT
  - model: indicates the non-zero entries in Omega[,colid]
  - mean_offdiag: posterior mean of the non-zero entries in Omega[,colid]
  - mean_diag: posterior mean of the diagonal entry Omega[colid,colid]
  - colid: index of the column to be updated (starts at 0)

  INPUT/OUTPUT
  - postmean: updated sum of posterior means of Omega[,colid]
  - postmeancount: updated number of times in which we added to postmean

*/
void save_postmean(arma::mat *postmean, arma::Mat<int> *postmeancount, arma::SpMat<short> *model, arma::mat *mean_offdiag, double *mean_diag, unsigned int colid) {
    int j = 0, p= (int) model->n_rows;
    arma::SpMat<short>::iterator itmodel;

    postmean->at(colid, colid) += *mean_diag;
    for (itmodel = model->begin(); itmodel != model->end(); ++itmodel) { 
        if (itmodel.row() == colid) continue;
        postmean->at(itmodel.row(), colid) += mean_offdiag->at(j, 0); 
        postmean->at(colid, itmodel.row()) = postmean->at(itmodel.row(), colid);
        j++;
    }

    for (j = 0; j < p; j++) { 
      postmeancount->at(j, colid)++; 
      postmeancount->at(colid, j) = postmeancount->at(j, colid); 
    }

}



/* Compute log-joint (log-marginal likelihood + log prior) for model specifying non-zero entries in column colid of Omega, given the entries selected by model of the inverse of Omega[-colid,-colid]

  INPUT
  - model: entries that are non-zero
  - colid: column id
  - ggm: object storing info about the Gaussian graphical model
  - Omegainv: inverse of Omega[-colid,-colid]
  - cholU_old: Cholesky decomposition of U for a previous model (modelold) that differs from model by adding/dropping at most 1 non-zero entry
  - modelold: entries that were non-zero in the previous model

  OUTPUT
  - logjoint: log-marginal + log-prior for model
  - m: mean of off-diagonal entries, i.e. Omega[-colid,colid] ~ multivariate Normal(m, Uinv)
  - cholUinv: Cholesky decomposition of Uinv, i.e. Uinv= t(cholUinv) * cholUinv
  - cholU: Cholesky decomposition of U

*/
void GGMrow_marg(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::mat *cholU, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv, arma::mat *cholU_old=nullptr, arma::SpMat<short> *modelold=nullptr) {

  unsigned int npar= model->n_nonzero -1;
  int colidint= (int) colid;
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
     
    //Obtain Cholesky decomposition of U, inverse, Cholesky decomposition and determinant of inverse of U
    bool fastupdate;
    double ct= ggm->prCoef_lambda + ggm->S.at(colid, colid);    
    double logdetUinv= 0;
    arma::mat Uinv(npar,npar);

    fastupdate= (npar > 1) && (cholU_old != nullptr);
   
    if (fastupdate) {

      int row_added, row_dropped, modelrow_dropped, modelrow_added;
      double *U_newcol;

      modelupdate_indexes(&row_dropped, &row_added, &modelrow_dropped, &modelrow_added, modelold, model);    
      if (modelrow_dropped >= colidint) row_dropped--;
      if (modelrow_added >= colidint) { row_added--; modelrow_added--; }

      if ((row_added == -1) && (row_dropped == -1)) fastupdate= false;

      if (fastupdate) {

        unsigned int idx;
        if (row_added != -1) { //if a row was added, copy the U entries onto U_newcol
          U_newcol= dvector(0, cholUinv->n_cols);
          for (it= model->begin(), i=0; it != model->end(); ++it) {
            idx= it.row();
            if (idx == colid) { continue; } else if (idx > colid) { idx--; }
            U_newcol[i]= ct * Omegainv->at(idx, modelrow_added);
            i++;
          }
          U_newcol[row_added] += tauinv;
        }

        if ((row_added == -1) && (row_dropped != -1)) {  //model drops a row from modelold

          choldcinv_det_droprow(&Uinv, cholUinv, &logdetUinv, cholU, cholU_old, row_dropped);
        
        } else if ((row_added != -1) && (row_dropped == -1)) {  //model adds a row from modelold

          choldcinv_det_addrow(&Uinv, cholUinv, &logdetUinv, cholU, cholU_old, U_newcol, row_added);

        } else { //model swaps rows (adds one row and removes another) relative to modelold

          arma::mat choltmp(npar-1, npar-1);
          choldcinv_det_droprow(nullptr, nullptr, nullptr, &choltmp, cholU_old, row_dropped);
          choldcinv_det_addrow(&Uinv, cholUinv, &logdetUinv, cholU, &choltmp, U_newcol, row_added);

        }

        if (row_added != -1) free_dvector(U_newcol, 0, cholUinv->n_cols);

      }

    }

    if (!fastupdate) {

      //Create U matrix
      arma::mat Omegainv_model(npar, npar); 
      get_Omegainv_model(&Omegainv_model, Omegainv, model, colid); //copy Omegainv[model,model] onto Omegainv_model (excluding colid)
      for (i=0; i<npar; i++) {
        U.at(i,i)= ct * Omegainv_model.at(i,i) + tauinv;
        for (j=i+1; j<npar; j++) U.at(i,j)= U.at(j,i)= ct * Omegainv_model.at(i,j);
      }

      choldcinv_det(&Uinv, cholUinv, &logdetUinv, cholU, &U);

    } 

    /* CHUNK TO DEBUG: PRINTS APPROXIMATION ERROR FOR FAST CHOLESKY UPDATE, IF EVER ABOVE 1.e-4 */
    /*if (fastupdate) {
      arma::mat Omegainv_model(npar, npar); 
      get_Omegainv_model(&Omegainv_model, Omegainv, model, colid); //copy Omegainv[model,model] onto Omegainv_model (excluding colid)
      arma::mat Uinv_slow(npar,npar);
      for (i=0; i<npar; i++) {     //Create U matrix
        U.at(i,i)= ct * Omegainv_model.at(i,i) + tauinv;
        for (j=i+1; j<npar; j++) U.at(i,j)= U.at(j,i)= ct * Omegainv_model.at(i,j);
      }
      choldcinv_det(&Uinv_slow, cholUinv, &logdetUinv, cholU, &U);
      double approxerror= norm(Uinv - Uinv_slow, "fro"); 
      //if (approxerror > 1.e-4) {
        Rprintf("GGMrow_marg approx error= %f\n", approxerror);
        modelold->print("modelold");
        model->print("model");
        Omegainv->print("Omegainv");
        cholU_old->print("cholUold");
        Uinv.print("Uinv (fast update)"); 
        Uinv_slow.print("Uinv (slow update)");
      //}
    }
    */
    
   
    (*m)= - (Uinv * s);
    (*logjoint) = 0.5 * (arma::as_scalar(s.t() * Uinv * s) - ((double) npar) * log( ggm->prCoef_tau ) + logdetUinv);
    //(*logjoint) = 0.5 * (arma::as_scalar(m->t() * U * (*m)) - ((double) npar) * log( ggm->prCoef_tau ) + logdetUinv); //same using U instead of Uinv

  } else { //if all off-diagonal entries are zero

    (*logjoint) = - ((double) npar) * log( ggm->prCoef_tau );

  }

  (*logjoint) += logprior_GGM(&model_offdiag, ggm);

}

// Return Omegainv[model,model], dropping column colid
void get_Omegainv_model(arma::mat *Omegainv_model, arma::mat *Omegainv, arma::SpMat<short> *model, unsigned int colid) {

  unsigned int npar= model->n_nonzero -1;
  if (Omegainv_model->n_cols != npar) Rf_error("Error in get_Omegainv_model: Omegainv_model has the wrong size");
  arma::SpMat<short>::iterator it;

  arma::SpMat<short> model_offdiag= *model;
  model_offdiag.shed_row(colid);  //remove row colid from model

  copy_submatrix(Omegainv_model, Omegainv, &model_offdiag);

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
  - cholXtX_old: Cholesky decomposition of X^T X for a previous model (modelold) that differs from model by adding/dropping at most 1 non-zero entry
  - modelold: entries that were non-zero in the previous model

  OUTPUT
  - logjoint: log-marginal + log-prior for model
  - m: ignored, kept for compatibility with GGMrow_marg
  - cholUinv: ignored, kept for compatibility with GGMrow_marg
  - cholXtX: Cholesky decomposition of X^T X, where X contains the columns of y selected by model

*/
void GGMrow_marg_regression(double *logjoint, arma::mat *m, arma::mat *cholUinv, arma::mat *cholXtX, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm, arma::mat *Omegainv_model, arma::mat *cholXtX_old=nullptr, arma::SpMat<short> *modelold= nullptr) {

  unsigned int npar= model->n_nonzero -1;
  arma::SpMat<short> model_offdiag(ggm->ncol -1, 1);
  double alphahalf= 0.5; //prior on diagonal entries is Gamma(alpha=1, lambda/2)
  double num, den, sumy2= ggm->S.at(colid, colid), ytXVXty;

  if (npar > 0) { //if some off-diagonal entry is non-zero

    unsigned int i;
    int colidint= (int) colid;
    double tauinv= 1.0/ ggm->prCoef_tau, logdetXtXinv, nuhalf, ss;
    arma::uvec covariate_indexes(npar);
    arma::vec Xty(npar);
    arma::mat XtX(npar,npar), cholXtXinv(npar,npar);
    arma::SpMat<short>::iterator it;

    //Create XtX and Xty sufficient statistic matrices for model specified by model_offdiag
    //Note: model_offdiag indicates non-zero off-diagonal elements
    for (it= model->begin(), i=0; it != model->end(); ++it) {
      if (it.row() == colid) continue;
      model_offdiag.at(i, 0)= model->at(it.row(), 0);
      Xty.at(i)= ggm->S.at(it.row(), colid);
      covariate_indexes[i]= it.row();
      i++;
    }

    //Obtain Cholesky decomposition of XtX, inverse, Cholesky decomposition and determinant of inverse of XtX
    bool fastupdate;
    fastupdate= (npar > 1) && (cholXtX_old != nullptr);

    if (fastupdate) {

      int row_added, row_dropped, modelrow_dropped, modelrow_added;
      double *XtX_newcol;

      modelupdate_indexes(&row_dropped, &row_added, &modelrow_dropped, &modelrow_added, modelold, model);
      if (modelrow_dropped >= colidint) row_dropped--;
      if (modelrow_added >= colidint) row_added--;

      if ((row_added == -1) && (row_dropped == -1)) fastupdate= false;

      if (fastupdate) {

        if (row_added != -1) { //if a row was added, copy the XtX entries
          XtX_newcol= dvector(0, cholXtXinv.n_cols);
          for (it= model->begin(), i=0; it != model->end(); ++it) {
            if (it.row() == colid) continue;
            XtX_newcol[i]= ggm->S.at(it.row(), modelrow_added);
            i++;
          }
          XtX_newcol[row_added] += tauinv;
        }

        if ((row_added == -1) && (row_dropped != -1)) {  //model drops a row from modelold

          choldcinv_det_droprow(nullptr, &cholXtXinv, &logdetXtXinv, cholXtX, cholXtX_old, row_dropped);
        
        } else if ((row_added != -1) && (row_dropped == -1)) {  //model adds a row from modelold

          choldcinv_det_addrow(nullptr, &cholXtXinv, &logdetXtXinv, cholXtX, cholXtX_old, XtX_newcol, row_added);

        } else { //model swaps rows (adds one row and removes another) relative to modelold

          arma::mat choltmp(npar-1, npar-1);
          choldcinv_det_droprow(nullptr, nullptr, nullptr, &choltmp, cholXtX_old, row_dropped);
          choldcinv_det_addrow(nullptr, &cholXtXinv, &logdetXtXinv, cholXtX, &choltmp, XtX_newcol, row_added);

        }

        if (row_added != -1) free_dvector(XtX_newcol, 0, cholXtXinv.n_cols);

      }

    }

    if (!fastupdate) {
     
      XtX = ggm->S.submat(covariate_indexes, covariate_indexes);
      for (i=0; i<npar; i++) XtX.at(i,i) += tauinv;
      choldcinv_det(nullptr, &cholXtXinv, &logdetXtXinv, cholXtX, &XtX);  //Cholesky decomposition and determinant of XtX^{-1} 

    }

    //Compute marginal likelihood    
    arma::vec r(npar);
    chol_times_vec(&cholXtXinv, &Xty, &r); //r= cholXtXinv Xty, so that r^T r= Xty^T XtXinv Xty
    for (i=0, ytXVXty=0; i < npar; i++) ytXVXty += r.at(i) * r.at(i);

    ss= ggm->prCoef_lambda + sumy2 - ytXVXty;
    //ss= ggm->prCoef_lambda + sumy2 - arma::as_scalar(Xty.t() * XtXinv * Xty);
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



/************************************************************************************************

 MATRIX MANIPULATION

************************************************************************************************/

/* Check if sparse matrices A and B have <= maxdif different non-zero entries

INPUT
- A: Pointer to the first sparse matrix.
- B: Pointer to the second sparse matrix.
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
