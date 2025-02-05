#include "modselIntegrals.h"
using namespace std;


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
  this->niter= as<int>(samplerPars["niter"]);
  this->updates_per_column= as<int>(samplerPars["updates_per_column"]);
  this->updates_per_iter= as<int>(samplerPars["updates_per_iter"]);
  this->tempering= as<double>(samplerPars["tempering"]);
  this->truncratio= as<double>(samplerPars["truncratio"]);
  this->pbirth= as<double>(samplerPars["pbirth"]);
  this->pdeath= as<double>(samplerPars["pdeath"]);
  this->pswap= 1 - this->pbirth - this->pdeath;
  this->log_pbirth= log(this->pbirth);
  this->log_pdeath= log(this->pdeath);
  this->log_pswap= log(this->pswap);
  this->use_tempering= use_tempering;

  this->lbound_death= as<double>(samplerPars["lbound_death"]);
  this->ubound_death= as<double>(samplerPars["ubound_death"]);
  this->lbound_birth= as<double>(samplerPars["lbound_birth"]);
  this->ubound_birth= as<double>(samplerPars["ubound_birth"]);

  this->prob_global= as<double>(samplerPars["prob_global"]);
  CharacterVector global_proposalR= samplerPars["global_proposal"];
  std::string global_proposalC = Rcpp::as<std::string>(global_proposalR);
  std::string regression("regression"), insample("in-sample");
  this->global_regression= (global_proposalC == regression);
  this->global_insample= (global_proposalC == insample);

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
  this->niter= ggm->niter;
  this->updates_per_column= ggm->updates_per_column;
  this->updates_per_iter= ggm->updates_per_iter;
  this->tempering= ggm->tempering;
  this->truncratio= ggm->truncratio;
  this->pbirth= ggm->pbirth;
  this->pdeath= ggm->pdeath;
  this->pswap= ggm->pswap;
  this->log_pbirth= ggm->log_pbirth;
  this->log_pdeath= ggm->log_pdeath;
  this->log_pswap= ggm->log_pswap;
  this->use_tempering= ggm->use_tempering;

  this->lbound_death= ggm->lbound_death;
  this->ubound_death= ggm->ubound_death;
  this->lbound_birth= ggm->lbound_birth;
  this->ubound_birth= ggm->ubound_birth;

  this->prob_global= ggm->prob_global;
  this->global_regression= ggm->global_regression;
  this->global_insample= ggm->global_insample;

  //Set print progress iteration to true/false
  this->verbose= ggm->verbose;

}

//Destructor
ggmObject::~ggmObject() {

}



/***********************************************************************************/
/* Integrated likelihoods for regression models                                    */
/***********************************************************************************/


modselIntegrals::modselIntegrals(pt2margFun marfun, pt2modelpriorFun priorfun, int nvars, lmObject *lm) {
  int i;

  this->maxVars= nvars;
  this->marginalFunction= marfun;
  this->priorFunction= priorfun;
  this->lm= lm;

  this->maxIntegral= -INFINITY;

  this->maxsave= 1000000000; //save first 10^9 models

  this->zerochar = (char *) calloc(nvars+1, sizeof(char));
  for (i=0; i<nvars; i++) this->zerochar[i]= '0';

}

modselIntegrals::~modselIntegrals() {

  free((char  *) this->zerochar);

  std::map<string, LM_logjoint>::iterator it;
  for (it= logjointSaved.begin(); it != logjointSaved.end(); ++it) {
    delete (it->second).mean;
    delete (it->second).cholVinv;
  }


}

//Return log(marginal likelihood) + log(prior). Uses logjointSaved if available, else adds result to logjointSaved. When maxVars>16, only models with a log-difference <=10 with the current mode are stored
// Input:
//
//   - sel: integer vector [0..maxVars-1], where 0's and 1's indicate covariates out/in the model (respectively)
//   - nsel: number of covariates in the model (i.e. sum(sel))
//   - lm: lmObject containing parameters needed to evaluate the marginal density of the data & prior on model space
//
// Output: evaluates log joint. It returns previously saved results in logjointSaved if available, else it performs the computation and saves the result in logjointSaved
double modselIntegrals::getJoint(arma::SpMat<short> *model, lmObject *lm, arma::SpMat<short> *modelold= nullptr) {
  bool delete_m_cholV= false;
  int nsel= model->n_nonzero, npar;
  double ans, d;
  arma::mat *m= nullptr, *cholVinv= nullptr, *cholVinv_old= nullptr;
  struct LM_logjoint newlogjoint;

  if ((*(lm->maxvars) >= 0) && (nsel > *lm->maxvars)) {
    ans= -INFINITY;
  } else {

    //Represent the model as a string
    std::string s= this->getModelid(model);
    if (*lm->family == 0) {
      int familycode= model->at(model->n_rows - 1, 0);
      s.back() = static_cast<char>(familycode + '0'); // note: assuming 'familycode' is a single digit number
      npar= nsel - 1;
    } else {
      npar= nsel;
    }

    if (logjointSaved.count(s) > 0) { //if model was previously saved

      ans= (logjointSaved[s]).logjoint;

    } else {

      //Retrieve Cholesky decomposition for modelold, if available
      if (modelold != nullptr) {
        std::string sold= this->getModelid(modelold);
        if (*lm->family == 0) {
          int familycode= modelold->at(modelold->n_rows - 1, 0);
          sold.back() = static_cast<char>(familycode + '0'); // note: assuming 'familycode' is a single digit number
        }
        if (logjointSaved.count(sold) > 0) {
          cholVinv_old= (logjointSaved[sold]).cholVinv;
        } else {
          cholVinv_old= nullptr;
        }
      }

      //Allocate memory for m and cholV
      if (npar > 0) {
        m= new arma::mat(npar, 1);
        cholVinv= new arma::mat(npar,npar);
      }

      ans= marginalFunction(model, lm, cholVinv_old, modelold, m, cholVinv);
      ans+= priorFunction(model, lm);
      d= maxIntegral - ans;
      if (d<10 || maxVars<=16 || logjointSaved.size() <= maxsave) {

        newlogjoint.logjoint= ans;
        newlogjoint.mean= m;
        newlogjoint.cholVinv= cholVinv;
        logjointSaved[s]= newlogjoint;

      } else {  //if not stored, free the allocated memory
        delete_m_cholV= true;
      }

      if (d<0) {
        maxIntegral= ans;
        maxModel= s;
      }

      //Free memory
      if (delete_m_cholV) {
        delete m;
        delete cholVinv;
      }

    }

  }

  return ans;
}


std::string modselIntegrals::getModelid(arma::SpMat<short> *model) {
  if (model->n_cols !=1) Rf_error("In getModelid, argument model must have 1 column");
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) this->zerochar[it.row()]= '1'; //Set zerochar to current model
  std::string s (this->zerochar);
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) this->zerochar[it.row()]= '0'; //Return zerochar to its original empty model status
  return s;
}



/************************************************************************************/
/* Integrated likelihoods for Gaussian graphical models with precision matrix Omega */
/*                                                                                  */
/* Models are defined by non-zero entries in column cold, and are conditional on    */
/* a given value of Omegainv, the inverse of Omega[-colid,-colid]                   */
/************************************************************************************/


//Class constructor
modselIntegrals_GGM::modselIntegrals_GGM(pt2GGM_rowmarg jointFunction, ggmObject *ggm, unsigned int colid, arma::mat *Omegainv) {

  int i;
  this->nvars= ggm->ncol - 1;

  this->jointFunction= jointFunction;

  this->ggm= ggm;

  this->colid= colid;

  this->Omegainv= Omegainv;

  this->maxIntegral= -1.0e250;

  this->maxsave= 1000000000; //save first 10^9 models

  this->zerochar = (char *) calloc(this->nvars + 2, sizeof(char));
  for (i=0; i <= this->nvars; i++) this->zerochar[i]= '0';

}

//Class destructor
modselIntegrals_GGM::~modselIntegrals_GGM() {

  free((char  *) this->zerochar);

  std::map<string, GGM_logjoint>::iterator it;
  for (it= logjointSaved.begin(); it != logjointSaved.end(); ++it) {
    delete (it->second).mean;
    delete (it->second).cholV;
    delete (it->second).cholVinv;
  }

}



/* The log-integrated likelihood + log-prior for model is returned in logjoint

If postSample=true, sample_offdiag returns a posterior sample for the off-diagonal parametes and sample_diag for the diagonal parameters

*/
void modselIntegrals_GGM::getJoint(double *logjoint, arma::mat *mean_offdiag, double *mean_diag, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, arma::SpMat<short> *modelold, bool postSample) {

  bool delete_m_cholV= false;  
  int npar= model->n_nonzero -1;
  arma::mat *m, *cholV, *cholVinv, *cholVinv_old= nullptr;
  arma::SpMat<short>::iterator it;
  arma::mat Omegainv_model(npar, npar); 
  struct GGM_logjoint *logjointptr, newlogjoint;

  std::string s= this->getModelid(model);

  if (logjointSaved.count(s) > 0) {  //if logjoint already computed in a previous call

    logjointptr= &(logjointSaved[s]);
    (*logjoint)= logjointptr->logjoint;
    m= logjointptr->mean;
    cholV= logjointptr->cholV;
     
  } else {

    //Retrieve Cholesky decomposition for modelold, if available
    if (modelold != nullptr) {
      std::string sold= this->getModelid(modelold);
      if (logjointSaved.count(sold) > 0) {
        cholVinv_old= (logjointSaved[sold]).cholVinv;
      } else {
        cholVinv_old= nullptr;
      }
    }

    //Allocate memory for m and cholV
    m= new arma::mat(npar, 1);
    cholV= new arma::mat(npar, npar);
    cholVinv= new arma::mat(npar, npar);

    jointFunction(logjoint, m, cholV, cholVinv, model, colid, this->ggm, Omegainv, cholVinv_old, modelold);

    if (ggm->use_tempering) (*logjoint) *= (ggm->tempering);

    //Store logjoint, m and cholV
    double d= maxIntegral - (*logjoint);
    if (d<15 || this->nvars<=16 || logjointSaved.size() <= maxsave) {
      newlogjoint.logjoint= *logjoint;
      newlogjoint.mean= m;
      newlogjoint.cholV= cholV;
      newlogjoint.cholVinv= cholVinv;
      logjointSaved[s]= newlogjoint;
      
    } else {  //if not stored, free the allocated memory
      delete_m_cholV= true;
    }

    //Update top model
    if (d<0) {
      maxIntegral= *logjoint;
      maxModel= s;
    }

  }

  if (postSample) {

    //Copy entries of Omegainv selected by model to Omegainv_model
    this->get_Omegainv_model(&Omegainv_model, model);

    if (mean_offdiag != nullptr) (*mean_offdiag)= (*m);

    //Sample off-diagonal elements
    rmvnormC(sample_offdiag, m, cholV);
     
    //Sample diagonal element
    double a= 0.5 * (double) ggm->n + 1.0;
    double b= 0.5 * (ggm->S).at(this->colid, this->colid) + 0.5 * ggm->prCoef_lambda; //Prior is Omega_{jj} ~ Exp(lambda)
    
    double ss= arma::as_scalar(sample_offdiag->t() * Omegainv_model  * (*sample_offdiag));
    (*sample_diag)= rgammaC(a, b) + ss;
    if (mean_diag != nullptr) (*mean_diag)= (a/b) + ss;

  }

  //Free memory
  if (delete_m_cholV) {
    delete m;
    delete cholV;
    delete cholVinv;
  }
  
}


/* The log-integrated likelihood + log-prior for model is returned in logjoint

 sample_offdiag returns the posterior mode for the off-diagonal parametes and sample_diag for the diagonal parameters

*/
void modselIntegrals_GGM::getMode(double *logjoint, arma::mat *mode_offdiag, double *mode_diag, arma::SpMat<short> *model, arma::SpMat<short> *modelold) {

  bool delete_m_cholV= false;  
  int npar= model->n_nonzero -1;
  arma::mat *m, *cholV, *cholVinv, *cholVinv_old= nullptr;
  arma::SpMat<short>::iterator it;
  arma::mat Omegainv_model(npar, npar); 
  struct GGM_logjoint *logjointptr, newlogjoint;

  std::string s= this->getModelid(model);

  if (logjointSaved.count(s) > 0) {  //if logjoint already computed in a previous call

    logjointptr= &(logjointSaved[s]);
    (*logjoint)= logjointptr->logjoint;
    m= logjointptr->mean;
    //(*logjoint)= logjointSaved[s];
    //m= meanSaved[s];
     
  } else {

    //Retrieve Cholesky decomposition for modelold, if available
    if (modelold != nullptr) {
      std::string sold= this->getModelid(modelold);
      if (logjointSaved.count(sold) > 0) {
        cholVinv_old= (logjointSaved[sold]).cholVinv;
      } else {
        cholVinv_old= nullptr;
      }
    }
    
    //Allocate memory for m and cholV
    m= new arma::mat(npar, 1);
    cholV= new arma::mat(npar, npar);
    cholVinv= new arma::mat(npar, npar);

    jointFunction(logjoint, m, cholV, cholVinv, model, colid, this->ggm, Omegainv, cholVinv_old, modelold);

    if (ggm->use_tempering) (*logjoint) *= (ggm->tempering);

    //Store logjoint, m and cholV
    double d= maxIntegral - (*logjoint);
    if (d<15 || this->nvars<=16 || logjointSaved.size() <= maxsave) {

      newlogjoint.logjoint= *logjoint;
      newlogjoint.mean= m;
      newlogjoint.cholV= cholV;
      newlogjoint.cholVinv= cholVinv;
      logjointSaved[s]= newlogjoint;

    } else {  //if not stored, free the allocated memory

      delete_m_cholV= true;

    }

    //Update top model
    if (d<0) {
      maxIntegral= *logjoint;
      maxModel= s;
    }

  }

  //Copy entries of Omegainv selected by model to Omegainv_model
  this->get_Omegainv_model(&Omegainv_model, model);

  //Posterior mode for off-diagonal elements
  (*mode_offdiag)= -(*m);
     
  //Posterior mode for diagonal element
  double a= 0.5 * (double) ggm->n + 1.0;
  double b= 0.5 * (ggm->S).at(this->colid, this->colid) + 0.5 * ggm->prCoef_lambda;
     
  (*mode_diag)= (a - 1)/b + arma::as_scalar(mode_offdiag->t() * Omegainv_model  * (*mode_offdiag));

  //Free memory
  if (delete_m_cholV) {
    delete m;
    delete cholV;
    delete cholVinv;
  }
  
}


std::string modselIntegrals_GGM::getModelid(arma::SpMat<short> *model) {
  if (model->n_cols !=1) Rf_error("In getModelid, argument model must have 1 column");
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) this->zerochar[it.row()]= '1'; //Set zerochar to current model
  std::string s (this->zerochar);
  for (arma::SpMat<short>::iterator it= model->begin(); it != model->end(); ++it) this->zerochar[it.row()]= '0'; //Return zerochar to its original empty model status
  return s;
}


// Return Omegainv[model,model], dropping column colid
void modselIntegrals_GGM::get_Omegainv_model(arma::mat *Omegainv_model, arma::SpMat<short> *model) {

  unsigned int npar= model->n_nonzero -1;
  if (Omegainv_model->n_cols != npar) Rf_error("Error in get_Omegainv_model: Omegainv_model has the wrong size");
  arma::SpMat<short>::iterator it;

  arma::SpMat<short> model_offdiag= *model;
  model_offdiag.shed_row(colid);  //remove row colid from model

  copy_submatrix(Omegainv_model, Omegainv, &model_offdiag);

}

