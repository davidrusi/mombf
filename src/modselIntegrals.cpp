#include "modselIntegrals.h"
using namespace std;


/***********************************************************************************/
/* Integrated likelihoods for regression models                                    */
/***********************************************************************************/


modselIntegrals::modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars) {
  int i;

  this->maxVars= nvars;
  this->marginalFunction= marfun;
  this->priorFunction= priorfun;

  this->maxIntegral= -1.0e250;

  this->maxsave= 1000000000; //save first 10^9 models

  this->zerochar = (char *) calloc(nvars+1, sizeof(char));
  for (i=0; i<nvars; i++) this->zerochar[i]= '0';

}

modselIntegrals::~modselIntegrals() {

  free((char  *) this->zerochar);

}

//Return log(marginal likelihood) + log(prior). Uses logjointSaved if available, else adds result to logjointSaved. When maxVars>16, only models with a log-difference <=10 with the current mode are stored
// Input:
//
//   - sel: integer vector [0..maxVars-1], where 0's and 1's indicate covariates out/in the model (respectively)
//   - nsel: number of covariates in the model (i.e. sum(sel))
//   - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
//
// Output: evaluates log joint. It returns previously saved results in logjointSaved if available, else it performs the computation and saves the result in logjointSaved
double modselIntegrals::getJoint(int *sel, int *nsel, struct marginalPars *pars) {
  int i;
  double ans;

  if (*nsel > *((*pars).maxvars)) {
    ans= -INFINITY;
  } else {

    for (i=0; i< *nsel; i++) zerochar[sel[i]]= '1';
    std::string s (zerochar);

    if (logjointSaved.count(s) > 0) {
      ans= logjointSaved[s];
    } else {
      ans= marginalFunction(sel,nsel,pars);
      ans+= priorFunction(sel,nsel,pars);
      double d= maxIntegral - ans;
      if (d<10 || maxVars<=16 || logjointSaved.size() <= maxsave) logjointSaved[s]= ans;
      if (d<0) {
        maxIntegral= ans;
        maxModel= s;
      }
    }

    for (i=0; i<= *nsel; i++) this->zerochar[sel[i]]= '0';
  }

  return ans;
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

  //std::map<string, arma::mat *>::iterator it;
  //for (it= meanSaved.begin(); it != meanSaved.end(); ++it) delete it->second;
  //for (it= cholVSaved.begin(); it != cholVSaved.end(); ++it) delete it->second;

}



/* The log-integrated likelihood + log-prior for model is returned in logjoint

If postSample=true, sample_offdiag returns a posterior sample for the off-diagonal parametes and sample_diag for the diagonal parameters

*/
void modselIntegrals_GGM::getJoint(double *logjoint, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, arma::SpMat<short> *modelold, bool postSample) {

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

    //Sample off-diagonal elements
    rmvnormC(sample_offdiag, m, cholV);
    (*sample_offdiag) *= -1.0;
     
    //Sample diagonal element
    double a= 0.5 * (double) ggm->n + 1.0;
    double b= 0.5 * (ggm->S).at(this->colid, this->colid) + 0.5 * ggm->prCoef_lambda; //Prior is Omega_{jj} ~ Exp(lambda)
     
    (*sample_diag)= rgammaC(a, b) + arma::as_scalar(sample_offdiag->t() * Omegainv_model  * (*sample_offdiag));

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

