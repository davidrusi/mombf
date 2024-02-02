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
  this->nvars= ggm->ncol() - 1;

  this->jointFunction= jointFunction;

  this->ggm= ggm;

  this->colid= colid;

  this->Omegainv= Omegainv;

  this->maxIntegral= -1.0e250;

  this->maxsave= 1000000000; //save first 10^9 models

  this->zerochar = (char *) calloc(this->nvars + 2, sizeof(char));
  for (i=0; i < this->nvars; i++) this->zerochar[i]= '0';

}

//Class destructor
modselIntegrals_GGM::~modselIntegrals_GGM() {

  free((char  *) this->zerochar);

  std::map<string, arma::mat *>::iterator it;
  for (it= meanSaved.begin(); it != meanSaved.end(); ++it) delete it->second;
  for (it= cholVSaved.begin(); it != cholVSaved.end(); ++it) delete it->second;

}


/* Log-integrated likelihood + log-prior for model. 

If postSample=true, sample_offdiag returns a posterior for the off-diagonal parametes and sample_diag for the diagonal parameters

*/
void modselIntegrals_GGM::getJoint(double *logjoint, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, bool postSample) {

  bool delete_m_cholV= false;  
  int npar= model->n_nonzero -1;
  arma::mat *m, *cholV;
  arma::SpMat<short>::iterator it;
  arma::mat Omegainv_model(npar, npar); 

  //Set zerochar to current model
  for (it= model->begin(); it != model->end(); ++it) zerochar[it.row()]= '1';
  std::string s (zerochar);

  //Copy entries of Omegainv selected by model to Omegainv_model
  //arma::mat Omegainv_model= this->get_Omegainv_model(model);

  if (logjointSaved.count(s) > 0) {  //if logjoint already computed in a previous call

    (*logjoint)= logjointSaved[s];
    m= meanSaved[s];
    cholV= cholVSaved[s];
     
  } else {

    //Copy entries of Omegainv selected by model to Omegainv_model
    this->get_Omegainv_model(&Omegainv_model, model);

    //Allocate memory for m and cholV
    m= new arma::mat(npar, 1);
    cholV= new arma::mat(npar, npar);

    jointFunction(logjoint, m, cholV, model, colid, this->ggm, &Omegainv_model);

    //Store logjoint, m and cholV
    double d= maxIntegral - (*logjoint);
    if (d<15 || this->nvars<=16 || logjointSaved.size() <= maxsave) {
      logjointSaved[s]= *logjoint;
      meanSaved[s]= m; 
      cholVSaved[s]= cholV;      
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
    arma::vec lambda= as<arma::vec>(ggm->prCoef["lambda"]); //Prior is Omega_{jj} ~ Exp(lambda)
    double a= 0.5 * (double) ggm->n() + 1.0;
    double b= 0.5 * (ggm->S).at(this->colid, this->colid) + 0.5 * lambda[0];
     
    (*sample_diag)= rgammaC(a, b) + arma::as_scalar(sample_offdiag->t() * Omegainv_model  * (*sample_offdiag));

  }

  //Free memory
  if (delete_m_cholV) {
    delete m;
    delete cholV;
  }
  
  //Return zerochar to its original empty model status
  for (it= model->begin(); it != model->end(); ++it) zerochar[it.row()]= '0';

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

