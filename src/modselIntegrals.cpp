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
    //Rprintf("marginal=%f, prior=%f\n",ans,priorFunction(sel,nsel,pars));
    ans+= priorFunction(sel,nsel,pars);
    double d= maxIntegral - ans;
    if (d<10 || maxVars<=16) logjointSaved[s]= ans;
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


modselIntegrals_GGM::modselIntegrals_GGM(pt2GGM_rowmarg marfun, ggmObject *ggm, unsigned int colid, arma::mat *Omegainv) {
  int i;
  this->nvars= ggm->ncol() - 1;

  this->jointFunction= marfun;

  this->Omegainv= Omegainv;

  this->maxIntegral= -1.0e250;

  this->zerochar = (char *) calloc(this->nvars + 1, sizeof(char));
  for (i=0; i < this->nvars; i++) this->zerochar[i]= '0';

}

modselIntegrals_GGM::~modselIntegrals_GGM() {

  free((char  *) this->zerochar);

  std::map<string, arma::mat *>::iterator it;
  for (it= meanSaved.begin(); it != meanSaved.end(); ++it) delete it->second;
  for (it= cholVSaved.begin(); it != cholVSaved.end(); ++it) delete it->second;

}



void modselIntegrals_GGM::getJoint(double *logjoint, arma::mat *sample_offdiag, double *sample_diag, arma::SpMat<short> *model, unsigned int colid, ggmObject *ggm) {
  
  arma::mat *m, *cholV;
  arma::SpMat<short>::iterator it;

  for (it= model.begin(); it != model.end; ++it) zerochar[it.nrow()]= '1';
  std::string s (zerochar);

  if (logjointSaved.count(s) > 0) {

    (*logjoint)= logjointSaved[s];
    m= meanSaved[s];
    cholV= cholVsaved[s];
     
  } else {

    //Allocate memory for m and cholV
    *m= new arma::mat(model->n_rows, 1);
    *cholV= new arma::mat(model->n_rows, model->n_rows);

    jointFunction(logjoint, m, cholV, model, colid, ggm);

    //TO DO: compute sample_offdiag, sample_diag


    //Store logjoint, m and cholV
    double d= maxIntegral - (*logjoint);
    if (d<15 || this->nvars<=16) {
      logjointSaved[s]= *logjoint;
      meanSaved[s]= m; 
      cholVSaved[s]= cholV;      
    } else {  //if not stored, free the allocated memory
      delete m;
      delete cholV;
    }

    //Update top model
    if (d<0) {
      maxIntegral= ans;
      maxModel= s;
    }
  }

  for (it= model.begin(); it != model.end; ++it) zerochar[it.nrow()]= '0';

  return ans;
}
