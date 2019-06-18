#include "modselIntegrals.h"
using namespace std;

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
