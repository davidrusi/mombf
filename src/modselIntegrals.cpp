#include "modselIntegrals.h"
#include "cstat.h"
using namespace std;

modselIntegrals::modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars) {

  this->maxVars= nvars;
  this->marginalFunction= marfun;
  this->priorFunction= priorfun;

  this->selchar = (char *) calloc(nvars, sizeof(char));
  //this->selchar= charvector(0,maxVars-1);
}

modselIntegrals::~modselIntegrals() {

  free((char  *) this->selchar);
  //free_charvector(this->selchar, 0, maxVars-1);

}

//Return log(marginal likelihood) + log(prior). Uses logjointSaved if available, else adds result to logjointSaved
// Input: 
//
//   - sel: integer vector [0..maxVars-1], where 0's and 1's indicate covariates out/in the model (respectively)
//   - nsel: number of covariates in the model (i.e. sum(sel))
//   - pars: struct of type marginalPars containing parameters needed to evaluate the marginal density of the data & prior on model space
//
// Output: evaluates log joint. It returns previously saved results in logjointSaved if available, else it performs the computation and saves the result in logjointSaved
double modselIntegrals::getJoint(int *sel, int *nsel, struct marginalPars *pars) {
  int i;
  for (i=0; i<= maxVars; i++) if (sel[i]==0) this->selchar[i]= '0'; else this->selchar[i]= '1';
  
  if (logjointSaved.count(selchar) > 0) return logjointSaved[selchar];

  double ans;
  ans= marginalFunction(sel,nsel,pars) + priorFunction(sel,nsel,pars);
  logjointSaved[selchar]= ans;
  return ans;
} 
