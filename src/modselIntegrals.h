#ifndef MODSELINTEGRALS
#define MODSELINTEGRALS 1

#include <map>
#include "modelSel.h"
#include "cstat.h"
using namespace std;


class modselIntegrals {

public:

  modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars);  //initialize logjoint to fun, maxVars to nvars
  
  ~modselIntegrals();

  double getJoint(int *sel, int *nsel, struct marginalPars *pars); //Return logjoint(). Uses logjointSaved if available, else adds result to logjointSaved

private:

  int maxVars; //Maximum number of covariates
  char *selchar;  //Store model id (vars in the model) in character format, e.g. "10100"
  pt2margFun marginalFunction;  //Function computing log(marginal likelihood)
  pt2margFun priorFunction;     //Function computing log(model prior)
  std::map<char*, double> logjointSaved; //Saves previously computed logjoint

};

#endif

