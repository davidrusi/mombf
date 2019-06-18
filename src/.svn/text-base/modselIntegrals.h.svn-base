#ifndef MODSELINTEGRALS
#define MODSELINTEGRALS 1

#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include <map>
#include <string>
#include "modelSel.h"
#include "cstat.h"
using namespace std;


class modselIntegrals {

public:

  modselIntegrals(pt2margFun marfun, pt2margFun priorfun, int nvars);  //initialize logjoint to fun, maxVars to nvars
  
  ~modselIntegrals();

  double getJoint(int *sel, int *nsel, struct marginalPars *pars); //Return logjoint(). Uses logjointSaved if available, else adds result to logjointSaved

  double maxIntegral; //Stores value of largest integral

  string maxModel; //Stores model with largest integral, e.g. "10001" 

private:

  int maxVars; //Maximum number of covariates
  char *zerochar;  //Store model id (vars in the model) in character format, e.g. "00000"
  pt2margFun marginalFunction;  //Function computing log(marginal likelihood)
  pt2margFun priorFunction;     //Function computing log(model prior)
  std::map<string, double> logjointSaved; //Saves previously computed logjoint

};

#endif

