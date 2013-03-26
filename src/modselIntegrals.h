//#ifndef MODSELINTEGRALS
//#define MODSELINTEGRALS 1
#include "modelSel.h"


class Seppel

{

public:

	Seppel();

	~Seppel();

	double calcIntegral(double * model);  //compute integral (requires computing mode)

private:

	double *aa;

};


class modselIntegrals 

{

public:

  modselIntegrals(pt2margFun fun, int nvars);  //initialize logjoint to fun, maxVars to nvars
  
  ~modselIntegrals();

  double getJoint(int *sel, int *nsel, struct marginalPars *pars); //Return logjoint(). Uses logjointSaved if available, else adds result to logjointSaved

private:

  int maxVars; //Maximum number of covariates
  char *selchar;  //Store model id (vars in the model) in character format, e.g. "10100"
  pt2margFun logjoint;  //Function computing log(marginal likelihood) + log(prior)
  std::map<char*, double> logjointSaved; //Saves previously computed logjoint

};

//#endif

