#ifndef NLPMARG_H
#define NLPMARG_H 1

//#include <R.h>
//#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include "crossprodmat.h"
#include "Polynomial.h"


/*
 * Function Prototypes
 */
double nlpMarginal(int *sel, int *nsel, struct marginalPars *pars);

#endif /* NLPMARG_H */
