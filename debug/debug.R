#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()

#library(Rcpp)
#compileAttributes("~/github/mombf")
library(mombf)
library(mvtnorm)
set.seed(1)
p= 5
Th= diag(p); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)
z= matrix(rnorm(1000*p), ncol=p)
y= z %*% chol(sigma)
 
#Obtain posterior samples
fit= modelSelectionGGM(y, scale=FALSE, almost_parallel="regression", sampler='birthdeath', niter=10^4)

#Compare marginal likelihood
#nlpMarginal(sel=c(TRUE,TRUE,TRUE,TRUE), y=y[,1], x=y[,-1], priorCoef=normalidprior(tau=1), priorVar=igprior(1, 0.5))