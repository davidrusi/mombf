#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()

#library(Rcpp)
#compileAttributes("~/github/mombf")
library(mombf)
library(mvtnorm)
set.seed(1)
Th= diag(4)
Th[abs(col(Th) - row(Th))==1]= 0.5
y= scale(rmvnorm(10^2, sigma=solve(Th)), center=TRUE, scale=FALSE)
fit.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=1000, burnin=0, scale=FALSE, almost_parallel=TRUE) 
