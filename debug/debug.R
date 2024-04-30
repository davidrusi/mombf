#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()

#library(Rcpp)
#compileAttributes("~/github/mombf")
library(mombf)
library(mvtnorm)
#set.seed(1)
#Th= diag(4)
#Th[abs(col(Th) - row(Th))==1]= 0.5
#y= scale(rmvnorm(10^2, sigma=solve(Th)), center=TRUE, scale=FALSE)
#fit.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=1000, burnin=0, scale=FALSE, almost_parallel=TRUE, tempering=0.5) 


#Valgrind issue example
p= 5
Th= diag(p); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)
z= matrix(rnorm(1000*p), ncol=p)
y= z %*% chol(sigma)
 
#Obtain posterior samples
fit= modelSelectionGGM(y, scale=FALSE, almost_parallel=TRUE, sampler='birthdeath', niter=10^4)
