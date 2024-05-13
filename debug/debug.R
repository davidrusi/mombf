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
#y= matrix(c(1,1,1,0.5),nrow=2)
#y= scale(y, center=TRUE, scale=FALSE)

#Obtain posterior samples
Omegaini= solve(cov(y))
fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini, almost_parallel="regression", sampler='birthdeath', niter=10^4)

