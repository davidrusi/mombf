#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#library(Rcpp)
#compileAttributes("~/github/mombf")

library(mombf)
library(mvtnorm)
library(survival)
set.seed(1)

#Cholesky update
#B= matrix(c(2,1,.5,.25, 1,1.5,.3,.4, .5,.3,1,.75, .25,.4,.75,2),nrow=4,byrow=TRUE)
#mombf:::testfunction(B, 3, 3)

## GGM valgrind example
#p= 5
#Th= diag(p); Th[1,2]= Th[2,1]= 0.5
#sigma= solve(Th)
#z= matrix(rnorm(1000*p), ncol=p)
#y= z %*% chol(sigma)
#Omegaini= Th
#fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini, almost_parallel="none", sampler='birthdeath', niter=100)

#Th= diag(20)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.9
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.37
#y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)

#niter= 10^4; burnin= 0.2 * niter
#fitr.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=solve(cov(y)), niter=niter, burnin=burnin, scale=FALSE, almost_parallel='regression', tempering=1, truncratio=100, save_proposal=TRUE, prob_parallel=1)

#Hard example p=5. Smallest eigenvalue very close to 0
Th= diag(5)
diag(Th)= 1.5
Th[abs(col(Th) - row(Th))==1]= 0.95
Th[abs(col(Th) - row(Th))==2]= 0.5
Th[abs(col(Th) - row(Th))==3]= 0.64
y= scale(rmvnorm(100, sigma=solve(Th)), center=TRUE, scale=FALSE)
Omegaini= Th
fit <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=Omegaini, niter=100, burnin=0, updates_per_iter=ncol(y), updates_per_column=ncol(Th), scale=FALSE, almost_parallel='none', tempering=1, truncratio=100)

