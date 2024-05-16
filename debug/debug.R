#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()

#library(Rcpp)
#compileAttributes("~/github/mombf")
library(mombf)
library(mvtnorm)
set.seed(1)
#p= 5
#Th= diag(p); Th[1,2]= Th[2,1]= 0.5
#sigma= solve(Th)
#z= matrix(rnorm(1000*p), ncol=p)
#y= z %*% chol(sigma)
#fit= modelSelectionGGM(y, scale=FALSE, almost_parallel="regression", sampler='birthdeath', niter=10)

#Hard example. Smallest eigenvalue very close to 0
Th= diag(50)
diag(Th)= 1.5
Th[abs(col(Th) - row(Th))==1]= 0.9
Th[abs(col(Th) - row(Th))==2]= 0.5
Th[abs(col(Th) - row(Th))==3]= 0.35
y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)
fit <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-bic', niter=5000, burnin=0, scale=FALSE, almost_parallel='regression', nbirth=1, tempering=1, save_proposal=TRUE)
#fitr.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-bic', niter=5000, burnin=0, scale=FALSE, almost_parallel='regression', nbirth=1, tempering=1, save_proposal=TRUE)
