#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#Rcpp::compileAttributes("~/github/mombf")

library(mombf)
library(mvtnorm)
set.seed(1)

## GGM examples
p= 5
Th= diag(p); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)
z= matrix(rnorm(1000*p), ncol=p)
y= z %*% chol(sigma)
Omegaini= Th
fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini, almost_parallel="none", sampler='LIT', niter=100)

#Th= diag(20)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.9
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.37
#y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)

#niter= 10^4; burnin= 0.2 * niter
#fitr.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=solve(cov(y)), niter=niter, burnin=burnin, scale=FALSE, almost_parallel='regression', tempering=1, truncratio=100, save_proposal=TRUE, prob_parallel=1)

#Hard example p=5. Smallest eigenvalue very close to 0
#Th= diag(5)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#y= scale(rmvnorm(100, sigma=solve(Th)), center=TRUE, scale=FALSE)
#Omegaini= Th
#fit <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=Omegaini, niter=100, burnin=0, updates_per_iter=ncol(y), updates_per_column=ncol(Th), scale=FALSE, almost_parallel='regression', tempering=1, truncratio=100)

