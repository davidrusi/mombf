#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#Rcpp::compileAttributes("~/github/mombf")

library(mombf)
library(mvtnorm)
set.seed(1)

## GGM examples
#p= 5
#Th= diag(p); Th[1,2]= Th[2,1]= 0.5
#sigma= solve(Th)
#z= matrix(rnorm(1000*p), ncol=p)
#y= z %*% chol(sigma)
#Omegaini= Th
#fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini, almost_parallel="none", sampler='Gibbs', niter=100)

p <- 3
n <- 5000
set.seed(1)
Omega= diag(p)
Omega[abs(col(Omega) - row(Omega))==1]= 0.5
y = rmvnorm(n = n, sigma = solve(Omega))

niter= 100
burnin= 0
updates_per_iter= p; updates_per_column= p
pbirth= 0.75; pdeath= 1-pbirth
#pbirth= 0.75; pdeath= 0.5*(1-pbirth)

fitbd= modelSelectionGGM(y, sampler='Gibbs', niter=niter, Omegaini=Omega, global_proposal="none", updates_per_iter=updates_per_iter, updates_per_column=updates_per_column, pbirth= pbirth, pdeath=pdeath)
