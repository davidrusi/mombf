#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#library(Rcpp)
#compileAttributes("~/github/mombf")

## GGM examples
library(mombf)
library(mvtnorm)
logprob2prob= function(x) { ans= exp(x - max(x)); ans/sum(ans) }
set.seed(1)
Th= diag(3)
Th[1,2]= Th[2,1]= 0.4
Th[2,3]= Th[3,2]= 0.9
y= scale(rmvnorm(10^3, sigma=solve(Th)), center=TRUE, scale=FALSE)

Omegaini= Th
fit <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=Omegaini, niter=10^4, updates_per_iter=3, updates_per_column=5, scale=FALSE, almost_parallel='none')


#Th= diag(5)
#diag(Th)= c(1.5, 1.6, 1.7, 1.8, 1.9)
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#Thnew= Th
#Thnew[,3]= Thnew[3,]= c(0.4,0.8,1.4,0.75,0.3)
#colid= 3
#mombf:::testfunction(Th, oldcol=1, newcol=2)

#Jack's simulated example
#load('~/Downloads/jack_ggm_simulation.RData')
#fit <- modelSelectionGGM(y, sampler='Gibbs', Omegaini='null', niter=10^4, burnin=0, scale=FALSE, almost_parallel='none')
#fit <- modelSelectionGGM(y, sampler='Gibbs', Omegaini=Omega_init, niter=10, burnin=0, scale=FALSE, almost_parallel='none', updates_per_iter= 5, updates_per_column=1, priorCoef=normalidprior(1), priorModel=modelbinomprior(0.1), priorDiag=exponentialprior(lambda=1))


#Hard example p=5. Smallest eigenvalue very close to 0
#Th= diag(5)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)
#fitr.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=5000, burnin=0, updates_per_iter=ncol(y), updates_per_column = 10, scale=FALSE, almost_parallel='regression', tempering=1, truncratio=100, save_proposal=TRUE, prob_parallel=0.5)


#Hard example. Smallest eigenvalue very close to 0
#Th= diag(200)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.9
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.35
#y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)
#fitr2.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=10^4, burnin=0, scale=FALSE, almost_parallel='none')
