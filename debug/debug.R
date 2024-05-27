#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#library(Rcpp)
#compileAttributes("~/github/mombf")

library(mvtnorm)
library(tidyverse)
library(mombf)
set.seed(784)
# Dimension
p <- 250
s <- 24
n <- 250
# Generate design matrix
A <- matrix(runif(p*p,-1,1),ncol=p)
Sigma <- t(A)%*%A
mu <- rep(runif(1,1,2),p)
X <- rmvnorm(n,mu,Sigma)
pcs <- princomp(X)$scores
norms.pcs<-diag(t(pcs) %*% pcs)/n
pcs.scaled <- pcs %*% diag(1/sqrt(norms.pcs))
X.design <- pcs %*% diag(1/sqrt(norms.pcs))
# True beta
betamin <- 0.35
beta_star <- c(betamin, runif(s-2,1.05*betamin,1.1*betamin),betamin,rep(0,p-s))
# Simulation
y <- X.design%*%(beta_star)+rnorm(n,0,1)

#debug
models= matrix(FALSE, nrow=1, ncol=ncol(X.design))
fit.null= bestBIC(y, X.design, models=models) #should be 1067.586



## GGM examples

library(mombf)
library(mvtnorm)
set.seed(1)
#p= 5
#Th= diag(p); Th[1,2]= Th[2,1]= 0.5
#sigma= solve(Th)
#z= matrix(rnorm(1000*p), ncol=p)
#y= z %*% chol(sigma)
#fit= modelSelectionGGM(y, scale=FALSE, almost_parallel="regression", sampler='birthdeath', niter=10)


#Th= diag(5)
#diag(Th)= c(1.5, 1.6, 1.7, 1.8, 1.9)
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#Thnew= Th
#Thnew[,3]= Thnew[3,]= c(0.4,0.8,1.4,0.75,0.3)
#colid= 3
#mombf:::testfunction(Th, oldcol=1, newcol=2)

#Hard example p=5. Smallest eigenvalue very close to 0
#Th= diag(5)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)
#fitr.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=5000, burnin=0, scale=FALSE, almost_parallel='regression', nbirth=1, tempering=1, truncratio=100, save_proposal=TRUE, fullscan=TRUE, prob_parallel=0.5)


#Hard example. Smallest eigenvalue very close to 0
Th= diag(50)
diag(Th)= 1.5
Th[abs(col(Th) - row(Th))==1]= 0.9
Th[abs(col(Th) - row(Th))==2]= 0.5
Th[abs(col(Th) - row(Th))==3]= 0.35
y= scale(rmvnorm(500, sigma=solve(Th)), center=TRUE, scale=FALSE)
fitr2.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=10^4, burnin=0, scale=FALSE, almost_parallel='regression', nbirth=1, tempering=1, truncratio=100, save_proposal=TRUE, fullscan=TRUE, prob_parallel=0.5)
