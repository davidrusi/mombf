#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#Rcpp::compileAttributes("~/github/mombf")

library(mombf)
library(mvtnorm)
set.seed(1)

#fast Cholesky updates
library(mombf)
set.seed(1234)
n <- 20
X3 <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
X3 <- cbind(matrix(1, nrow = n, ncol = 1), X3) # add intercept
theta3_truth <- matrix(c(1, 0, 1, 0), ncol = 1)
theta3_truth_bool <- as.logical(theta3_truth)
theta3_truth_idx <- which(theta3_truth_bool)
y3 <- X3 %*% theta3_truth + rnorm(n)

n <- 200
X6 <- matrix(rnorm(n * 6), nrow = n, ncol = 6)
X6 <- cbind(matrix(1, nrow = n, ncol = 1), X6) # add intercept
theta6_truth <- matrix(c(0, 0, 1, 1, 0, 1, 1), ncol = 1)
theta6_truth_bool <- as.logical(theta6_truth)
theta6_truth_idx <- which(theta6_truth_bool)
y6 <- X6 %*% theta6_truth + rnorm(n)

n <- 150
X9 <- matrix(rnorm(n * 9), nrow = n, ncol = 9)
X9 <- cbind(matrix(1, nrow = n, ncol = 1), X9) # add intercept
theta9_truth <- matrix(c(0, 1, 1, 0, 1, 1, 1, 0, 0, 0), ncol = 1)
theta9_truth_bool <- as.logical(theta9_truth)
theta9_truth_idx <- which(theta9_truth_bool)
y9 <- X9 %*% theta9_truth + rnorm(n)
groups9 <- c(1,2,3,4,5,5,5,6,6,6)


#modelSelection
family="normal"; pCoef=normalidprior(); pGroup=groupzellnerprior()
pDelta = modelunifprior()
groups <- c(1, 1, 2, 2, 3, 4, 4)
#nlpMarginal(sel=rep(TRUE,7), y=y6, x=X6, priorCoef=pCoef, priorGroup=pGroup, groups=groups) #full model
fit <- modelSelection(y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, XtXprecomp=FALSE, family=family, priorSkew=pCoef, priorGroup=pGroup, groups=groups, center=FALSE, scale=FALSE)


#icarplus prior
#p= 6; n=200
#x= matrix(rnorm(n*p), nrow=n, ncol=p)
#b= matrix(c(1,.9,.8,.7,rep(0,p-4)), ncol=1)
#y= x %*% b + rnorm(n)
#neighbours= mombf:::neighbours_1d(1:p)
#nlpMarginal(sel=c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE), y=y, x=x, priorCoef=icarplusprior(taustd=1), neighbours=neighbours)


#Cholesky update
#B= matrix(c(2,1,.5,.25, 1,1.5,.3,.4, .5,.3,1,.75, .25,.4,.75,2),nrow=4,byrow=TRUE)
#mombf:::testfunction(B, 3, 3)


## GGM examples
#library(mombf)
#library(mvtnorm)
#library(survival)
#set.seed(1)

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
#Th= diag(5)
#diag(Th)= 1.5
#Th[abs(col(Th) - row(Th))==1]= 0.95
#Th[abs(col(Th) - row(Th))==2]= 0.5
#Th[abs(col(Th) - row(Th))==3]= 0.64
#y= scale(rmvnorm(100, sigma=solve(Th)), center=TRUE, scale=FALSE)
#Omegaini= Th
#fit <- modelSelectionGGM(y, sampler='birthdeath', Omegaini=Omegaini, niter=100, burnin=0, updates_per_iter=ncol(y), updates_per_column=ncol(Th), scale=FALSE, almost_parallel='regression', tempering=1, truncratio=100)

