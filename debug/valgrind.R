library(mombf)
library(mvtnorm)

Th= diag(3); Th[1,2]= Th[2,1]= 0.5
sigma= solve(Th)

z= matrix(rnorm(1000*3), ncol=3)
y= z %*% chol(sigma)

#Obtain posterior samples
fit= modelSelectionGGM(y, scale=FALSE)

#Parameter estimates, intervals, prob of non-zero
coef(fit)

#Estimated inverse covariance
icov(fit)

#Estimated inverse covariance, entries set to 0
icov(fit, threshold=0.95)

#Shows first posterior samples
head(fit$postSample)

