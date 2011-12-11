###
### zellnerbf.unknownsig.R
###

zellnerbf.unknownsig <- function(theta1hat,V1,n,nuisance.theta,g=1,theta0,ssr,logbf=FALSE) {
#Bayes factor based on Zellner's g-prior for linear models (unknown sigma^2 case). 
# - theta1hat: vector with estimated value of the coefficients that are to be tested
# - V1: submatrix of covariance corresponding to elements in theta1hat
# - n: sample size used to fit the model
# - g: prior parameter
# - theta0: hypothesized value for theta1hat (defaults to 0)
# - ssr: sum of squared residuals
# - logbf: if log==TRUE, the natural logarithm of the Bayes factor is returned
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat); p <- p1 + nuisance.theta
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1)  
sigma2hat <- (ssr + l/(1+n*g))/(n-nuisance.theta)
muk <- p1+ l* n*g/((1+n*g)*sigma2hat)
bf <- (-(n-nuisance.theta)/2)*log(1+n*g*ssr/(ssr+l)) + ((n-p)/2)*log(1+n*g)
if (!logbf) bf <- exp(bf)
return(bf)
}

