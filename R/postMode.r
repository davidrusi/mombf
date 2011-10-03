postMode <- function(y,x,priorCoef) {
  f2min <- function(th) {
    phi <- exp(th[length(th)])
    logl <- sum(dnorm(y,x %*% matrix(th[-length(th)],ncol=1),sd=sqrt(phi),log=TRUE))
    logpr <- priorFun(th)
    return(-logl-logpr)
  }
  #Set up priorFun
  tau <- as.double(priorCoef@priorPars['tau'])
  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r'])
    priorFun <- function(th) dmom(matrix(th[-length(th)],nrow=1),tau=tau,phi=exp(th[length(th)]),r=r,logscale=TRUE,penalty="product")
  } else if (priorCoef@priorDistr=='piMOM') {
    priorFun <- function(th) dimom(matrix(th[-length(th)],nrow=1),tau=tau,phi=exp(th[length(th)]),logscale=TRUE,penalty="product") 
  } else if (priorCoef@priorDistr=='peMOM') {
    priorFun <- function(th) demom(matrix(th[-length(th)],nrow=1),tau=tau,phi=exp(th[length(th)]),logscale=TRUE)
  } else {
    stop('Prior specified in priorDistr not recognized')
  }
  #Optimize
  n <- nrow(x); p <- ncol(x)
  lm1 <- lm(y~-1+x)
  thini <- c(coef(lm1),log(sum(residuals(lm1)^2/(n-p))))
  opt <- nlminb(thini, objective=f2min)$par
  return(list(coef=opt[-length(opt)],sigma2=exp(opt[length(opt)])))
}
