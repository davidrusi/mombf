# Estimate GLM parameters under the full model
# - y: response variable
# - x: design matrix
# - family: 'normal', 'binomial', 'poisson'. For families 'laplace', 'twopiecenormal', 'twopiecelaplace', initialization is based on the normal model
# - initpar: method to initialize the parameters. It can be 'MLE', 'L1', 'MLE-aisgd', 'L2-aisgd'. If missing, it uses MLE when p <= n/2 and L1 when p > n/2. If p>200 or n>10,000 then ai-sgd is used (averaged intrinsic stochatic gradient descent)
initParameters <- function(y, x, family, initpar) {
    n= length(y); p= ncol(x)
    if (!(initpar %in% c('none','auto','MLE','L1'))) stop("initpar must be 'none', 'auto', 'MLE'")
    #Disabled sgd since sgd package was discontinued
    #if (!(initpar %in% c('none','auto','MLE','MLE-aisgd','L1','L2-aisgd'))) stop("initpar must be 'none', 'auto', 'MLE', 'MLE-aisgd', 'L1' or 'L2-aisgd'")
    if (initpar =='auto') {
        initpar= ifelse(p <= n/2, 'MLE', 'L1')
        #large= (n > 10000) || (p > 200)
        #if ((p<=n/2) && (!large)) {
        #    initpar= 'MLE'
        #} else if ((p<=n/2) && large) {
        #    initpar= 'MLE-aisgd'
        #} else if ((p>n/2) && (!large)) {
        #    initpar= 'L1'
        #} else {
        #    initpar= 'L2-aisgd'
        #}
    }
    if (!(family %in% c('binomial','poisson'))) family= 'gaussian'
    if (initpar == 'none') {
        ans= rep(0,p)
    } else if (initpar == 'MLE') {
        fit= glm(y ~ x -1, family=family)
        ans= coef(fit)
    #} else if (initpar == 'MLE-aisgd') {
    #    fit= sgd(x=x, y=y, model="glm", model.control=list(family=family), sgd.control=list(method="ai-sgd"))
    #    ans= coef(fit)
    }  else if (initpar == 'L1') {
        int= (colSums(x==1) == nrow(x))
        fit= glmnet::glmnet(x=x[,!int], y=y, intercept=TRUE, family=family, nlambda=50, alpha=1)
        b= matrix(NA, nrow=ncol(x), ncol=length(fit$lambda))
        b[!int,]= as.matrix(coef(fit)[-1,])
        b[int,]= coef(fit)[1,]
        bicvals= glmBIC(b, y=y, x=x, family=family)
        ans= coef(fit)[,which.min(bicvals)]
    #} else if (initpar == 'L1-aisgd') {
        #Function sgd does not appear to return correct results for aisgd with L1
      #} else if (initpar == 'L2-aisgd') {
      #  #Ridge penalty
      #  fit= sgd(x=x, y=y, model="glm", model.control=list(family=family, lambda2=0.001), sgd.control=list(method="ai-sgd"))
      #  ans= coef(fit)
    #} else stop("Invalid 'initpar'. Allowed values are 'MLE', 'L1', 'MLE-aisgd', 'L2-aisgd'")
    } else stop("Invalid 'initpar'. Allowed values are 'MLE', 'L1'")
    return(ans)    
}


#Wrapper to initParameters, returning 0 when initpar=='none'
getthinit = function(y, x, family, initpar, enumerate) {
  if ((!missing(initpar)) && is.numeric(initpar)) {
      thinit= as.double(initpar)
      if (length(thinit) != ncol(x)) stop(paste("x has",ncol(x),"columns but initpar has",length(thinit),"elements"))
  } else {
      if ((!missing(initpar)) && (initpar=='none')) {
          usethinit= as.integer(ifelse(enumerate, 0, 1)) #tell C code to not init or do its own init
      } else {
          usethinit= as.integer(3) #tell C code to use thinit
      }
      thinit= as.double(initParameters(y=y, x=x, family=family, initpar=initpar))
  }
  ans= list(thinit=thinit, usethinit=usethinit)
  return(ans)
}


#Obtain BIC for GLMs
# - beta: regression coefficients. Either a vector with length equal to ncol(x), or a matrix with nrow equal to ncol(x) and different estimates in the columns of beta
# - y: response variable
# - x: design matrix
# - family: 'binomial', 'poisson' or otherwise the normal model is used
glmBIC = function(beta, y, x, family) {
    n= length(y); p= ncol(x)
    if (is.vector(beta)) {
        beta= matrix(beta,ncol=1)
        npars= sum(beta!=0)
    } else {
        npars= colSums(beta!=0)
    }
    lp= x %*% beta
    if (family=='binomial') {
        logl= double(ncol(lp))
        for (i in 1:ncol(lp)) logl[i]= logL_logreg(beta[,i], ytX= t(y) %*% x, Xbeta=lp[,i], n=n, X=x, logscale=TRUE)
        ans = -2 * logl + log(n)*p
    } else if (family=='poisson') {
        logl= double(ncol(lp))
        for (i in 1:ncol(lp)) logl[i]= logL_poisreg(beta[,i], ytX= t(y) %*% x, Xbeta=lp[,i], n=n, X=x, sumlfactorialy=sum(lfactorial(y)), logscale=TRUE)
        ans = -2 * logl + log(n)*p
    } else {
        ans= n * log(colSums((y-lp)^2)/n) + n*(log(2*pi)+1) + log(n)*p
    }
    return(ans)
}


###########################################################################################
# AUXILIARY FUNCTIONS FOR LOGISTIC REGRESSION
###########################################################################################

logit= function(z) log(z/(1-z))  #logit function
expit= function(z) 1/(1+exp(-z)) #inverse of logit function

#Negative logistic regression log-likelihood
logL_logreg= function(beta, ytX, Xbeta, n, X, logscale=TRUE) {
    if (any(beta != 0)) {
        if (missing(Xbeta)) Xbeta= as.vector(X %*% matrix(beta,ncol=1))
        ans= sum(ytX * beta) - sum(log(1+exp(Xbeta)))
    } else {
        if (missing(n)) stop("If beta==0, n must be specified")
        ans= -n * log(2)
    }
    if (!logscale) ans= exp(ans)
    return(ans)
}



###########################################################################################
# AUXILIARY FUNCTIONS FOR POISSON REGRESSION
###########################################################################################

#Negative Poisson regression log-likelihood
logL_poisreg= function(beta, ytX, Xbeta, n, X, sumlfactorialy, logscale=TRUE) {
    #if (missing(sumlfactorialy)) sumlfactorialy= sum(lfactorial(y))
    if (any(beta != 0)) {
        if (missing(Xbeta)) Xbeta= as.vector(X %*% matrix(beta,ncol=1))
        ans= sum(ytX * beta) + sum(exp(Xbeta)) - sumlfactorialy
    } else {
        if (missing(n)) stop("If beta==0, n must be specified")
        ans= - n - sumlfactorialy
    }
    if (!logscale) ans= exp(ans)
    return(ans)
}

