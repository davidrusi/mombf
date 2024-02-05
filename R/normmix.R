######################################################################################
## FUNCTIONS TO PERFORM MODEL SELECTION ON NORMAL MIXTURES
######################################################################################


setMethod("show", signature(object='mixturebf'), function(object) {
  cat('mixturebf object with',object@p,'variables\n\n')
  cat("Use draw() to obtain posterior samples. postProb() returns posterior probabilities as given below\n\n")
  print(object@postprob)
}
)


setMethod("postProb", signature(object='mixturebf'), function(object, nmax, method='norm') {
  object@postprob
}
)

setMethod("postSamples", signature(object='mixturebf'), function(object) {
    if (!is.null(object@mcmc)) { ans= object@mcmc } else { ans= NULL }
    return(ans)
}
)

#setMethod("coef", signature(object='mixturebf'), function(object, ...) {})
coef.mixturebf <- function(object,...) {
    if (length(object@mcmc)==0) stop('No MCMC samples are available, posterior means cannot be computed')
    ans= vector("list",length(object@mcmc))
    names(ans)= names(object@mcmc)
    for (i in 1:length(ans)) {
        k= as.integer(sub('k=','',names(object@mcmc)[i]))
        if ((k>1) & (!is.null(object@mcmc[[i]]))) {
            sel= which(names(object@mcmc[[i]])=='cholSigmainv')
            sel2= which(names(object@mcmc[[i]])=='momiw.weight')
            ans[[i]]= lapply(object@mcmc[[i]][-c(sel,sel2)],colMeans)
            cholSigmainv= object@mcmc[[i]][[sel]]
            ans[[i]]$Sigmainv= vector("list",k)
            for (j in 1:k) { ans[[i]]$Sigmainv[[j]]= Reduce("+",lapply(1:nrow(cholSigmainv),function(z) getSigmainv(j,p=object@p,cholSigmainv[z,])))/nrow(cholSigmainv) }
        }
        if ((k==1) & (length(grep('Normal',object@postprob[i,'model']))==1)) {
            ans[[i]]= list(eta=1, mu= object@postpars$mu1, Sigmainv=solve(object@postpars$S1)*object@postpars$nu1)
        }
    }
    return(ans)
}


getSigmainv= function(j,p,cholSigmainv) {
    ans= matrix(0,nrow=p,ncol=p)
    ans[lower.tri(ans,diag=TRUE)]= cholSigmainv[(1+(j-1)*p*(p+1)/2):(j*p*(p+1)/2)]
    return(ans %*% t(ans))
}


######################################################################################
## BAYES FACTORS
######################################################################################

bfnormmix <- function(x, k=1:2, mu0=rep(0,ncol(x)), g, nu0, S0, q=3, q.niw=1, B=10^4, burnin= round(B/10), logscale=TRUE, returndraws=TRUE, verbose=TRUE) {
    #Bayes factor comparing a Normal mixture with k-1 versus k components (different component-specific covariances)
    #
    #
    # Likelihood p(x[i,] | mu,Sigma,eta)= sum_j eta_j N(x[i,]; mu_j,Sigma_j)
    # Prior: p(mu_j, Sigma_j)= N(mu_j; mu0, g Sigma) IW(Sigma_j; nu0, S0) indep j=1,...,k
    #        p(eta)= Dir(eta; q)
    #
    # Input arguments
    # - x: n x p input data matrix
    # - k: number of components
    # - mu0, g, S0, nu0, q: prior parameters
    # - B: number of MCMC iterations
    # - burnin: number of burn-in iterations
    # - logscale: if set to TRUE the log-BF is returned
    # Output
    # - k: number of components
    # - pp.momiw: posterior prob of k components under a MOM-IW-Dir(q) prior
    # - pp.niw: posterior prob of k components under a Normal-IW-Dir(q.niw) prior
    # - probempty: probability that any one cluster is empty under a Normal-IW-Dir(q.niw) prior
    # - bf.momiw: Bayes factor comparing 1 vs k components under a MOM-IW-Dir(q) prior
    # - logpen: log of the posterior mean of the MOM-IW-Dir(q) penalty term
    # - logbf.niw: Bayes factor comparing 1 vs k components under a Normal-IW-Dir(q.niw) prior
    if (length(k)<2) stop('k must have length >=2 specifying the number of components to be compared')
    if (any(k)<1) stop('All values of k must be >=1')
    if (!is.matrix(x)) x= as.matrix(x)
    n= nrow(x); p= ncol(x)
    if (missing(g)) {
        f= function(g) (pgamma(4,p/2+1,1/(4*g)) - 0.05)^2
        g= optimize(f,c(0,5.7))$minimum
    }
    if (missing(nu0)) nu0= p+4
    if (missing(S0)) { if (p==1) { S0= matrix(1/nu0) } else { S0= diag(p)/nu0 } }
    logprobempty= logpen= double(length(k))
    logbf.niw= double(length(k)-1)
    mcmcout= vector("list",length(k)); names(mcmcout)= paste('k=',k,sep='')
    for (i in length(k):2) {
        ak= lgamma(k[i]*q.niw) - lgamma(k[i]*q.niw-q.niw) + lgamma(n+k[i]*q.niw-q.niw) - lgamma(n+k[i]*q.niw)
        #Initialize clusters via MLE
        em= Mclust(data=x, G=k[i], modelNames=ifelse(p==1,'V','VVV'),verbose=0)
        #if MLE failed, set prior to prevent 0 variances
        if ("NULL" %in% class(em)) em= try(Mclust(data=x,G=k[i],modelNames=ifelse(p==1,'V','VVV'),prior=priorControl(functionName="defaultPrior"),verbose=0))
        if ('try-error' %in% class(em)) { z= as.integer(kmeans(x, centers=k[i])$cluster) } else { z= as.integer(em$classification) }
        #Gibbs sampling (version using mombf)
        mcmcfit= "normalmixGibbsCI"(as.double(x),as.integer(n),as.integer(p),as.integer(k[i]),z,as.double(mu0),as.double(g),as.integer(nu0),as.double(S0),as.double(q.niw),as.integer(B),as.integer(burnin),as.integer(verbose))
        #mcmcfit= .Call("normalmixGibbsCI",as.double(x),as.integer(n),as.integer(p),as.integer(k[i]),z,as.double(mu0),as.double(g),as.integer(nu0),as.double(S0),as.double(q.niw),as.integer(B),as.integer(burnin),as.integer(verbose))
        names(mcmcfit)= c('logprobempty','logpen','eta','mu','cholSigmainv')
        logprobempty[i]= mcmcfit$logprobempty
        eta= t(matrix(mcmcfit$eta,ncol=(B-burnin)))
        mu= t(matrix(mcmcfit$mu,ncol=(B-burnin)))
        cholSigmainv= t(matrix(mcmcfit$cholSigmainv,ncol=(B-burnin))) #Cholesky decomp. Sigma^{-1}= cholSigmainv %*% t(cholSigmainv)
        mcmcout[[i]]= list(eta=eta,mu=mu,cholSigmainv=cholSigmainv)
        #probone= ppOneEmpty(x=x,g=g,q=q,q.niw=q.niw,mcmcfit=mcmcout[[i]],logscale=TRUE,verbose=TRUE) #to debug only, can be ignored
        qdif= q-q.niw; constddir= lgamma(k[i]*q) - lgamma(k[i]*q.niw) + k[i]*(lgamma(q.niw)-lgamma(q))
        logpen.mcmc= mcmcfit$logpen - normctNMix(k=k[i],p=p) + constddir + qdif * rowSums(log(eta))
        mcmcout[[i]]$momiw.weight= exp(logpen.mcmc-max(logpen.mcmc))
        logpen[i]= max(logpen.mcmc) + log(mean(mcmcout[[i]]$momiw.weight)) #log MOM-IW penalty
        mcmcout[[i]]$momiw.weight= mcmcout[[i]]$momiw.weight / sum(mcmcout[[i]]$momiw.weight)
        #Bayes factors
        logbf.niw[i-1]= mcmcfit$logprobempty - ak
    }
    logbf.niw= -c(0,cumsum(logbf.niw))  #BF of all models vs k=1 under Normal-IW-Dir(q.niw) prior
    logbf.momiw= logbf.niw + logpen    #BF of all models vs k=1 under MOM-IW-Dir(q) prior
    logprobempty[1]= -Inf
    pp.momiw= exp(logbf.momiw - max(logbf.momiw)); pp.momiw= pp.momiw/sum(pp.momiw)
    pp.niw= exp(logbf.niw - max(logbf.niw)); pp.niw= pp.niw/sum(pp.niw)
    ans= cbind(k,pp.momiw,pp.niw,logprobempty,logbf.momiw,logpen,logbf.niw)
    if (!logscale) { ans[,-1:-3]= exp(ans[,-1:-3]); names(ans)= sub('log','',names(ans)) }
    ans= data.frame(ans, model='Normal, VVV')
    #For 1 component, (mu,Sigma)|x ~ N(mu; mu1, Sigma/prec) IWishart(Sigma; nu=nu1, S=S1)
    xbar = colMeans(x)
    nu1 = nu0 + n
    S1 = S0 + cov(x)*(n-1) + n/(1+n*g) * (matrix(xbar-mu0,ncol=1) %*% matrix(xbar-mu0,nrow=1))
    w = n/(n + 1/g)
    m1 = w * xbar + (1 - w) * mu0
    postpars= list(mu1=m1, prec= n+1/g, nu1=nu1, S1=S1)
    #Return answer
    ans= new("mixturebf", postprob=ans, p=p, n=n, priorpars=list(mu0=mu0, g=g, nu0=nu0, S0=S0, q=q, q.niw=q.niw), postpars=postpars)
    if (returndraws) ans@mcmc= mcmcout
    return(ans)
}








######################################################################################
## POSTERIOR PROBABILITY OF EMPTY CLUSTERS
######################################################################################

ppOneEmpty <- function(x,g,q,q.niw,mcmcfit,logscale=TRUE,verbose=TRUE) {
    # Posterior probability that one cluster is empty and posterior mean of MOM-IW penalty
    # Input arguments
    # - x: n * p data matrix
    # - mcmcfit: MCMC output fitting a Normal mixture. A list with elements list(eta=eta,mu=mu,cholSigmainv=cholSigmainv), where Sigma^{-1}= cholSigmainv %*% t(cholSigmainv) is the precision matrix
    # - logscale: if TRUE log-prob is returned
    # - verbose: set to TRUE to print iteration progress
    # Output
    # - logprobempty: mean P(n_j=0|y,k) under a Dir(q.niw) prior, where n_j is the number of individuals in cluster j
    # - logpen: posterior mean of the MOM penalty d(theta) under a product Normal-IW-Dirichlet prior, that is E(d(theta) | y,k)
    #   where d(theta)= (1/C_k) [ Dir(eta; q) / Dir(eta; q.niw)] [ prod_{i<j} (mu_i-mu_j)' A^{-1} (mu_i-mu_j)/g ] [ prod_j N(mu_j; 0, g A) IW(Sigma_j; nu0, S0) ]
    #         C_k: MOM-IW prior norm constant
    #         g: prior dispersion parameter
    #         A^{-1}= sum_j Sigma_j^{-1} / k
    niter= nrow(mcmcfit$eta) #number of MCMC draws
    p= ncol(x); k= ncol(mcmcfit$eta)
    probOneEmpty= matrix(NA,nrow=niter,ncol=k)
    logpen= double(niter)
    if (verbose) cat("Post-processing MCMC output")
    qdif= q-q.niw; constddir= lgamma(k*q) - lgamma(k*q.niw) + k*(lgamma(q.niw)-lgamma(q))
    for (i in 1:niter) {
        #Extract (eta,mu,Sigma)
        eta= mcmcfit$eta[i,]
        mu= lapply(1:k, function(j) mcmcfit$mu[i,(1+(j-1)*p):(j*p)])
        Sigmainv= lapply(1:k, function(j) { getSigmainv(j,p,mcmcfit[[3]][i,]) })
        Sigma= lapply(Sigmainv, solve)
        #Cluster allocation probabilities (n x k matrix)
        pp= sapply(1:k, function(j) dmvnorm(x,mu[[j]],sigma=Sigma[[j]],log=TRUE) + log(eta[j]))
        pp= exp(pp-apply(pp,1,'max')); pp= pp/rowSums(pp)
        pp[pp> 1-1e-13]= 1-1e-13
        #log-probability of each cluster being empty
        probOneEmpty[i,]= colSums(log(1-pp))
        #MOM-IW penalty
        Ainv= Reduce("+",Sigmainv) / k
        mumat= do.call(rbind,mu)
        logprior= sum(dmvnorm(mumat,sigma=g*solve(Ainv),log=TRUE)) + constddir + qdif * sum(log(eta))
        #same as logprior= sum(dmvnorm(mumat,sigma=g*solve(Ainv),log=TRUE)) + ddir(eta,q=q,logscale=TRUE) - ddir(eta,q=q.niw,logscale=TRUE)
        for (jj in 1:k) logprior= logprior - dmvnorm(mu[[jj]],sigma=g*Sigma[[jj]],log=TRUE)
        d= as.vector(dist(mumat %*% t(chol(Ainv))))^2 #Pairwise Mahalanobis distances between mu's
        logpen[i]= logprior + sum(log(d)) - k*(k-1)/2*log(g) - normctNMix(k=k,p=p)
        if (verbose & ((i %% (niter/10))==0)) cat(".")
    }
    if (verbose) cat("\n")
    #Compute mean log-prob, in a way that helps prevent numerical overflow
    logprob= double(k)
    for (i in 1:k) {
        m= max(probOneEmpty[,i])
        logprob[i]= m + log(mean(exp(probOneEmpty[,i] - m)))
    }
    logprob= max(logprob) + log(mean(exp(logprob - max(logprob))))
    logpen= max(logpen) + log(mean(exp(logpen-max(logpen))))
    ans= c(logprobempty=logprob, logpen=logpen)
    if (!logscale) ans= exp(ans)
    return(ans)
}





######################################################################################
## PRIOR NORMALIZATION CONSTANT
######################################################################################

normctNMix= function(k,p) {
#Prior normalization constant for MOM-IW prior
    if (p==1) {
        ans= sum(lgamma(2:(k+1)))
    } else if (k==2) {
        ans= log(2*p)
    } else {
        if ((p<=ncol(Cktab)) & (k<=nrow(Cktab)+1)) { ans= Cktab[k-1,p] } else { ans= normctNMixEstimate(p=p,k=k,B=20000) }
    }
    return(ans)
}

#Pre-tabulated log norm constant for MOM-IW
Cktab= matrix(NA,nrow=19,ncol=20)
rownames(Cktab)= paste('k=',2:(nrow(Cktab)+1),sep=''); colnames(Cktab)= paste('p=',1:ncol(Cktab),sep='')
Cktab[1,]= c(0.693,1.386,1.792,2.079,2.303,2.485,2.639,2.773,2.890,2.996,3.091,3.178,3.258,3.332,3.401,3.466,3.526,3.584,3.638,3.689)
Cktab[2,]= c(2.485,4.561,5.709,6.507,7.133,7.655,8.095,8.481,8.817,9.118,9.386,9.647,9.878,10.092,10.295,10.484,10.663,10.827,10.984,11.140)
Cktab[3,]= c(5.663,9.810,11.976,13.504,14.713,15.658,16.526,17.253,17.900,18.478,19.024,19.498,19.962,20.374,20.759,21.132,21.479,21.812,22.111,22.412)
Cktab[4,]= c(10.450,17.050,20.791,23.242,25.145,26.677,28.074,29.204,30.284,31.207,32.042,32.816,33.576,34.252,34.892,35.487,36.039,36.586,37.097,37.565)
Cktab[5,]= c(17.030,27.187,32.119,35.893,38.505,40.914,42.929,44.680,46.042,47.294,48.619,49.744,50.860,51.782,52.796,53.619,54.470,55.258,55.968,56.681)
Cktab[6,]= c(25.555,38.463,45.795,51.076,54.848,58.037,60.621,62.954,65.149,67.026,68.901,70.553,71.715,73.189,74.438,75.502,76.727,77.797,78.748,79.829)
Cktab[7,]= c(36.159,52.974,61.186,69.241,73.942,79.432,81.386,85.021,87.891,90.849,92.470,94.466,98.722,98.475,100.097,101.348,103.371,104.420,105.544,107.230)
Cktab[8,]= c(48.961,66.164,80.042,88.514,96.708,101.179,109.560,110.135,113.193,117.850,119.457,122.496,124.638,127.080,128.651,131.149,132.972,134.709,136.385,137.994)
Cktab[9,]= c(64.066,82.847,102.055,113.733,120.787,126.917,135.411,138.371,143.527,146.264,151.659,152.687,156.542,159.355,161.416,164.206,166.930,168.785,171.930,173.261)
Cktab[10,]= c(81.568,97.246,124.129,136.741,144.772,157.996,161.171,169.065,174.822,179.265,183.386,187.168,191.950,195.910,199.041,203.199,203.853,207.579,211.241,212.259)
Cktab[11,]= c(101.555,122.362,146.287,162.528,180.388,187.270,192.937,203.938,212.889,217.734,221.541,228.209,231.441,236.452,238.341,243.116,246.672,249.123,252.716,254.760)
Cktab[12,]= c(124.107,141.233,168.917,196.092,214.977,221.311,234.298,240.367,248.477,260.534,260.992,265.749,272.312,278.072,283.667,287.527,291.691,296.176,301.515,301.772)
Cktab[13,]= c(149.299,156.957,206.438,225.427,243.774,259.805,267.883,285.804,291.234,299.880,309.264,311.079,316.570,323.971,331.077,334.229,342.464,344.108,348.812,354.336)
Cktab[14,]= c(177.198,179.557,237.431,253.155,274.733,298.454,315.369,323.422,336.471,341.744,353.027,359.312,367.176,377.951,378.114,386.572,393.042,400.376,404.191,406.108)
Cktab[15,]= c(207.870,202.500,250.975,288.691,316.324,340.346,351.251,370.607,382.489,395.351,402.748,411.758,422.648,431.317,438.742,441.203,447.649,451.989,461.505,467.732)
Cktab[16,]= c(241.375,229.919,285.081,338.441,357.618,381.691,401.358,415.652,436.766,442.482,465.513,470.246,477.489,480.464,493.833,501.363,506.500,513.625,519.481,527.751)
Cktab[17,]= c(277.770,254.477,321.324,370.048,395.518,425.777,449.232,463.807,483.476,503.058,510.875,527.850,533.347,547.950,555.441,564.276,575.912,578.512,585.708,593.824)
Cktab[18,]= c(317.110,274.729,354.350,418.900,445.156,478.831,500.711,524.482,547.148,557.461,578.541,584.673,594.984,611.518,618.742,628.901,638.195,648.375,656.415,666.447)
Cktab[19,]= c(359.446,317.930,397.421,452.398,504.381,536.868,556.158,578.169,600.491,627.437,633.916,657.388,661.210,671.309,683.771,700.666,706.412,725.335,732.537,739.385)


normctNMixEstimate= function(k,p,B=10^4) {
#Monte Carlo estimate of MOM-IW normalization constant (B is the number of MC samples)
    d= double(B)
    for (i in 1:B) {
        mu= matrix(rnorm(p*k),ncol=p)
        d[i]= sum(log(as.vector(dist(mu))^2)) #Pairwise Mahalanobis distances between mu's
    }
    max(d) + log(mean(exp(d-max(d))))
}



