modelsearchBlockDiag <- function(y, x, priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(1,1), priorVar=igprior(0.01,0.01), blocksize=10, maxiter=10, maxvars=100, maxlogmargdrop=20, maxenum=10, verbose=TRUE) {
    n <- length(y); p <- ncol(x)
    #Check and format input parameters
    if (priorDelta@priorDistr=='binomial' & ('p' %in% names(priorDelta@priorPars))) {
        rho <- priorDelta@priorPars[['p']]
        if (length(rho)>1) stop("modelsearchBlockDiag not implemented for vector p in modelbinomprior")
        priorModel <- function(nvar) nvar*log(rho) + (p-nvar)*log(1-rho)
        priorModelBlock <- function(nvar,blocksize) nvar*log(rho) + (blocksize-nvar)*log(1-rho)
    } else if (priorDelta@priorDistr=='binomial' & !('p' %in% names(priorDelta@priorPars))) {
        alpha=priorDelta@priorPars['alpha.p']; beta=priorDelta@priorPars['beta.p']
        priorModel <- function(nvar) lbeta(nvar + alpha, p - nvar + beta) - lbeta(alpha, beta)
        priorModelBlock <- function(nvar,blocksize) rep(-blocksize*log(2),length(nvar)) #set coolblock prior to unif. post prob based on correct priorModel
    } else if (priorDelta@priorDistr=='uniform') {
        rho <- 0.5
        priorModel <- function(nvar) rep(-p*log(2),length(nvar))
        priorModelBlock <- function(nvar,blocksize) rep(-blocksize*log(2),length(nvar))
    } else { stop("Prior on model space not recognized. Use modelbbprior(), modelunifprior() or modelbinomprior()") }
    if (priorCoef@priorDistr == 'zellner') {
        g <- as.double(priorCoef@priorPars['tau'])
    } else if (priorCoef@priorDistr == 'pMOM') {
        tau <- as.double(priorCoef@priorPars['tau'])
        g <- n * tau
    } else {
        stop("priorCoef must be pMOM or Zellner. Use momprior() or zellnerprior()")
    }
    if (priorVar@priorDistr != 'invgamma') stop('priorVar must be an inverse gamma. Use igprior()')
    a.phi <- as.double(priorVar@priorPars['alpha'])
    l.phi <- as.double(priorVar@priorPars['lambda'])
    #Pre-compute useful quantities
    sumy2 <- l.phi + sum(y^2)
    apos <- a.phi+n
    shrinkage <- g/(1+g)
    Xty <- t(x) %*% matrix(y,ncol=1)
    #Block search forward
    bestlogpp <- bestmarg <- -Inf; iter <- 0; improved <- TRUE
    e <- y; selfixed <- integer(0)
    ans <- data.frame(nvars=integer(0), modelid=character(0), logpp=double(0), marglhood=double(0))
    while ((iter< maxiter) & improved) {
        improved <- FALSE
        #Forward pass
        ansfwd <- blocksearch(y=y,e=e,x=x,selfixed=selfixed,blocksize=blocksize,maxvars=maxvars,priorCoef=priorCoef,priorVar=priorVar,priorModel=priorModel,priorModelBlock=priorModelBlock,g=g,a.phi=a.phi,l.phi=l.phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,maxlogmargdrop=maxlogmargdrop)
        ansfwd <- ansfwd[!is.na(ansfwd$logpp),]
        selfwd.maxlogpp <- which.max(ansfwd$logpp)
        selfwd.bestmarg <- which.max(ansfwd$marglhood)
        if (ansfwd$logpp[selfwd.maxlogpp] > bestlogpp) { bestlogpp <- ansfwd$logpp[selfwd.maxlogpp]; sel.maxlogpp <- nrow(ans)+selfwd.maxlogpp; improved <- TRUE }
        if (ansfwd$marglhood[selfwd.bestmarg] > bestmarg) { bestmarg <- ansfwd$marglhood[selfwd.bestmarg]; sel.bestmarg <- nrow(ans)+selfwd.bestmarg }
        sel <- as.numeric(strsplit(ansfwd$modelid[max(selfwd.maxlogpp,selfwd.bestmarg)],split=',')[[1]])  #choose post mode or maximal integrated likelihood model, whichever is largest
        ans <- rbind(ans,ansfwd[!(ansfwd$modelid %in% ans$modelid),])
        #Backward pass
        if (length(sel)>1) {
            ansbackward <- blocksearch(y=y,x=x[,sel],blocksize=blocksize,maxvars=maxvars,priorCoef=priorCoef,priorVar=priorVar,priorModel=priorModel,priorModelBlock=priorModelBlock,g=g,a.phi=a.phi,l.phi=l.phi,sumy2=sumy2,apos=apos,shrinkage=shrinkage,maxlogmargdrop=maxlogmargdrop)
            ansbackward$modelid <- sapply(lapply(strsplit(ansbackward$modelid,','), function(z) sel[as.numeric(z)]), paste, collapse=',')
            ansbackward <- ansbackward[!(ansbackward$modelid %in% ans$modelid),]
            if (nrow(ansbackward)>0) {
                selback.maxlogpp <- which.max(ansbackward$logpp)
                selback.bestmarg <- which.max(ansbackward$marglhood)
                if (ansbackward$logpp[selback.maxlogpp] > bestlogpp) { bestlogpp <- ansbackward$logpp[selback.maxlogpp]; sel.maxlogpp <- nrow(ans)+selback.maxlogpp; improved <- TRUE }
                if (ansbackward$marglhood[selback.bestmarg] > bestmarg) { bestmarg <- ansbackward$marglhood[selback.bestmarg]; sel.bestmarg <- nrow(ans)+selback.bestmarg }
                ans <- rbind(ans,ansbackward)
            }
        }
        if (improved) {
            selfixed <- as.numeric(strsplit(ans[sel.bestmarg,'modelid'],split=',')[[1]])
            e <- y - predict(lm(y ~ x[,selfixed]))
        }
        if (verbose) {
            cat("Iteration",iter,"\n")
            cat("  Posterior mode",ans[sel.maxlogpp,'modelid'],ans[sel.maxlogpp,'logpp'],"\n")
            cat("  Mode marg likelihood",ans[sel.bestmarg,'modelid'],ans[sel.bestmarg,'marglhood'],"\n")
        }
        iter <- iter+1
    }
    #Exhaustive enumeration for submodels of highest marginal likelihood model (or highest post prob model, if the former has >maxenum variables)
    sel2enum <- which(max(ans$marglhood,na.rm=TRUE) - ans$marglhood <= maxlogmargdrop)
    sel2enum <- sel2enum[length(sel2enum)]
    if (ans$nvars[sel2enum] > maxenum) { sel2enum <- which(ans$nvars<= maxenum); sel2enum <- sel2enum[length(sel2enum)] }
    sel <- as.numeric(strsplit(ans$modelid[sel2enum],split=',')[[1]])
    sel <- sel[order(sel)]
    xsel <- x[,sel,drop=FALSE]
    if (is.null(colnames(xsel))) colnames(xsel) <- paste("x",sel,sep='')
    ms <- modelSelection(y=y,x=xsel,center=FALSE,scale=FALSE,enumerate=TRUE,maxvars=min(c(ncol(xsel),10)),priorCoef=priorCoef,priorDelta=modelunifprior(),priorVar=priorVar,verbose=FALSE)
    nvarsenum <- sapply(strsplit(as.character(ms$models$modelid),split=','),length)
    modelid <- lapply(strsplit(as.character(ms$models$modelid),split=','), function(i) sel[as.numeric(i)])
    ct <- nlpMarginal(sel=integer(0),family='normal',priorCoef=priorCoef,priorVar=priorVar,y=y,x=x,logscale=TRUE) - log(ms$models[ms$models$modelid=='','pp'])
    ansenum <- data.frame(nvars=sapply(modelid,length), modelid=sapply(modelid,paste,collapse=','), logpp=log(ms$models$pp) + priorModel(nvarsenum) + ct)
    ans <- unique(rbind(ans[,c('nvars','modelid','logpp')],ansenum[!(ansenum$modelid %in% ans$modelid),]))
    return(ans[order(ans$logpp,decreasing=TRUE),])
}


blocksearch <- function(y,e=y,x,selfixed=integer(0),blocksize,maxvars,priorCoef,priorVar,priorModel,priorModelBlock,g,a.phi,l.phi,sumy2,apos,shrinkage,maxlogmargdrop=20) {
    #Use spectral clustering followed by coolblock to define a sequence of models gamma of increasing size. Return actual p(y|gamma)p(gamma) for that sequence.
    notfixed <- setdiff(1:ncol(x), selfixed)
    nfixed <- length(selfixed)
    if (nfixed>0) {
        xn <- x[,notfixed,drop=FALSE]
    } else {
        xn <- x
    }
    #Cluster variables
    blocks <- spectralClus(xn, blocksize=blocksize)
    blocks <- as.integer(factor(blocks))
    nblocks <- max(blocks)
    blocksize <- table(blocks)
    nn <- 1:ncol(xn)
    varidx <- lapply(as.integer(names(blocksize)), function(k) nn[blocks==k])
    #Compute u-scores, least squares estimates and covariances within blocks
    Xty <- t(xn) %*% matrix(e,ncol=1)
    XtX <- lapply(1:nblocks,function(i) { sel <- which(blocks==i); t(xn[,sel]) %*% xn[,sel] })
    ubyblocks <- blockuscores(blocks, blocksize=blocksize, y=e, Xty=Xty, XtX=XtX)
    u <- ubyblocks$u; umv <- ubyblocks$umv; ubest <- ubyblocks$ubest
    #Enumerate models with coolblock and compute p(y,model)
    models <- coolblock(ubest,g=g,priorModelBlock=priorModelBlock,varidx=varidx,maxvars=maxvars)
    models <- models[!is.na(models$u),]
    modelid <- strsplit(models$modelid,split=',')
    models$modelid <- sapply(lapply(modelid,function(z) z[order(as.numeric(z))]), paste, collapse=',') #put variables ids in increasing order
    priorp <- priorModel(nfixed+models$nvars)
    marglhood <- pp <- rep(NA,nrow(models))
    i <- 1; bestpp <- -Inf; bestmarg <- -Inf; stopearly <- FALSE
    modelnames <- character(nrow(models))
    while ((i<=nrow(models)) & (!stopearly)) {
        sel <- c(selfixed,notfixed[as.numeric(modelid[[i]])])
        sel <- sel[order(sel)]
        modelnames[i] <- paste(sel,collapse=',')
        marglhood[i] <- nlpMarginal(sel=sel,family='normal',priorCoef=priorCoef,priorVar=priorVar,y=y,x=x,logscale=TRUE)
        pp[i] <- priorp[i] + marglhood[i]
        if (pp[i] > bestpp) { bestpp <- pp[i] }
        if (marglhood[i] > bestmarg) { bestmarg <- marglhood[i] } else if (marglhood[i] + maxlogmargdrop < bestmarg) { stopearly <- TRUE }
        i <- i+1
    }
    ans <- data.frame(nvars=models[,'nvars'],modelid=modelnames,logpp=pp,marglhood=marglhood,stringsAsFactors=FALSE)
    return(ans)
}

initKmeans <- function(coord, nclus) {
    #Initialize clusters based on quantile grid on first component
    clini= as.numeric(cut(coord[,1], breaks=c(-Inf,quantile(coord[,1],probs=seq(1/nclus,1-1/nclus,length=nclus-1)),Inf)))
    centers= as.matrix(aggregate(coord,by=list(clini),FUN='mean')[,-1])
    centers= sapply(1:nclus, function(i) { sel= which(clini==i); coord[sel[which.min(colSums((t(coord[sel,])-centers[i,])^2))],] })
    if (ncol(coord)>1) centers= t(centers) else centers= matrix(centers,ncol=1)
    return(centers)
}

spectralClus <- function(x, blocksize=10, scale=FALSE, ndim= min(c(ncol(x)/blocksize,nrow(x)))) {
    #Spectral clustering of columns in x via k-means on largest eigenvalues of cor(x)^2. A deterministic initialization using quantiles in the first eigenvector is used. Clusters are iteratively subdivided until all clusters have size<=blocksize.
    # - x: data matrix whose columns we want to cluster
    # - blocksize: maximum cluster size
    # - scale: if TRUE
    # - ndim: number of eigenvectors to use in the cluster
    p <- ncol(x)
    if (p<blocksize) {
        cl <- rep(1,p)
    } else {
        xstd= scale(x,center=TRUE,scale=TRUE)
        S= cov(xstd)^2
        #if (!inverse) { S= cov(xstd)^2 } else { S= pseudoinverse(cov(xstd))^2 } #pseudoinv from package corpcor
        D= diag(p); diag(D)= 1/sqrt(rowSums(S))
        S= D %*% S %*% D
        eig= eigen(S)
        ndim <- min(c(ndim,nrow(x),p))
        if (!scale) {
            coord= t(t(eig$vectors[,1:ndim,drop=FALSE]) * sqrt(eig$values[1:ndim]))
        } else {
            coord= eig$vectors[,1:ndim,drop=FALSE]
        }
        #Initialize clusters based on quantile grid on first component
        centers= initKmeans(coord, nclus=max(round(p/blocksize),2))
        km= kmeans(coord, centers=centers)
        cl= km$cluster
        maxcl= max(cl)
        tab= table(cl)
        clussel= as.numeric(names(tab[tab>blocksize]))
        while (length(clussel)>0) {
            for (i in 1:length(clussel)) {
                sel= cl==clussel[i]
                centers= initKmeans(coord[sel,,drop=FALSE], nclus=2)
                cl[sel]= maxcl + kmeans(coord[sel,,drop=FALSE],centers=centers)$cluster
                tab[clussel[i]]= 0
                tab[maxcl+1:2]= table(cl[sel])
                names(tab)[maxcl+1:2]= as.character(maxcl+1:2)
                maxcl= maxcl+2
            }
            clussel= as.numeric(names(tab[tab>blocksize]))
        }
    }
    return(cl)
}



