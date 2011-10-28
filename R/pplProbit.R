## Evaluate Posterior Predictive Loss under a probit model 
## Input:
## - fit: probit model fit, e.g. as returned by pmomPM
## - x: covariates used to fit the model
## - xadj: adjustment covariates
## - y: response variables (e.g. 0/1 vector or TRUE/FALSE)
## - kPen: Loss is Dev(yp,a) + kPen*Dev(yobs,a), where yp: draw from post predictive, yobs: observed data and a is E(yp|yobs).
##         i.e. kPen is a penalty term specifying the relative importance of deviations from the observed data. 
## Ouput
## - D: P + G
## - P: P_k(m) - (Penalty), i.e. sum(hm - h(mu))
## - G: G_k(m) - (Fit)
pplProbit <- function(fit, x, xadj, y, kPen=1){
  m <- nrow(fit$postModel); n <- nrow(x);
# generate predictive distribution -------------------------------------
  yp <- matrix(0, nrow=m, ncol=n);
  for(i in 1:m){
    th1 <- fit$postCoef1[i,]; th2 <- fit$postCoef2[i,];
    lpred <- as.vector(x %*% matrix(th1,ncol=1) + xadj %*% matrix(th2,ncol=1));
    p <- pnorm(lpred);
    yp[i,] <- rbinom(n, 1, p);
  }
# Compute ppl ----------------------------------------------------------
  h <- function(z) (z + 0.5)*log(z + 0.5) + (1.5-z)*log(1.5 - z);
  msize <- mean(apply(fit$postModel, 1, sum)) + ncol(xadj);
  if (kPen=='msize') kPen <- msize
# P_k(m) - (Penalty) ---------------------------------------------------
  mu    <- apply(yp, 2, mean, na.rm=TRUE);
  hi    <- h(yp);
  hm    <- apply(hi, 2, mean, na.rm=TRUE);
  P     <- sum(hm - h(mu));
# G_k(m) - (Fit) -------------------------------------------------------
  Gm    <- (h(mu) + kPen*h(y))/(kPen+1) - h((mu + kPen*y)/(kPen+1));
  G     <- (kPen+1)*sum(Gm, na.rm=TRUE); 
# D_k(m) ---------------------------------------------------------------
  D = P + G;
# Return Output -------------------------------------------------------- 
  return(list(d=D, g=G, p=P, msize=msize));
}           
# END ppl() ----------------------------------------------------------------------------------- # 
