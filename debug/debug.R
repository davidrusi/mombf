#Uncomment the 2 lines below to re-compile all C code each time
#devtools::clean_dll()
#devtools::load_all()
#Rcpp::compileAttributes("~/github/mombf")

library(mombf)
library(mvtnorm)
set.seed(1)

## GGM examples
#p= 5
#Th= diag(p); Th[1,2]= Th[2,1]= 0.5
#sigma= solve(Th)
#z= matrix(rnorm(1000*p), ncol=p)
#y= z %*% chol(sigma)
#Omegaini= Th
#fit= modelSelectionGGM(y, scale=FALSE, Omegaini=Omegaini, almost_parallel="none", sampler='Gibbs', niter=100)

## Jack's example
library(mnormt)
p <- 5
n <- 500
N_rep <- 1
set.seed(3)
size <- p
#
G <- diag(1, p)
diag(G[,-1][1:size-1, 1:size-1]) <- rep(1,size-1)
diag(G[-1,][1:size-1, 1:size-1]) <- rep(1,size-1)
#
min_rho <- 0.1
max_rho <- 0.5
#
Rho <- diag(0.5, p)
Rho[lower.tri(Rho, diag = FALSE)] <- sample(c(-1, 1), p*(p-1)/2, replace = TRUE)*
  sample(seq(min_rho, max_rho, by = 0.1), p*(p-1)/2, replace = TRUE)*
G[lower.tri(G, diag = FALSE)]/2
Rho <- Rho + t(Rho)
sqrt_omega_jj <- sqrt(rgamma(p, shape = 3, rate = 1))
Theta <- outer(sqrt_omega_jj, sqrt_omega_jj)*Rho
#  
y = rmvnorm(n = n, sigma = solve(Theta))
Theta_sim <- Theta
data_y <- y
#
prior_params_1p50 <- list("s_slab" = 1.78, "w_slab" = 2/(p-1), "lambda" = 0.01)
warmup <- 0
iterations <- 100
### Serial Gibbs
updates_per_iter <- p
updates_per_column <- p
## mombf_SerialGibbs
fit <- modelSelectionGGM(data_y, sampler='Gibbs', Omegaini=Theta, niter=iterations+warmup, burnin=warmup, scale=FALSE, almost_parallel='none', updates_per_iter = updates_per_iter, updates_per_column = updates_per_column, priorCoef=normalidprior(prior_params_1p50$s_slab^2), priorModel=modelbinomprior(prior_params_1p50$w_slab), priorDiag=exponentialprior(lambda=prior_params_1p50$lambda))
