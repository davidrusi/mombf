library(mombf)
set.seed(1)
Th = diag(4)
Th[abs(col(Th) - row(Th))==1] = 0.5 # nolint
y = scale(rmvnorm(10^2, sigma=solve(Th)), center=TRUE, scale=FALSE)
fit.ap <- modelSelectionGGM(y, sampler='birthdeath', Omegaini='glasso-ebic', niter=100, burnin=0, scale=FALSE, almost_parallel=TRUE)

