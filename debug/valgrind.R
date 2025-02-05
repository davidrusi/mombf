library(mombf)
library(mvtnorm)

p <- 3
n <- 5000
set.seed(1)
Omega= diag(p)
Omega[abs(col(Omega) - row(Omega))==1]= 0.5
y = rmvnorm(n = n, sigma = solve(Omega))

niter= 100
burnin= 0
updates_per_iter= p; updates_per_column= p
pbirth= 0.75; pdeath= 1-pbirth

fitbd= modelSelectionGGM(y, sampler='Gibbs', niter=niter, Omegaini=Omega, global_proposal="none", updates_per_iter=updates_per_iter, updates_per_column=updates_per_column, pbirth= pbirth, pdeath=pdeath)

