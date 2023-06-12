# mombf

Model Selection with Bayesian Methods and Information Criteria

## Installation

``` r
# Install mombf from CRAN
install.packages("mombf")

# from GitHub:
# install.packages("devtools")
devtools::install_github("davidrusi/mombf")
```

## Quick start

The main Bayesian model selection (BMS) function is `modelSelection`. For information criteria
consider `bestBIC`, `bestEBIC`, `bestAIC`, `bestIC`. 
Bayesian model averaging (BMA) is also available for some models,
mainly linear and generalized linear models.
Local variable selection is implemented in `localnulltest` and `localnulltest_fda`.
Details are in [`mombf`'s vignette](https://CRAN.R-project.org/package=mombf/vignettes/mombf.pdf),
here we illustrate quickly how to get posterior model probabilities,
marginal posterior inclusion probabilities, BMA point estimates and posterior
intervals for the regression coefficients and predicted outcomes.

```r
library(mombf)
set.seed(1234)
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)

priorCoef <- momprior(tau=0.348)  # Default MOM prior on parameters
priorDelta <- modelbbprior(1,1)   # Beta-Binomial prior for model space
fit1 <- modelSelection(y ~ x[,1]+x[,2]+x[,3], priorCoef=priorCoef, priorDelta=priorDelta)
# Output
# Enumerating models...
# Computing posterior probabilities................ Done.
```

from here, we can also get the posterior model probabilities:

```r
postProb(fit1)
# Output
#    modelid family           pp
# 7      2,3 normal 9.854873e-01
# 8    2,3,4 normal 7.597369e-03
# 15   1,2,3 normal 6.771575e-03
# 16 1,2,3,4 normal 1.437990e-04
# 3        3 normal 3.240602e-17
# 5        2 normal 7.292230e-18
# 4      3,4 normal 2.150174e-19
# 11     1,3 normal 9.892869e-20
# 6      2,4 normal 5.615517e-20
# 13     1,2 normal 2.226164e-20
# 12   1,3,4 normal 1.477780e-21
# 14   1,2,4 normal 3.859388e-22
# 1          normal 2.409908e-25
# 2        4 normal 1.300748e-27
# 9        1 normal 2.757778e-28
# 10     1,4 normal 3.971521e-30
```

also the BMA estimates, 95% intervals, marginal posterior probability

```r
coef(fit1)
# Output
#              estimate        2.5%      97.5%      margpp
# (Intercept) 0.007230966 -0.02624289 0.04085951 0.006915374
# x[, 1]      1.134700387  0.93487948 1.33599873 1.000000000
# x[, 2]      1.135810652  0.94075622 1.33621298 1.000000000
# x[, 3]      0.000263446  0.00000000 0.00000000 0.007741168
# phi         1.100749637  0.83969879 1.44198567 1.000000000
```

and BMA predictions for y, 95% intervals

```r
ypred <- predict(fit1)
head(ypred)
# Output
#         mean       2.5%       97.5%
# 1 -0.8936883 -1.1165154 -0.67003262
# 2 -0.2162846 -0.3509188 -0.08331286
# 3  1.3152329  1.0673711  1.56348261
# 4 -3.2299241 -3.6826696 -2.77728625
# 5 -0.4431820 -0.6501280 -0.23919345
# 6  0.7727824  0.6348189  0.90977798
cor(y, ypred[,1])
# Output
#           [,1]
# [1,] 0.8468436
```

## Bug report

Please submit bug reports to the [issue tracker](https://github.com/davidrusi/mombf/issues).
