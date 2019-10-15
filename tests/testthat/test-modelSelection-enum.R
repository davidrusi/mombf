context("Test modelSelection with enumerate=TRUE")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-5

patrick::with_parameters_test_that(
  "modelSelection with pmom works for", {
    pCoef <- momprior(tau=0.348)  # Default MOM prior on parameters
    pDelta <- modelbbprior(1,1)   # Beta-Binomial prior for model space
    fit1 <- modelSelection(y=y3, x=X3, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family=family, priorSkew=pCoef)
  },
  family=c("normal", "twopiecenormal", "laplace", "twopiecelaplace"),
  test_name=family
)
