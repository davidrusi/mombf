context("nlpMarginal")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-5

patrick::with_parameters_test_that(
  "nlpMarginal is correctly impemented for normal family and",
  {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    true_covs <- which(theta3 != 0)
    ans_max <- nlpMarginal(true_covs, y3, X3, priorCoef=pCoef, priorVar=pVar)
    ans_all <- nlpMarginal(
      1:length(theta3), y3, X3, priorCoef=pCoef, priorVar=pVar
    )
    expect_true(ans_max > ans_all)
    expect_equal(ans_max, expected_max, tolerance=tolerance)
    expect_equal(ans_all, expected_all, tolerance=tolerance)
  },
  patrick::cases(
    momprior=list(pCoef=momprior(tau=0.328, r=1), expected_max=-37.47804, expected_all=-41.73017),
    imomprior=list(pCoef=imomprior(tau=0.328), expected_max=-38.35765, expected_all=-43.9903),
    emomprior=list(pCoef=emomprior(tau=0.328), expected_max=-37.31882, expected_all=-42.475),
    zellnerprior=list(pCoef=zellnerprior(tau=0.328), expected_max=-39.81047, expected_all=-39.97039)
  )
)
