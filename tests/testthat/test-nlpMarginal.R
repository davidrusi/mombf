context("Test nlpMarginal without groups")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-5

patrick::with_parameters_test_that(
  "nlpMarginal is correctly impemented for normal family and", {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    ans_max <- nlpMarginal(theta3_truth_idx, y3, X3, priorCoef=pCoef, priorVar=pVar)
    ans_all <- nlpMarginal(
      seq_along(theta3_truth), y3, X3, priorCoef=pCoef, priorVar=pVar
    )
    expect_true(ans_max > ans_all)
    expect_equal(ans_max, expected_max, tolerance=tolerance)
    expect_equal(ans_all, expected_all, tolerance=tolerance)
  },
  patrick::cases(
    momprior=list(pCoef=momprior(tau=0.328, r=1), expected_max=-37.47804, expected_all=-41.73017),
    imomprior=list(pCoef=imomprior(tau=0.328), expected_max=-38.35765, expected_all=-43.9903),
    emomprior=list(pCoef=emomprior(tau=0.328), expected_max=-37.31882, expected_all=-42.475),
    zellnerprior=list(pCoef=zellnerprior(tau=0.328), expected_max=-39.81047, expected_all=-39.97039),
    normalidprior=list(pCoef=normalidprior(tau=0.328), expected_max=-38.35193, expected_all=-39.4416)
  )
)

patrick::with_parameters_test_that(
  "nlpMarginal is correctly impemented for twopiecenormal family and", {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    ans_max <- nlpMarginal(
      theta3_truth_idx, y3, X3, family="twopiecenormal", priorCoef=pCoef,
      priorVar=pVar, priorSkew=pSkew
    )
    ans_all <- nlpMarginal(
      seq_along(theta3_truth), y3, X3, family="twopiecenormal", priorCoef=pCoef,
      priorVar=pVar, priorSkew=pSkew
    )
    expect_true(ans_max > ans_all)
    expect_equal(ans_max, expected_max, tolerance=tolerance)
    expect_equal(ans_all, expected_all, tolerance=tolerance)
  },
  patrick::cases(
    momprior=list(pCoef=momprior(tau=0.328, r=1), pSkew=momprior(tau=0.3, r=1), expected_max=-38.64731, expected_all=-44.05285),
    imomprior=list(pCoef=imomprior(tau=0.328), pSkew=imomprior(tau=0.3), expected_max=-39.69978, expected_all=-46.27032),
    emomprior=list(pCoef=emomprior(tau=0.328), pSkew=emomprior(tau=0.3), expected_max=-38.33845, expected_all=-44.76512)
  )
)

patrick::with_parameters_test_that(
  "nlpMarginal is correctly impemented for laplace family and",{
    pVar <- igprior(alpha=0.01, lambda=0.01)
    ans_max <- nlpMarginal(
      theta3_truth_idx, y3, X3, family="laplace", priorCoef=pCoef, priorVar=pVar
    )
    ans_all <- nlpMarginal(
      seq_along(theta3_truth), y3, X3, family="laplace", priorCoef=pCoef, priorVar=pVar
    )
    expect_true(ans_max > ans_all)
    expect_equal(ans_max, expected_max, tolerance=tolerance)
    expect_equal(ans_all, expected_all, tolerance=tolerance)
  },
  patrick::cases(
    momprior=list(pCoef=momprior(tau=0.328, r=1), expected_max=-38.25163, expected_all=-45.1899),
    imomprior=list(pCoef=imomprior(tau=0.328), expected_max=-39.13097, expected_all=-47.51524),
    emomprior=list(pCoef=emomprior(tau=0.328), expected_max=-37.92642, expected_all=-46.37662)
  )
)

patrick::with_parameters_test_that(
  "nlpMarginal is correctly impemented for twopiecelaplace family and", {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    ans_max <- nlpMarginal(
      theta3_truth_idx, y3, X3, family="twopiecelaplace", priorCoef=pCoef,
      priorVar=pVar, priorSkew=pSkew
    )
    ans_all <- nlpMarginal(
      seq_along(theta3_truth), y3, X3, family="twopiecelaplace", priorCoef=pCoef,
      priorVar=pVar, priorSkew=pSkew
    )
    expect_true(ans_max > ans_all)
    expect_equal(ans_max, expected_max, tolerance=tolerance)
    expect_equal(ans_all, expected_all, tolerance=tolerance)
  },
  patrick::cases(
    momprior=list(pCoef=momprior(tau=0.328, r=1), pSkew=momprior(tau=0.3, r=1), expected_max=-39.32618, expected_all=-45.58486),
    imomprior=list(pCoef=imomprior(tau=0.328), pSkew=imomprior(tau=0.3), expected_max=-40.31074, expected_all=-48.32312),
    emomprior=list(pCoef=emomprior(tau=0.328), pSkew=emomprior(tau=0.3), expected_max=-38.85072, expected_all=-47.06661)
  )
)

test_that(
  "logical or integer sel give equal results", {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    pCoef=momprior(tau=0.328, r=1)
    ans_int <- nlpMarginal(
      theta3_truth_idx, y3, X3, family="normal", priorCoef=pCoef, priorVar=pVar
    )
    ans_bool <- nlpMarginal(
      theta3_truth_bool, y3, X3, family="normal", priorCoef=pCoef, priorVar=pVar
    )
    expect_equal(ans_int, ans_bool)
  }
)
