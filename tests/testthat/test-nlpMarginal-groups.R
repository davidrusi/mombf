context("Test nlpMarginal with groups")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-5

patrick::with_parameters_test_that(
  "check_sel_groups works", {
    groups <- c(1,2,2,2,3,3,4,5)
    expect_error(check_sel_groups(sel, groups), error_msg)
  },
  patrick::cases(
    sel_outofbounds=list(sel=c(1,7,20), error_msg="sel larger than"),
    sel_outofbounds=list(sel=c(1,2,7), error_msg="incompatible with .+ groups"),
    sel_outofbounds=list(sel=c(2,3,4,7,8), error_msg=NA),
    sel_outofbounds=list(sel=c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE, FALSE), error_msg=NA)
  )
)

test_that(
  "nlpMarginal ignores sel argument when its a formula", {
    pCoef <- momprior(tau=0.348)
    ans1 <- nlpMarginal(seq(4), y3, X3, priorCoef=pCoef)
    formula <- y3~X3[,2]+X3[,3]+X3[,4]
    expect_warning(ans2 <- nlpMarginal(c(1,2), y=formula, priorCoef=pCoef), "ignoring sel")
    expect_equal(ans1, ans2)
    expect_warning(ans3 <- nlpMarginal(y=formula, priorCoef=pCoef), NA)
    expect_equal(ans1, ans3)
  }
)

patrick::with_parameters_test_that(
  "nlpMarginal with same kind of priorCoef and Group is correctly impemented for normal family:", {
    pVar <- igprior(alpha=0.01, lambda=0.01)
    groups <- c(1,2,3,4,5,5,5,6,6,6)
    ans_max <- nlpMarginal(theta9_truth_idx, y9, X9, family="normal", priorCoef=pCoef, priorVar=pVar)
    ans_max_group <- nlpMarginal(
      theta9_truth_idx, y9, X9, groups=groups, family="normal",
      priorCoef=pCoef, priorGroup=pCoef, priorVar=pVar
    )

    ans_all <- nlpMarginal(seq_along(theta9_truth), y9, X9, family="normal",priorCoef=pCoef, priorVar=pVar)
    ans_all_group <- nlpMarginal(
      seq_along(theta9_truth), y9, X9, groups=groups, family="normal",
      priorCoef=pCoef, priorGroup=pCoef, priorVar=pVar
    )
    expect_true(ans_max_group > ans_all_group)
    expect_equal(ans_max, ans_max_group)
    expect_equal(ans_all, ans_all_group)
  },
  test_name=c("mom", "imom", "emom", "zellner", "normalid"),
  pCoef=c(momprior(tau=0.35), imomprior(tau=0.35), emomprior(tau=0.35), zellnerprior(tau=0.35), normalidprior(tau=0.35))
)
