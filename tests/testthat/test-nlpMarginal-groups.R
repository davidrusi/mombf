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
