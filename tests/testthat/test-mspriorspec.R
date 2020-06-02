context("Test msPriorSpec initializers")
library("mombf")

source(test_path("data-for-tests.R"))

patrick::with_parameters_test_that(
  "msPriorSpec can be initialized with tau:", {
    pr_default <- priorfun()
    pr_taustd <- priorfun(taustd=2)
    pr_tau <- priorfun(tau=.3)
    expect_equal(as.double(pr_default@priorPars['taustd']), 1)
    expect_equal(as.double(pr_taustd@priorPars['taustd']), 2)
    expect_equal(as.double(pr_tau@priorPars['tau']), .3)
  },
  test_name=c("momprior", "zellnerprior", "normalidprior", "groupmomprior", "groupzellnerprior"),
  priorfun=c(momprior, zellnerprior, normalidprior, groupmomprior, groupzellnerprior)
)
