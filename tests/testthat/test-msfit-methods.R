context("Test msfit methods")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-6

test_that(
  "msfit show method works", {
    log <- capture.output(fit <- modelSelection(y=y3, x=X3))
    expect_output(show(fit), "msfit object")
  }
)

test_that(
  "msfit coef method works", {
    log <- capture.output(fit <- modelSelection(y=y3, x=X3))
    fit_coef <- coef(fit)
    expect_equal(names(fit_coef[1,])[4], "margpp")
    expect_true(all(fit_coef[,4]>=0))
    expect_true(all(fit_coef[,4]<=1))
  }
)

test_that(
  "msfit predict method works", {
    log <- capture.output(
      fit <- modelSelection(y3~X3[,2]+X3[,3]+X3[,4]),
      ypred <- predict(fit)
    )
    expect_equal(names(ypred[1,])[1], "mean")
    expect_lt(mean((predict(fit)[,1]-y3)^2), 2)
  }
)

patrick::with_parameters_test_that(
  "msfit postProb method works", {
    log <- capture.output(fit <- modelSelection(y=y3, x=X3, enumerate=FALSE))
    pprobs <- postProb(fit, method=method)
    pprobs_cut <- postProb(fit, nmax=5, method=method)
    expect_equal(nrow(pprobs_cut), 5)
    expect_equal(pprobs[1:5,3], pprobs_cut[,3], tolerance=tolerance)
  },
  method=c("norm", "exact"),
  test_name=method
)
