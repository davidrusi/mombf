context("Test modelSelection with enumerate=TRUE")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 1e-5

patrick::with_parameters_test_that(
  "modelSelection without groups works for", {
    pDelta <- modelbbprior(1,1)
    log <- capture.output(
      fit1 <- modelSelection(y=y3, x=X3, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family=family, priorSkew=pCoef),
      fit2 <- modelSelection(y3~X3[,2]+X3[,3]+X3[,4], priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family=family, priorSkew=pCoef),
      fit3 <- modelSelection(as.formula("y~X2+X3+X4"), data=data.frame(X3, y=y3), priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family=family, priorSkew=pCoef)
    )
    pp1 <- postProb(fit1)
    pp2 <- postProb(fit2)
    pp3 <- postProb(fit3)
    expect_equal(pp1$modelid, pp2$modelid)
    expect_equal(pp1$modelid, pp3$modelid)
    expect_equal(pp1$pp, pp2$pp, tolerance=tolerance)
    expect_equal(pp1$pp, pp3$pp, tolerance=tolerance)
  },
  patrick::cases(
    mom_auto=list(family="auto", pCoef=momprior(tau=0.348)),
    mom_normal=list(family="normal", pCoef=momprior(tau=0.348)),
    mom_twopiecenormal=list(family="twopiecenormal", pCoef=momprior(tau=0.348)),
    mom_laplace=list(family="laplace", pCoef=momprior(tau=0.348)),
    mom_twopiecelaplace=list(family="twopiecelaplace", pCoef=momprior(tau=0.348)),
    imom_auto=list(family="auto", pCoef=imomprior(tau=0.348)),
    imom_normal=list(family="normal", pCoef=imomprior(tau=0.348)),
    imom_twopiecenormal=list(family="twopiecenormal", pCoef=imomprior(tau=0.348)),
    imom_laplace=list(family="laplace", pCoef=imomprior(tau=0.348)),
    imom_twopiecelaplace=list(family="twopiecelaplace", pCoef=imomprior(tau=0.348)),
    emom_auto=list(family="auto", pCoef=emomprior(tau=0.348)),
    emom_normal=list(family="normal", pCoef=emomprior(tau=0.348)),
    emom_twopiecenormal=list(family="twopiecenormal", pCoef=emomprior(tau=0.348)),
    emom_laplace=list(family="laplace", pCoef=emomprior(tau=0.348)),
    emom_twopiecelaplace=list(family="twopiecelaplace", pCoef=emomprior(tau=0.348))
  )
)

patrick::with_parameters_test_that(
  "modelSelection with groups works for", {
    pDelta <- modelbbprior(1,1)
    groups <- c(1, 1, 2, 2, 3, 4, 4)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE,
        family=family, priorSkew=pCoef, priorGroup=pCoef, groups=groups
      )
    )
    pprobs <- postProb(fit)
    expect_equal(length(pprobs$modelid), 16)
    expect_equal(as.character(pprobs$modelid[1]), "3,4,6,7")
  },
  patrick::cases(
    mom_normal=list(family="normal", pCoef=momprior(tau=0.348)),
    mom_twopiecenormal=list(family="twopiecenormal", pCoef=momprior(tau=0.348)),
    mom_laplace=list(family="laplace", pCoef=momprior(tau=0.348)),
    mom_twopiecelaplace=list(family="twopiecelaplace", pCoef=momprior(tau=0.348)),
    imom_normal=list(family="normal", pCoef=imomprior(tau=0.348)),
    imom_twopiecenormal=list(family="twopiecenormal", pCoef=imomprior(tau=0.348)),
    imom_laplace=list(family="laplace", pCoef=imomprior(tau=0.348)),
    imom_twopiecelaplace=list(family="twopiecelaplace", pCoef=imomprior(tau=0.348)),
    emom_normal=list(family="normal", pCoef=emomprior(tau=0.348)),
    emom_twopiecenormal=list(family="twopiecenormal", pCoef=emomprior(tau=0.348)),
    emom_laplace=list(family="laplace", pCoef=emomprior(tau=0.348)),
    emom_twopiecelaplace=list(family="twopiecelaplace", pCoef=emomprior(tau=0.348))
  )
)

patrick::with_parameters_test_that(
  "model space prior work in modelSelection:", {
    pCoef <- momprior(tau=0.348)
    log <- capture.output(
      fit <- modelSelection(y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family="normal")
    )
    expect_output(show(fit))
    pprobs <- postProb(fit)
    expect_true(any(pprobs$modelid[1:5] == "3,4,6,7"))
  },
  test_name=c("uniform", "binomial", "betabinomial", "complex"),
  pDelta=c(modelunifprior(), modelbinomprior(p=0.5), modelbbprior(alpha.p=1, beta.p=1), modelcomplexprior(c=1))
)

patrick::with_parameters_test_that(
  "asymetric binomial prior works in modelSelection (pDelta_pConstr):", {
    i0 <- integer(0)
    constraints <- list(i0, 1, i0, 3, i0, i0, 6)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorDelta=pDelta, priorConstraints=pConstr,
        enumerate=TRUE, constraints=constraints
      )
    )
    log <- capture.output(
      fit_asym <- modelSelection(
        y=y6, x=X6, priorDelta=modelbinomprior(p=c(0.7, 0.6, 0.2, 0.4)),
        priorConstraints=pConstr, enumerate=TRUE, constraints=constraints
      )
    )
    pprobs <- postProb(fit)
    pprobs_asym <- postProb(fit_asym)
    expect_true(pprobs[1, "modelid"] == "3,4,6,7")
    expect_true(pprobs_asym[1, "modelid"] == "3,4,6,7")
    expect_equal(pprobs[1, "pp"], .9957818, tolerance=tolerance)
    expect_true(abs(pprobs[1, "pp"] - pprobs_asym[1, "pp"]) > tolerance)
  },
  patrick::cases(
    vect_vect=list(pDelta=modelbinomprior(p=c(.5, .5, .5, .5)), pConstr=modelbinomprior(p=c(0.5, .5, .5))),
    vect_scalar=list(pDelta=modelbinomprior(p=c(.5, .5, .5, .5)), pConstr=modelbinomprior(p=0.5)),
    scalar_vect=list(pDelta=modelbinomprior(p=.5), pConstr=modelbinomprior(p=c(0.5, .5, .5))),
    scalar_scalar=list(pDelta=modelbinomprior(p=.5), pConstr=modelbinomprior(p=0.5))
  )
)

patrick::with_parameters_test_that(
  "modelSelection methods in normal family work:", {
    if (method == "Hybrid") {pCoef <- imomprior(tau=0.348)} else {pCoef <- momprior(tau=0.348)}
    pDelta <- modelbbprior(1,1)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family="normal",
        method=method, B=200, optimMethod=optimMethod, hess=hess
      )
    )
    expect_output(show(fit))
    pprobs <- postProb(fit)
    expect_true(any(pprobs$modelid[1:5] == "3,4,6,7"))
  },
  patrick::cases(
    auto=list(method="auto", optimMethod="CDA", hess="asymp"),
    laplace_cda=list(method="Laplace", optimMethod="CDA", hess="asymp"),
    laplace_cda_asympdiag=list(method="Laplace", optimMethod="CDA", hess="asympDiagAdj"),
    laplace_lma=list(method="Laplace", optimMethod="LMA", hess="asymp"),
    hybrid=list(method="Hybrid", optimMethod="CDA", hess="asymp"),
    mc=list(method="MC", optimMethod="CDA", hess="asymp"),
    bic=list(method="plugin", optimMethod="CDA", hess="asymp")
  )
)

patrick::with_parameters_test_that(
  "modelSelection methods in laplace family work:", {
    if (method == "Hybrid") {pCoef <- imomprior(tau=0.348)} else {pCoef <- momprior(tau=0.348)}
    pDelta <- modelbbprior(1,1)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=TRUE, family="laplace",
        method=method, B=200, optimMethod=optimMethod, hess=hess
      )
    )
    expect_output(show(fit))
    pprobs <- postProb(fit)
    expect_true(any(pprobs$modelid[1:5] == "3,4,6,7"))
  },
  patrick::cases(
    auto=list(method="auto", optimMethod="CDA", hess="asymp"),
    laplace_cda=list(method="Laplace", optimMethod="CDA", hess="asymp"),
    laplace_cda_asympdiag=list(method="Laplace", optimMethod="CDA", hess="asympDiagAdj"),
    laplace_lma=list(method="Laplace", optimMethod="LMA", hess="asymp"),
    hybrid=list(method="Hybrid", optimMethod="CDA", hess="asymp"),
    mc=list(method="MC", optimMethod="CDA", hess="asymp"),
    bic=list(method="plugin", optimMethod="CDA", hess="asymp")
  )
)
