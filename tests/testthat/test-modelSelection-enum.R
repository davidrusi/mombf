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
    # mom_auto=list(family="auto", pCoef=momprior(tau=0.348)),
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
