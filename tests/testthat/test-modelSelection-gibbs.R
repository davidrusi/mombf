context("Test modelSelection with Gibbs")
library("mombf")

source(test_path("data-for-tests.R"))
tolerance <- 5e-3

patrick::with_parameters_test_that(
  "modelSelection without groups works for", {
    pDelta <- modelbbprior(1,1)
    log <- capture.output(
      fit1 <- modelSelection(y=y3, x=X3, priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE, family=family, priorSkew=pCoef),
      fit2 <- modelSelection(y3~X3[,2]+X3[,3]+X3[,4], priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE, family=family, priorSkew=pCoef),
      fit3 <- modelSelection(as.formula("y~X2+X3+X4"), data=data.frame(X3, y=y3), priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE, family=family, priorSkew=pCoef)
    )
    pp1 <- postProb(fit1)
    pp2 <- postProb(fit2)
    pp3 <- postProb(fit3)
    expect_equal(as.character(pp1$modelid[1]), as.character(pp2$modelid[1]))
    expect_equal(as.character(pp1$modelid[1]), as.character(pp3$modelid[1]))
    expect_equal(as.numeric(pp1$pp[1]), as.numeric(pp2$pp[1]), tolerance=tolerance)
    expect_equal(as.numeric(pp1$pp[1]), as.numeric(pp3$pp[1]), tolerance=tolerance)
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
  "modelSelection with groups (pCoef=pGroup) works for", {
    pDelta <- modelbbprior(1,1)
    groups <- c(1, 1, 2, 2, 3, 4, 4)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE,
        family=family, priorSkew=pCoef, priorGroup=pCoef, groups=groups
      )
    )
    pprobs <- postProb(fit)
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
#    emom_twopiecenormal=list(family="twopiecenormal", pCoef=emomprior(tau=0.348)),
    emom_laplace=list(family="laplace", pCoef=emomprior(tau=0.348)),
    emom_twopiecelaplace=list(family="twopiecelaplace", pCoef=emomprior(tau=0.348))
  )
)

patrick::with_parameters_test_that(
  "modelSelection with groups (pCoef!=pGroup), crossprodmat and covariancemat work for", {
    pDelta <- modelbbprior(1,1)
    groups <- c(1, 1, 2, 2, 3, 4, 4)
    log <- capture.output(
      fit <- modelSelection(
        y=y6, x=X6, priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE, XtXprecomp=FALSE,
        family=family, priorSkew=pCoef, priorGroup=pGroup, groups=groups
      )
    )
    pprobs <- postProb(fit)
    expect_true("3,4,6,7" %in% pprobs$modelid[1:4])
  },
  patrick::cases(
    normid_gzell=list(family="normal", pCoef=normalidprior(tau=0.348), pGroup=groupzellnerprior(tau=0.4)),
    zell_gzell=list(family="normal", pCoef=zellnerprior(tau=0.348), pGroup=groupzellnerprior(tau=0.4))
  )
)


test_that(
  "modelSelection with smoothterms and no groups works", {
    pCoef=momprior(tau=0.348)
    pDelta <- modelbbprior(1,1)
    log <- capture.output(
      fit <- modelSelection(
        y6~X6[,2:7], priorCoef=pCoef, priorDelta=pDelta, enumerate=FALSE,
        smoothterms= ~ X6[,2:7], family="normal", priorSkew=pCoef, priorGroup=pCoef
      )
    )
    pprobs <- postProb(fit)
    expect_equal(names(pprobs)[3], "pp")

    #expect_true(any(pprobs$modelid[1:5] == "3,4,6,7"))
  }
)

test_that(
  "modelSelection with includevars works", {
    log <- capture.output(
      fit <- modelSelection(
        y9, X9, enumerate=FALSE, includevars=c(rep(TRUE, 3), rep(FALSE, 7))
      )
    )
    models <- postProb(fit)$modelid
    expect_true(all(grepl("1,2,3", models)))
  }
)
