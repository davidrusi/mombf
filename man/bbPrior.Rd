\name{bbPrior}
\alias{bbPrior}
\alias{unifPrior}
\alias{binomPrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Priors on model space for variable selection problems }
\description{
  \code{unifPrior}
  implements a uniform prior (equal a priori probability for all
  models). \code{binomPrior} implements a Binomial prior.
  \code{bbPrior} implements a Beta-Binomial prior. 
}
\usage{
unifPrior(sel, logscale=TRUE, groups=1:length(sel),
constraints=lapply(1:length(unique(groups)), function(z) integer(0)))

binomPrior(sel, prob=.5, logscale=TRUE, probconstr=prob, groups=1:length(sel),
constraints=lapply(1:length(unique(groups)), function(z) integer(0)))

bbPrior(sel, alpha=1, beta=1, logscale=TRUE, alphaconstr=alpha,
betaconstr=beta, groups=1:length(sel),
constraints=lapply(1:length(unique(groups)), function(z) integer(0)))
}
\arguments{
  \item{sel}{ Logical vector indicating which variables are included in
    the model}
  \item{logscale}{ Set to \code{TRUE} to return the log-prior
    probability.}
  \item{groups}{Group that each variable belongs to (e.g. dummy
    indicators for categorical variables with >2 categories). The idea
    is that all variables in a group are jointly added/removed from the
    model. By default all variables are assumed to be in separate groups}
  \item{constraints}{List with length equal to the number of groups
    (distinct elements in \code{groups}). Element j in the list should
    indicate any hierarchical constraints on the group, for instance
    constraints[[3]]==c(1,2) indicates that group 3 can only be included
  in the model if groups 1 and 2 are also in the model. This can be used
  to enforce that an interaction can only be in the model if the main
  effects are also in the model.
  }

  \item{prob}{ Success probability for the Binomial prior}
  \item{probconstr}{Success probability for the Binomial prior for
  groups that are subject to constraints}
  
  \item{alpha}{ First parameter of the Beta-Binomial prior, which is equivalent
    to specifying a Beta(alpha,beta) prior on \code{prob}.}
  \item{beta}{ First parameter of the Beta-Binomial prior, which is equivalent
    to specifying a Beta(alpha,beta) prior on \code{prob}.}

  \item{alphaconstr}{Same as alpha for the groups that are subject to constraints}
  \item{betaconstr}{Same as beta for the groups that are subject to constraints}
}
\value{
  Prior probability of the specified model
}
\author{ David Rossell }
\examples{
library(mombf)
sel <- c(TRUE,TRUE,FALSE,FALSE)
unifPrior(sel,logscale=FALSE)
binomPrior(sel,prob=.5,logscale=FALSE)
bbPrior(sel,alpha=1,beta=1,logscale=FALSE)
}
\keyword{ distribution }
