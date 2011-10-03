\name{hald}
\alias{hald}
\alias{x.hald}
\alias{y.hald}
\title{Hald Data}
\usage{
data(hald)
}
\description{
Montgomery and Peck (1982) illustrated variable selection techniques on the Hald cement data and gave several references to other analysis. The response variable \emph{y} is the \emph{heat evolved} in a cement mix. The four explanatory variables are ingredients of the mix, i.e., x1: \emph{tricalcium aluminate}, x2: \emph{tricalcium silicate}, x3: \emph{tetracalcium alumino ferrite}, x4: \emph{dicalcium silicate}. An important feature of these data is that the variables x1 and x3 are highly correlated (corr(x1,x3)=-0.824), as well as the variables x2 and x4 (with corr(x2,x4)=-0.975). Thus we should expect any subset of (x1,x2,x3,x4) that includes one variable from highly correlated pair to do as any subset that also includes the other member.
}

\format{
  \code{hald} is a matrix with 13 observations (rows) and 5 variables (columns), the first column is the dependent variable. \code{y.hald} and \code{x.hald} are also availables.
}

\source{
  Montgomery, D.C., Peck, E.A. (1982) 
  \emph{Introduction to linear regression analysis,}
  John Wiley, New York.
}

%\seealso{}
%\examples{}

\keyword{datasets}
