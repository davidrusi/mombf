\name{priorp2g}
\alias{priorp2g}
\title{ Moment and inverse moment prior elicitation }
\description{
  \code{priorp2g} finds the \code{g} value giving \code{priorp} prior
  probability to the interval (\code{-q},\code{q}).
}
\usage{
priorp2g(priorp, q, nu=1, prior=c("iMom", "normalMom", "tMom"))
}
\arguments{

  \item{prior}{\code{prior=='normalMom'} does computations for the normal
    moment prior, \code{prior=='tMom'} for the T moment prior,
    \code{prior=='iMom'} does computations for the inverse moment prior.
    Currently \code{prior=='tMom'} is not implemented in \code{priorp2g}.}
  \item{q}{ \code{priorp2g} returns g giving \code{priorp} prior
    probability to the interval \code{(-q,q)}.}
  \item{nu}{ Prior degrees of freedom for the T moment prior or the iMom
    prior (ignored if \code{prior=='normalMom'}).}
  \item{priorp}{ \code{priorp2g} returns g giving \code{priorp} prior
    probability to the interval \code{(-q,q)}}

}
\details{
  See \code{pmom} and \code{pimom} for the MOM/iMOM cumulative
  distribution functions.
}
\value{
  \code{priorp2g} returns g giving \code{priorp} prior probability to the
  interval \code{(-q,q)}.
}
\references{
 See http://rosselldavid.googlepages.com for technical reports.
}
\author{
  David Rossell \email{rosselldavid@gmail.com}
}
\seealso{
  \code{\link{pmom}},
  \code{\link{pimom}}
}
\examples{
data(hald)
lm1 <- lm(hald[, 1] ~ hald[, 2] + hald[, 3] + hald[, 4] + hald[, 5])

#find g value giving 0.05 probability to interval (-.2,.2)
priorp <- .05; q <- .2
gmom <- priorp2g(priorp=priorp, q=q, prior='normalMom')
gimom <- priorp2g(priorp=priorp, q=q, prior='iMom')
gmom
gimom
}
\keyword{distribution}
\keyword{htest}

