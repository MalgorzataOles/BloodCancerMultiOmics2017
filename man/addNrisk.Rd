\name{addNrisk}
\alias{addNrisk}
\title{Add number-at-risk annotations to a plot}
\description{
Add number-at-risk (NAR) annotations to an existing survival plot, underneath the X-axis.
}
\usage{
addNrisk(x, at = axTicks(1), 
         line = 4, hadj = 0.5, 
         title = "Number at risk", title.adj = 0, 
         labels, hoff = 5, col = 1)
}
\arguments{
  \item{x}{A list as returned by \code{\link{survfit}}.}
  \item{at}{Time points at which the NAR values are calculated and placed.}
  \item{line}{Number of lines into the margin to start displaying the NAR. }
  \item{hadj}{Horizontal adjustment for the NAR values.}
  \item{title}{Optional title above the NAR.}
  \item{title.adj}{Text adjustment for the title}
  \item{labels}{Labels for each stratum.}
  \item{hoff}{Horizontal offset for the labels}
  \item{col}{Color for each stratum.}
}
\details{
This function was written and documented by Aron Charles Eklund in his package
survplot version 0.0.7.
}
\value{Invisibly, a matrix containing the number-at-risk values}
\author{Aron Charles Eklund (survplot version 0.0.7)}
\seealso{
See \code{\link{nrisk}} to retrieve number-at-risk values without plotting them.
See also \code{\link{survplot}}.
}
\examples{
library(survival)
s <- Surv(colon$time / 365, colon$status)

## Need to increase margins a bit
par(mar = c(10,6,2,1))

## no stratification
fit1 <- survfit(s ~ 1)
plot(fit1)
addNrisk(fit1)

## with stratification
fit2 <- survfit(s ~ rx, data = colon)
plot(fit2, xlab = 'Time (years)', ylab = 'Survival')
addNrisk(fit2)

}
\keyword{ survival }
\keyword{ aplot }
