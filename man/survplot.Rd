\name{survplot}
\alias{survplot}
\title{ Draw augmented K-M survival curves  }
\description{ Plot Kaplan-Meier survival curves, automatically generate a key for each strata, and calculate and display hazard ratio if there are exactly two strata.  Optionally, indicate the number-at-risk below the main plot.}
\usage{
survplot(x, data = NULL, subset = NULL, 
         snames, stitle, 
         col, lty, lwd,
         show.nrisk = TRUE, color.nrisk = TRUE,
         hr.pos = 'topright', legend.pos = 'bottomleft', ...)
}
\arguments{
  \item{x}{ A formula, as would be appropriate for \code{\link[survival]{survfit}} and \code{\link[survival]{coxph}}. }
  \item{data, subset}{ Arguments passed to \code{\link[survival]{survfit}} and \code{\link[survival]{coxph}}. }
  \item{snames}{ Names for each stratum, to be used in the legend. If missing, these are inferred from the data. }
  \item{stitle}{ Title for the strata legend. If missing, this is inferred from \code{x}. }
  \item{col, lty, lwd}{ Colors, line type, and line width for each stratum (optional).  }
  \item{show.nrisk}{ Indicate the number-at-risk for each stratum below the plot? }
  \item{color.nrisk}{ Color the number-at-risk to match the plot? }
  \item{hr.pos}{ Where to put the hazard ratio information, or NA to omit (see \code{\link{legend}}) }
  \item{legend.pos}{ Where to put the legend, or NA to omit (see \code{\link{legend}}) }
  \item{\dots}{ Further parameters sent to \code{\link[survival]{plot.survfit}}. }
}
\details{
This function was written and documented by Aron Charles Eklund in his package
survplot version 0.0.7.

Hazard ratio (and 95\% confidence intervals) and logrank P are calculated and displayed if there are exactly two groups. 

If there is exactly one group (no stratification), the legend is omitted.
}
\value{ If there are exactly two groups, a character vector with the HR and P value is returned invisibly. }
\author{Aron Charles Eklund (survplot version 0.0.7)}
\note{ The lower figure margin is increased if the number-at-risk is displayed. }
\seealso{ \code{\link{nrisk}} }
\examples{ 
library(survival)
surv <- Surv(colon$time / 365, colon$status)

survplot(surv ~ rx,
  data = colon,
  lty = 1:3, 
  main = 'Patients stratified by treatment', 
  xlab = 'Time (Years)')

survplot(surv ~ colon$sex, 
  main = 'Patients stratified by sex', 
  xlab = 'Time (Years)',
  snames = c('F', 'M'),
  stitle = 'Gender')

survplot(surv ~ sex, 
  data = colon, 
  subset = colon$surg == 1)

## Example without stratification
survplot(surv ~ 1, data = colon)
  
}
\keyword{ survival }
\keyword{ hplot }
