\name{whichInGrob}
\alias{whichInGrob}
\title{
Return indices of layers of interest from the grob object
}
\description{
The function matches the supplied vector of grob's layer names to the grob object and returns the indices of those layer names.
}
\usage{
whichInGrob(grob, layer)
}
\arguments{
  \item{grob}{grob}
  \item{layer}{character vector}
}
\details{
If the layer doesn't exist the function returns NA.
}
\value{
numeric vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
library("ggplot2")
gg = ggplotGrob(qplot(1,1))
whichInGrob(gg, "xlab-b")
whichInGrob(gg, c("xlab-b","panel"))
}
