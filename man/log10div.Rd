\name{log10div}
\alias{log10div}
\title{
log10 with sign
}
\description{
The function calculates the log10 of the given value and returns
it together with the sign of the input value.
}
\usage{
log10div(x)
}
\arguments{
  \item{x}{numeric vector}
}
\details{
This function is useful when coloring p-values stratified by two
possible effect directions (sensitivity and resistance in our case).

}
\value{
numeric vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
log10div(x=c(-10,10))
}
