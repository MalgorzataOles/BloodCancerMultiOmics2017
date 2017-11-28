\name{moround}
\alias{moround}
\title{
Round numbers to the ceiling of a given base
}
\description{
The function rounds the value up (either numeric or a numeric vector) to the
multiplication of the specified base.
}
\usage{
moround(x,base)
}
\arguments{
  \item{x}{numeric vector}
  \item{base}{numeric vector}
}
\details{
Both arguments could be either single numeric or numeric vectors.
Base argument should be either of length 1 or the divisible of the length of
argument x.
}
\value{
numeric vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
moround(x=c(1.23, 5, 5.1, 8), base=5)
moround(x=c(1.23, 5, 5.1, 8), base=c(2, 5))
}
