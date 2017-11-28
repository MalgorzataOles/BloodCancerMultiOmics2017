\name{col2hex}
\alias{col2hex}
\title{
Converts color names with alpha to hex
}
\description{
The function takes the color names as specified in colors() together
with alpha levels and transforms it to hex representation. Optionally
it can also name the returned vector by names provided by the user.
}
\usage{
col2hex(cols, alpha=1, names=NA)
}
\arguments{
  \item{cols}{character vector}
  \item{alpha}{numeric, ranged 0-1}
  \item{names}{character vector, default NA}
}
\value{
numeric vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
col2hex(cols=c("hotpink","skyblue"), alpha=0.5)
col2hex(cols=c("hotpink","skyblue"), alpha=0.5, names = c("A","B"))
}
