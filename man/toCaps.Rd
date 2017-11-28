\name{toCaps}
\alias{toCaps}
\title{
Capitalize first character
}
\description{
The function capitalizes the first character of the given string or every
element of the character vector.
}
\usage{
toCaps(word)
}
\arguments{
  \item{word}{character vector}
}
\value{
character vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
toCaps("abc")
toCaps(c("abc", "Abc", "aBC", "ABC", "4you"))
}