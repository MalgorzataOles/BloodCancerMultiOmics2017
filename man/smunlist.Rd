\name{smunlist}
\alias{smunlist}
\title{
Unlist with name preservation
}
\description{
Collapses list to a named vector with keeping the names as they were in the
lowest leaves in a list.
}
\usage{
smunlist(li)
}
\arguments{
  \item{li}{list}
}
\details{
The function works for the lists of multiple levels. These levels can be named,
unnamed, or mixture of both.
The names of the returned vector are preserved exactly as they were in a lowest
leaves of the list, which means that they can be duplicated.
}
\value{
named character vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
mylist = list(A=setNames(1:3, nm=letters[1:3]), B=list(D=3:4, setNames("a", nm=2)))
smunlist(mylist)
}