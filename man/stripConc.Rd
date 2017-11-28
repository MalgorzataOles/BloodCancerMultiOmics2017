\name{stripConc}
\alias{stripConc}
\title{
Convert extended drug IDs to drug names
}
\description{
Out of drug IDs like D_001_1, it extracts the concentration step '_1'. 
}
\usage{
stripConc(x)
}
\arguments{
  \item{x}{character vector}
}
\details{
x has to be present in row names of drugs object.
}
\value{
character vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
data("drugs")
stripConc(c("D_001_1"))
}