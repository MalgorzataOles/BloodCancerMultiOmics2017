\name{meltWholeDF}
\alias{meltWholeDF}
\title{
Wide format to long format data conversion
}
\description{
The function converts wide format data which is either a data.frame or a matrix
(with dimnames present) to a long format structure. The output data.frame have
three columns: X, Y, and Measure. These are: column names, row names and
values of the input object, respectively.
}
\usage{
meltWholeDF(df)
}
\arguments{
  \item{df}{data.frame}
}
\details{
This function is particularly useful to prepare data for plotting with
ggplot2 package.
}
\value{
data.frame
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
df = data.frame(A=1:4, B=4:7, row.names=letters[1:4])
meltWholeDF(df)
}
