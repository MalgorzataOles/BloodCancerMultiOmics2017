\name{deckel}
\alias{deckel}
\title{
treshold an array from below and above
}
\description{
treshold an array from below and above
}
\usage{
deckel(x, lower = -Inf, upper = +Inf)
}
\arguments{
  \item{x}{numeric matrix}
  \item{lower}{numeric}
  \item{upper}{numeric}
}
\details{
The function takes the matrix and censores the values from below (if lower param is set) or from above (if upper param is set), or from both of them. If neither lower nor upper param is set or if none of the values meet the criteria for thresholding, then function returns unmodified object.
}
\value{
matrix
}
\author{
Wolfgang Huber <wolfgang.huber@embl.de>
}
\examples{
mat = matrix(1:40, nrow=5)

# threshold values below 5
deckel(mat, lower=5)

# threshold values above 15
deckel(mat, upper=15)

# threshold values below 5 and above 15
deckel(mat, lower=5, upper=15)

# threshold values below 0 and above 50 -> no thresholding will be done!
identical(mat, deckel(mat, lower=0, upper=50))
}