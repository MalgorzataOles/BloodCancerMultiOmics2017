\name{exprTreat}
\alias{exprTreat}
\docType{data}
\title{
Gene expression before and after drug treatment
}
\description{
This \code{"\linkS4class{ExpressionSet}"} object contains microarray data for 12 patient samples
before and 12 h after treatment with everolimus, ibrutinib,
selumetinib, idelalisib and a negative control (chaetoglobosin A).
For annotation, please refer to \code{pData(exprTreat)}. The data underwent
variance stabilization (\code{\link[vsn]{vsn2}} function from vsn R/Bioconductor package)
and quantile normalization (\code{\link[limma]{normalizeQuantiles}} function from limma R/Bioconductor package).
}
\author{Sascha Dietrich}
\usage{
exprTreat
}
\format{
\code{"\linkS4class{ExpressionSet}"} object with 12 patient samples and 48107 features.
}
\keyword{datasets}
