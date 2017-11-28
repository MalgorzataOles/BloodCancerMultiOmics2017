\name{methData}
\alias{methData}
\docType{data}
\title{
DNA methylation data
}
\description{
The data set includes methylation data for the 5000 most variable CpG sites of
the CLL samples. The data was produced with the use of either 450k or 850k
methylation arrays.
Preprocessing of raw IDAT files was made using minfi R/Bioconductor package
 version 1.19.16. Intensities were normalized using the functional
normalization algorithm. CpG sites containing SNPs inside the probe body 
were removed.
}
\author{Andreas Mock, Malgorzata Oles}
\references{
Oakes CC, Seifert M, Assenov Y, Gu L, Przekopowitz M, Ruppert AS, Wang Q,
Imbusch CD, Serva A, Koser SD, et al. DNA methylation dynamics during
B cell maturation underlie a continuum of disease phenotypes in chronic
lymphocytic leukemia. Nat Genet. 2016;48(3):253-64
}
\usage{
methData
}
\format{
\code{"\linkS4class{RangedSummarizedExperiment}"} object with Features 5000 and Samples 196.
}
\keyword{datasets}
