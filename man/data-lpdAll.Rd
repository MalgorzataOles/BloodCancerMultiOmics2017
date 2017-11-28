\name{lpdAll}
\alias{lpdAll}
\docType{data}
\title{
An assembly of drug viability data, methylation clusters, important mutations and copy number variants 
}
\description{
Columns indicate patients and rows different omics features.

\code{pData()} contains some basic patient characteristics. 

\code{fData()$type} contains information to which omics data type each feature belongs to:


1) viab (viability values for n=448 data points): 'D_001' stands for the drug as coded in the object \code{\link{drugs}} and '_01' indicates the concentration step. '_1:5' corresponds to the average across all five concentration steps and '_4:5' corresponds to all the concentration steps .
2) gen (n=89): important gene mutations and copy number variants derived from WES, SNP arrays, FISH and targeted sequencing.       
3) Methylation_Cluster: The association of each CLL patient with one of the three Methylation Cluster was determined as described in the methods section.
4) IGHV mutation status for CLL patients was determined as described in the methods section.

}
\author{Wolfgang Huber}
\usage{
lpdAll
}
\format{
\code{"\linkS4class{ExpressionSet}"} with Features 539 and Samples 249.
}
\keyword{datasets}
