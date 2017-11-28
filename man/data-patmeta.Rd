\name{patmeta}
\alias{patmeta}
\docType{data}
\title{
Meta data of the patient samples
}
\description{
This data set contains basic clinical information of patients who donated the
samples. Row names code for Patient IDs. The data.frame contains such
information as diagnosis ('Diagnosis'), sex ('Gender'), IGHV status ('IGHV'),
methyation cluster assignment ('ConsClust'), age of patient at which the
sample was taken ('Age4Main'). Moreover, the binary columns:
'treatedAfter' - \code{TRUE} if the patient
was treated after the sample was taken, 'died' - \code{TRUE} if the patient died,
'IC50beforeTreatment' - \code{TRUE} if the patient was treated before the sample was
taken. Column 'T5' includes time (in years) which passed from taking the sample
to the next treatment. Column 'T6' includes time (in years) which passed from
taking the sample to patients' death. 
}
\author{Malgorzata Oles}
\usage{
patmeta
}
\format{
data.frame with 265 rows and 10 columns
}
\keyword{datasets}
