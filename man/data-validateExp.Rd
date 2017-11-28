\name{validateExp}
\alias{validateExp}
\docType{data}
\title{
Data of the validation drug sensitivity screen using five additional drugs
}
\description{

To validate some of the associations observed in the main screen,
including the associations between IGHV status and HSP90 inhibitors
and the associations between trisomy 12 and MEK/ERK pathway inhibitors,
the effect of five additional drugs, cobimetinib (MEK inhibitor),
trametinib (MEK inhibitor), SCH772984 (ERK inhibitor), Ganetespib
(HSP90 inhibitor) and Onalespib (HSP90 inhibitor) were tested on 128 CLL
samples that were also used in the main screen. 

The results were stored in a tidy table (tibble) with four columns:

1) 'patientID' - The patient identifiers.

2) 'Drug' - The names of the drugs used in this screen.  

3) 'Concentration' - The concentrations of the drugs in the unit of uM.  

4) 'viab' - Viabilities of the samples after drug treatment, normalized by negative controls (DMSO).  
}
\author{Junyan Lu}
\usage{
validateExp
}
\format{
Tidy table with 3200 rows and 4 columns.
}
\keyword{datasets}
