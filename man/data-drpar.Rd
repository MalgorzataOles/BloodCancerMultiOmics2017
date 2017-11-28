\name{drpar}
\alias{drpar}
\docType{data}
\title{
Cell viability data from the high-throughput drug screen
}
\description{
This \code{"\linkS4class{NChannelSet}"} object contains normalized (to the control wells) viability
data for 48 h incubation period within the drug screen. Patient samples are
annotated in columns and drugs are annotated in rows. Seven channels are
available: 'viaraw.1', 'viaraw.2', 'viaraw.3', 'viaraw.4', 'viaraw.5' - containing viability
information for drug concentrations from c1 (highest) to c5 (lowest)
respectively (see also \code{\link{conctab}}), and 'viaraw.1_5', 'viaraw.4_5' - 
containing the mean viability of all five concentrations and the two lowest
concentrations used, respectively. \code{pData} contains two columns: 'PatientID'
and 'ExpDate'. The second one contains the date at which the ATP content of the
wells after 48 h of incubation was measured.
}
\author{Malgorzata Oles}
\usage{
drpar
}
\format{
\code{"\linkS4class{NChannelSet}"} object with 249 patient samples (columns) and 64 drugs (rows).
}
\keyword{datasets}
