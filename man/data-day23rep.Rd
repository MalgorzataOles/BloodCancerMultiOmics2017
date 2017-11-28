\name{day23rep}
\alias{day23rep}
\docType{data}
\title{
Cell viability data for 3 replicated samples
}
\description{
This \code{"\linkS4class{NChannelSet}"} object contains normalized (to the negative control wells)
viability data for 48 h ('day2' channel) and 72 h ('day3' channel) incubation
period for the replicated experiment comprising 3 patient samples.
Patient samples are annotated in columns and compounds are annotated in rows.
The screen was performed for 67 drugs in 1-2 different drug concentrations
(16 drugs in 1 and 51 drugs in 2 concentration steps; see \code{fData(day23rep)}).  
}
\author{Malgorzata Oles}
\usage{
day23rep
}
\format{
\code{"\linkS4class{NChannelSet}"} object with 4 channels, 3 patient samples (columns) and 118 drugs (rows).
}
\keyword{datasets}
