\name{mutCOM}
\alias{mutCOM}
\docType{data}
\title{
Genetic information of patient samples
}
\description{
This \code{"\linkS4class{NChannelSet}"} object contains genetic data for samples investigated in any
of the three experiments: whole exome sequencing, targeted sequencing or
fluorescent in situ hybridization. Object consists of one channel called
binary, with values: 0 if the mutation was absent, 1 if mutation was present
or NA if the mutation was not investigated. Feature data of the object
contains detailed information about mutation in TP53 and BRAF genes - the
variant(s) detected ('*_CDS' and '*_AA' columns) and the percentage at which
each variant was detected ('*_%' columns).
For TP53, BRAF, KRAS, del17p13, UMODL1, CREBBP, PRPF8 and trisomy12 mutation
an additional column 'cs' summarizes the clone size of
the mutated population. This value is a fraction at which the most abundant
variant is present in a sample.
}
\author{Malgorzata Oles}
\usage{
mutCOM
}
\format{
\code{"\linkS4class{NChannelSet}"} object with 89 genes (columns) and 265 patient samples (rows).
}
\keyword{datasets}
