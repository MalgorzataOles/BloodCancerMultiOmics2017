\name{cytokineViab}
\alias{cytokineViab}
\docType{data}
\title{
Response of CLL to exposure to cytokines
}
\description{
The data set include the response measurements of 18 CLL patient samples exposed to six cytokines: IL-2, IL-10, IL-4, IL-21 (c1=0.001, c2=0.1, c3=10 ng/ul), LPS (c1=1, c2=10, c3=100 ng/ul) and IgM (c1=10 nM, c2=1, c3 = 10 uM) for 48 hours. Viability was measured using a CellTitre Glo assay, and luminescence was normalized to unstimulated controls. The results were stored in a tidy table (tibble) with 11 columns:
'Patient' is a patient sample ID,
'Timepoint' is a screening timepoint (48 h),
'Recording_date' is a date when the measurements were collected,
'Seeding_date' is a date when the experiment was started,
'Stimulation' is a name of cytokine used,
'Cytokine_Concentration' is a concentration of cytokine,
'Duplicate' is an information about the duplicates,
'Normalized_DMSO' is a drug response value after normalization by untreated control,
'mtor' is an information on whether the sample was classified into mtor group by our study,
'Edge' is an information of the position of the well respective to the whole screening plate,
'Cytokine_Concentration2 is again the concentration of the cytokine but in a different format.
}
\author{Sascha Dietrich, Malgorzata Oles}
\usage{
cytokineViab
}
\format{
tibble with 324 rows and 11 columns. 
}
\keyword{datasets}
