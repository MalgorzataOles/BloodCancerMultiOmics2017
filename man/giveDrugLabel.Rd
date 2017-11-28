\name{giveDrugLabel}
\alias{giveDrugLabel}
\title{
Convert extended drug IDs to drug names
}
\description{
The function converts drug IDs given in a format X1_X2_Y1 or X1_X2-Y2,
where X1_X2 is a drug id, Y1 is a number of drug concentration step and Y2
is a drug concentration, to format "Z Y2 µM", where Z is a drug name.
}
\usage{
giveDrugLabel(drid, ctab, dtab)
}
\arguments{
  \item{drid}{character vector}
  \item{ctab}{data frame}
  \item{dtab}{data frame}
}
\details{
The drug ID (X) has to be present in row names of dtab object.
ctab is a data frame with drug concentrations 
(columns are concentrations and rows are the drugs).
dtab is a data frame with drugs in the rows and at least one 
column with drug characteristics. Here the column "name" with the
name of the drug is needed.
}
\value{
character vector
}
\author{
Małgorzata Oleś <malgorzata.oles@embl.de>
}
\examples{
data("drugs","conctab")
giveDrugLabel(c("D_001-4", "D_002-0.02", "D_001_4", "D_002_1"), 
conctab, drugs)
}