\name{defineResponseGroups}
\alias{defineResponseGroups}

\title{divides patients into response groups: BTK, MEK, mTOR, non-responders}
\usage{
defineResponseGroups(lpd)
}
\arguments{
  \item{lpd}{lpd object with comprehensive data}
}
\description{
The function divides patients into 4 groups depending on their mean response to the two lowest concentrations of BTK inhibitor (ibrutinib), mTOR inhibitor (everolimus) and MEK inhibitor (selumetinib). Division is done by looking at the distribution of viabilities for the three drugs mentioned above and using the mirror method to derive, first, a measure of the background variation of the values for these drugs (`ssd') and then define a cutoff as multiple (`z_factor') of that. The mirror method assumes that the observed values are a mixture of two components:

\itemize{
\item  a null distribution, which is symmetric about 1, and
\item  responder distribution, which has negligible mass above 1.
}

The choice of `z_factor' is a crucial step, because it determines the trade-off between falsely called responders (false positives) versus falsely called non-responders (false negatives). Under normality assumption, it is related to the false positive rate (FPR) by

\code{FPR = 1 - pnorm(z)}

An FPR of 0.05 thus corresponds to 

\code{z_factor <- qnorm(0.05, lower.tail = FALSE)}

The threshold is then calculated by:
\code{1 - z_factor * ssd}

Each patient is then assigned to a group as follows. If the response to ibrutinib was lower than the calculated threshold, we assign patient to BTK group. If not, we check the drug response value to everlimus in the same fashion. If still the value is not lower than the threshold, the procedure is repeated for selumetinib. If none of the responses mentioned above is below the threshold, we assign patient to the non-responder group.
}
\value{
data.frame
}
\author{
Wolfgang Huber <wolfgang.huber@embl.de>, Małgorzata Oleś <malgorzata.oles@embl.de>
}
