\name{dds}
\alias{dds}
\docType{data}
\title{
Gene expression data
}
\description{
The object contains the gene expression data after differential
gene expression analysis performed with DESeq2 R/Bioconductor package.
The preprocessing of the RNA-Seq data included
read alignment to the human reference genome (GRCh 37.1 / hg 19;
STAR version 2.3.0), and read counting done with htseq-count
(default mode union).
}
\author{Sascha Dietrich}
\references{
Love MI, Huber W, and Anders S. Moderated estimation of fold change
and dispersion for RNA-seq data with DESeq2.
Genome Biol. 2014;15(12):550

Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C,
Jha S, Batut P, Chaisson M, and Gingeras TR. STAR: ultrafast universal
RNA-seq aligner. Bioinformatics. 2013;29(1):15-21

Anders S, Pyl PT, and Huber W. HTSeq--a Python framework to work
with high-throughput sequencing data. Bioinformatics. 2015;31(2):166-9
}
\usage{
dds
}
\format{
\code{"\linkS4class{DESeqDataSet}"} object with 136 CLL samples and 63677 features.
}
\keyword{datasets}
