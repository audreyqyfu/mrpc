\name{Data_GEUVADIS}
\alias{Data_GEUVADIS}

\title{
Data from GEUVADIS 
}
\description{
GEUVADIS (Lappalainen et al., 2013) data (i.e., gene expression measured in LCLs of a subset of individuals genotyped in the 1000 Genomes Project) from 373 Europeans and 89 Africans samples.
}

\details{
Among the most stringent set of eQTLs, ~70 have more than one target gene. Finally, we found 62 unique eQTLs discovered at the most stringent criteria exhibit pleiotropy. We extracted the genotypes of these 62 eQTLs and the expression of the target genes in the 373 Europeans and 89 Africans.
}
\value{
Matrix
}
\references{
1. Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506-511 (2013).
}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
#Data for 373 Europeans of eQTL 1
Data_GEUVADIS$Data_Q1$Data_EUR
#Data for 89 Africans of eQTL 1
Data_GEUVADIS$Data_Q1$Data_AFR
}

