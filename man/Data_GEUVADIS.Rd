\name{data_GEUVADIS}
\alias{data_GEUVADIS}

\title{
GEUVADIS data with slelected 62 eQTL-gene sets 
}
\description{
The GEUVADIS (Lappalainen et al., 2013) data (i.e., gene expression) measured in Lymphoblastoid Cell Lines (LCLs) on a subset of individuals from the 1,000 Genomes Project including 373 Europeans and 89 Africans.
}

\details{
The GEUVADIS (Genetic European Variation in Disease) project identified eQTLs across the human genome. Among the most stringent set of eQTLs, ~70 have more than one target gene. Finally, we found 62 unique eQTLs discovered at the most stringent criteria exhibit pleiotropy. We extracted the genotypes of these 62 eQTLs and the expression of the target genes in the 373 Europeans and 89 Africans (see Badsha and Fu, 2019).
}
\value{
An object of \link{class} that contains slelected 62 eQTL-gene sets data for 373 Europeans and 89 Africans.
}

\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Lappalainen T, et al. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. Nature, 501, 506-511.

}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
# Data for 373 Europeans of eQTL #1
data_GEUVADIS$Data_Q1$Data_EUR

# Data for 89 Africans of eQTL #1
data_GEUVADIS$Data_Q1$Data_AFR
}
