\name{data_GEUVADIS}
\alias{data_GEUVADIS}

\title{
GEUVADIS data with 62 eQTL-gene sets 
}
\description{
The GEUVADIS (Lappalainen et al., 2013) data (i.e., gene expression) measured in Lymphoblastoid Cell Lines (LCLs) on a subset of individuals from the 1,000 Genomes Project including 373 Europeans and 89 Africans.
}

\details{
The GEUVADIS (Genetic European Variation in Disease) project identified eQTLs across the human genome. Among these eQTLs, ~70 have more than one target gene. Additionally, we found 62 unique eQTLs which exhibit pleiotropy. We extracted the genotypes of these 62 eQTLs and the expression of the target genes for 373 Europeans and 89 Africans (see Badsha and Fu, 2019; Badsha et al., 2021).
}
\value{
A list that contains 62 eQTL-gene sets data for 373 Europeans and 89 Africans.
}

\references{
1. Lappalainen T, et al. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. Nature, 501, 506-511.

2. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10:460.

3. Badsha MB, Martin EA and Fu AQ (2021). MRPC: An R package for inference of causal graphs. Frontiers in Genetics, 10:651812.
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\examples{
# Data for 373 Europeans of eQTL #1
data_GEUVADIS$Data_Q1$Data_EUR

# Data for 89 Africans of eQTL #1
data_GEUVADIS$Data_Q1$Data_AFR
}
