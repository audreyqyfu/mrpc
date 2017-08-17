\name{simu.data_starshaped }
\alias{simu.data_starshaped}
\title{
Data for Star Model
}
\description{
Generating simulated data for start model, where one gene has more than two child.
}

\details{
The column order of the input matrix is the genotype from  Single-nucleotide polymorphism (SNPs), indels,copy number variation (CNV) or Expression quantitative trait loci (eQTL) to genes. If we consider one SNPs/CNV/indels/CNV/eQTL, then the first column of the input matrix is SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data.
}
\value{
Matrix
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
## Star Model
simu.data_starshaped=SimulatedData(N = 10^3,p = 0.45,'starshaped',b0.1 = 0,b1.1 = 1,b1.2 = 1,
b1.3=1,sd.1 = 1)
}
