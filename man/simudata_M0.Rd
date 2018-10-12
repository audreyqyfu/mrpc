\name{simu.data_M0}
\alias{simu.data_M0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data for Model 0
}
\description{
Simulate Data for Model 0.
}

\details{
The column order of the input matrix is the genotype from Single-nucleotide polymorphism (SNPs), indels, copy number variation (CNV), or Expression quantitative trait loci (eQTL) to genes. If we consider one SNPs/CNV/indels/CNV/eQTL then the first column of the input matrix is the SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{SimulatedData} for simulated data generating function.
}

\examples{
## Model 0
simu.data_M0 <- SimulatedData(N = 10^3, p = 0.45,
                              'model0', b0.1 = 0,
                              b1.1 = 1, b1.2 = 1,
                              b1.3 = 1, sd.1 = 1)
}
