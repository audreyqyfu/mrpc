\name{simu.data_layered}
\alias{simu.data_layered}
\title{
Data for the Layered Model
}

\description{
Simulate Data for the Layered Model.
}

\details{
The column order of the input matrix is the genotype from Single-nucleotide polymorphism (SNPs), indels, copy number variation (CNV), or Expression quantitative trait loci (eQTL) to genes. If we consider one SNPs/CNV/indels/CNV/eQTL then the first column of the input matrix is the SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>, Evan A Martin <mart9986@vandals.uidaho.edu>, and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}

\seealso{
\link{SimulatedData} for simulated data generating function.
}
\examples{
## Layered Model
simu.data_layered <- SimulatedData(N = 10^3, p = 0.45,
                                   'layered', b0.1 = 0,
                                   b1.1 = 1, b1.2 = 1,
                                   b1.3 = 1, sd.1 = 1)
}
