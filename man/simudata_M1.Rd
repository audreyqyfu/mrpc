\name{simu.data_M1}
\alias{simu.data_M1}

\title{
Data for Model 1
}
\description{
Generating Simulated Data for Model 1
}
\details{
The column order of the input matrix is the genotype from  Single-nucleotide polymorphism (SNPs), indels,copy number variation (CNV) or Expression quantitative trait loci (eQTL) to genes. If we consider one SNPs/CNV/indels/CNV/eQTL, then the first column of the input matrix is SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data.
}
\value{
Matrix
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}
\seealso{
\link{SimulatedData} for simulated data generating function.
}

\examples{
## Model 1
simu.data_M1<-SimulatedData(N = 10^3,p = 0.45, 'model1',b0.1 = 0,b1.1
=1,b1.2 = 1,b1.3=1,sd.1 = 1)
}
