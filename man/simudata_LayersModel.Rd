\name{simu.data_LayersModel}
\alias{simu.data_LayersModel}
\title{
Data for Layers Model
}
\description{
The column order of the input matrix is the genotype from  Single-nucleotide polymorphism (SNPs), indels,copy number variation (CNV) or Expression quantitative trait loci (eQTL) to genes.
}

\details{
If we consider one SNPs/CNV/indels/CNV/eQTL, then the first column of the input matrix is SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data.}
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
## Layers Model
simu.data_LayersModel=SimulatedData(N = 10^3,p = 0.45, seed = 5,'layersmodel'
,b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3=1,sd.1 = 1)
}
