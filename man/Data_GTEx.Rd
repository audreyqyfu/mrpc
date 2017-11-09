\name{Data_GTEx}
\alias{Data_GTEx}

\title{
Data from GETx Portal for the GEUVADIS SNPs
}
\description{
Genotypes and gene expression data for the GEUVADIS SNPs among GTEx (The GTEx Consortium. 2013) samples.
}
\details{
First column is the genotype and the remaning are the gene expression. 

}
\value{
Matrix
}
\references{
1. The GTEx Consortium. (2013). The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45(6), 580–585.
}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{Data_GEUVADIS}
}
\examples{
#Data from GETx Portal for the GEUVADIS 
#SNPs 32 for 369 wholw blood sample
Data_GTEx$Q32$WB
#Data from GETx Portal for the GEUVADIS 
#SNPs 32 for 117 LCL sample
Data_GTEx$Q32$LCL
}
