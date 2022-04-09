\name{IdentifyAssociatedPCs}
\alias{IdentifyAssociatedPCs}

\title{
Identifyprincipal components (PCs) that are significantly associated with eQTLs and genes
}
\description{
This function identifies PCs that are significantly associated the eQTLs or genes, and merge the associated PCs with the data on the eQTL and genes.  PCs may be derived from Principal Component Analysis (PCA) of the entire gene expression matrix, and may be viewed as potential confounders in the sequent causal network analysis on the eQTLs and genes.  See details in  Badsha and Fu (2019) and Badsha et al. (2018).
}
\usage{
IdentifyAssociatedPCs(PCs.matrix,no.PCs,data,fdr.level,corr.threshold,corr.value)
}

\arguments{
\item{PCs.matrix}{A matrix of PCs.}
\item{no.PCs}{Number of top PCs to test for association. The default is 10.}
\item{data}{Data of the eQTLs and genes, containing the genotypes of the eQTLs and the expression of the genes.}
\item{fdr.level}{(optional) The false discover rate (FDR) for association tests. Must be in (0,1]. The default is "0.05".}
\item{corr.threshold}{(optional). The default is "FALSE". If "TRUE" then a constraint on the correlation between a PC and an eQTL or a gene is applied in addition to the FDR control.}
\item{corr.value}{The threshold for the Pearson correlation between a PC and an eQTL or a gene when \code{corr.threshold} is "TRUE". The default is 0.3.}    
 
}

\value{
A \link{list} of object that containing the following:
\itemize{
\item \code{AssociatedPCs}: All the PCs that are significantly associated with the eQTLs and genes.
\item \code{data.withPC}: The data matrix that contains eQTLs, gene expression, and associated PCs.
\item \code{corr.PCs}: The matrix of correlations between PCs and eQTLs/genes.
\item \code{PCs.asso.list}: List of all associated PCs for each of the eQTLs and genes.
\item \code{qobj}: The output from applying the qvalue function.
}

}

\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899.
}
\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\seealso{
\link{data_GEUVADIS_combined}
}
\examples{

\dontrun{

# Load genomewide gene expression data in GEUVADIS 
# 373 individuals
# 23722 genes
data_githubURL <- "https://github.com/audreyqyfu/mrpc_data/raw/master/data_GEUVADIS_allgenes.RData"
load(url(data_githubURL))
PCs <- prcomp(data_GEUVADIS_allgenes,scale=TRUE)
# Extract the PCs matrix 
PCs.matrix <- PCs$x

# The eQTL-gene set contains eQTL rs7124238 and genes SBF1-AS1 and SWAP70

data_GEU_Q50 <- data_GEUVADIS$Data_Q50$Data_EUR
colnames(data_GEU_Q50) <- c("rs7124238","SBF2-AS1","SWAP70")
data <- data_GEU_Q50

# Identify associated PCs for this eQTL-gene set

Output <- IdentifyAssociatedPCs(PCs.matrix,no.PCs=10,data,fdr.level=0.05,corr.threshold=TRUE
,corr.value = 0.3)

# Gene SBF2-AS1 is significantly associated with PC2 
# Data with PC2 as a potential confounder

data_withPC <- Output$data.withPC

n <- nrow(data_withPC)         # Number of rows
V <- colnames(data_withPC)     # Column names

# Calculate Pearson correlation for MRPC analysis

suffStat <- list(C = cor(data_withPC,use = "complete.obs"),
                 n = n)

# Infer the graph by MRPC 

MRPC.fit_FDR<- MRPC(data_withPC,
                    suffStat,
                    GV = 1,
                    FDR = 0.05,
                    indepTest = 'gaussCItest',
                    labels = V,
                    FDRcontrol = 'LOND',
                    verbose = TRUE)
                    
plot(MRPC.fit_FDR, main="MRPC with PCs (potential confounders)")

}
}
