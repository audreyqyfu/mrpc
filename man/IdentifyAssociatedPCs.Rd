\name{IdentifyAssociatedPCs}
\alias{IdentifyAssociatedPCs}

\title{
Identify principal components (PCs) that are significantly associated with eQTL-gene sets
}
\description{
We performed Principal Component Analysis (PCA) on the entire gene expression matrix, and used the top PCs as potential confounders. For a confounder to affect the causal inference, it should be associated with the genes of interest. We examined the association between each of the top PCs and the eQTL-gene sets, identified statistically significant associations, and then applied MRPC inference jointly for the eQTL-gene set and the associated PCs. See details in  Badsha and Fu, 2019 and Badsha et al., 2018.
}
\usage{
IdentifyAssociatedPCs(PCs.matrix,no.PCs,data,fdr.level,corr.threshold,corr.value)
}

\arguments{
\item{PCs.matrix}{PCA matrix on the entire gene expression.}
\item{no.PCs}{Number of top PCs as potential confounders. The default is 10.}
\item{data}{Data of interest used to identify associated PCs.}
\item{fdr.level}{(optional) A level at which to control the FDR. Must be in (0,1]. The defult "0.05", if this is selected, a vector of TRUE and FALSE is returned that specifies whether each q-value is less than fdr.level or not.}
\item{corr.threshold}{(optional). The default is "FALSE". If "TRUE" then the most significat impact of the PCs on eQTL-gene sets will be selected based on given correlation value. }
\item{corr.value}{Pearson correlation value when correlation threshold is "TRUE". The default is 0.3.}    
 
}

\value{
A \link{list} of object that containing the following:
\itemize{
\item \code{AssociatedPCs}: The associated PCs are significantly associated with eQTL-gene sets.
\item \code{data.withPC}: Data with the associated PCs (associated PCs are the additional nodes in the data) .
\item \code{corr.PCs}: The matrix of correlations with PCs  
\item \code{PCs.asso.list}: List of all associated PCs in each nodes.
\item \code{qobj}: A list of all objects after applied qvalue function.
}

}

\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899.
}
\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{data_GEUVADIS_combined}
}
\examples{

# Load genomewide gene expression data in GEUVADIS 
# 373 individuals
# 23722 genes
data_githubURL <- "https://github.com/audreyqyfu/mrpc_data/raw/master/data_GEUVADIS_allgenes.RData"
load(url(data_githubURL))
PCs <- prcomp(data_GEUVADIS_allgenes,scale=TRUE)
# Extract the PCs matrix 
PCs.matrix <- PCs$x

# Gene SBF2-AS1 is significantly associated with PC2 

data_GEU_Q50 <- data_GEUVADIS$Data_Q50$Data_EUR
colnames(data_GEU_Q50) <- c("rs7124238","SBF2-AS1","SWAP70")
data <- data_GEU_Q50

# Identify associated PCs

Output <- IdentifyAssociatedPCs(PCs.matrix,no.PCs=10,data,fdr.level=0.05,corr.threshold=TRUE
,corr.value = 0.3)

# Data with associated PCs as potential confounders

data_withPC <- Output$data.withPC

n <- nrow(data_withPC)         # Number of rows
V <- colnames(data_withPC)     # Column names

# Calculate Pearson correlation

suffStat <- list(C = cor(data_withPC,use = "complete.obs"),
                 n = n)

# Infer the graph by MRPC 

MRPC.fit_FDR<- MRPC(data_withPC,
                    suffStat,
                    GV = 1,
                    FDR = 0.05,
                    alpha = 0.05,
                    indepTest = 'gaussCItest',
                    labels = V,
                    FDRcontrol = TRUE,
                    verbose = TRUE)
                    
plot(MRPC.fit_FDR, main="MRPC with PCs (potential confounders)")

}

