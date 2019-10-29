\name{data_GEUVADIS_wPCs}
\alias{data_GEUVADIS_wPCs}
\docType{data}
\title{
GEUVADIS data with Principal Components (PCs)
}
\description{
We performed Principal Component Analysis (PCA) on the entire gene expression
matrix in GEUVADIS (Lappalainen et al., 2013), and used the top PCs as potential confounders.
}


\details{
A confounder to affect the causal inference, it should be associated with the genes of interest. We examined the association between each of the top 10 PCs and the eQTL-gene sets, identified statistically significant associations, and then applied MRPC jointly for the eQTL-gene set and the associated PC. See details in Badsha and Fu, 2019.
}

\references{
1. Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506-511 (2013).

2. Md Bahadur Badsha and Audrey Qiuyan Fu. Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics 2019, 10(460).
}
\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{

#Examining top PCs as potential confounders in GEUVADIS data

#load data with all genes in GEUVADIS (sample=373 and genes=23722)
data_githubURL <- "https://github.com/audreyqyfu/mrpc_data/raw/master/data_GEUVADIS_allgenes.RData"

load(url(data_githubURL))

#Apply pca
library(stats) #library for prcomp
PCs <- prcomp(data_GEUVADIS_allgenes,scale=TRUE)
#extract pca matrix 
PCs_matrix <- PCs$x

#load GEUVADIS 62 eQTL-gene sets data (sample=373 and variables=194 (eQTLs=62 and
#genes=132))
data("data_GEUVADIS_combined")  

library(psych) #library for corr.test
#correlations and pvalues using top 10 PCs
no_PCs <- 10
corr_PCs <- corr.test(PCs_matrix[,1:no_PCs],data_GEUVADIS_combined)
#The matrix of correlations
corr_matrix <- corr_PCs$r
# All p-values
Pvalues <- corr_PCs$p

#Statistically significant associations between each of the top 10 PCs and the eQTL-gene sets
library(WGCNA) #library for qvalue
qobj <- qvalue(Pvalues, fdr.level=0.05,robust = TRUE) 

#Significant associations list
Significant_asso=qobj$significant
List_significant_asso <- which(Significant_asso, arr.ind = TRUE, useNames = TRUE)
#1st columns are the PCs and 2nd columns are the corresponding associated eQTLs/genes.
List_significant_asso[1:10,]

#Example list of eQTLs or genes names with statistically associated with PC2
eqtl.genes_PC2 <- colnames(data_GEUVADIS_combined)[List_significant_asso
                            [which(List_significant_asso[,1]=="2"),2]]
#For example: gene SBF2-AS1 statistically associated with PC2 
#Note, gene SBF2-AS1 is a member of eQTL 50 with snp rs7124238 and gene SWAP70 in
#GEUVADIS. 
eqtl.genes_PC2[24]

#Then we applied MRPC jointly for the eQTL-gene set and the associated PC2 (rs7124238,
#SWAP70, SBF2-AS1, PC2)
}

