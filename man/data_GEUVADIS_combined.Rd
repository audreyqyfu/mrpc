\name{data_GEUVADIS_combined}
\alias{data_GEUVADIS_combined}
\docType{data}
\title{
Combined genotype and gene expression data from 62 eQTL-gene sets in 373 Europeans from GEUVADIS
}
\description{
The genotype and gene expression data of 62 eQTL-gene sets in 373 Europeans from the GEUVADIS consortium (Lappalainen et al., 2013) are combined into one data matrix.  Each of these eQTLs has been identified to be associated with more than one gene (see details in Badsha and Fu, 2019).
}

\details{
The data set contains 373 samples in rows and 194 variables (62 eQTLs and 132 genes) in columns. Specifically, the columns are: eQTL1, gene1 for eQTL1, gene2 for eQTL1, eQTL2, gene1 for eQTL2, gene2 for eQTL2 and so on. 
  
For analysis, we account for potential confounding variables as additional nodes in the graph.  To do so, we first perform Principal Component Analysis (PCA) on the entire gene expression matrix from the European samples in GEUVADIS, and extract the top 10 PCs as potential confounding variables.  We next examine the statistical association between each of the top PCs and the eQTL-gene sets, and identify statistically significant associations (accounting for multiple testing with the q vlaue method).  We then apply MRPC to each eQTL-gene set with its associated PCs. See details in the examples below. Also see Badsha and Fu (2019) and Badsha et al. (2018).
}

\value{
Matrix
}

\references{
1. Lappalainen T, et al. (2013). Transcriptome and genome sequencing uncovers functional variation in humans. Nature, 501, 506-511.
  
2. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

3. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899.
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\seealso{
\link{data_GEUVADIS} 
}

\examples{

\donttest{
# Examining principal components (PCs) as potential confounders in analysis of the GEUVADIS data
  
library(MRPC) # MRPC

# Load genomewide gene expression data in GEUVADIS 
# 373 individuals
# 23722 genes
data_githubURL <- "https://github.com/audreyqyfu/mrpc_data/raw/master/data_GEUVADIS_allgenes.RData"
load(url(data_githubURL))

# Run PCA
library(stats) # prcomp
PCs <- prcomp(data_GEUVADIS_allgenes,scale=TRUE)
# Extract the PCs 
PCs_matrix <- PCs$x

# Load the 62 eQTL-gene sets
# 373 individuals
# 194 variables (eQTLs=62 and genes=132)
data("data_GEUVADIS_combined")  

# Identify PCs that are significantly associated with eQTL-gene sets
# Compute the correlation and corresponding p values between the top PCs and the eQTLs and genes
library(psych) # to use corr.test
no_PCs <- 10   
corr_PCs <- corr.test(PCs_matrix[,1:no_PCs],data_GEUVADIS_combined)
# The correlation matrix
corr_matrix <- corr_PCs$r
# The p values
Pvalues <- corr_PCs$p
# Apply the q value method at FDR of 0.05
library(WGCNA) # qvalue
qobj <- qvalue(Pvalues, fdr.level=0.05,robust = TRUE) 

# Significant associations
Significant_asso <- qobj$significant
List_significant_asso <- which(Significant_asso, arr.ind = TRUE, useNames = TRUE)
# 1st column contains the PCs
# 2nd column contains the associated eQTLs or genes
List_significant_asso[1:10,]

# Examples of eQTLs or genes that are significantly associated with selected PCs
# PC1
eqtl.genes_PC1 <- colnames(data_GEUVADIS_combined)[List_significant_asso 
                           [which(List_significant_asso[,1]=="1"),2]]
print(eqtl.genes_PC1)
# PC2
eqtl.genes_PC2 <- colnames(data_GEUVADIS_combined)[List_significant_asso 
                           [which(List_significant_asso[,1]=="2"),2]]
print(eqtl.genes_PC2)
# PC3
eqtl.genes_PC3 <- colnames(data_GEUVADIS_combined)[List_significant_asso 
                           [which(List_significant_asso[,1]=="3"),2]]
print(eqtl.genes_PC3)


#------------- 
# Example 1 
# Gene SBF2-AS1 is significantly associated with PC2 
print(eqtl.genes_PC2[24])

# Gene SBF2-AS1 is in the eQTL-gene set #50 with snp rs7124238 and gene SWAP70
data_GEU_Q50 <- data_GEUVADIS$Data_Q50$Data_EUR
colnames(data_GEU_Q50) <- c("rs7124238","SBF2-AS1","SWAP70")

# Analyze the eQTL-gene set without PC2
n <- nrow (data_GEU_Q50)        # Number of rows
V <- colnames(data_GEU_Q50)     # Column names

# Calculate Pearson correlation
suffStat_C_Q50 <- list(C = cor(data_GEU_Q50, use = 'pairwise.complete.obs'),
                       n = n)

# Infer the graph by MRPC 
MRPC.fit_withoutPC_GEU_Q50 <- MRPC(data_GEU_Q50, 
                                  suffStat = suffStat_C_Q50, 
                                  GV = 1, 
                                  FDR = 0.05, 
                                  indepTest = 'gaussCItest', 
                                  labels = V,  
                                  FDRcontrol = TRUE, 
                                  verbose = FALSE)

# Analyze the eQTL-gene set with PC2
data_withPC_Q50 <- cbind(data_GEU_Q50,PCs_matrix[,2])
colnames(data_withPC_Q50)[4] <- "PC2"

n <- nrow (data_withPC_Q50)        # Number of rows
V <- colnames(data_withPC_Q50)     # Column names

# Calculate Pearson correlation
suffStat_C_withPC_Q50 <- list(C = cor(data_withPC_Q50, use = 'pairwise.complete.obs'),  
                              n = n)

# Infer the graph by MRPC
MRPC.fit_withPC_GEU_Q50 <- MRPC(data_withPC_Q50, 
                                suffStat = suffStat_C_withPC_Q50, 
                                GV = 1, 
                                FDR = 0.05, 
                                indepTest = 'gaussCItest', 
                                labels = V,  
                                FDRcontrol = TRUE, 
                                verbose = FALSE)

# Plot inferred graphs 
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q50,
     main = "Without PC" )
plot(MRPC.fit_withPC_GEU_Q50,
     main = "Without PC")

#------------- 
# Example 2
# Gene LCMT2 is significantly associated with PC1 
print(eqtl.genes_PC1[8])

# Gene LCMT2 is in the eQTL-gene set #29 with snp rs2278858 and gene ADAL 
data_GEU_Q29 <- data_GEUVADIS$Data_Q29$Data_EUR
colnames(data_GEU_Q29) <- c("rs2278858", "LCMT2", "ADAL")

# Analyze the eQTL-gene set without PC1
n <- nrow (data_GEU_Q29)        # Number of rows
V <- colnames(data_GEU_Q29)     # Column names

# Calculate Pearson correlation
suffStat_C_Q29 <- list(C = cor(data_GEU_Q29, use = 'pairwise.complete.obs'),
                       n = n)

# Infer the graph by MRPC
MRPC.fit_withoutPC_GEU_Q29 <- MRPC(data_GEU_Q29, 
                                   suffStat = suffStat_C_Q29, 
                                   GV = 1, 
                                   FDR = 0.05, 
                                   indepTest = 'gaussCItest', 
                                   labels = V,  
                                   FDRcontrol = TRUE, 
                                   verbose = FALSE)
 
# Analyze the eQTL-gene set with PC1
data_withPC_Q29 <- cbind(data_GEU_Q29,PCs_matrix[,1])
colnames(data_withPC_Q29)[4] <- "PC1"

n <- nrow (data_withPC_Q29)        # Number of rows
V <- colnames(data_withPC_Q29)     # Column names

# Calculate Pearson correlation
suffStat_C_withPC_Q29 <- list(C = cor(data_withPC_Q29, use = 'pairwise.complete.obs'),
                              n = n)

# Infer graph by MRPC
MRPC.fit_withPC_GEU_Q29 <- MRPC(data_withPC_Q29, 
                                suffStat = suffStat_C_withPC_Q29, 
                                GV = 1, 
                                FDR = 0.05, 
                                indepTest = 'gaussCItest', 
                                labels = V,  
                                FDRcontrol = TRUE, 
                                verbose = FALSE)
 
# Plot inferred graphs 
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q29,
     main = "Without PC" )
plot(MRPC.fit_withPC_GEU_Q29,
     main = "With PC")

#------------- 
# Example 3
# Genes SERPINB8 and HMSD are significantly associated with PC2 
print(eqtl.genes_PC2[c(20,21)])

# Genes SERPINB8 and HMSD are in the eQTL-gene set #43 with snp rs55928920 
data_GEU_Q43 <- data_GEUVADIS$Data_Q43$Data_EUR
colnames(data_GEU_Q43) <- c("rs55928920",	"SERPINB8",	"HMSD")

# Analyze the eQTL-gene set without PC2
n <- nrow (data_GEU_Q43)        # Number of rows
V <- colnames(data_GEU_Q43)     # Column names

# Calculate Pearson correlation
suffStat_C_Q43 <- list(C = cor(data_GEU_Q43, use = 'pairwise.complete.obs'), 
                       n = n)

# Infer the graph by MRPC
MRPC.fit_withoutPC_GEU_Q43 <- MRPC(data_GEU_Q43, 
                                   suffStat = suffStat_C_Q43, 
                                   GV = 1, 
                                   FDR = 0.05, 
                                   indepTest = 'gaussCItest', 
                                   labels = V,  
                                   FDRcontrol = TRUE, 
                                   verbose = FALSE)
 
# Analyze the eQTL-gene set with PC2
data_withPC_Q43 <- cbind(data_GEU_Q43,PCs_matrix[,2])
colnames(data_withPC_Q43)[4] <- "PC2"

n <- nrow (data_withPC_Q43)        # Number of rows
V <- colnames(data_withPC_Q43)     # Column names

# Calculate Pearson correlation
suffStat_C_withPC_Q43 <- list(C = cor(data_withPC_Q43, use = 'pairwise.complete.obs'), 
                              n = n)

# Infer the graph by MRPC
MRPC.fit_withPC_GEU_Q43 <- MRPC(data_withPC_Q43, 
                                suffStat = suffStat_C_withPC_Q43, 
                                GV = 1, 
                                FDR = 0.05, 
                                indepTest = 'gaussCItest', 
                                labels = V,  
                                FDRcontrol = TRUE, 
                                verbose = FALSE)

# Plot inferred graphs
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q43,
     main = "Without PC" )
plot(MRPC.fit_withPC_GEU_Q43,
     main = "With PC")

#------------- 
# Example 4
# Gene PLAC8 is significantly associated with PC2 and PC3
print(eqtl.genes_PC2[17])
print(eqtl.genes_PC3[12])

# Gene PLAC8 is in the eQTL-gene set  #34 with snp rs28718968 and gene COQ2 
data_GEU_Q34 <- data_GEUVADIS$Data_Q34$Data_EUR
colnames(data_GEU_Q34) <- c("rs28718968",	"COQ2", "PLAC8")

# Analyze the eQTL-gene set without PC2 and PC3
n <- nrow (data_GEU_Q34)        # Number of rows
V <- colnames(data_GEU_Q34)     # Column names

# Calculate Pearson correlation
suffStat_C_Q34 <- list(C = cor(data_GEU_Q34, use = 'pairwise.complete.obs'), 
                       n = n)

# Infer the graph by MRPC
MRPC.fit_withoutPC_GEU_Q34 <- MRPC(data_GEU_Q34, 
                                   suffStat = suffStat_C_Q34, 
                                   GV = 1, 
                                   FDR = 0.05, 
                                   indepTest = 'gaussCItest', 
                                   labels = V,  
                                   FDRcontrol = TRUE, 
                                   verbose = FALSE)
 
# Analyze the eQTL-gene set with PC2 and PC3
data_withPC_Q34 <- cbind(data_GEU_Q34,PCs_matrix[,c(2,3)])
colnames(data_withPC_Q34)[4:5] <- c("PC2", "PC3")

n <- nrow (data_withPC_Q34)        # Number of rows
V <- colnames(data_withPC_Q34)     # Column names

# Calculate Pearson correlation
suffStat_C_withPC_Q34 <- list(C = cor(data_withPC_Q34, use = 'pairwise.complete.obs'), 
                              n = n)

# Infer the graph by MRPC
MRPC.fit_withPC_GEU_Q34 <- MRPC(data_withPC_Q34, 
                                suffStat = suffStat_C_withPC_Q34, 
                                GV = 1, 
                                FDR = 0.05, 
                                indepTest = 'gaussCItest', 
                                labels = V,  
                                FDRcontrol = TRUE, 
                                verbose = FALSE)

# Plot inferred graphs 
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q34,
     main = "Without PC" )
plot(MRPC.fit_withPC_GEU_Q34,
     main = "With PC")

#-------------
# Example 5
# Genes PIP4P1 and PNP are significantly associated with PC1 and PC3, respectively.
print(eqtl.genes_PC1[1])
print(eqtl.genes_PC3[7])

# Genes PIP4P1 and PNP are in the eQTL-gene set #8 with snp rs11305802 and gene  AL355075.3 
data_GEU_Q8 <- data_GEUVADIS$Data_Q8$Data_EUR
colnames(data_GEU_Q8) <- c("rs11305802","PIP4P1", "AL355075.3", "PNP")

# Analyze the eQTL-gene set without PC1 and PC3
n <- nrow (data_GEU_Q8)        # Number of rows
V <- colnames(data_GEU_Q8)     # Column names

# Calculate Pearson correlation
suffStat_C_Q8 <- list(C = cor(data_GEU_Q8, use = 'pairwise.complete.obs'), 
                      n = n)

# Infer the graph by MRPC
MRPC.fit_withoutPC_GEU_Q8 <- MRPC(data_GEU_Q8, 
                                  suffStat = suffStat_C_Q8, 
                                  GV = 1, 
                                  FDR = 0.05, 
                                  indepTest = 'gaussCItest', 
                                  labels = V,  
                                  FDRcontrol = TRUE, 
                                  verbose = FALSE)
 
# Analyze the eQTL-gene set with PC1 and PC3
data_withPC_Q8 <- cbind(data_GEU_Q8,PCs_matrix[,c(1,3)])
colnames(data_withPC_Q8)[5:6] <- c("PC1","PC3")


n <- nrow (data_withPC_Q8)        # Number of rows
V <- colnames(data_withPC_Q8)     # Column names

# Calculate Pearson correlation
suffStat_C_withPC_Q8 <- list(C = cor(data_withPC_Q8, use = 'pairwise.complete.obs'), 
                             n = n)

# Infer the graph by MRPC
MRPC.fit_withPC_GEU_Q8 <- MRPC(data_withPC_Q8, 
                               suffStat = suffStat_C_withPC_Q8, 
                               GV = 1, 
                               FDR = 0.05, 
                               indepTest = 'gaussCItest', 
                               labels = V,  
                               FDRcontrol = TRUE, 
                               verbose = FALSE)
# Plot inferred graphs 
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q8,
     main = "Without PC" )
plot(MRPC.fit_withPC_GEU_Q8,
     main = "With PC")

}
}

