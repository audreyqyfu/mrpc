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


#Example list of eQTLs or genes names with statistically associated with PC1
eqtl.genes_PC1 <- colnames(data_GEUVADIS_combined)[List_significant_asso
                                                   [which(List_significant_asso[,1]=="1"),2]]
print(eqtl.genes_PC1)
#Example list of eQTLs or genes names with statistically associated with PC2
eqtl.genes_PC2 <- colnames(data_GEUVADIS_combined)[List_significant_asso
                                                   [which(List_significant_asso[,1]=="2"),2]]
print(eqtl.genes_PC2)
#Example list of eQTLs or genes names with statistically associated with PC3
eqtl.genes_PC3 <- colnames(data_GEUVADIS_combined)[List_significant_asso
                                                   [which(List_significant_asso[,1]=="3"),2]]
print(eqtl.genes_PC3)

#For example: 
#gene SBF2-AS1 statistically associated with PC2 
print(eqtl.genes_PC2[24])

#Note, gene SBF2-AS1 is a member of eQTL 50 with snp rs7124238 and gene SWAP70 in
#GEUVADIS. Then we applied MRPC jointly for the eQTL-gene set and the associated PC2 
#(rs7124238, SWAP70, SBF2-AS1, PC2)

#eQTL 50 data without PC2 in GEUVADIS
data_GEU_Q50 <- data_GEUVADIS$Data_Q50$Data_EUR
colnames(data_GEU_Q50) <- c("rs7124238","SBF2-AS1","SWAP70")

n <- nrow (data_GEU_Q50)        # Number of rows
V <- colnames(data_GEU_Q50)     # Column names

# Classical correlation
suffStat_C_Q50 <- list(C = cor(data_GEU_Q50, use = 'pairwise.complete.obs'),
                            n = n)

# Estimated graph by MRPC 
MRPC.fit_withoutPC_GEU_Q50 <- MRPC(D1_GEU_Q50,
                           suffStat_C_Q50,
                           GV = 1,
                           FDR = 0.05,
                           alpha = 0.05,
                           FDRcontrol = TRUE,
                           indepTest = 'gaussCItest',
                           labels = V,
                           verbose = TRUE)

#eQTL 50 data with PC2 in GEUVADIS
data_withPC_Q50 <- cbind(data_GEU_Q50,PCs_matrix[,2])
colnames(data_withPC_Q50)[4] <- "PC2"

n <- nrow (data_withPC_Q50)        # Number of rows
V <- colnames(data_withPC_Q50)     # Column names

# Classical correlation
suffStat_C_withPC_Q50 <- list(C = cor(data_withPC_Q50, use = 'pairwise.complete.obs'),
                                 n = n)

# Estimated graph by MRPC 
MRPC.fit_witPC_GEU_Q50 <- MRPC(Data_withPC_Q50,
                                suffStat_C_withPC_Q50,
                                GV = 1,
                                FDR = 0.05,
                                alpha=0.05,
                                FDRcontrol = TRUE,
                                indepTest = 'gaussCItest',
                                labels = V,
                                verbose = TRUE)


#Plot
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q50,
     main = "Without PC" )
plot(MRPC.fit_witPC_GEU_Q50,
     main = "Without PC")

#gene LCMT2 statistically associated with PC1 
print(eqtl.genes_PC1[8])

#Note, gene LCMT2 is a member of eQTL 29 with snp rs2278858 and gene ADAL in
#GEUVADIS. Then we applied MRPC jointly for the eQTL-gene set and the associated PC1 
#(rs2278858, LCMT2, ADAL, PC1)

#eQTL 29 data without PC1 in GEUVADIS
data_GEU_Q29 <- data_GEUVADIS$Data_Q29$Data_EUR
colnames(data_GEU_Q29) <- c("rs2278858", "LCMT2", "ADAL")

n <- nrow (data_GEU_Q29)        # Number of rows
V <- colnames(data_GEU_Q29)     # Column names

# Classical correlation
suffStat_C_Q29 <- list(C = cor(data_GEU_Q29, use = 'pairwise.complete.obs'),
                       n = n)

# Estimated graph by MRPC 
MRPC.fit_withoutPC_GEU_Q29 <- MRPC(D1_GEU_Q29,
                                   suffStat_C_Q29,
                                   GV = 1,
                                   FDR = 0.05,
                                   alpha = 0.05,
                                   FDRcontrol = TRUE,
                                   indepTest = 'gaussCItest',
                                   labels = V,
                                   verbose = TRUE)

#eQTL 29 data with PC1 in GEUVADIS
data_withPC_Q29 <- cbind(data_GEU_Q29,PCs_matrix[,1])
colnames(data_withPC_Q29)[4] <- "PC1"

n <- nrow (data_withPC_Q29)        # Number of rows
V <- colnames(data_withPC_Q29)     # Column names

# Classical correlation
suffStat_C_withPC_Q29 <- list(C = cor(data_withPC_Q29, use = 'pairwise.complete.obs'),
                              n = n)

# Estimated graph by MRPC
MRPC.fit_witPC_GEU_Q29 <- MRPC(Data_withPC_Q29,
                               suffStat_C_withPC_Q29,
                               GV = 1,
                               FDR = 0.05,
                               alpha=0.05,
                               FDRcontrol = TRUE,
                               indepTest = 'gaussCItest',
                               labels = V,
                               verbose = TRUE)

#Plot
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q29,
     main = "Without PC" )
plot(MRPC.fit_witPC_GEU_Q29,
     main = "With PC")

#genes SERPINB8 and HMSD are statistically associated with PC2 
print(eqtl.genes_PC2[c(20,21)])

#Note, genes SERPINB8 and HMSD are a member of eQTL 43 with snp rs55928920 in GEUVADIS.
#Then we applied MRPC jointly for the eQTL-gene set and the associated PC2 
#(rs55928920, SERPINB8, HMSD, PC2)

#eQTL 43 data without PC2 in GEUVADIS
data_GEU_Q43 <- data_GEUVADIS$Data_Q43$Data_EUR
colnames(data_GEU_Q43) <- c("rs55928920",	"SERPINB8",	"HMSD")

n <- nrow (data_GEU_Q43)        # Number of rows
V <- colnames(data_GEU_Q43)     # Column names

# Classical correlation
suffStat_C_Q43 <- list(C = cor(data_GEU_Q43, use = 'pairwise.complete.obs'),
                       n = n)

# Estimated graph by MRPC 
MRPC.fit_withoutPC_GEU_Q43 <- MRPC(D1_GEU_Q43,
                                   suffStat_C_Q43,
                                   GV = 1,
                                   FDR = 0.05,
                                   alpha = 0.05,
                                   FDRcontrol = TRUE,
                                   indepTest = 'gaussCItest',
                                   labels = V,
                                   verbose = TRUE)

#eQTL 43 data with PC2 in GEUVADIS
data_withPC_Q43 <- cbind(data_GEU_Q43,PCs_matrix[,2])
colnames(data_withPC_Q43)[4] <- "PC2"

n <- nrow (data_withPC_Q43)        # Number of rows
V <- colnames(data_withPC_Q43)     # Column names

# Classical correlation
suffStat_C_withPC_Q43 <- list(C = cor(data_withPC_Q43, use = 'pairwise.complete.obs'),
                              n = n)

# Estimated graph by MRPC 
MRPC.fit_witPC_GEU_Q43 <- MRPC(Data_withPC_Q43,
                               suffStat_C_withPC_Q43,
                               GV = 1,
                               FDR = 0.05,
                               alpha=0.05,
                               FDRcontrol = TRUE,
                               indepTest = 'gaussCItest',
                               labels = V,
                               verbose = TRUE)


#Plot
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q43,
     main = "Without PC" )
plot(MRPC.fit_witPC_GEU_Q43,
     main = "With PC")

#gene PLAC8 is statistically associated with PC3
print(eqtl.genes_PC3[12])
#Note, gene PLAC8 is a member of eQTL 34 with snp rs28718968 and gene COQ2 in
#GEUVADIS.Then we applied MRPC jointly for the eQTL-gene set and the associated PC3
#(rs28718968, COQ2, PLAC8, PC3)

#eQTL 34 data without PC3 in GEUVADIS
data_GEU_Q34 <- data_GEUVADIS$Data_Q34$Data_EUR
colnames(data_GEU_Q34) <- c("rs28718968",	"COQ2",	"PLAC8")

n <- nrow (data_GEU_Q34)        # Number of rows
V <- colnames(data_GEU_Q34)     # Column names

# Classical correlation
suffStat_C_Q34 <- list(C = cor(data_GEU_Q34, use = 'pairwise.complete.obs'),
                       n = n)

# Estimated graph by MRPC 
MRPC.fit_withoutPC_GEU_Q34 <- MRPC(D1_GEU_Q34,
                                   suffStat_C_Q34,
                                   GV = 1,
                                   FDR = 0.05,
                                   alpha = 0.05,
                                   FDRcontrol = TRUE,
                                   indepTest = 'gaussCItest',
                                   labels = V,
                                   verbose = TRUE)

#eQTL 34 data with PC3 in GEUVADIS
data_withPC_Q34 <- cbind(data_GEU_Q34,PCs_matrix[,3])
colnames(data_withPC_Q34)[4] <- "PC3"

n <- nrow (data_withPC_Q34)        # Number of rows
V <- colnames(data_withPC_Q34)     # Column names

# Classical correlation
suffStat_C_withPC_Q34 <- list(C = cor(data_withPC_Q34, use = 'pairwise.complete.obs'),
                              n = n)

# Estimated graph by MRPC using gaussCItest and beta = 0.005
MRPC.fit_witPC_GEU_Q34 <- MRPC(Data_withPC_Q34,
                               suffStat_C_withPC_Q34,
                               GV = 1,
                               FDR = 0.05,
                               alpha=0.05,
                               FDRcontrol = TRUE,
                               indepTest = 'gaussCItest',
                               labels = V,
                               verbose = TRUE)

par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q34,
     main = "Without PC" )
plot(MRPC.fit_witPC_GEU_Q34,
     main = "With PC")

#gene PNP is statistically associated with PC3
print(eqtl.genes_PC3[7])
#Note, gene PNP is a member of eQTL 8 with snp rs11305802 and genes PIP4P1 and AL355075.3 in
#GEUVADIS. Then we applied MRPC jointly for the eQTL-gene set and the associated PC3 
#(rs11305802,	PIP4P1, AL355075.3,PNPPC3)

#eQTL 8 data without PC3 in GEUVADIS
data_GEU_Q8 <- data_GEUVADIS$Data_Q8$Data_EUR
colnames(data_GEU_Q8) <- c("rs11305802","PIP4P1", "AL355075.3", "PNP")

n <- nrow (data_GEU_Q8)        # Number of rows
V <- colnames(data_GEU_Q8)     # Column names

# Classical correlation
suffStat_C_Q8 <- list(C = cor(data_GEU_Q8, use = 'pairwise.complete.obs'),
                      n = n)

# Estimated graph by MRPC 
MRPC.fit_withoutPC_GEU_Q8 <- MRPC(D1_GEU_Q8,
                                  suffStat_C_Q8,
                                  GV = 1,
                                  FDR = 0.05,
                                  alpha = 0.05,
                                  FDRcontrol = TRUE,
                                  indepTest = 'gaussCItest',
                                  labels = V,
                                  verbose = TRUE)

#eQTL 8 data with PC3 in GEUVADIS
data_withPC_Q8 <- cbind(data_GEU_Q8,PCs_matrix[,3])
colnames(data_withPC_Q8)[5] <- "PC3"

n <- nrow (data_withPC_Q8)        # Number of rows
V <- colnames(data_withPC_Q8)     # Column names

# Classical correlation
suffStat_C_withPC_Q8 <- list(C = cor(data_withPC_Q8, use = 'pairwise.complete.obs'),
                             n = n)

# Estimated graph by MRPC using gaussCItest and beta = 0.005
MRPC.fit_witPC_GEU_Q8 <- MRPC(Data_withPC_Q8,
                              suffStat_C_withPC_Q8,
                              GV = 1,
                              FDR = 0.05,
                              alpha=0.05,
                              FDRcontrol = TRUE,
                              indepTest = 'gaussCItest',
                              labels = V,
                              verbose = TRUE)
#Plot
par(mfrow=c(1,2))
plot(MRPC.fit_withoutPC_GEU_Q8,
     main = "Without PC" )
plot(MRPC.fit_witPC_GEU_Q8,
     main = "With PC")

}

