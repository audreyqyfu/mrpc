\name{data_without_outliers}
\alias{data_without_outliers}

\title{
Example data without outliers (noises)
}
\description{
The data contain two genotype nodes, V1 and V2, and three phenotype nodes, T1, T2 and T3. The code below compares the performance of \link{MRPC}, \link{pc} \link{pc.stable}, \link{mmpc}, \link{mmhc}, and \link{hc} on this data set.

}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\examples{
# Load packages

library(MRPC)     # MRPC
library(pcalg)    # pc
library(bnlearn)  # pc.stable, mmpc, mmhc, and hc

# Truth without outlier
tarmat <- matrix(0,
                 nrow = ncol(data_with_outliers),
                 ncol = ncol(data_with_outliers))
                
colnames(tarmat) <- colnames(data_with_outliers)
rownames(tarmat) <- colnames(data_with_outliers)

tarmat[1,2] <- 1
tarmat[2,1] <- 1
tarmat[1,3] <- 1
tarmat[4,3] <- 1
tarmat[4,5] <- 1

Truth <- as(tarmat,
            "graphNEL")

# Data without outliers
n <- nrow(data_without_outliers)     # Number of rows
V <- colnames(data_without_outliers) # Column names

# Calculate Pearson correlation
suffStat_C1 <- list(C = cor(data_without_outliers),
                    n = n)

# Infer the graph by MRPC
MRPC.fit_withoutoutliers <- MRPC (data_without_outliers, 
                                  suffStat = suffStat_C1, 
                                  GV = 2, 
                                  FDR = 0.05, 
                                  indepTest ='gaussCItest', 
                                  labels = V, 
                                  FDRcontrol = 'LOND', 
                                  verbose = FALSE)

# Infer the graph by pc with Pearson correlation
pc.fit_withoutoutliers <- pc(suffStat = suffStat_C1,
                            indepTest = gaussCItest,
                            alpha = 0.05, 
                            labels = V,
                            verbose = FALSE)

# arcs not to be included from gene expression to genotype for blacklist argument 
# in pc.stable and mmpc

GV <- 2
to <- rep (colnames (data_without_outliers)[1:GV], each = (ncol (data_without_outliers) - GV))
from <- rep (colnames (data_without_outliers)[(GV + 1):ncol (data_without_outliers)], GV)
bl <- cbind (from, to)

# Infer the graph by pc.stable
pc.stable_withoutoutliers <- pc.stable (data.frame (data_without_outliers), 
                                        blacklist = bl, 
                                        alpha = 0.05, 
                                        debug = FALSE, 
                                        undirected = FALSE)
# Infer the graph by mmpc
mmpc_withoutoutliers <- mmpc (data.frame (data_without_outliers), 
                              blacklist = bl, 
                              alpha = 0.05, 
                              debug = FALSE, 
                              undirected = FALSE)
# Infer the graph by mmhc
mmhc_withoutoutliers <- mmhc (data.frame (data_without_outliers), 
                              blacklist = bl, 
                              debug = FALSE)
# Infer the graph by hc
hc_withoutoutliers <- hc (data.frame (data_without_outliers), 
                          blacklist = bl, 
                          debug = FALSE)
# True graph
plot (Truth, main = "Truth")

#-------------
# Plot inferred graphs
par (mfrow = c (2,3))

# Data without outliers
# Inference with Pearson correlation
plot (MRPC.fit_withoutoutliers, main = "MRPC")
plot (pc.fit_withoutoutliers, main = "pc")
graphviz.plot (pc.stable_withoutoutliers, main = "pc.stable")
graphviz.plot (mmpc_withoutoutliers, main = "mmpc")
graphviz.plot (mmhc_withoutoutliers, main = "mmhc")
graphviz.plot (hc_withoutoutliers, main = "hc")
}
