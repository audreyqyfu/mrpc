\name{data_without_outliers}
\alias{data_without_outliers}

\title{
Example data without outliers (noises)
}
\description{
The data contain two genotype nodes, V1 and V2, and three phenotype nodes, T1, T2 and T3. The code below compares the performance of \link{MRPC}, \link{mmhc} and \link{pc} on this data set.

}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
# Load packages

library(MRPC)     # MRPC
library(pcalg)    # pc
library(bnlearn)  # mmhc

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

# Calculate robust correlation (Beta = 0.005)
Rcor_R1 <- RobustCor(data_without_outliers,
                     Beta=0.005) 
suffStat_R1 <- list(C = Rcor_R1$RR,
                    n = n)

# Infer the graph by MRPC with robust correlation 
MRPC.fit_withoutoutlier <- MRPC(data_without_outliers, 
                                suffStat = suffStat_R1,
                                GV = 2, 
                                FDR = 0.05,
                                indepTest ='gaussCItest',
                                labels = V, 
                                FDRcontrol = TRUE,
                                verbose = TRUE)

# Infer the by pc with Pearson correlation
pc.fit_withoutoutlier <- pc(suffStat = suffStat_C1,
                            indepTest = gaussCItest,
                            alpha = 0.05, 
                            labels = V,
                            verbose = TRUE)

# Infer the graph by mmhc
data <- data.frame(data_without_outliers)
mmhc_withoutoutlier <- mmhc(data)

# Data with outliers
n <- nrow (data_with_outliers)    # Number of rows
V <- colnames(data_with_outliers) # Column names

# Calculate Pearson correlation
suffStat_C1 <- list(C = cor(data_with_outliers),
                    n = n)

# Calculate robust correlation (Beta = 0.005)
Rcor_R1 <- RobustCor(data_with_outliers,
                     Beta=0.005) 
suffStat_R1 <- list(C = Rcor_R1$RR,
                    n = n)

# Infer the graph by MRPC with robust correlation
MRPC.fit_withoutlier <- MRPC(data_with_outliers, 
                             suffStat = suffStat_R1,
                             GV = 2, 
                             FDR = 0.05,
                             indepTest ='gaussCItest',
                             labels = V, 
                             FDRcontrol = TRUE,
                             verbose = TRUE)

# Infer the graph by pc with Pearson correlation
pc.fit_withoutlier <- pc(suffStat = suffStat_C1,
                         indepTest = gaussCItest,
                         alpha = 0.05, 
                         labels = V,
                         verbose = TRUE)

# Infer the graph by mmhc
data <- data.frame(data_with_outliers)
mmhc_withoutlier <- mmhc(data)

# Plot the inferred graphs
par(mfrow = c(2, 4))
plot(Truth,
     main = "Truth")
plot(MRPC.fit_withoutoutlier,
     main = "MRPC")
plot(pc.fit_withoutoutlier,
     main = "pc")
graphviz.plot(mmhc_withoutoutlier,
              main = "mmhc")
plot(Truth,
     main = " ")
plot(MRPC.fit_withoutlier,
     main = " ")
plot(pc.fit_withoutlier,
     main = " ")
graphviz.plot(mmhc_withoutlier,
              main = " ")
}
