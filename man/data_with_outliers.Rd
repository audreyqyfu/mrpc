\name{data_with_outliers}
\alias{data_with_outliers}

\title{
Example data with outliers
}
\description{
The data contain two genotype nodes, V1 and V2, and 3 phenotype nodes, T1, T2 and T3. The genotype nodes are discrete, whereas the phenotype nodes are continuous. The data matrix includes 10 outliers generated from a uniform distribution. The code below compares the performance of \link{MRPC}, \link{mmhc} and \link{pc} on this data set.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
library(pcalg) #for pc
library(bnlearn) #for mmhc

#Truth without outlier
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

# Classical correlation
suffStat_C1 <- list(C = cor(data_without_outliers),
                    n = n)

# Robust correlation (Beta = 0.005)
Rcor_R1 <- RobustCor(data_without_outliers,
                     0.005) 
                     
suffStat_R1 <- list(C = Rcor_R1$RR,
                    n = n)

# Estimated graph by MRPC using gaussCItest and beta = 0.005
MRPC.fit_withoutoutlier <- MRPC(data_without_outliers, suffStat_R1,
                                GV = 2, FDR = 0.05,
                                indepTest ='gaussCItest',
                                labels = V, FDRcontrol = TRUE,verbose = TRUE)

# Estimated graph by mmhc
data <- data.frame(data_without_outliers)
mmhc_withoutoutlier <- mmhc(data)

# Estimated graph by pc with alpha = 0.05
pc.fit_withoutoutlier <- pc(suffStat_C1,
                            indepTest = gaussCItest,
                            alpha = 0.05, labels = V,
                            verbose = TRUE)

# Data with outliers
n <- nrow (data_with_outliers)    # Number of rows
V <- colnames(data_with_outliers) # Column names

# Classical correlation
suffStat_C1 <- list(C = cor(data_with_outliers),
                    n = n)

# Robust correlation (Beta = 0.005)
Rcor_R1 <- RobustCor(data_with_outliers,
                     0.005) 
                     
suffStat_R1 <- list(C = Rcor_R1$RR,
                    n = n)

# Estimated graph by MRPC using gaussCItest and beta = 0.005
MRPC.fit_withoutlier <- MRPC(data_with_outliers, suffStat_R1,
                             GV = 2, FDR = 0.05,
                             indepTest ='gaussCItest',
                             labels = V, FDRcontrol = TRUE,verbose = TRUE)

# Estimated graph by mmhc
data <- data.frame(data_with_outliers)
mmhc_withoutlier <- mmhc(data)

# Estimated graph by pc with alpha = 0.05
pc.fit_withoutlier <- pc(suffStat_C1,
                         indepTest = gaussCItest,
                         alpha = 0.05, labels = V,
                         verbose = TRUE)

# Plot the results
# Show the estimated graph
par(mfrow = c(2, 4))
plot(Truth,
     main = "Truth")
plot(MRPC.fit_withoutoutlier,
     main = "MRPC")
graphviz.plot(mmhc_withoutoutlier,
              main = "mmhc")
plot(pc.fit_withoutoutlier,
     main = "pc")
plot(Truth,
     main = " ")
plot(MRPC.fit_withoutlier,
     main = " ")
graphviz.plot(mmhc_withoutlier,
              main = " ")
plot(pc.fit_withoutlier,
     main = " ")
}
