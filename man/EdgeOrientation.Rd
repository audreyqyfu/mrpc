\name{EdgeOrientation}
\alias{EdgeOrientation}

\title{
Edge orientation rules for the MRPC algorithm
}
\description{
This function performs the last step of the \link{MRPC} algorithm where it determines the edge direction of the undirected graph. The function first determines the edges between genetic variants and gene expression nodes based on MR. Then it orients the v-structures followed by the remaining edges weather MR is applicable or not. MR is a new way for edge direction determination based on four different cases. See below for the details.
}
\usage{
EdgeOrientation(gInput, GV = GV, suffStat, FDR, indepTest = indepTest, verbose = FALSE)

}

\arguments{
\item{gInput}{Object containing skeleton, marginal and conditional independence information.
}
\item{GV}{The number of genetic variants (SNPs/indels/CNV/eQTL) in the input data matrix. For example, if the data has one genetic variant, first column, then GV = 1, if 2, 1st and 2nd Column, then GV = 2, and so on.
}
\item{suffStat}{A \link{list} of sufficient statistics containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistics consist of the correlation matrix of the data and the sample size.
}

\item{FDR}{
Need to specify pre-assigned level. If FDR = 0.05, that ensures FDR and mFDR remains below 0.05.
}
\item{indepTest}{A \link{function} for testing conditional independence. It is used to test the conditional independence of x and y given S, called as indepTest(x, y, S, suffStat). Where, x and y are variables, and S is a vector, possibly empty, of variables. suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. The different indepTest is used for different data types, for example, Gaussian data = gaussCItest, Discrete data = disCItest and Binary data = binCItest. See help(gaussCItest)}
\item{verbose}{(optional) 1: detailed output is provided; 0: No output is provided}
}
\details{
The orientation of the edge directions based on Mendelian randomization using the four different cases. Here, we consider x is a genetic variant, y and z are the gene expression data. The 1st column of the input matrix will be the genetic variant and the remaining columns are the gene expression data.

Four different cases are as follows:

\bold{Case-1}: Relation between x, genetic variant, and the other nodes. Then genetic variant will regulate the other node, genes, and direction will be genetic variant --> other node. Note that if the data has more than one genetic variant and two genetic variant have edges, then direction will be genetic variant <--> genetic variant, which indicates that there is evidence that the two genetic variants are not independent, but we do not have enough information to determine which genetic variant is the regulator and which is the target.

\bold{Case-2}: If y and z are adjacent and, x and z are conditionally independent given y, then gene y will regulate the expression of gene z and the edge direction will be y --> z.


\bold{Case-3}: If y and z are adjacent and, x and z are conditionally dependent given y, then gene z will regulate the expression of gene y and the edge direction will be z --> y.

\bold{Case-4}: If y and z are adjacent with x and y conditionally dependent given z and x and z conditionally dependent given y, then the edge direction will be y <--> z.

}

\value{
An object of \link{class} that contains an estimate of the equivalence class of the underlying DAG.
\describe{
    
    \item{\code{call}:}{A \link{call} object: the original function call.}
    
    \item{\code{n}:}{The sample 
    size used to estimate the graph.}
        \item{\code{max.ord}:}{The 
    maximum size of the conditioning set used 
    in the conditional independence tests of   
    the first part of the algorithm.}
    \item{\code{n.edgetests}:}{The number of 
    conditional independence tests performed by
    the first part of the algorithm.}
    \item{\code{sepset}:}{Separation sets.}
    \item{\code{pMax}:}{A square matrix
    , where the (i, j)th entry contains the 
    maximal p-value of all conditional
    independence tests for edge i--j.}
    \item{\code{graph}:}{An object of class \code{"\linkS4class{graph}"}:
      The undirected or partially directed graph that was estimated.}
          \item{\code{zMin}:}{Deprecated.}
    \item{\code{test}:}{The number of tests that have been performed.}
    \item{\code{alpha}:}{The level of significance 
    for the current test.}
    \item{\code{R}:}{A vector of all the decisions made so far from the tests that have been performed.}
    
  }
}


\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based (PC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{SimulateData} for simulated data generating function.

}

\examples{
# Load packages
library(pcalg)  #library for existing pc

# Load predefined data
# Data pre-processing

# The 1st column of the input matrix will be the
# genetic variant and the remaining columns are the gene expression data.

# Model 1
Truth <- MRPCtruth$M1   #Truth for model 1
data <- simu_data_M1    #data load for model 1
n <- nrow (data)        #Number of row
V <- colnames(data)     #Column names

Rcor_R <- RobustCor(data,
                    0.005) #Robust correlation (Beta = 0.005)
                    
suffStat_R <- list(C = Rcor_R$RR,
                   n = n)

# Estimate skeleton
Skel.fit <- ModiSkeleton(data, suffStat_R, FDR = 0.05,
                         indepTest = 'gaussCItest',
                         labels = V, verbose = TRUE)

# Edge Orientation
Edge_orientation <- EdgeOrientation(Skel.fit, GV = 1,
                                    suffStat_R, FDR = 0.05,
                                    indepTest = 'gaussCItest', verbose = 1)

# Plot the results
par(mfrow = c(1, 2))
plot(Truth,
     main = "(A) Truth")
plot(Edge_orientation,
     main = "(B) MRPC ")

# Other models are available and may be called as follows:
# Model 0
# Truth <- MRPCtruth$M0
# data <- simu.data_M0

# Model 2
# Truth <- MRPCtruth$M2
# data <- simu_data_M2

# Model 3
# Truth <- MRPCtruth$M3
# data <- simu_data_M3

# Model 4
# Truth <- MRPCtruth$M4
# data <- simu_data_M4

# Model Multiparent
# Truth <- MRPCtruth$Multiparent
# data <- simu_data_multiparent

# Model Star
# Truth <- MRPCtruth$Star
# data <- simu_data_starshaped

# Model Layered
# Truth <- MRPCtruth$Layered
# data <- simu_data_layered

}

