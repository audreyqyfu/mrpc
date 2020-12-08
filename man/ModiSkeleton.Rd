\name{ModiSkeleton}
\alias{ModiSkeleton}

\title{
Infer a graph skeleton (undirected graph)
}
\description{
This function infers a graph skeleton (i.e., an undirected graph). It is based on the function \link{skeleton} from the \code{pcalg} package. Both functions perform marginal and conditional indpenendence tests. However, \code{ModiSkeleton} implements an online false discovery rate (FDR) control method in order to control the overall FDR, whereas \link{skeleton} controls only the type I error rate for each individual test. See details below.
}
\usage{
ModiSkeleton(data, suffStat, FDR, alpha, indepTest, labels, p,
             method = c("stable", "original", "stable.fast"),
             m.max = Inf, fixedGaps = NULL, fixedEdges = NULL,
             NAdelete = TRUE, FDRcontrol = c("LOND", "ADDIS"),
             tau, lambda, verbose = FALSE)
}

\arguments{
Many arguments are similar to those in \link{skeleton} and \link{pc} in the \code{pcalg} package. Several arguments here are also arguments for the function \link{MRPC}.

\item{data}{
Data matrix, where the rows are samples and the columns are features (e.g., genetic variants (GVs) and phenotypes). Columns are for GVs, if available, appear before other columns for phenotypes (e.g., gene expression). For example, if there is one GV, then the first column of the data matrix is the GV and the remaining columns are the gene expression data.
}

\item{suffStat}{A list consisting of the correlation matrix of the data and the sample size.  
}

\item{FDR}{
Desired overall FDR level.
}
\item{alpha}{significance level (number in (0,1) for the individual tests.}

\item{indepTest}{Name of the statistical test. It is used to test the independence of x and y given S, where x and y are variables and S is a vector, possibly empty, of variables. The return value of indepTest is the p-value of the test for conditional independence.  Different tests may used for different data types. For example, \code{indepTest='gaussCItest'} for Gaussian data, \code{indepTest='disCItest'} for discrete data, and \code{indepTest='binCItest'} for binary data. See additional details in help(gaussCItest).

\link{ci.test} in the \code{bnlearn} package (Marco Scutari, 2010) may also be used for testing conditional independence and return a p-value. The default test statistic is the mutual information for categorical variables, the Jonckheere-Terpstra test for ordered factors and the linear correlation for continuous variables. See help(ci.test).
}

\item{labels}{A character vector of names of variables (nodes). These are typically the column names of the data matrix.}

\item{p}{(optional) The number of variables (nodes). Need to be specified if the labels are not provided, in which case the labels are set to 1:p.}

  \item{method}{
(optional) Character string specifying method. The default, "stable" provides an order-independent skeleton.
}

  \item{m.max}{
(optional) Maximum size of the conditioning sets that are considered in the conditional independence tests.
}

\item{fixedGaps}{(optional) A logical matrix of dimension p*p. If entry [x, y], [y, x], or both are TRUE, the edge x---y is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.}

\item{fixedEdges}{ (optional) A logical matrix of dimension p*p. If entry [x, y], [y, x], or both are TRUE, the edge x---y is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.}

\item{NAdelete}{(optional) If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.}

\item{FDRcontrol}{(optional) A character string specifying the method used for online FDR control. The two options are either LOND (Javanmard and Montanari, 2015) or ADDIS (Tian and Ramdas, 2019). If this is left blank a fixed significance level, alpha, will be used for each test.}

\item{tau}{A number between 0 and 1. This value is used to determine if a p-value will be considered for testing. For example, if a p-value is greater than tau then it is discarded and no test will be performed.}

\item{lambda}{A number between 0 and tau. This value is used to determine if a p-value is a candidate for rejection. For example, if a p-value is smaller than lambda then it can be rejected when testing the hypothesis (if the p-value is smaller than alphai).}

\item{verbose}{
(optional) If TRUE, detailed output is provided. Default is FALSE for no output details
}

}
\details{
The \link{ModiSkeleton} function incorporates sequential hypothesis testing to infer the graph skeleton. This function starts with a complete graph (all nodes are connected with undirected edges) and performs a series of marginal and conditional independence tests, removing the corresponding edge if the test is not rejected.

\bold{First}, all pairs of nodes are tested for marginal independence. If two nodes x and y are judged to be marginally independent at a type I error rate alpha, the edge between them is deleted and the empty set is saved as separation sets S[x, y] and S[y, x]. After all pairs have been tested for marginal independence, some edges may be removed.

\bold{Second}, nodes (x, y) with an edge are tested for conditional independence given all subsets of the neighboring nodes. If there is any node z such that x and y are conditionally independent given z, the edge between x and y is removed and node z is saved as separation set, sepset, S[x, y] and S[y, x]. The algorithm continues in this way by increasing the size of the conditioning set step by step. The algorithm stops if all adjacency sets in the current graph are smaller than the size of the conditioning set. The result is the skeleton in which every edge is still undirected.

Unlike existing algorithms, which control only the type I error rate for each individual test, MRPC implements the LOND (Level On the Number of Discoveries) method (Javanmard and Montanari, 2015), which is a sequential hypothesis testing procedure and sets value of alpha for each test based on the number of discoveries (i.e., rejections), to control the overall false discovery rate.
}

\value{
An object containing an estimate of the skeleton of the underlying DAG as follow:
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
    maximum p-value of all conditional
    independence tests for edge i--j.}
    \item{\code{graph}:}{Object of class \code{"\linkS4class{graph}"}:
      The undirected or partially directed graph that was estimated.}
          \item{\code{zMin}:}{Deprecated.}
    \item{\code{test}:}{The number of tests that have been performed.}
    \item{\code{alpha}:}{The level of significance 
    for the current test.}
    \item{\code{R}:}{All of the decisions made from tests that have been performed. A 1 indicates a rejected null hypothesis and 0 represents a null hypothesis that was not rejected.}
    \item{\code{K}:}{The total number of rejections.}
    \item{\code{pval}:}{A vector of p-values calculated for each test.}
    \item{\code{normalizer}:}{The value that ensures the vector gammai sums to one.}
    \item{\code{exponent}:}{The exponent of the p-series used to calculate each value of the gammai vector.}
    \item{\code{alphai}:}{A vector containing the alpha value calculated for each test.}
    \item{\code{kappai}:}{A vector containing the iteration at which each rejected test occurs.}
    \item{\code{kappai_star}:}{Each element of this vector is the sum of the Si vector up to the iteration at which each rejection occurs.}
    \item{\code{Ci}:}{A vector indicating whether or not a p-value is a candidate for being rejected.}
    \item{\code{Si}:}{A vector indicating whether or not a p-value was discarded.}
    \item{\code{Ci_plus}:}{Each element of this vector represents the number of times each kappai value was counted when calculating each alphai value.}
    \item{\code{gammai}:}{The elements of this vector are the values of the p-series 0.4374901658/(m^(1.6)), where m is the iteration at which each test is performed.}
    
  }
}
\references{
1. Badsha MB, Mollah MN, Jahan N and Kurata H (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Benjamini Y and Hochberg Y (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing, J. R. Statist. Soc. B, B, 57, 289-300.

3. Javanmard A and Montanari A (2015). On Online Control of False Discovery Rate. arXiv:150206197 [statME].

4. Kalisch M and Buhlmann P (2007). Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

5. Scutari M (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

6. Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.

7.Tsamardinos I, Brown LE and Aliferis CF (2006). The Max-Min Hill-Climbing Bayesian Network Structure Learning Algorithm. Machine Learning, 65 (1), 31-78.

8. Tian J and Ramdas A (2019). ADDIS: an adaptive discarding algorithm for online FDR control with conservative nulls. In Advances in Neural Information Processing Systems (pp. 9388-9396).

}
\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\seealso{
 \link{MRPC}; \link{EdgeOrientation}; \link{SimulateData}.
}
\examples{

# Load predefined simulated data
# Data pre-processing

# The 1st column of the input matrix will be the
# genotype of the expression quantitative trait loci
# (eQTL)/Copy number variation (CNVs) and the remaining
# columns are the gene expression data.
# We used pre-assigned level alpha = 0.05 that ensures
# FDR and mFDR remains below 0.05.

# Model 1

data <- simu_data_M1 # load data for model 1
n <- nrow(data)      # Number of row
V <- colnames(data)  # Column names

# Calculate Pearson correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Infer a graph skeleton
Skel.fit <- ModiSkeleton(data, 
                         suffStat = suffStat_C,
                         FDR = 0.05, 
                         alpha = 0.05, 
                         indepTest = 'gaussCItest',
                         labels = V, 
                         FDRcontrol = 'LOND', 
                         verbose = FALSE)

# Plot the results

plot(Skel.fit@graph,
     main ="Estimated Skeleton")

# Other models are available and may be called as follows:
# Model 0
# data <- simu_data_M0

# Model 2
# data <- simu_data_M2

# Model 3
# data <- simu_data_M3

# Model 4
# data <- simu_data_M4

# Model Multiparent
# data <- simu_data_multiparent

# Model Star
# data <- simu_data_starshaped

# Model Layered
# data <- simu_data_layered

}

