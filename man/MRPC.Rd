\name{MRPC}
\alias{MRPC}
\title{
Estimate a causal network using the MRPC algorithm
}
\description{
The MRPC () function is used to estimate a causal network (or causal graph) with directed and undirected edges from observational data, by modifying the existing \link{pc} algorithm (Spirtes et al., 2000). Existing \link{pc} algorithm is used for estimating the DAG from observational data and several variants are implemented in the R package called pcalg. Several challenges remain when the \link{pc} algorithm is applied to noisy genomic data. To tackle these challenges, we developed a novel machine learning algorithm called Mendelian Randomization (MR) based PC (MRPC) to efficiently learn a causal gene regulatory network (or a causal graph of genes) using genotype and gene expression data. Unlike existing \link{pc} algorithm, the four major updates in the MRPC algorithm are that (i) we incorporate a sequential testing method to control the False Discovery Rate (FDR), (ii) take a robust approach to reduce the impact of outliers, (iii) improve v-structure identification and (iv) implement a new way for edge direction determination based on the Principle of Mendelian Randomization (PMR) (Badsha and Fu, 2018 and Badsha et al., 2018). See details below. 

}
\usage{
MRPC(data, suffStat, GV, FDR = 0.05, indepTest, labels, p,
    fixedGaps = NULL, fixedEdges = NULL,
    NAdelete = TRUE, m.max = Inf,
    u2pd = c("relaxed", "rand", "retry"),
    skel.method = c("stable", "original", "stable.fast"),
    conservative = FALSE,
    maj.rule = FALSE, solve.confl = FALSE,
    verbose = FALSE)
}

\arguments{
The MRPC algorithm is a modification of the existing \link{pc} algorithm. Therefore, most of the terms and sentences are taken from \link{pc} (Spirtes et al., 2000).

\item{data}{
Data matrix, where the rows are the observations and the columns are the genes. The columns are ordered from  single-nucleotide polymorphism (SNPs), indels, copy number variation (CNV), or expression quantitative trait loci (eQTL) to genes. If we consider one genetic variant, then the first column of the input matrix is the genetic variant and the remaining columns are the gene expression data and so on.
}

\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistics consist of the correlation matrix of the data and the sample size.
}

\item{GV}{The number of genetic variants (SNPs/indels/CNV/eQTL) in the input data matrix. For example, if the data has one SNPs/indels/CNV/eQTL, first column, then GV = 1, if 2 SNPs/indels/CNV/eQTL, 1st and 2nd Column, then GV = 2, if no GV then GV = 0, and so on.
}

\item{FDR}{
Need to specify pre-assigned level. If FDR = 0.05 this ensures FDR and mFDR remains below 0.05.
}
\item{indepTest}{A \link{function} for testing conditional independence. It is used to test the conditional independence of x and y given S, called as indepTest(x, y, S, suffStat). Where, x and y are variables, and S is a vector, possibly empty, of variables. suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. The different indepTest is used for different data types, for example, Gaussian data = gaussCItest, Discrete data = disCItest and Binary data = binCItest. See help(gaussCItest)

The \link{ci.test} (Marco Scutari, 2010) is also used for testing conditional independence and return value of indepTest is the p-value. If none is specified, the default test statistic is the mutual information for categorical variables, the Jonckheere-Terpstra test for ordered factors and the linear correlation for continuous variables.See help(ci.test)

Remember that need to specify the which indepTest would like for independence testing. For example, if you would like to use gaussCItest you would type indepTest = 'gaussCItest' into the function otherwise indepTest = 'citest'. Note that, we used \link{gaussCItest} to compare our \link{MRPC} with the existing \link{pc}, because of \link{ci.test} is not robust. See details in example.
}

\item{labels}{A character vector of variable, or node, names. All variables are denoted in column in the input matrix.}

\item{p}{(optional) The number of variables, or nodes. May be specified if the labels are not provided, in which case the labels are set to 1:p.}

\item{fixedGaps}{(optional) A logical matrix of dimension p*p. If entry [x, y], [y, x], or both are TRUE, the edge x---y is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.}

\item{fixedEdges}{ (optional) A logical matrix of dimension p*p. If entry [x, y], [y, x], or both are TRUE, the edge x---y is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.}

\item{NAdelete}{(optional) If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.}

\item{m.max}{ (optional) Maximum size of the conditioning sets that are considered in the conditional independence tests.}

\item{u2pd}{ (optional) Character string specifying the method for dealing with conflicting information.}

\item{skel.method}{ (optional) Character string specifying method; the default, "stable" provides an order-independent skeleton.}

\item{conservative}{(optional) Logical. Indicates if the conservative PC algorithm is used. In this case, only option u2pd = "relaxed" is supported. Note that the resulting object might not be extendable to a DAG. See details for more information.}

\item{maj.rule}{(optional) Logical. Indicates that the triplets will be checked for ambiguity using a majority rule idea, which is less strict than the conservative PC algorithm. For more information, see details.}

\item{solve.confl}{(optional) If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. In this case, only option u2pd = "relaxed"" is supported. Note that the resulting object might not be a CPDAG because bi-directed edges might be present. See details for more information.}

\item{verbose}{(optional) If TRUE, detailed output is provided. The default is FALSE which does not print output details.}

}
\details{
Nodes are used as major building blocks of the models, which represent random variables and edges, which encode conditional dependence relations of the enclosing vertices (Kalisch and Buhlmann, 2007). The structure of conditional independence among the random variables can be explored using the Markov properties. The directed edges show the presence and direction of direct causal effects. A bidirected edge means that the edge orientation should be forward (--->) or bacward (<---). The PC algorithm is computationally efficient for learning the underlying DAG (Kalisch and Buhlmann, 2007). The MPRC algorithm improves on the existing PC algorithm in the following aspects:

(i) the genotype data at genetic variants (e.g., SNPs and copy number variation) provide additional information that helps distinguish the casual direction between two genes; this is the rationale behind Mendelian Randomization (MR). MR can greatly reduce the space of possible graphs and increase the inference efficiency for inference of genomic data but PC does not rely on MR and determination of the edge direction can be tricky when the graph is large.

(ii) the number of statistical tests is unknown beforehand in all existing variants of the PC algorithm. It is unclear how these algorithms control the false discovery rate (FDR) or marginal FDR (mFDR), because it is used alpha as a fixed pre-assigned significance level for testing the independence test, but results are sensitive to the choice of alpha, which is an important issue for multiple hypotheses testing.

(iii) the gene expression data often have outliers, even after normalization, which may drastically alter the topology of the inferred graph. 

(iv) to find the v-structure the PC algorithm assumes that if there is no evidence that two nodes are conditionally independent, then the nodes are conditionally dependent, although we don't have any evidence about their conditional dependence. In that case the v-structure is worng.

To overcome the above challenges, we propose a novel algorithm, MRPC (MR-based PC), which can be applied to genotype and gene expression data and efficiently learn a causal graph of genes. \bold{First}, we take a robust approach (Badsha et al., 2013) and calculate the robust correlation matrix on which the series of hypothesis testing is performed. \bold{Second}, we adopt a sequential method LOND (significance Levels based On Number of Discoveries) (Javanmard and Montanari, 2015) that controls the FDR or marginal FDR (mFDR) in an online manner, where level of significance is used a function of the previous decision made so far (we adjusted at each step in the test) and that ensures that FDR and mFDR remains below alpha. \bold{Third}, we improve the v-structure identification using the addtional conditional test. \bold{Finally}, we implement a new way for edge direction determination based on Mendelian randomization.

The two main steps in the MRPC algorithm are as follows:

\bold{Step-1}: We incorporated sequential hypothesis testing to draw the undirected graph (skeleton) by \link{ModiSkeleton} function  (similar as the pc algorithm by \link{skeleton} function but, difference is that the \link{ModiSkeleton} function is used alpha (significance level) as a function of the prior outcomes (reject or accept hypothesis) but the \link{skeleton} function used alpha as a fixed pre-assigned significance level and \link{ModiSkeleton} take a robust correlation approach reduces the impact of outliers but the \link{skeleton} used classical correlation which is highly influnced by outliers. See details in in \link{ModiSkeleton} and \link{skeleton}).

\bold{Step-2}: We implemented a new way for edge direction determination based on Mendelian randomization (MR). We consider the first column of the input matrix will be the genotype of the SNPs/indels/CNV/eQTL and the remaining column are the gene expression data. See details in \link{EdgeOrientation}.

All statistical inference is done in \bold{Step 1}, while \bold{Step-2} is just application of deterministic rules on the results of \bold{Step-1}. If the first part is done correctly, the second part will never fail. If, however, there are errors in \bold{Step-1}, \bold{Step-2} will be more sensitive since it depends on the inferential results of \bold{Step-1}.
}

\value{
An object of \link{class} that contains an estimate of the equivalence class of the underlying DAG.
\describe{
    
    \item{\code{call}:}{a \link{call} object: the original function call.}
    
    \item{\code{n}:}{The sample 
    size used to estimate the graph.}
        \item{\code{max.ord}:}{The 
    maximum size of the conditioning set used 
    in the conditional independence tests in   
    the first part of the algorithm.}
    \item{\code{n.edgetests}:}{The number of 
    conditional independence tests performed by
    the first part of the algorithm.}
    \item{\code{sepset}:}{Separation sets.}
    \item{\code{pMax}:}{A numeric square matrix
    , where the (i, j)th entry contains the 
    maximal p-value of all conditional
    independence tests for edge i--j.}
    \item{\code{graph}:}{Object of class \code{"\linkS4class{graph}"}:
      the undirected or partially directed graph that was estimated.}
          \item{\code{zMin}:}{Deprecated.}
    \item{\code{test}:}{The number of tests that have been performed.}
    \item{\code{alpha}:}{The level of significance 
    for the current test.}
    \item{\code{R}:}{All of the decisions made so far from tests that have been performed.}
    
  }
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\references{
1.Badsha, M. B. Mollah, M. N. Jahan, N. and Kurata, H. (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Javanmard and Montanari (March 5, 2015) On Online Control of False Discovery Rate. arXiv:150206197 [statME].

3.Kalisch, M., Machler, M., Colombo, D., Maathuis, M.H. & Buhlmann, P. Causal Inference Using Graphical Models with the R Package pcalg. J. Stat. Softw. 47, 26 (2012).

4. Kalisch, M. and Buhlmann, P. (2007) Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

5.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

6. Md Bahadur Badsha and Audrey Qiuyan Fu. Learning causal biological networks with the 
principle of Mendelian randomization. bioRxiv, 2018. doi:10.1101/171348.

7. Md Bahadur Badsha, Evan A Martin, and Audrey Qiuyan Fu. MRPC: An R package for accurate inference of causal graphs. arXiv, 2018. arXiv:1806.01899

8. Spirtes,P.  Glymour, C. and Scheines, R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.

}
\seealso{
\link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{EdgeOrientation} for orientation rules determination for edges in MRPC; \link{SimulateData} for simulated data generating function.
}
\examples{
# Load packages
# We compare different simulated data across three methods: MRPC, 
# PC in pcalg (Kalisch et al., 2012), and mmhc in bnlearn (Marco Scutari, 2010)
library(pcalg)    # library for pc algorithm
library(bnlearn)  # library for mmhc

# Load predefined simulated data
# Data pre-processing

# The 1st column of the input matrix will be the genotype of the
# expression quantitative trait loci (eQTL)/Copy number variation (CNV)
# and the remaining columns are the gene expression data.

# We used pre-assigned level alpha = 0.05 that ensures FDR and mFDR
# will remain below 0.05.

# Model 1
Truth <- MRPCtruth$M1   # Truth for model 1
data <- simu_data_M1    # load data for model 1
n <- nrow (data)        # Number of rows
V <- colnames(data)     # Column names

# Classical correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Robust correlation (Beta = 0.005)
Rcor_R <- RobustCor(data, 0.005) 
suffStat_R <- list(C = Rcor_R$RR,
                   n = n)

# Estimated graph by MRPC using gaussCItest and beta = 0.005
MRPC.fit<- MRPC(data,
                suffStat_R,
                GV = 1,
                FDR = 0.05,
                indepTest = 'gaussCItest',
                labels = V,
                verbose = TRUE)

# Estimated graph by mmhc
mmhc.fit <- mmhc(data)

# Estimated graph by PC with alpha = 0.05
pc.fit <- pc(suffStat_C,
             indepTest = gaussCItest,
             alpha = 0.05,
             labels = V,
             verbose = TRUE)

# Plot the results
# Show estimated graph
par(mfrow = c(2, 2))
plot(Truth,
     main = "(A) Truth")
plot(MRPC.fit,
     main = "(B) MRPC")
graphviz.plot(mmhc.fit,
              main = "(C) mmhc")
plot(pc.fit,
     main ="(D) pc")

# Another option for plot of the results. First fig is the nodes
# dendrogram with colored modules. Second fig is the plot of the graph
# with color based on modules.
# To idendify modules and complex graph (Suitable if you have many nodes)
# Adjacency matrix from directed graph
Adj_directed <- as(MRPC.fit@graph,
                   "matrix")

# Plot of the graph with colored modules.
DendroModuleGraph(Adj_directed,
                  minModuleSize = 1,
                  GV = 1)

# Other models are available and may be called as follows:
# Model 0
# Truth <- MRPCtruth$M0
# data <- simu_data_M0

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
