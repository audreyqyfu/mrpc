\name{MRPC}
\alias{MRPC}
\title{
Infer a causal network using the MRPC algorithm
}
\description{
This function is used to infer a causal network (or a causal graph) with directed and undirected edges from observational data. It implements the MRPC (PC with the principle of Mendelian randomization) algorithm described in Badsha and Fu, 2019 and Badsha et al., 2018, and the implementation is based on the \link{pc} algorithm in the \code{pcalg} package. The MRPC algorithm contains four major updates over the \link{pc} algorithm: (i) incorporating a sequential testing method to control the False Discovery Rate (FDR), (ii) improved v-structure identification; (iii) allowing for calculation of robust correlation to reduce the impact of outliers, and (iv) a new procedure for edge orientation based on the principle of Mendelian randomization (PMR) (Badsha and Fu, 2019 and Badsha et al., 2018). See details below. 

}
\usage{
MRPC(data, suffStat, GV, FDR = 0.05, alpha = 0.05, indepTest, labels, p,
    fixedGaps = NULL, fixedEdges = NULL,
    NAdelete = TRUE, m.max = Inf,
    u2pd = c("relaxed", "rand", "retry"),
    skel.method = c("stable", "original", "stable.fast"),
    conservative = FALSE,
    maj.rule = FALSE, solve.confl = FALSE, FDRcontrol = TRUE,
    verbose = FALSE)
}

\arguments{
This function is based on the \link{pc} function in the \code{pcalg} package. Therefore, many arguments are similar to those in \link{pc}.

\item{data}{
Data matrix, where the rows are observations and the columns are features (i.e., variables, or nodes). If genetic variants (GVs) are included, then the columns start from GVs (e.g., single-nucleotide polymorphisms, or SNPs; insertions and deletions, or indels; copy number variation, or CNVs; and expression quantitative trait loci, or eQTLs to genes), and followed by phenotypes (e.g., gene expression). For example, if the data contains one GV, then the first column of the input matrix is the GV and the remaining columns are the gene expression data.
}

\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistics consist of the correlation matrix of the data and the sample size.
}

\item{GV}{The number of genetic variants (SNPs/indels/CNV/eQTL) in the input data matrix. For example, if the data has one SNPs/indels/CNV/eQTL, first column, then GV = 1, if 2 SNPs/indels/CNV/eQTL, 1st and 2nd Column, then GV = 2, if no GV then GV = 0, and so on.
}

\item{FDR}{
Need to specify the desired level of the overall false discovery rate. 
}
\item{alpha}{significance level (number in (0,1) for the individual tests.}
    
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

\item{solve.confl}{(optional) If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. In this case, only option u2pd = "relaxed" is supported. Note that the resulting object might not be a CPDAG because bi-directed edges might be present. See details for more information.}

\item{FDRcontrol}{(optional) The default is TRUE which is used sequential FDR control method, otherwise used fixed significance level for the individual tests.}

\item{verbose}{(optional) If TRUE, detailed output is provided. The default is FALSE which does not print output details.}

}
\details{
The PC algorithm is computationally efficient for learning a directed acyclic graph (Spirtes et al., 2000).  Several variants of the original PC algorithms are available (Kalisch and Buhlmann, 2007; Kalisch et al., 2012). Similar to these PC-like algorithms, our MRPC algorithm also contains two main steps:
    
\bold{Step-1}: Inference of the graph skeleton.  A graph skeleton is an undirected graph with edges that are supported by the data.  Similar to existing PC-like algorithms, we perform statistical tests for marginal and conditional independence tests.  If the null hypothesis of independence is not rejected, then the corresponding edge is removed and never tested again.  
    
However, unlike existing algorithms, which control only the type I error rate for each individual test, MRPC implements the LOND (Level On the Number of Discoveries) method (Javanmard and Montanari, 2015), which is a sequential hypothesis testing procedure and sets the significance level for each test based on the number of discoveries (i.e., rejections), to control the overall false discovery rate (FDR).  See \link{ModiSkeleton}.
    
Genome data may have outliers that drastically alter the topology of the inferred graph. MRPC allows for the estimate of robust correlation, which may be the substitute of the Pearson correlation as the input to graph inference (Badsha et al., 2013).    

\bold{Step-2}: Edge orientation.  With the graph skeleton inferred from Step 1, we orient each edge that is present in the graph.  MRPC is fundamentally different from algorithms in the \code{pcalg} (Kalisch and Buhlmann, 2007; Kalisch et al., 2012) and \code{bnlearn} (Scutari, 2010) packages in the following ways (see \link{EdgeOrientation}):
    
(i) When analyzing genomic data, genetic variants provide additional information that helps distinguish the casual direction between two genes.  Our MRPC algorithm incorporates the principle of Mendelian randomization in graph inference, which greatly reduces the space of possible graphs and increases the inference efficiency.
    
(ii) Next or if the input is not genomic data, we search for possible triplets 
that may form a v-structure (e.g.,X-->Y<--Z). We check conditional test results
from step I to see whether X and Z are independent given Y. If they are, then 
this is not a v-structure; alternative models for the triplet may be any of the 
following three Markov equivalent graphs: X-->Y-->Z, X<--Y<--Z, and 
X<--Y-->Z. If this test is not performed in the first step, we conduct it in 
this step. This step improves the accuracy of the v-structure identification over
existing methods. 
    
(iii) If there are undirected edges after steps (i) and (ii), we examine iteratively triplets of nodes with at least one directed edge and no more than one undirected edge. We check the marginal and conditional test results from Step I to determine which of the basic models is consistent with the test results. It is plausible that some undirected edges cannot be oriented, and we leave them as undirected.
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
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899.

3. Badsha MB, Mollah MN, Jahan N and Kurata H (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng, 116(3): 397-407.

4. Javanmard A and Montanari A (2015). On Online Control of False Discovery Rate. arXiv:150206197 [statME].

5. Kalisch M and Buhlmann P (2007). Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

6. Kalisch M, Machler M, Colombo D, Maathuis MH and Buhlmann P (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47, 26.

7. Scutari M (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

8. Spirtes P, Glymour C and Scheines R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
}
                    
\seealso{
\link{ModiSkeleton} for inferring a graph skeleton (i.e., an undirected graph); \link{EdgeOrientation} for edge orientation in the inferred graph skeleton; \link{SimulateData} for generating data under a topology.
}
\examples{
# Load packages
# We compare different simulated data across six methods: MRPC, 
# PC in pcalg (Kalisch et al., 2012), and pc.stable, mmpc mmhc and hc in
# bnlearn (Marco Scutari, 2010)

library(MRPC)     # MRPC
library(pcalg)    # pc
library(bnlearn)  # pc.stable, mmpc, mmhc and hc

# Data pre-processing
# The 1st column of the input matrix will be the genotype of the
# expression quantitative trait loci (eQTL)/Copy number variation (CNV)
# and the remaining columns are the gene expression data.

# We used pre-assigned level alpha = 0.05 that ensures FDR and mFDR
# will remain below 0.05.

# Load predefined simulated data
# Model 1
Truth <- MRPCtruth$M1   # Truth for model 1
data <- simu_data_M1    # load data for model 1
n <- nrow (data)        # Number of rows
V <- colnames(data)     # Column names

# Calculate Pearson correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Infer the graph by MRPC
MRPC.fit <- MRPC(data,
                suffStat = suffStat_C,
                GV = 1,
                FDR = 0.05,
                indepTest = 'gaussCItest',
                labels = V,
                FDRcontrol = TRUE,
                verbose = FALSE)

# Infer the graph by PC
pc.fit <- pc(suffStat = suffStat_C,
             indepTest = gaussCItest,
             alpha = 0.05,
             labels = V,
             verbose = FALSE)

# arcs not to be included from gene expression to genotype used in pc.stable, mmpc 
bl <- data.frame (from = colnames (data)[-1], 
                  to = 'V1')

# Infer the graph by pc.stable
pc.stable.fit <- pc.stable(data.frame(data),
                           blacklist = bl,
                           undirected = FALSE) 

# Infer the graph by mmpc
mmpc.fit <- mmpc(data.frame(data),
                 blacklist = bl,
                 undirected = FALSE,
                 debug = FALSE) 

# Infer the graph by mmhc
mmhc.fit <- mmhc(data.frame(data),
                 blacklist = bl,
                 debug = FALSE)

# Infer the graph by hc
hc.fit <- hc (data.frame (data), 
                          blacklist = bl, 
                          debug = FALSE)

# Plot the inferred graphs
par(mfrow = c(1, 7))
plot(Truth,
     main = "(A) Truth")
plot(MRPC.fit,
     main = "(B) MRPC")
plot(pc.fit,
    main ="(C) pc")
graphviz.plot(pc.stable.fit,
    main = "(D) pc.stable")
graphviz.plot(mmpc.fit,
    main = "(E) mmpc")
graphviz.plot(mmhc.fit,
    main = "(F) mmhc")
graphviz.plot(hc.fit,
    main = "(G) hc")
    
# Another option for plot of the results. First fig is the nodes
# dendrogram with colored modules. Second fig is the plot of the graph
# with color based on modules.
# To idendify modules and complex graph (Suitable if you have many nodes)
# Adjacency matrix from directed graph
Adj_directed <- as(MRPC.fit@graph,
                   "matrix")

# Plot of dendrogram with modules colors of nodes
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                                    minModuleSize = 2)
                  
# Visualization of inferred graph with modules colors
PlotGraphWithModulesObj <- PlotGraphWithModules(Adj_directed,
                                                PlotDendrogramObj,
                                                GV = 1,
                                                node.size = 8,
                                                arrow.size = 5,
                                                label.size = 3,
                                                alpha = 1) 

# Plot 
plot(PlotGraphWithModulesObj) 


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
