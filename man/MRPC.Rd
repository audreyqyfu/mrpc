\name{MRPC}
\alias{MRPC}
\title{
Estimate the DAG using Principle of Mendelian Randomization (PMR) based PC (MRPC) Algorithm
}
\description{
The MRPC () function is used to estimate a directed acyclic graph (DAG) from observational data, by modifying of the existing \link{pc} algorithm (Spirtes et al., 2000). Existing \link{pc} algorithms is used for estimating the DAG from observational data and several variants are implemented in the R package called pcalg. Several challenges remain when \link{pc} algorithms are applied to noisy genomic data. To tackle these challenges, we develop a novel machine learning algorithm called Mendelian Randomization (MR) based PC (MRPC) to efficiently learn a causal gene regulatory network (or a causal graph of genes) using genotype and gene expression data. Unlike existing \link{pc} algorithms, the four-major updated in the MRPC algorithm is that (i) we are incorporating sequential testing method to controls the False Discovery Rate (FDR), (ii) take a robust approach to reduces the impact of outliers, (iii)improved v-structure identification and (iv) implemented a new way for edge direction determination based on PMR (Badsha and Fu, 2018 and Badsha et al., 2018). See details below. 

}
\usage{
MRPC(data,suffStat,GV,FDR=0.05,indepTest, labels, p,
    fixedGaps = NULL, fixedEdges = NULL,
    NAdelete = TRUE, m.max = Inf,u2pd =
    c("relaxed", "rand", "retry"),
    skel.method = c("stable", "original",
    "stable.fast"),conservative = FALSE,
    maj.rule= FALSE, solve.confl = FALSE,
    verbose = FALSE)
}

\arguments{
The MRPC is the modification of the existing \link{pc}. Therefore, most of the terms and sentences are taken from \link{pc} (Spirtes et al., 2000).

\item{data}{
Data matrix, where row is the samples and column is the genes. The column order from  Single-nucleotide polymorphism (SNPs), indels,copy number variation (CNV) or Expression quantitative trait loci (eQTL) to genes.
}
\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistic consists of the correlation matrix of the data and the sample size.
}
\item{GV}{Need to specify for how many genetic variant (SNPs/indels/CNV/eQTL) in your input data matrix. For example, if dataset have one SNPs/indels/CNV/eQTL (first column), then GV=1, if 2 SNPs/indels/CNV/eQTL (1st and 2nd Column), then GV=2, if no GV then GV=0 and so on.
}
\item{FDR}{
Need to specify pre-assigned level.If FDR=0.05, that ensures FDR and mFDR remains below 0.05.
}
\item{indepTest}{A \link{function} for testing conditional independence. It is used to tests the conditional independence of x and y given S, called as indepTest(x,y,S,suffStat). Where, x and y are variables, and S is a (possibly empty) vector of variables. suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. The different indepTest is used for different data types, for example, Gaussian data = gaussCItest, Discrete data = disCItest and Binary data = binCItest. See help(gaussCItest)

The \link{ci.test} (Marco Scutari, 2010) is also used for testing conditional independence and return value of indepTest is the p-value. If none is specified, the default test statistic is the mutual information for categorical variables, the Jonckheere-Terpstra test for ordered factors and the linear correlation for continuous variables.See help(ci.test)

Remember that need to specify the which indepTest would like for independence testing. For example, if you would like to use gaussCItest you would type indepTest='gaussCItest' into the function otherwise indepTest='citest'. Note that, we used \link{gaussCItest} to compare our \link{MRPC} with the existing \link{pc}, because of \link{ci.test} is not robust. See details in example.
}

\item{labels}{character vector of variable (or node) names (all variables are denoted in column in the input matrix).}

\item{p}{(optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.}

\item{fixedGaps}{(optional) A logical matrix of dimension p*p. If entry [x,y] or [y,x] (or both) are TRUE, the edge x---y is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.}

\item{fixedEdges}{ (optional) A logical matrix of dimension p*p. If entry [x,y] or [y,x] (or both) are TRUE, the edge x---y is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.}

\item{NAdelete}{(optional) If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.}

\item{m.max}{ (optional) Maximal size of the conditioning sets that are considered in the conditional independence tests.}

\item{u2pd}{ (optional) String specifying the method for dealing with conflicting information.}

\item{skel.method}{ (optional) Character string specifying method; the default, "stable" provides an order-independent skeleton.}

\item{conservative}{(optional) Logical indicating if the conservative PC is used. In this case, only option u2pd = "relaxed" is supported. Note that therefore the resulting object might not be extendable to a DAG. See details for more information.}

\item{maj.rule}{(optional) Logical indicating that the triplets shall be checked for ambiguity using a majority rule idea, which is less strict than the conservative PC algorithm. For more information, see details.}

\item{solve.confl}{(optional) If TRUE, the orientation of the v-structures and the orientation rules work with lists for candidate sets and allow bi-directed edges to resolve conflicting edge orientations. In this case, only option u2pd = relaxed is supported. Note that therefore the resulting object might not be a CPDAG because bi-directed edges might be present. See details for more information.}

\item{verbose}{(optional) If TRUE, detailed output is provided. Default is FALSE for no output details.}

}
\details{
Nodes are used as a major building blocks of the models, which represent random variables and edges, which encode conditional dependence relations of the enclosing vertices (Kalisch and Buhlmann, 2007). The structure of conditional independence among the random variables can be explored using the Markov properties. The directed edges show the presence and direction of direct causal effects. A bidirected edges means that the edge orientation should be forward (--->) or bacward (<---). PC algorithms are computationally efficient for learning the underlying DAG (Kalisch and Buhlmann, 2007). Our MPRC algorithm improves on the existing PC algorithms in the following aspects:

(i) the genotype data at genetic variants (e.g., SNPs and copy number variation) provide additional information that helps to distinguish the casual direction between two genes; this is the rationale behind the Mendelian Randomization (MR). MR can greatly reduce the space of possible graphs and increase the inference efficiency for inference of genomic data but PC does not rely on MR and determination of the edge direction can be tricky when graph is large.

(ii) the number of statistical tests is unknown beforehand in all existing variants of PC algorithms. It is unclear how these algorithms control the false discovery rate (FDR) or marginal FDR (mFDR), because it is used alpha as a fixed pre-assigned significance level for testing the independence test, but results are sensitive to the choice of alpha, which is an important issue for multiple hypotheses testing.

(iii) the gene expression data often have outliers, even after normalization, which may drastically alter the topology of the inferred graph. 

(iv) to find the v-structure PC-algorithm assumed that if no evidence of two nodes are conditionally independent,then the nodes are conditionally dependent, although we don't have any evidence about their conditional dependence. In that case the orient v-structure is worng.

To overcome the above challenges, we proposed a novel algorithm MRPC (MR-based PC), which can be applied to genotype and gene expression data and efficiently learn a causal graph of genes. \bold{Firstly}, we take a robust approach (Badsha et al., 2013) and calculate the robust correlation matrix on which the series of hypothesis testing is performed. \bold{Secondly}, we adopt a sequential method LOND (significance Levels based On Number of Discoveries) (Javanmard and Montanari, 2015) that controls the FDR or marginal FDR (mFDR) in an online manner, where level of significance is used a function of the previous decision made so far (we adjusted  at each step in the test) and that ensures that FDR and mFDR remains below the total alpha.
\bold{Thirdly}, we improved v-structure identification using the addtional conditional test. \bold{Finally}, we implemented a new way for edge direction determination based on Mendelian randomization.

The two main step in MRPC algorithm as follows:

\bold{Step-1}: We incorporated sequential hypothesis testing to draw the undirected graph (skeleton) by \link{ModiSkeleton} function  (similar as the pc algorithm by \link{skeleton} function but, difference is that the \link{ModiSkeleton} function is used alpha (significance level) as a function of the prior outcomes (reject or accept hypothesis) but the \link{skeleton} function used alpha as a fixed pre-assigned significance level and \link{ModiSkeleton} take a robust correlation approach reduces the impact of outliers but the \link{skeleton} used classical correlation which is highly influnced by outliers. See details in in \link{ModiSkeleton} and \link{skeleton}).

\bold{Step-2}: We implemented a new way for edge direction determination based on Mendelian randomization (MR). We consider the first column of the input matrix will be the genotype of the SNPs/indels/CNV/eQTL and the remaining column are the gene expression data. See details in \link{EdgeOrientation}.

All statistical inference is done in the \bold{Step-1}, while the \bold{Step-2} is just application of deterministic rules on the results of the \bold{Step-1}. If the first part is done correctly, the second part will never fail. If, however, there occur errors in the \bold{Step-1}, the \bold{Step-2} will be more sensitive to it, since it depends on the inferential results of \bold{Step-1}.
}

\value{
An object of \link{class} that containing an estimate of the equivalence class of the underlying DAG.
\describe{
    
    \item{\code{call}:}{a \link{call} object: the original function call.}
    
    \item{\code{n}:}{an "integer", the sample 
    size used to estimate the graph.}
        \item{\code{max.ord}:}{an integer, the 
    maximum size of the conditioning set used 
    in the conditional independence tests of   
    the (first part of the algorithm).}
    \item{\code{n.edgetests}:}{the number of 
    conditional independence tests performed by
    the (first part of the) algorithm.}
    \item{\code{sepset}:}{ separation sets.}
    \item{\code{pMax}:}{a numeric square matrix
    , where the (i, j)th entry contains the 
    maximal p-value of all conditional
    independence tests for edge i--j.}
    \item{\code{graph}:}{Object of class \code{"\linkS4class{graph}"}:
      the undirected or partially directed graph that was estimated.}
          \item{\code{zMin}:}{Deprecated.}
    \item{\code{test}:}{an "integer", no. of tests have been performed.}
    \item{\code{alpha}:}{"Non integer", level of significance 
    for the current test.}
    \item{\code{R}:}{a "vector", all decisions made so far from tests have been performed.}
    
  }
}

\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>,Evan A Martin <mart9986@vandals.uidaho.edu> and Audrey Qiuyan Fu <audreyf@uidaho.edu>
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
\link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{EdgeOrientation} for orientation rules determination for edges in MRPC; \link{SimulatedData} for simulated data generating function.
}
\examples{
## Load packages
#We compare different simulated data across three methods: MRPC, 
#the pc method in pcalg (Kalisch et al., 2012) and the mmhc method 
#in bnlearn (Marco Scutari, 2010)
library(pcalg)    #library for existing pc
library(bnlearn)  #library for mmhc

## Load predefined simulated data
## Data pre-processing

#The 1st column of the input matrix will be the genotype of the
#expression quantitative trait loci (eQTL)/Copy number variation (CNVs)
#and the remaining columns are the gene expression data.

# We used pre-assigned level alpha=0.05 that ensures FDR and mFDR
#remains below 0.05.

## Model 1
Truth<-MRPCtruth$M1   #Truth for model 1
data<-simu.data_M1    #data load for model 1
n<-nrow (data)        #Number of row
V<-colnames(data)     #Column names

#Classical correlation
suffStat_C<- list(C = cor(data), n = n)

#Robust correlation (Beta=0.005)
Rcor_R<- RobustCor(data, 0.005) 
suffStat_R<- list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit<- MRPC(data,suffStat_R,GV=1,FDR=0.05,
indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc.fit<- mmhc(data)

## Estimated graph by existing pc with alpha=0.05
pc.fit<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

## Plot of the results
## Show estimated graph
#if (require(Rgraphviz)){
par(mfrow=c(2,2))
plot(Truth,main ="(A) Truth")
plot(MRPC.fit, main ="(B) MRPC")
graphviz.plot(mmhc.fit, main ="(C) mmhc")
plot(pc.fit,  main ="(D) pc")
#}
#Another option for plot of the results. First fig is the nodes
#dendrogram with modules color. second fig is the plot of the graph
#with color based on modules.
#To idendify modules and complex graph (Suitable if you have many
#nodes)
#Adjacency matrix from directed graph
Adj_directed<-as(MRPC.fit@graph,"matrix")
#Plot of the graph with modules color
DendroModuleGraph(Adj_directed,minModuleSize = 1,GV=1)

## Other models are available and may be called as follows:
# Model 0
# Truth<-MRPCtruth$M0
# data<-simu.data_M0,
# Model 2
# Truth<-MRPCtruth$M2
# data<-simu.data_M2
# Model 3
# Truth<-MRPCtruth$M3
# data<-simu.data_M3
# Model 4
# Truth<-MRPCtruth$M4
# data<-simu.data_M4
# Model Multiparent
# Truth<-MRPCtruth$Multiparent
# data<-simu.data_multiparent
# Model Star
# Truth<-MRPCtruth$Star
# data<-simu.data_starshaped
# Model Layered
# Truth<-MRPCtruth$Layered
# data<-simu.data_layered
}
