\name{MRPC}
\alias{MRPC}
\title{
Estimate the DAG using the Mendelian Randomization (MR) based PC (MRPC) Algorithm
}
\description{
The MRPC () function is used to estimate a directed acyclic graph (DAG) from observational data, by modifying of the existing \link{pc} algorithm (Spirtes et al., 2000). Existing \link{pc} algorithms is used for estimating the DAG from observational data and several variants are implemented in the R package called pcalg. Several challenges remain when \link{pc} algorithms are applied to noisy genomic data. To tackle these challenges, we develop a novel machine learning algorithm called Mendelian Randomization (MR) based PC (MRPC) to efficiently learn a causal gene regulatory network (or a causal graph of genes) using genotype and gene expression data. Unlike existing \link{pc} algorithms, the four-major updated in the MRPC algorithm is that (i) we are incorporating sequential testing method to controls the False Discovery Rate (FDR), (ii) take a robust approach to reduces the impact of outliers, (iii)improved v-structure identification and (iv) implemented a new way for edge direction determination based on MR. See details below.

}
\usage{
MRPC(data,suffStat,NQ,FDR,indepTest, labels, p,
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
\item{NQ}{Need to specify for how many SNPs/indels/CNV/eQTL in your input data matrix. For example, If your dataset have one SNPs/indels/CNV/eQTL (first column), then NQ=1, if 2 SNPs/indels/CNV/eQTL (1st and 2nd Column), then NQ=2..and so on.
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

\bold{Step-2}: We implemented a new way for edge direction determination based on Mendelian randomization (MR). We consider the first column of the input matrix will be the genotype of the SNPs/indels/CNV/eQTL and the remaining column are the gene expression data. See details in \link{EdgesOrientation}.

All statistical inference is done in the \bold{Step-1}, while the \bold{Step-2} is just application of deterministic rules on the results of the \bold{Step-1}. If the first part is done correctly, the second part will never fail. If, however, there occur errors in the \bold{Step-1}, the \bold{Step-2} will be more sensitive to it, since it depends on the inferential results of \bold{Step-1}.
}

\value{
An object of \link{class} that are containing an estimate of the underlying DAG.

}
\references{
1.Badsha, M. B. Mollah, M. N. Jahan, N. and Kurata, H. (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Javanmard and Montanari (March 5, 2015) On Online Control of False Discovery Rate. arXiv:150206197 [statME].

3. Kalisch, M. and Buhlmann, P. (2007) Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

4.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

5. Spirtes,P.  Glymour, C. and Scheines, R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.

}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{EdgesOrientation} for orientation rules determination for edges in MRPC; \link{SimulatedData} for simulated data generating function.
}
\examples{
## Load packages
library(pcalg)  #library for existing pc
library(bnlearn)#library for ci.test

## Load predefined simulated data
## Data pre-processing

#The 1st column of the input matrix will be the genotype of the
#expression quantitative trait loci (eQTL)/Copy number variation (CNVs)
#and the remaining columns are the gene expression data.

# We used pre-assigned level alpha=0.05 that ensures FDR and mFDR
#remains below 0.05.

## Model 0
n<-nrow (simu.data_M0)    #Number of row
V<-colnames(simu.data_M0) #Column names

Rcor_C=RobustCor(simu.data_M0, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M0, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)


## Estimated graph by existing pc
#We used apha= 0.05 and 0.01, because, two different
#value of alpha produces two different topologies for

#Model 0 by existing PC.
pc.fit_01 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =
V,verbose = TRUE)
pc.fit_02 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.01,labels =
V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_0_GT <- MRPC(simu.data_M0,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_0_CT <- MRPC(simu.data_M0,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels=V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_01,  main ="(A) Model 0(by pc (alpha=0.05))")
plot(pc.fit_02,  main ="(B) Model 0(by pc(alpha=0.01))")
plot(MRPC.fit_0_GT, main ="(C) Model 0(by MRPC_GT)")
plot(MRPC.fit_0_CT, main ="(C) Model 0(by MRPC_CT)")

#Another option for plot of the results. First fig is the nodes
#dendrogram with modules color. second fig is the plot of the graph
#with color based on modules.
#To idendify modules and complex graph (Suitable if you have many
#nodes)
#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_0_GT@graph,"matrix")
#Plot of the graph with modules color
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)

## Model 1
n<-nrow (simu.data_M1)     #Number of row
V<-colnames(simu.data_M1)  #Column names

Rcor_C=RobustCor(simu.data_M1, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate graph by existing pc
pc.fit_1 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V
,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_1_GT <- MRPC(simu.data_M1,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_1_CT <- MRPC(simu.data_M1,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels=V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_1,  main ="(A) Model 1(by pc)")
plot(MRPC.fit_1_GT, main ="(B) Model 1(by MRPC_GT)")
plot(MRPC.fit_1_CT, main ="(C) Model 1(by MRPC_CT)")

#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_1_GT@graph,"matrix")
#Plot of the graph with modules color
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)

## Model 2
n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M2) #Column names

Rcor_C=RobustCor(simu.data_M2, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate graph by existing pc
pc.fit_2 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V
,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_2_GT <- MRPC(simu.data_M2,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_2_CT <- MRPC(simu.data_M2,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels=V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_2,  main ="(A) Model 2(by pc)")
plot(MRPC.fit_2_GT, main ="(B) Model 2(by MRPC_GT)")
plot(MRPC.fit_2_CT, main ="(C) Model 2(by MRPC_CT)")

#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_2_GT@graph,"matrix")
#Plot of the graph with modules color
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)

## Model 3
n<-nrow (simu.data_M3)    #Number of row
V<-colnames(simu.data_M3) #Column names

Rcor_C=RobustCor(simu.data_M3, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M3, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate graph by existing pc
#We used apha= 0.05 and 0.01, because, two different
#value of alpha produces two different topologies for

#model 3 by existing PC.
pc.fit_31 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =
V,verbose = TRUE)
pc.fit_32 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.01,labels =
V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_3_GT <- MRPC(simu.data_M3,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_3_CT <- MRPC(simu.data_M3,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels=V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_31,  main ="(A) Model 3(by pc(alpha=0.05))")
plot(pc.fit_32,  main ="(B) Model 3(by pc(alpha=0.01))")
plot(MRPC.fit_3_GT, main ="(C) Model 3(by MRPC_GT)")
plot(MRPC.fit_3_CT, main ="(D) Model 3(by MRPC_CT)")

#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_3_GT@graph,"matrix")
#Plot of the graph with modules color
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)

## Model 4
n<-nrow (simu.data_M4)    #Number of row
V<-colnames(simu.data_M4) #Column names

Rcor_C=RobustCor(simu.data_M4, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M4, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate graph by existing pc
pc.fit_4 <- pc(suffStat_R,indepTest = gaussCItest,alpha=0.05,labels = V
,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_4_GT <- MRPC(simu.data_M4,suffStat_R,NQ=1,FDR=0.05,indepTest =
'gaussCItest', labels= V, verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_4_CT <- MRPC(simu.data_M4,suffStat_R,NQ=1,FDR=0.05,indepTest =
'citest',labels=V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_4, main ="(A) Model 4(by pc)")
plot(MRPC.fit_4_GT,main ="(B) Model 4(by MRPC_GT)")
plot(MRPC.fit_4_CT,main ="(C) Model 4(by MRPC_CT)")

#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_4_GT@graph,"matrix")
#Plot of the graph with modules color
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)

##Start shape topology
n<-nrow (simu.data_StarModel)    #Number of row
V<-colnames(simu.data_StarModel) #Column names

Rcor_C=RobustCor(simu.data_StarModel, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_StarModel, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_StarModel_05<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05
,labels = V,verbose = TRUE)
pc.fit_StarModel_01<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.01
,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_StarModel_GT <- MRPC(simu.data_StarModel,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_StarModel_CT <- MRPC(simu.data_StarModel,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'citest',labels=V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_StarModel_05,  main ="(A) Star topology (by pc (alpha=0.05))")
plot(pc.fit_StarModel_01,  main ="(B) Star topology (by pc (alpha=0.01))")
plot(MRPC.fit_StarModel_GT , main ="(C) Star topology (by MRPC_GT)")
plot(MRPC.fit_StarModel_CT , main ="(D) Star topology (by MRPC_CT)")

##Multiple Parents Model
n<-nrow (simu.data_MultipleParentsModel)    #Number of row
V<-colnames(simu.data_MultipleParentsModel) #Column names

Rcor_C=RobustCor(simu.data_MultipleParentsModel, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_MultipleParentsModel, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)


### Estimated graph by pc
pc.fit_MultipleParentsModel_05 <- pc(suffStat_C,indepTest = gaussCItest
,alpha=0.05,labels = V,verbose = TRUE)
pc.fit_MultipleParentsModel_01 <- pc(suffStat_C,indepTest = gaussCItest
,alpha=0.01,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_MultipleParentsModel_GT <- MRPC(simu.data_MultipleParentsModel
,suffStat_R,NQ=1,FDR=0.05,indepTest ='gaussCItest',labels = V,verbose =TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_MultipleParentsModel_CT <- MRPC(simu.data_MultipleParentsModel
,suffStat_R,NQ=1,FDR=0.05,indepTest = 'citest',labels=V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_MultipleParentsModel_05,  main ="(A) Multiple Parents topology
(by pc (alpha=0.05))")
plot(pc.fit_MultipleParentsModel_01 ,  main ="(B) Multiple Parents
topology (by pc (alpha=0.01))")
plot(MRPC.fit_MultipleParentsModel_GT,  main ="(C) Multiple Parents (by MRPC_GT)")
plot(MRPC.fit_MultipleParentsModel_CT,  main ="(D) Multiple Parents (by MRPC_CT)")


##Layer shape topology
n<-nrow (simu.data_LayersModel)    #Number of row
V<-colnames(simu.data_LayersModel) #Column names

Rcor_C=RobustCor(simu.data_LayersModel, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_LayersModel, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_LayersModel_05 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05
,labels = V,verbose = TRUE)
pc.fit_LayersModel_01 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.01
,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_LayersModel_GT <- MRPC(simu.data_LayersModel,suffStat_R,NQ=1,
FDR=0.05,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_LayersMode_CT <- MRPC(simu.data_LayersModel,suffStat_R,NQ=1,
FDR=0.05,indepTest = 'citest',labels=V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_LayersModel_05 ,  main ="(A) Layer topology (by pc (alpha=0.05))")
plot(pc.fit_LayersModel_01 ,  main ="(B) Layer topology (by pc (alpha=0.01))")
plot(MRPC.fit_LayersModel_GT,  main ="(C) Layer topology (by MRPC_GT)")
plot(MRPC.fit_LayersMode_CT,  main ="(D) Layer topology (by MRPC_CT)")


#In SimulatedData generating function, we generate the basic 5
#topologies for three nodes and star shape topology, Layer shape
#topology. In this example, we consider four nodes and we generated
#another 36 topologies. We generate those 36 topologies based on the
#previous 5 basic topologies.
#We also,showed example model from pcalg and multiple SNPs/indels/CNV/eQTL topology

#Model 5 (extention of model 0) based on model 0
#Data

V=simu.data_M0[ ,1];
T1=simu.data_M0[ ,2];
T3=simu.data_M0[ ,3];
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M5=cbind(V,T1,T2,T3)
colnames(simu.data_M5)=c("V","T1.5","T2.5","T3.5")

# Estimated graph
n<-nrow (simu.data_M5)    #Number of row
V<-colnames(simu.data_M5) #Column names

Rcor_C=RobustCor(simu.data_M5, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M5, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_5<- pc(suffStat_C,indepTest=gaussCItest,alpha=0.05,labels=V
,verbose= TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_5_GT <- MRPC(simu.data_M5,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_5_CT <- MRPC(simu.data_M5,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_5,  main ="(A) Model 5(by pc (alpha=0.05))")
plot(MRPC.fit_5_GT, main ="(B) Model 5(by MRPC_GT)")
plot(MRPC.fit_5_CT, main ="(C) Model 5(by MRPC_CT)")


#Model 6 (extention of model 0) based on model 1

#Data
V=simu.data_M1[ ,1];
T1=simu.data_M1[ ,2];
T3=simu.data_M1[ ,3];
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M6=cbind(V,T1,T2,T3)
colnames(simu.data_M6)=c("V","T1.6","T2.6","T3.6")

#Estimated graph
n<-nrow (simu.data_M6)    #Number of row
V<-colnames(simu.data_M6) #Column names

Rcor_C=RobustCor(simu.data_M6, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M6, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by pc
pc.fit_6<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_6_GT <- MRPC(simu.data_M6,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_6_CT <- MRPC(simu.data_M6,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_6,  main ="(A) Model 6(by pc (alpha=0.05))")
plot(MRPC.fit_6_GT, main ="(B) Model 6(by MRPC_GT)")
plot(MRPC.fit_6_CT, main ="(C) Model 6(by MRPC_CT)")


#Model 7 (extention of model 0) based on model 2

#Data
V=simu.data_M2[ ,1];
T1=simu.data_M2[ ,2];
T3=simu.data_M2[ ,3];
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M7=cbind(V,T1,T2,T3)
colnames(simu.data_M7)=c("V","T1.7","T2.7","T3.7")

#Estimated graph
n<-nrow (simu.data_M7)    #Number of row
V<-colnames(simu.data_M7) #Column names

Rcor_C=RobustCor(simu.data_M7, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M7, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_7<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_7_GT <- MRPC(simu.data_M7,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_7_CT <- MRPC(simu.data_M7,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_7,  main ="(A) Model 7(by pc (alpha=0.05))")
plot(MRPC.fit_7_GT, main ="(B) Model 7(by MRPC_GT)")
plot(MRPC.fit_7_CT, main ="(C) Model 7(by MRPC_CT)")

#Model 8 (extention of model 0) based on model 3

#Data
V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T3=simu.data_M3[ ,3]
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M8=cbind(V,T1,T2,T3)
colnames(simu.data_M8)=c("V","T1.8","T2.8","T3.8")

#Estimated graph
n<-nrow (simu.data_M8)    #Number of row
V<-colnames(simu.data_M8) #Column names

Rcor_C=RobustCor(simu.data_M8, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M8, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_8<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_8_GT <- MRPC(simu.data_M8,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_8_CT <- MRPC(simu.data_M8,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_8,  main ="(A) Model 8(by pc (alpha=0.05))")
plot(MRPC.fit_8_GT, main ="(B) Model 8(by MRPC_GT)")
plot(MRPC.fit_8_CT, main ="(C) Model 8(by MRPC_CT)")

#Model 9 (extention of model 0) based on model 4

#Data
V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T3=simu.data_M4[ ,3]
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M9=cbind(V,T1,T2,T3)
colnames(simu.data_M9)=c("V","T1.9","T2.9","T3.9")

#Estimated graph
n<-nrow (simu.data_M9)    #Number of row
V<-colnames(simu.data_M9) #Column names

Rcor_C=RobustCor(simu.data_M9, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M9, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_9<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V
,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_9_GT <- MRPC(simu.data_M9,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_9_CT <- MRPC(simu.data_M9,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_9,  main ="(A) Model 9(by pc (alpha=0.05))")
plot(MRPC.fit_9_GT, main ="(B) Model 9(by MRPC_GT)")
plot(MRPC.fit_9_CT, main ="(C) Model 9(by MRPC_CT)")

#Model 10 (extention of model 1) based on model 0 (triplet=V,T1,T2)
#Data

V=simu.data_M0[ ,1]
T1=simu.data_M0[ ,2]
T3=simu.data_M0[ ,3]
T2=Case_1P(N=10^3,P1=T1, b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M10=cbind(V,T1,T2,T3)
colnames(simu.data_M10)=c("V","T1.10","T2.10","T3.10")

#Estimated graph
n<-nrow (simu.data_M10)    #Number of row
V<-colnames(simu.data_M10) #Column names

Rcor_C=RobustCor(simu.data_M10, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M10, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_10<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_10_GT <- MRPC(simu.data_M10,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_10_CT <- MRPC(simu.data_M10,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_10,  main ="(A) Model 10(by pc (alpha=0.05))")
plot(MRPC.fit_10_GT, main ="(B) Model 10(by MRPC_GT)")
plot(MRPC.fit_10_CT, main ="(C) Model 10(by MRPC_CT)")


#Model 11 (extention of model 1) based on model 1 (triplet=V,T1,T2)
#Data

V=simu.data_M1[ ,1]
T1=simu.data_M1[ ,2]
T2=simu.data_M1[ ,3]
T3=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M11=cbind(V,T1,T2,T3)
colnames(simu.data_M11)=c("V","T1.11","T2.11","T3.11")

#Estimated graph
n<-nrow (simu.data_M11)    #Number of row
V<-colnames(simu.data_M11) #Column names

Rcor_C=RobustCor(simu.data_M11, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M11, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_11<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_11_GT <- MRPC(simu.data_M11,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_11_CT <- MRPC(simu.data_M11,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_11,  main ="(A) Model 11(by pc (alpha=0.05))")
plot(MRPC.fit_11_GT, main ="(B) Model 11(by MRPC_GT)")
plot(MRPC.fit_11_CT, main ="(C) Model 11(by MRPC_CT)")

#Model 12 (extention of model 1) based on model 2 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T3=simu.data_M2[ ,3]
T2=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M12=cbind(V,T1,T2,T3)
colnames(simu.data_M12)=c("V","T1.12","T2.12","T3.12")

#Estimated graph
n<-nrow (simu.data_M12)    #Number of row
V<-colnames(simu.data_M12) #Column names

Rcor_C=RobustCor(simu.data_M12, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M12, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_12<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_12_GT <- MRPC(simu.data_M12,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_12_CT <- MRPC(simu.data_M12,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_12,  main ="(A) Model 12(by pc (alpha=0.05))")
plot(MRPC.fit_12_GT, main ="(B) Model 12(by MRPC_GT)")
plot(MRPC.fit_12_CT, main ="(C) Model 12(by MRPC_CT)")

#Model 13 (extention of model 1) based on model 3 (triplet=V,T1,T2)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T3=simu.data_M3[ ,3]
T2=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M13=cbind(V,T1,T2,T3)
colnames(simu.data_M13)=c("V","T1.13","T2.13","T3.13")

#Estimated graph
n<-nrow (simu.data_M13)    #Number of row
V<-colnames(simu.data_M13) #Column names

Rcor_C=RobustCor(simu.data_M13, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M13, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_13<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_13_GT <- MRPC(simu.data_M13,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_13_CT <- MRPC(simu.data_M13,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_13,  main ="(A) Model 13(by pc (alpha=0.05))")
plot(MRPC.fit_13_GT, main ="(B) Model 13(by MRPC_GT)")
plot(MRPC.fit_13_CT, main ="(C) Model 13(by MRPC_CT)")

#Model 14 (extention of model 1) based on model 4 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T3=simu.data_M4[ ,3]
T2=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M14=cbind(V,T1,T2,T3)
colnames(simu.data_M14)=c("V","T1.14","T2.14","T3.14")

#Estimated graph
n<-nrow (simu.data_M14)    #Number of row
V<-colnames(simu.data_M14) #Column names

Rcor_C=RobustCor(simu.data_M14, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M14, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_14<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_14_GT <- MRPC(simu.data_M14,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_14_CT <- MRPC(simu.data_M14,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)


## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_14,  main ="(A) Model 14(by pc (alpha=0.05))")
plot(MRPC.fit_14_GT, main ="(B) Model 14(by MRPC_GT)")
plot(MRPC.fit_14_CT, main ="(C) Model 14(by MRPC_CT)")

#Model 15 (extention of model 1) based on model 1 (triplet=T1,T2,T3)
#Data

V=simu.data_M1[ ,1]
T1=simu.data_M1[ ,2]
T2=simu.data_M1[ ,3]
T3=Case_1P(N=10^3,P1=T2,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M15=cbind(V,T1,T2,T3)
colnames(simu.data_M15)=c("V","T1.15","T2.15","T3.15")

#Estimated graph
n<-nrow (simu.data_M15)    #Number of row
V<-colnames(simu.data_M15) #Column names

Rcor_C=RobustCor(simu.data_M15, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M15, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_15<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_15_GT<- MRPC(simu.data_M15,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_15_CT<- MRPC(simu.data_M15,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_15,  main ="(A) Model 15(by pc (alpha=0.05))")
plot(MRPC.fit_15_GT, main ="(B) Model 15(by MRPC_GT)")
plot(MRPC.fit_15_CT, main ="(C) Model 15(by MRPC_CT)")

#Model 16 (extention of model 1) based on model 2 (triplet=T1,T2,T3)
#Data

V=simu.data_M1[ ,1]
T1=simu.data_M1[ ,2]
T3=simu.data_M0[ ,3]
T2=Case_2P(N=10^3,P1=T1,P2=T3,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
simu.data_M16=cbind(V,T1,T2,T3)
colnames(simu.data_M16)=c("V","T1.16","T2.16","T3.16")

#Estimated graph
n<-nrow (simu.data_M16)    #Number of row
V<-colnames(simu.data_M16) #Column names

Rcor_C=RobustCor(simu.data_M16, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M16, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_16<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_16_GT <- MRPC(simu.data_M16,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_16_CT <- MRPC(simu.data_M16,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_16,  main ="(A) Model 16(by pc (alpha=0.05))")
plot(MRPC.fit_16_GT, main ="(B) Model 16(by MRPC_GT)")
plot(MRPC.fit_16_CT, main ="(C) Model 16(by MRPC_CT)")

#Model 17 (extention of model 1) based on model 4 (triplet=T1,T2,T3)
#Data

V=simu.data_M1[ ,1]
T1=simu.data_M1[ ,2]

T2.a <- Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T3.a <- Case_2P(N=10^3,P1=T1,P2=T2.a,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
T3.b <- Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T2.b <- Case_2P(N=10^3,P1=T1,P2=T3.a,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T2 <- rep(0, 10^3)
T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]

T3 <- rep(0, 10^3)
T3[which(coinToss == 0)] <- T3.a[which(coinToss == 0)]
T3[which(coinToss == 1)] <- T3.b[which(coinToss == 1)]

simu.data_M17=cbind(V,T1,T2,T3)
colnames(simu.data_M17)=c("V","T1.17","T2.17","T3.17")

#Estimated graph
n<-nrow (simu.data_M17)    #Number of row
V<-colnames(simu.data_M17) #Column names

Rcor_C=RobustCor(simu.data_M17, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M17, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_17<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_17_GT <- MRPC(simu.data_M17,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_17_CT <- MRPC(simu.data_M17,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_17,  main ="(A) Model 17(by pc (alpha=0.05))")
plot(MRPC.fit_17_GT, main ="(B) Model 17(by MRPC_GT)")
plot(MRPC.fit_17_CT, main ="(C) Model 17(by MRPC_CT)")

#Model 18 (extention of model 2) based on model 0 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T2=simu.data_M2[ ,3]
T3=Case_NP(N=10^3, b0.1 = 0,sd.1 = 1)
simu.data_M18=cbind(V,T1,T2,T3)
colnames(simu.data_M18)=c("V","T1.18","T2.18","T3.18")

#Estimated graph
n<-nrow (simu.data_M18)    #Number of row
V<-colnames(simu.data_M18) #Column names

Rcor_C=RobustCor(simu.data_M18, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M18, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_18<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_18_GT <- MRPC(simu.data_M18,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_18_CT <- MRPC(simu.data_M18,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_18,  main ="(A) Model 18(by pc (alpha=0.05))")
plot(MRPC.fit_18_GT, main ="(B) Model 18(by MRPC_GT)")
plot(MRPC.fit_18_CT, main ="(C) Model 18(by MRPC_CT)")

#Model 19 (extention of model 2) based on model 1 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T2=simu.data_M2[ ,3]
T3=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M19=cbind(V,T1,T2,T3)
colnames(simu.data_M19)=c("V","T1.19","T2.19","T3.19")

#Estimated graph
n<-nrow (simu.data_M19)    #Number of row
V<-colnames(simu.data_M19) #Column names

Rcor_C=RobustCor(simu.data_M19, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M19, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_19<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_19_GT <- MRPC(simu.data_M19,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_19_CT <- MRPC(simu.data_M19,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_19,  main ="(A) Model 19(by pc (alpha=0.05))")
plot(MRPC.fit_19_GT, main ="(B) Model 19(by MRPC_GT)")
plot(MRPC.fit_19_CT, main ="(C) Model 19(by MRPC_CT)")

#Model 20 (extention of model 2) based on model 2 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T2=simu.data_M2[ ,3]
T3=Case_NP(N=10^3, b0.1 = 0,sd.1 = 1)
T1=Case_3P(N=10^3,P1=V,P2=T2,P3=T3, b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)
simu.data_M20=cbind(V,T1,T2,T3)
colnames(simu.data_M20)=c("V","T1.20","T2.20","T3.20")

#Estimated graph
n<-nrow (simu.data_M20)    #Number of row
V<-colnames(simu.data_M20) #Column names

Rcor_C=RobustCor(simu.data_M20, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M20, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_20<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_20_GT <- MRPC(simu.data_M20,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_20_CT <- MRPC(simu.data_M20,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_20,  main ="(A) Model 20 (by pc (alpha=0.05))")
plot(MRPC.fit_20_GT, main ="(B) Model 20 (by MRPC_GT)")
plot(MRPC.fit_20_CT, main ="(C) Model 20 (by MRPC_CT)")

#Model 21 (extention of model 2) based on model 3 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T2=simu.data_M2[ ,3]
T3=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M21=cbind(V,T1,T2,T3)
colnames(simu.data_M21)=c("V","T1.21","T2.21","T3.21")

#Estimated graph
n<-nrow (simu.data_M21)    #Number of row
V<-colnames(simu.data_M21) #Column names

Rcor_C=RobustCor(simu.data_M21, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M21, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_21<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_21_GT <- MRPC(simu.data_M21,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_21_CT <- MRPC(simu.data_M21,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_21,  main ="(A) Model 21 (by pc (alpha=0.05))")
plot(MRPC.fit_21_GT, main ="(B) Model 21 (by MRPC_GT)")
plot(MRPC.fit_21_CT, main ="(C) Model 21 (by MRPC_CT)")

#Model 22 (extention of model 2) based on model 4 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T2=simu.data_M2[ ,3]

T1.a <- Case_2P(N=10^3,P1=V,P2=T2,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
T3.a <- Case_2P(N=10^3,P1=V,P2=T1.a,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

T3.b <- Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T1.b <- Case_3P(N=10^3,P1=V,P2=T2,P3=T3.b,b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T1 <- rep(0, 10^3)
T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]

T3 <- rep(0, 10^3)
T3[which(coinToss == 0)] <- T3.a[which(coinToss == 0)]
T3[which(coinToss == 1)] <- T3.b[which(coinToss == 1)]
simu.data_M22=cbind(V,T1,T2,T3)
colnames(simu.data_M22)=c("V","T1.22","T2.22","T3.22")

#Estimated graph
n<-nrow (simu.data_M22)    #Number of row
V<-colnames(simu.data_M22) #Column names

Rcor_C=RobustCor(simu.data_M22, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M22, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_22<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_22_GT <- MRPC(simu.data_M22,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_22_CT <- MRPC(simu.data_M22,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_22,  main ="(A) Model 22 (by pc (alpha=0.05))")
plot(MRPC.fit_22_GT, main ="(B) Model 22 (by MRPC_GT)")
plot(MRPC.fit_22_CT, main ="(C) Model 22 (by MRPC_CT)")

#Model 23 (extention of model 2) based on model 3 (triplet=T2,T1,T3)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T2=simu.data_M2[ ,3]
T3=Case_1P(N=10^3,P1=T2,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M23=cbind(V,T1,T2,T3)
colnames(simu.data_M23)=c("V","T1.23","T2.23","T3.23")

#Estimated graph
n<-nrow (simu.data_M23)    #Number of row
V<-colnames(simu.data_M23) #Column names

Rcor_C=RobustCor(simu.data_M23, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M23, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_23<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_23_GT <- MRPC(simu.data_M23,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_23_CT <- MRPC(simu.data_M23,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_23,  main ="(A) Model 23 (by pc (alpha=0.05))")
plot(MRPC.fit_23_GT, main ="(B) Model 23 (by MRPC_GT)")
plot(MRPC.fit_23_CT, main ="(C) Model 23 (by MRPC_CT)")

#Model 24 (extention of model 2) based on model 4 (triplet=T2,T1,T3)
#Data

V=simu.data_M2[ ,1]
T2=simu.data_M2[ ,3]

T1.a <- Case_2P(N=10^3,P1=V,P2=T2,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
T3.a <- Case_2P(N=10^3,P1=T2,P2=T1.a,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

T3.b <- Case_1P(N=10^3,P1=T2,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T1.b <- Case_3P(N=10^3,P1=V,P2=T2,P3=T3.b,b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T1 <- rep(0, 10^3)
T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]

T3 <- rep(0, 10^3)
T3[which(coinToss == 0)] <- T3.a[which(coinToss == 0)]
T3[which(coinToss == 1)] <- T3.b[which(coinToss == 1)]

simu.data_M24=cbind(V,T1,T2,T3)
colnames(simu.data_M24)=c("V","T1.24","T2.24","T3.24")

#Estimated graph
n<-nrow (simu.data_M24)    #Number of row
V<-colnames(simu.data_M24) #Column names

Rcor_C=RobustCor(simu.data_M24, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M24, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_24<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_24_GT <- MRPC(simu.data_M24,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_24_CT <- MRPC(simu.data_M24,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_24,  main ="(A) Model 24 (by pc (alpha=0.05))")
plot(MRPC.fit_24_GT, main ="(B) Model 24 (by MRPC_GT)")
plot(MRPC.fit_24_CT, main ="(C) Model 24 (by MRPC_CT)")

#Model 25 (extention of model 3) based on model 0 (triplet=V,T1,T2)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T2=simu.data_M3[ ,3]
T3=Case_NP(N=10^3, b0.1 = 0,sd.1 = 1)
simu.data_M25=cbind(V,T1,T2,T3)
colnames(simu.data_M25)=c("V","T1.25","T2.25","T3.25")

#Estimated graph
n<-nrow (simu.data_M25)    #Number of row
V<-colnames(simu.data_M25) #Column names

Rcor_C=RobustCor(simu.data_M25, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M25, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_25<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_25_GT <- MRPC(simu.data_M25,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_25_CT <- MRPC(simu.data_M25,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_25,  main ="(A) Model 25 (by pc (alpha=0.05))")
plot(MRPC.fit_25_GT, main ="(B) Model 25 (by MRPC_GT)")
plot(MRPC.fit_25_CT, main ="(C) Model 25 (by MRPC_CT)")

#Model 26 (extention of model 3) based on model 1 (triplet=V,T1,T2)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T2=simu.data_M3[ ,3]
T3=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M26=cbind(V,T1,T2,T3)
colnames(simu.data_M26)=c("V","T1.26","T2.26","T3.26")

#Estimated graph
n<-nrow (simu.data_M26)    #Number of row
V<-colnames(simu.data_M26) #Column names

Rcor_C=RobustCor(simu.data_M26, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M26, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_26<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_26_GT <- MRPC(simu.data_M26,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_26_CT <- MRPC(simu.data_M26,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_26,  main ="(A) Model 26 (by pc (alpha=0.05))")
plot(MRPC.fit_26_GT, main ="(B) Model 26 (by MRPC_GT)")
plot(MRPC.fit_26_CT, main ="(C) Model 26 (by MRPC_CT)")

#Model 27 (extention of model 3) based on model 2 (triplet=V,T1,T2)
#Data

V=simu.data_M2[ ,1]
T1=simu.data_M2[ ,2]
T3=simu.data_M2[ ,3]
T2=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M27=cbind(V,T1,T2,T3)
colnames(simu.data_M27)=c("V","T1.27","T2.27","T3.27")

#Estimated graph
n<-nrow (simu.data_M27)    #Number of row
V<-colnames(simu.data_M27) #Column names

Rcor_C=RobustCor(simu.data_M27, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M27, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_27<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_27_GT <- MRPC(simu.data_M27,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_27_CT <- MRPC(simu.data_M27,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_27,  main ="(A) Model 27 (by pc (alpha=0.05))")
plot(MRPC.fit_27_GT, main ="(B) Model 27 (by MRPC_GT)")
plot(MRPC.fit_27_CT, main ="(C) Model 27 (by MRPC_CT)")

#Model 28 (extention of model 3) based on model 3 (triplet=V,T1,T2)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T3=simu.data_M3[ ,3]
T2=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M28=cbind(V,T1,T2,T3)
colnames(simu.data_M28)=c("V","T1.28","T2.28","T3.28")

#Estimated graph
n<-nrow (simu.data_M28)    #Number of row
V<-colnames(simu.data_M28) #Column names

Rcor_C=RobustCor(simu.data_M28, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M28, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_28<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_28_GT <- MRPC(simu.data_M28,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_28_CT <- MRPC(simu.data_M28,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_28,  main ="(A) Model 28 (by pc (alpha=0.05))")
plot(MRPC.fit_28_GT, main ="(B) Model 28 (by MRPC_GT)")
plot(MRPC.fit_28_CT, main ="(C) Model 28 (by MRPC_CT)")

#Model 29 (extention of model 3) based on model 4 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T3=simu.data_M4[ ,3]
T2=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M29=cbind(V,T1,T2,T3)
colnames(simu.data_M29)=c("V","T1.29","T2.29","T3.29")

#Estimated graph
n<-nrow (simu.data_M29)    #Number of row
V<-colnames(simu.data_M29) #Column names

Rcor_C=RobustCor(simu.data_M29, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M29, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_29<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_29_GT <- MRPC(simu.data_M29,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_29_CT <- MRPC(simu.data_M29,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_29,  main ="(A) Model 29 (by pc (alpha=0.05))")
plot(MRPC.fit_29_GT, main ="(B) Model 29 (by MRPC_GT)")
plot(MRPC.fit_29_CT, main ="(C) Model 29 (by MRPC_CT)")

#Model 30 (extention of model 3) based on model 1 (triplet=V,T1,T3)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T2=simu.data_M3[ ,3]
T3=Case_1P(N=10^3,P1=T2,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M30=cbind(V,T1,T2,T3)
colnames(simu.data_M30)=c("V","T1.30","T2.30","T3.30")

#Estimated graph
n<-nrow (simu.data_M30)    #Number of row
V<-colnames(simu.data_M30) #Column names

Rcor_C=RobustCor(simu.data_M30, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M30, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_30<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_30_GT <- MRPC(simu.data_M30,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_30_CT <- MRPC(simu.data_M30,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_30,  main ="(A) Model 30 (by pc (alpha=0.05))")
plot(MRPC.fit_30_GT, main ="(B) Model 30 (by MRPC_GT)")
plot(MRPC.fit_30_CT, main ="(C) Model 30 (by MRPC_CT)")

#Model 31 (extention of model 3) based on model 2 (triplet=V,T1,T3)
#Data

V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T3=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
T2=Case_2P(N=10^3,P1=V,P2=T3,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
simu.data_M31=cbind(V,T1,T2,T3)
colnames(simu.data_M31)=c("V","T1.31","T2.31","T3.31")

#Estimated graph
n<-nrow (simu.data_M31)    #Number of row
V<-colnames(simu.data_M31) #Column names

Rcor_C=RobustCor(simu.data_M31, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M31, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_31<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_31_GT <- MRPC(simu.data_M31,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_31_CT <- MRPC(simu.data_M31,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_31,  main ="(A) Model 31 (by pc (alpha=0.05))")
plot(MRPC.fit_31_GT, main ="(B) Model 31 (by MRPC_GT)")
plot(MRPC.fit_31_CT, main ="(C) Model 31 (by MRPC_CT)")

#Model 32 (extention of model 3) based on model 4 (triplet=V,T1,T3)
#Data

V=simu.data_M4[ ,1]
T2=simu.data_M4[ ,2]
T3=simu.data_M4[ ,3]
T1=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
simu.data_M32=cbind(V,T1,T2,T3)
colnames(simu.data_M32)=c("V","T1.32","T2.32","T3.32")

#Estimated graph
n<-nrow (simu.data_M32)    #Number of row
V<-colnames(simu.data_M32) #Column names

Rcor_C=RobustCor(simu.data_M32, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M32, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_32<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_32_GT <- MRPC(simu.data_M32,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_32_CT <- MRPC(simu.data_M32,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_32,  main ="(A) Model 32 (by pc (alpha=0.05))")
plot(MRPC.fit_32_GT, main ="(B) Model 32 (by MRPC_GT)")
plot(MRPC.fit_32_CT, main ="(C) Model 32 (by MRPC_CT)")

#Model 33 (extention of model 3) based on model 2 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]
T3=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M33=cbind(V,T1,T2,T3)
colnames(simu.data_M33)=c("V","T1.33","T2.33","T3.33")

#Estimated graph
n<-nrow (simu.data_M33)    #Number of row
V<-colnames(simu.data_M33) #Column names

Rcor_C=RobustCor(simu.data_M33, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M33, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_33<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_33_GT <- MRPC(simu.data_M33,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_33_CT <- MRPC(simu.data_M33,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_33,  main ="(A) Model 33 (by pc (alpha=0.05))")
plot(MRPC.fit_33_GT, main ="(B) Model 33 (by MRPC_GT)")
plot(MRPC.fit_33_CT, main ="(C) Model 33 (by MRPC_CT)")

#Model 34 (extention of model 4) based on model 1 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]
T3=Case_1P(N=10^3,P1=T1,b0.1 = 0,b1.1=1,sd.1 = 1)
simu.data_M34=cbind(V,T1,T2,T3)
colnames(simu.data_M34)=c("V","T1.34","T2.34","T3.34")

#Estimated graph
n<-nrow (simu.data_M34)    #Number of row
V<-colnames(simu.data_M34) #Column names

Rcor_C=RobustCor(simu.data_M34, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M34, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_34<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_34_GT <- MRPC(simu.data_M34,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_34_CT <- MRPC(simu.data_M34,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_34,  main ="(A) Model 34 (by pc (alpha=0.05))")
plot(MRPC.fit_34_GT, main ="(B) Model 34 (by MRPC_GT)")
plot(MRPC.fit_34_CT, main ="(C) Model 34 (by MRPC_CT)")

#Model 35 (extention of model 4) based on model 2 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T3=simu.data_M2[ ,3]

T1.a <- Case_2P(N=10^3,P1=V,P2=T3,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
T2.a <- Case_2P(N=10^3,P1=V,P2=T1.a,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

T2.b <- Case_1P(N=10^3,P1=V,b0.1=0,b1.1=1,sd.1=1)
T1.b <- Case_3P(N=10^3,P1=V,P2=T3,P3=T2.b,b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)


coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T1 <- rep(0, 10^3)
T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]

T2 <- rep(0, 10^3)
T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]


simu.data_M35=cbind(V,T1,T2,T3)
colnames(simu.data_M35)=c("V","T1.35","T2.35","T3.35")

#Estimated graph
n<-nrow (simu.data_M35)    #Number of row
V<-colnames(simu.data_M35) #Column names

Rcor_C=RobustCor(simu.data_M35, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M35, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_35<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_35_GT<- MRPC(simu.data_M35,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_35_CT<- MRPC(simu.data_M35,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_35,  main ="(A) Model 35 (by pc (alpha=0.05))")
plot(MRPC.fit_35_GT, main ="(B) Model 35 (by MRPC_GT)")
plot(MRPC.fit_35_CT, main ="(C) Model 35 (by MRPC_CT)")

#Model 36 (extention of model 4) based on model 3 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]
T3=Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1=1,sd.1 = 1)
simu.data_M36=cbind(V,T1,T2,T3)
colnames(simu.data_M36)=c("V","T1.36","T2.36","T3.36")

#Estimated graph
n<-nrow (simu.data_M36)    #Number of row
V<-colnames(simu.data_M36) #Column names

Rcor_C=RobustCor(simu.data_M36, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M36, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_36<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_36_GT <- MRPC(simu.data_M36,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_36_CT <- MRPC(simu.data_M36,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_36,  main ="(A) Model 36 (by pc (alpha=0.05))")
plot(MRPC.fit_36_GT, main ="(B) Model 36 (by MRPC_GT)")
plot(MRPC.fit_36_CT, main ="(C) Model 36 (by MRPC_CT)")

#Model 37 (extention of model 4) based on model 4 (triplet=V,T1,T2)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]

T3.a <-Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T3.b=Case_2P(N=10^3,P1=V,P2=T1,b0.1 = 0,b1.1 = 1,1,sd.1 = 1)

T3 <- rep(0, 10^3)
T3[which(coinToss == 0)] <- T3.a[which(coinToss == 0)]
T3[which(coinToss == 1)] <- T3.b[which(coinToss == 1)]

simu.data_M37=cbind(V,T1,T2,T3)
colnames(simu.data_M37)=c("V","T1.37","T2.37","T3.37")

#Estimated graph
n<-nrow (simu.data_M37)    #Number of row
V<-colnames(simu.data_M37) #Column names

Rcor_C=RobustCor(simu.data_M37, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M37, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_37<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_37_GT <- MRPC(simu.data_M37,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_37_CT <- MRPC(simu.data_M37,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_37,  main ="(A) Model 37 (by pc (alpha=0.05))")
plot(MRPC.fit_37_GT, main ="(B) Model 37 (by MRPC_GT)")
plot(MRPC.fit_37_CT, main ="(C) Model 37 (by MRPC_CT)")

#Model 38 (extention of model 4) based on model 1 (triplet=V,T1,T3)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]
T3=Case_1P(N=10^3,P1=T2,b0.1 = 0,b1.1=1,sd.1 = 1)
simu.data_M38=cbind(V,T1,T2,T3)
colnames(simu.data_M38)=c("V","T1.38","T2.38","T3.38")

#Estimated graph
n<-nrow (simu.data_M38)    #Number of row
V<-colnames(simu.data_M38) #Column names

Rcor_C=RobustCor(simu.data_M38, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M38, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_38<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_38_GT <- MRPC(simu.data_M38,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_38_CT <- MRPC(simu.data_M38,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_38,  main ="(A) Model 38 (by pc (alpha=0.05))")
plot(MRPC.fit_38_GT, main ="(B) Model 38 (by MRPC_GT)")
plot(MRPC.fit_38_CT, main ="(C) Model 38 (by MRPC_CT)")

#Model 39 (extention of model 4) based on model 2 (triplet=V,T1,T3)
#Data

V=simu.data_M4[ ,1]
T3=simu.data_M2[ ,3]

T1.a <- Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T2.a <- Case_3P(N=10^3,P1=V,P2=T1.a,P3=T3,b0.1 = 0,b1.1 = 1,b1.2 = 1,b1.3 = 1,sd.1 = 1)

T2.b <- Case_2P(N=10^3,P1=V,P2=T3,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)
T1.b <- Case_2P(N=10^3,P1=V,P2=T2.b,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T1 <- rep(0, 10^3)
T1[which(coinToss == 0)] <- T1.a[which(coinToss == 0)]
T1[which(coinToss == 1)] <- T1.b[which(coinToss == 1)]

T2 <- rep(0, 10^3)
T2[which(coinToss == 0)] <- T2.a[which(coinToss == 0)]
T2[which(coinToss == 1)] <- T2.b[which(coinToss == 1)]

simu.data_M39=cbind(V,T1,T2,T3)
colnames(simu.data_M39)=c("V","T1.39","T2.39","T3.39")

#Estimated graph
n<-nrow (simu.data_M39)    #Number of row
V<-colnames(simu.data_M39) #Column names

Rcor_C=RobustCor(simu.data_M39, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M39, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_39<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_39_GT <- MRPC(simu.data_M39,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_39_CT <- MRPC(simu.data_M39,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_39,  main ="(A) Model 39 (by pc (alpha=0.05))")
plot(MRPC.fit_39_GT, main ="(B) Model 39 (by MRPC_GT)")
plot(MRPC.fit_39_CT, main ="(C) Model 39 (by MRPC_CT)")

#Model 40 (extention of model 4) based on model 4 (triplet=V,T1,T3)
#Data

V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T2=simu.data_M4[ ,3]

T3.a <-Case_1P(N=10^3,P1=V,b0.1 = 0,b1.1 = 1,sd.1 = 1)
T3.b=Case_2P(N=10^3,P1=V,P2=T2,b0.1 = 0,b1.1 = 1,b1.2 = 1,sd.1 = 1)

coinToss <- rbinom(n = 10^3,
                   size = 1,
                   prob = 0.5)

T3 <- rep(0, 10^3)
T3[which(coinToss == 0)] <- T3.a[which(coinToss == 0)]
T3[which(coinToss == 1)] <- T3.b[which(coinToss == 1)]

simu.data_M40=cbind(V,T1,T2,T3)
colnames(simu.data_M40)=c("V","T1.40","T2.40","T3.40")

#Estimated graph
n<-nrow (simu.data_M40)    #Number of row
V<-colnames(simu.data_M40) #Column names

Rcor_C=RobustCor(simu.data_M40, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_M40, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_40<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_40_GT <- MRPC(simu.data_M40,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_40_CT <- MRPC(simu.data_M40,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_40,  main ="(A) Model 40 (by pc (alpha=0.05))")
plot(MRPC.fit_40_GT, main ="(B) Model 40 (by MRPC_GT)")
plot(MRPC.fit_40_CT, main ="(C) Model 40 (by MRPC_CT)")



#To examine the robustness (Model 0-Model 4) of the MRPC method in
#comparision of the PC method,we added some artificial outliers to
#our data matrix randomly.First, we generated uniform distribution
#for outlying expression and add with the original data.

## Model 0
Outdata=matrix(c(runif(n = 30, min = 80, max = 100)),ncol = 3)
simu.data_M0[c(991:1000),c(2,3)]=Outdata[,1] #Adding outliers
simu.outdata_M0=simu.data_M0 #New data with outliers

n<-nrow (simu.outdata_M0)    #Number of row
V<-colnames(simu.outdata_M0) #Column names

Rcor_C=RobustCor(simu.outdata_M0, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.outdata_M0, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by existing PC.
pc.fit_0<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by classical correlation
MRPC.fit_0_C_GT <- MRPC(simu.outdata_M0,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test by classical correlation
MRPC.fit_0_C_CT <- MRPC(simu.outdata_M0,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'citest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by robust correlation
MRPC.fit_0_R_GT <- MRPC(simu.outdata_M0,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_0,  main ="(A) Model 0(by pc (alpha=0.05))")
plot(MRPC.fit_0_C_GT,  main ="(B) Model 0(MRPC_GT(Classical Correlation))")
plot(MRPC.fit_0_C_CT, main ="(C) Model 0(by MRPC_CT)")
plot(MRPC.fit_0_R_GT, main ="(D) Model 0(by MRPC_GT(Robust Correlation))")


## Model 1
Outdata=matrix(c(runif(n = 15, min = 80, max = 100)),ncol = 3)
simu.data_M1[c(996:1000),2]=Outdata[,1] #Adding outliers
simu.outdata_M1=simu.data_M1 #New data with outliers

n<-nrow (simu.outdata_M1)    #Number of row
V<-colnames(simu.outdata_M1) #Column names

Rcor_C=RobustCor(simu.outdata_M1, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.outdata_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)


## Estimated graph by existing PC.
pc.fit_1<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by classical correlation
MRPC.fit_1_C_GT <- MRPC(simu.outdata_M1,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test by classical correlation
MRPC.fit_1_C_CT <- MRPC(simu.outdata_M1,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'citest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by robust correlation
MRPC.fit_1_R_GT <- MRPC(simu.outdata_M1,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_1,  main ="(A) Model 1(by pc (alpha=0.05))")
plot(MRPC.fit_1_C_GT,  main ="(B) Model 1(MRPC_GT(Classical Correlation))")
plot(MRPC.fit_1_C_CT, main ="(C) Model 1(by MRPC_CT)")
plot(MRPC.fit_1_R_GT, main ="(D) Model 1(by MRPC_GT(Robust Correlation))")


## Model 2
Outdata=matrix(c(runif(n = 15, min = 80, max = 100)),ncol = 3)
simu.data_M2[c(996:1000),2]=Outdata[,1] #Adding outliers
simu.outdata_M2=simu.data_M2 #New data with outliers

n<-nrow (simu.outdata_M2)    #Number of row
V<-colnames(simu.outdata_M2) #Column names

Rcor_C=RobustCor(simu.outdata_M2, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.outdata_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by existing PC.
pc.fit_2<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by classical correlation
MRPC.fit_2_C_GT <- MRPC(simu.outdata_M2,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test by classical correlation
MRPC.fit_2_C_CT <- MRPC(simu.outdata_M2,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'citest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by robust correlation
MRPC.fit_2_R_GT <- MRPC(simu.outdata_M2,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_2,  main ="(A) Model 0(by pc (alpha=0.05))")
plot(MRPC.fit_2_C_GT,  main ="(B) Model 0(MRPC_GT(Classical Correlation))")
plot(MRPC.fit_2_C_CT, main ="(C) Model 0(by MRPC_CT)")
plot(MRPC.fit_2_R_GT, main ="(D) Model 0(by MRPC_GT(Robust Correlation))")


## Model 3
Outdata=matrix(c(runif(n = 15, min = 80, max = 100)),ncol = 3)
simu.data_M3[c(996:1000),2]=Outdata[,1] #Adding outliers
simu.outdata_M3=simu.data_M3 #New data with outliers

n<-nrow (simu.outdata_M3)    #Number of row
V<-colnames(simu.outdata_M3) #Column names

Rcor_C=RobustCor(simu.outdata_M3, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.outdata_M3, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by existing PC.
pc.fit_3<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by classical correlation
MRPC.fit_3_C_GT <- MRPC(simu.outdata_M3,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test by classical correlation
MRPC.fit_3_C_CT <- MRPC(simu.outdata_M3,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'citest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by robust correlation
MRPC.fit_3_R_GT <- MRPC(simu.outdata_M3,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_3,  main ="(A) Model 0(by pc (alpha=0.05))")
plot(MRPC.fit_3_C_GT,  main ="(B) Model 0(MRPC_GT(Classical Correlation))")
plot(MRPC.fit_3_C_CT, main ="(C) Model 0(by MRPC_CT)")
plot(MRPC.fit_3_R_GT, main ="(D) Model 0(by MRPC_GT(Robust Correlation))")

## Model 4
Outdata=matrix(c(runif(n = 300, min = 80, max = 100)),ncol = 3)
simu.data_M4[c(901:1000),2]=Outdata[,1] #Adding outliers
simu.outdata_M4=simu.data_M4 #New data with outliers

n<-nrow (simu.outdata_M4)    #Number of row
V<-colnames(simu.outdata_M4) #Column names

Rcor_C=RobustCor(simu.outdata_M4, 0)  #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.outdata_M4, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by existing PC.
pc.fit_4<- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by classical correlation
MRPC.fit_4_C_GT <- MRPC(simu.outdata_M3,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test by classical correlation
MRPC.fit_4_C_CT <- MRPC(simu.outdata_M3,suffStat_C,NQ=1,FDR=0.05
,indepTest = 'citest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using GaussCItest by robust correlation
MRPC.fit_4_R_GT <- MRPC(simu.outdata_M3,suffStat_R,NQ=1,FDR=0.05
,indepTest = 'gaussCItest',labels = V,verbose = TRUE)

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(pc.fit_4,  main ="(A) Model 4(by pc (alpha=0.05))")
plot(MRPC.fit_4_C_GT,  main ="(B) Model 4(MRPC_GT(Classical Correlation))")
plot(MRPC.fit_4_C_CT, main ="(C) Model 4(by MRPC_CT)")
plot(MRPC.fit_4_R_GT, main ="(D) Model 4(by MRPC_GT(Robust Correlation))")

#Example from pclag
data(gmG)
n<-nrow (gmG$x)    #Number of row
#V<-ncol(gmG$ x) #Column names
V=colnames(gmG$ x)=c("Author","Bar","Ctrl","Goal","V5","V6","V7","V8")

Rcor_C=RobustCor(gmG$ x, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(gmG$ x, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)



### Estimated graph by pc
pc.fit_Ex <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using gaussCItest
MRPC.fit_Ex_GT <- MRPC(gmG$ x,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)

## Estimated graph by MRPC using ci.test
MRPC.fit_Ex_CT <- MRPC(gmG$ x,suffStat_R,NQ=1,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_Ex,  main ="(A) Model Ex(by pc (alpha=0.05))")
plot(MRPC.fit_Ex_GT, main ="(B) Model Ex(by MRPC_GT)")
plot(MRPC.fit_Ex_CT, main ="(C) Model Ex(by MRPC_CT)")

##Data for the multiple SNPs/indels/CNV/eQTL
V1 <- c(sample(c(0, 1, 2),
               size = 10^3,
               replace = TRUE,
               prob = c((1 - 0.45)^2,
                        2*0.45*(1 - 0.45),
                        0.45^2)))

V2<- Case_1P(N=10^3,P1=V1,b0.1 = 1,b1.1=1,sd.1 = 1) #2nd eQTL
T1<- Case_1P(N=10^3,P1=V1,b0.1 = 1,b1.1=1,sd.1 = 1)
T2<- Case_1P(N=10^3,P1=V2,b0.1 = 1,b1.1=1,sd.1 = 1)

T3<- Case_2P(N=10^3,P1=T2,P2=T1,b0.1 = 1,b1.1=1,b1.2=1,sd.1 = 1)
T4<- Case_2P(N=10^3,P1=T1,P2=T2,b0.1 = 1,b1.1=1,b1.2=1,sd.1 = 1)

simu.data_MQ=data.frame(V1=V1,
                V2=V2,
                T1= T1,
                T2 = T2,
                T3 = T3,
                T4 = T4)

n<-nrow (simu.data_MQ)    #Number of row
V<-colnames(simu.data_MQ) #Column names

Rcor_C=RobustCor(simu.data_MQ, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(simu.data_MQ, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)
### Estimated graph by pc
pc.fit_MQ <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using gaussCItest
MRPC.fit_MQ_GT <- MRPC(simu.data_MQ,suffStat_R,NQ=2,FDR=0.05,
indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_MQ_CT<- MRPC(simu.data_MQ,suffStat_R,NQ=2,FDR=0.05,
indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_MQ ,  main ="(A) Multiple eQTL topology (by pc (alpha=0.05))")
plot(MRPC.fit_MQ_GT, main ="(B) Multiple eQTL topology (by MRPC_GT)")
plot(MRPC.fit_MQ_CT, main ="(C) Multiple eQTL topology (by MRPC_CT)")


#Three new models are not including above

##v-->T1 and T2, T1-->T3, T2-->T3

N=1000
p=.45
b0.1=0
b1.1=1
b1.2=1
sd.1=1

V <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))

T1<- Case_1P(N=N,P1=V,b0.1 = b0.1,b1.1 = b1.1,sd.1 = sd.1)
T2<- Case_1P(N=N,P1=V,b0.1 = b0.1,b1.1 = b1.1,sd.1 = sd.1)
T3<- Case_2P(N=N,P1=T1,P2=T2,b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,sd.1 = sd.1)

Sim_New1=data.frame(V=V,
                    T1= T1,
                    T2 = T2,
                    T3=T3)

n<-nrow (Sim_New1)    #Number of row
V<-colnames(Sim_New1) #Column names

Rcor_C=RobustCor(Sim_New1, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(Sim_New1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

### Estimated graph by pc
pc.fit_New1 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using gaussCItest
MRPC.fit_New1_GT <- MRPC(Sim_New1,suffStat_R,NQ=1,FDR=0.05,
                       indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_New1_CT<- MRPC(Sim_New1,suffStat_R,NQ=1,FDR=0.05,
                      indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_New1,  main ="(A) New model 1 (by pc (alpha=0.05))")
plot(MRPC.fit_New1_GT, main ="(B) New model 1 (by MRPC_GT)")
plot(MRPC.fit_New1_CT, main ="(C) New model 1 (by MRPC_CT)")

##v-->T1 and T2, T1-->T3-->T3

N=1000
p=.45
b0.1=0
b1.1=1
b1.2=1
sd.1=1

V <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))

T1<- Case_1P(N=N,P1=V,b0.1 = b0.1,b1.1 = b1.1,sd.1 = sd.1)
T3<- Case_1P(N=N,P1=T1,b0.1 = b0.1,b1.1 = b1.1,sd.1 = sd.1)
T2<- Case_2P(N=N,P1=V,P2=T3,b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,sd.1 = sd.1)


Sim_New2=data.frame(V=V,
                    T1= T1,
                    T2 = T2,
                    T3=T3)

n<-nrow (Sim_New2)    #Number of row
V<-colnames(Sim_New2) #Column names

Rcor_C=RobustCor(Sim_New2, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(Sim_New2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)
### Estimated graph by pc
pc.fit_New2 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using gaussCItest
MRPC.fit_New2_GT <- MRPC(Sim_New2,suffStat_R,NQ=1,FDR=0.05,
                         indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_New2_CT<- MRPC(Sim_New2,suffStat_R,NQ=1,FDR=0.05,
                        indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_New2,  main ="(A) New model 2 (by pc (alpha=0.05))")
plot(MRPC.fit_New2_GT, main ="(B) New model 2 (by MRPC_GT)")
plot(MRPC.fit_New2_CT, main ="(C) New model 2 (by MRPC_CT)")


##v-->T1 and T2, T2-->T1,T2

N=1000
p=.45
b0.1=0
b1.1=1
b1.2=1
sd.1=1

V <- c(sample(c(0, 1, 2),size = N,replace = TRUE,prob = c((1 - p)^2,2*p*(1 - p),p^2)))

T3<- Case_NP(N=N,b0.1 = b0.1,sd.1 = sd.1)
T1<- Case_2P(N=N,P1=V,P2=T3,b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,sd.1 = sd.1)
T2<- Case_2P(N=N,P1=V,P2=T3,b0.1 = b0.1,b1.1 = b1.1,b1.2 = b1.2,sd.1 = sd.1)


Sim_New3=data.frame(V=V,
                    T1= T1,
                    T2 = T2,
                    T3=T3)

n<-nrow (Sim_New3)    #Number of row
V<-colnames(Sim_New3) #Column names

Rcor_C=RobustCor(Sim_New3, 0) #Classical correlation (Beta=0)
suffStat_C= list(C = Rcor_C$RR, n = n)

Rcor_R=RobustCor(Sim_New3, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)
### Estimated graph by pc
pc.fit_New3 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels = V,verbose = TRUE)
## Estimated graph by MRPC using gaussCItest
MRPC.fit_New3_GT <- MRPC(Sim_New3,suffStat_R,NQ=1,FDR=0.05,
                         indepTest ='gaussCItest',labels = V,verbose = TRUE)
## Estimated graph by MRPC using ci.test
MRPC.fit_New3_CT<- MRPC(Sim_New3,suffStat_R,NQ=1,FDR=0.05,
                        indepTest ='citest',labels = V,verbose = TRUE)

## Show estimated graph
par(mfrow=c(1,3))
plot(pc.fit_New3, main ="(A) New model 3 (by pc (alpha=0.05))")
plot(MRPC.fit_New3_GT, main ="(B) New model 3 (by MRPC_GT)")
plot(MRPC.fit_New3_CT, main ="(C) New model 3 (by MRPC_CT)")

}
