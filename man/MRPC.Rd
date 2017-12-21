\name{MRPC}
\alias{MRPC}
\title{
Estimate the DAG using Generalized Mendelian Randomization (gMR) based PC (MRPC) Algorithm
}
\description{
The MRPC () function is used to estimate a directed acyclic graph (DAG) from observational data, by modifying of the existing \link{pc} algorithm (Spirtes et al., 2000). Existing \link{pc} algorithms is used for estimating the DAG from observational data and several variants are implemented in the R package called pcalg. Several challenges remain when \link{pc} algorithms are applied to noisy genomic data. To tackle these challenges, we develop a novel machine learning algorithm called Mendelian Randomization (MR) based PC (MRPC) to efficiently learn a causal gene regulatory network (or a causal graph of genes) using genotype and gene expression data. Unlike existing \link{pc} algorithms, the four-major updated in the MRPC algorithm is that (i) we are incorporating sequential testing method to controls the False Discovery Rate (FDR), (ii) take a robust approach to reduces the impact of outliers, (iii)improved v-structure identification and (iv) implemented a new way for edge direction determination based on gMR. See details below.

}
\usage{
MRPC(data,suffStat,GV,FDR,indepTest, labels, p,
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
\itemize{
\item \code{obj}: object of class "graph" of the directed and undirected graph that was estimated.
\item \code{test}: current test number
\item \code{alpha}: level of significance that updated in every step of sequential FDR test.
\item \code{R}: rejection number.
}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\references{
1.Badsha, M. B. Mollah, M. N. Jahan, N. and Kurata, H. (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Javanmard and Montanari (March 5, 2015) On Online Control of False Discovery Rate. arXiv:150206197 [statME].

3. Kalisch, M. and Buhlmann, P. (2007) Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

4.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

5. Spirtes,P.  Glymour, C. and Scheines, R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.

6.Kalisch, M., Machler, M., Colombo, D., Maathuis, M.H. & Buhlmann, P. Causal Inference Using Graphical Models with the R Package pcalg. J. Stat. Softw. 47, 26 (2012).

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

## Model 0
n<-nrow (simu.data_M0)    #Number of row
V<-colnames(simu.data_M0) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_M0), n = n)

#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M0, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_0<- MRPC(simu.data_M0,suffStat_R,GV=1,FDR=0.05,
indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_0<- mmhc(simu.data_M0)

## Estimated graph by existing pc with alpha=0.05
pc.fit_0 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for Model 0
tarmat=matrix(0,nrow=ncol(simu.data_M0),ncol = ncol(simu.data_M0))
colnames(tarmat)=colnames(simu.data_M0)
rownames(tarmat)=colnames(simu.data_M0)
tarmat[1,2]=1
Truth_0<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_0,main ="(A) Turth Model 0")
plot(MRPC.fit_0$graph, main ="(B) MRPC")
graphviz.plot(mmhc_0, main ="(C) mmhc")
plot(pc.fit_0,  main ="(D) pc")


#Another option for plot of the results. First fig is the nodes
#dendrogram with modules color. second fig is the plot of the graph
#with color based on modules.
#To idendify modules and complex graph (Suitable if you have many
#nodes)
#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_0$graph,"matrix")
#Plot of the graph with modules color
DendroModuleGraph(Adj_directed,minModuleSize = 1,GV=1)


## Model 1
n<-nrow (simu.data_M1)    #Number of row
V<-colnames(simu.data_M1) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_M1), n = n)

#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M1, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_1 <- MRPC(simu.data_M1,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_1<- mmhc(simu.data_M1)

## Estimated graph by existing pc with alpha=0.05
pc.fit_1 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for Model 1
tarmat=matrix(0,nrow=ncol(simu.data_M1),ncol = ncol(simu.data_M1))
colnames(tarmat)=colnames(simu.data_M1)
rownames(tarmat)=colnames(simu.data_M1)
tarmat[1,2]=1
tarmat[2,3]=1
Truth_1<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_1,main ="(A) Turth Model 1")
plot(MRPC.fit_1$graph, main ="(B) MRPC")
graphviz.plot(mmhc_1, main ="(C) mmhc")
plot(pc.fit_1,main ="(D) pc")

## Model 2
n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M2) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_M2), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M2, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_2 <- MRPC(simu.data_M2,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_2<- mmhc(simu.data_M2)

## Estimated graph by existing pc with alpha=0.05
pc.fit_2 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for Model 2
tarmat=matrix(0,nrow=ncol(simu.data_M2),ncol = ncol(simu.data_M2))
colnames(tarmat)=colnames(simu.data_M2)
rownames(tarmat)=colnames(simu.data_M2)
tarmat[1,2]=1
tarmat[3,2]=1
Truth_2<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_2,main ="(A) Turth Model 2")
plot(MRPC.fit_2$graph, main ="(B) MRPC")
graphviz.plot(mmhc_2, main ="(C) mmhc")
plot(pc.fit_2,  main ="(D) pc")

## Model 3
n<-nrow (simu.data_M3)    #Number of row
V<-colnames(simu.data_M3) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_M3), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M3, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_3<- MRPC(simu.data_M3,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_3<- mmhc(simu.data_M3)

## Estimated graph by existing pc with alpha=0.05
pc.fit_3 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for Model 3
tarmat=matrix(0,nrow=ncol(simu.data_M3),ncol = ncol(simu.data_M3))
colnames(tarmat)=colnames(simu.data_M3)
rownames(tarmat)=colnames(simu.data_M3)
tarmat[1,2]=1
tarmat[1,3]=1

Truth_3<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_3,main ="(A) Turth Model 3")
plot(MRPC.fit_3$graph, main ="(B) MRPC")
graphviz.plot(mmhc_3, main ="(C) mmhc")
plot(pc.fit_3,  main ="(D) pc")

## Model 4
n<-nrow (simu.data_M4)    #Number of row
V<-colnames(simu.data_M4) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_M4), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M4, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_4<- MRPC(simu.data_M4,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_4<- mmhc(simu.data_M4)

## Estimated graph by existing pc with alpha=0.05
pc.fit_4 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for Model 4
tarmat=matrix(0,nrow=ncol(simu.data_M4),ncol = ncol(simu.data_M4))
colnames(tarmat)=colnames(simu.data_M4)
rownames(tarmat)=colnames(simu.data_M4)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1

Truth_4<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_4,main ="(A) Turth Model 4")
plot(MRPC.fit_4$graph, main ="(B) MRPC")
graphviz.plot(mmhc_4, main ="(C) mmhc")
plot(pc.fit_4,  main ="(D) pc")

## Multiple Parents Model
n<-nrow (simu.data_multiparent)    #Number of row
V<-colnames(simu.data_multiparent) #Column names

#Classical correlation
suffStat_C= list(C = cor(simu.data_multiparent), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_multiparent, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_Multiparents<- MRPC(simu.data_multiparent,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_Multiparents<- mmhc(simu.data_multiparent)

## Estimated graph by existing pc with alpha=0.05
pc.fit_Multiparents <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for multiparent
tarmat=matrix(0,nrow=ncol(simu.data_multiparent),ncol = ncol(simu.data_multiparent))
colnames(tarmat)=colnames(simu.data_multiparent)
rownames(tarmat)=colnames(simu.data_multiparent)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[4,2]=1

Truth_Multiparents<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_Multiparents,main ="(A) Turth Model Multiparents")
plot(MRPC.fit_Multiparents$graph, main ="(B) MRPC")
graphviz.plot(mmhc_Multiparents, main ="(C) mmhc")
plot(pc.fit_Multiparents,  main ="(D) pc")


## Star Model
n<-nrow (simu.data_starshaped)    #Number of row
V<-colnames(simu.data_starshaped) #Column names

#Classical correlation 
suffStat_C= list(C=cor(simu.data_starshaped), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_starshaped , 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_StarModel<- MRPC(simu.data_starshaped ,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_StarModel<- mmhc(simu.data_starshaped )

## Estimated graph by existing pc with alpha=0.05
pc.fit_StarModel <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth Star Model
tarmat=matrix(0,nrow=ncol(simu.data_starshaped),ncol = ncol(simu.data_starshaped))
colnames(tarmat)=colnames(simu.data_starshaped)
rownames(tarmat)=colnames(simu.data_starshaped)

tarmat[1,2]=1
tarmat[2,3]=1
tarmat[2,4]=1
tarmat[2,5]=1
tarmat[2,6]=1

Truth_starshaped<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_starshaped, main ="(A) Truth Star Model")
plot(MRPC.fit_StarModel$graph, main ="(B) MRPC")
graphviz.plot(mmhc_StarModel, main ="(C) mmhc")
plot(pc.fit_StarModel,  main ="(D) pc")

## Layer shape topology
n<-nrow (simu.data_layered)    #Number of row
V<-colnames(simu.data_layered) #Column names

#Classical correlation (Beta=0)
suffStat_C= list(C=cor(simu.data_layered), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_layered, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_Layered<- MRPC(simu.data_layered,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_Layered<- mmhc(simu.data_layered)

## Estimated graph by existing pc with alpha=0.05
pc.fit_Layered <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth layer shape topology
tarmat=matrix(0,nrow=ncol(simu.data_layered),ncol = ncol(simu.data_layered))
colnames(tarmat)=colnames(simu.data_layered)
rownames(tarmat)=colnames(simu.data_layered)

tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[2,5]=1
tarmat[2,6]=1
tarmat[3,6]=1
tarmat[3,7]=1
tarmat[4,8]=1

Truth_layered<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_layered, main ="(A) Truth Layered Model")
plot(MRPC.fit_Layered$graph, main ="(B) MRPC")
graphviz.plot(mmhc_Layered, main ="(C) mmhc")
plot(pc.fit_Layered,  main ="(D) pc")


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

## Model 5
n<-nrow (simu.data_M5)    #Number of row
V<-colnames(simu.data_M5) #Column names

#Classical correlation (Beta=0)
suffStat_C= list(C = cor(simu.data_M5), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M5, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_5<- MRPC(simu.data_M5,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_5<- mmhc(data.frame(simu.data_M5))

## Estimated graph by existing pc with alpha=0.05
pc.fit_5 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 5
tarmat=matrix(0,nrow=ncol(simu.data_M5),ncol = ncol(simu.data_M5))
colnames(tarmat)=colnames(simu.data_M5)
rownames(tarmat)=colnames(simu.data_M5)
tarmat[1,2]=1

Truth_5<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_5,main ="(A) Turth Model 5")
plot(MRPC.fit_5$graph, main ="(A) MRPC")
graphviz.plot(mmhc_5, main ="(C) mmhc")
plot(pc.fit_5,  main ="(D) pc")

#Model 6 (extention of model 0) based on model 1
#Data
V=simu.data_M1[ ,1];
T1=simu.data_M1[ ,2];
T3=simu.data_M1[ ,3];
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M6=cbind(V,T1,T2,T3)
colnames(simu.data_M6)=c("V","T1.6","T2.6","T3.6")

#Estimated graph
## Model 6
n<-nrow (simu.data_M6)    #Number of row
V<-colnames(simu.data_M6) #Column names

#Classical correlation (Beta=0)
suffStat_C= list(C = cor(simu.data_M6), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M6, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_6<- MRPC(simu.data_M6,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_6<- mmhc(data.frame(simu.data_M6))

## Estimated graph by existing pc with alpha=0.05
pc.fit_6 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 6
tarmat=matrix(0,nrow=ncol(simu.data_M6),ncol = ncol(simu.data_M6))
colnames(tarmat)=colnames(simu.data_M6)
rownames(tarmat)=colnames(simu.data_M6)
tarmat[1,2]=1
tarmat[2,4]=1

Truth_6<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_6,main ="(A) Turth Model 6")
plot(MRPC.fit_6$graph, main ="(A) MRPC")
graphviz.plot(mmhc_6, main ="(C) mmhc")
plot(pc.fit_6,  main ="(D) pc")

#Model 7 (extention of model 0) based on model 2
#Data
V=simu.data_M2[ ,1];
T1=simu.data_M2[ ,2];
T3=simu.data_M2[ ,3];
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M7=cbind(V,T1,T2,T3)
colnames(simu.data_M7)=c("V","T1.7","T2.7","T3.7")

#Estimated graph
## Model 7
n<-nrow (simu.data_M7)    #Number of row
V<-colnames(simu.data_M7) #Column names

#Classical correlation (Beta=0)
suffStat_C= list(C = cor(simu.data_M7), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M7, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_7<- MRPC(simu.data_M7,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_7<- mmhc(data.frame(simu.data_M7))

## Estimated graph by existing pc with alpha=0.05
pc.fit_7 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 7
tarmat=matrix(0,nrow=ncol(simu.data_M7),ncol = ncol(simu.data_M7))
colnames(tarmat)=colnames(simu.data_M7)
rownames(tarmat)=colnames(simu.data_M7)
tarmat[1,2]=1
tarmat[4,2]=1

Truth_7<-as(tarmat, "graphNEL")
## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_7,main ="(A) Turth Model 7")
plot(MRPC.fit_7$graph, main ="(A) MRPC")
graphviz.plot(mmhc_7, main ="(C) mmhc")
plot(pc.fit_7,  main ="(D) pc")

#Model 8 (extention of model 0) based on model 3
#Data
V=simu.data_M3[ ,1]
T1=simu.data_M3[ ,2]
T3=simu.data_M3[ ,3]
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M8=cbind(V,T1,T2,T3)
colnames(simu.data_M8)=c("V","T1.8","T2.8","T3.8")

#Estimated graph
## Model 8
n<-nrow (simu.data_M8)    #Number of row
V<-colnames(simu.data_M8) #Column names

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M8), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M8, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_8<- MRPC(simu.data_M8,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_8<- mmhc(data.frame(simu.data_M8))

## Estimated graph by existing pc with alpha=0.05
pc.fit_8 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 8
tarmat=matrix(0,nrow=ncol(simu.data_M8),ncol = ncol(simu.data_M8))
colnames(tarmat)=colnames(simu.data_M8)
rownames(tarmat)=colnames(simu.data_M8)
tarmat[1,2]=1
tarmat[1,4]=1

Truth_8<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_8,main ="(A) Turth Model 8")
plot(MRPC.fit_8$graph, main ="(A) MRPC")
graphviz.plot(mmhc_8, main ="(C) mmhc")
plot(pc.fit_8,  main ="(D) pc")

#Model 9 (extention of model 0) based on model 4
#Data
V=simu.data_M4[ ,1]
T1=simu.data_M4[ ,2]
T3=simu.data_M4[ ,3]
T2=Case_NP(N=10^3,b0.1 = 0,sd.1 = 1)
simu.data_M9=cbind(V,T1,T2,T3)
colnames(simu.data_M9)=c("V","T1.9","T2.9","T3.9")

#Estimated graph
## Model 9
n<-nrow (simu.data_M9)    #Number of row
V<-colnames(simu.data_M9) #Column names

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M9), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M9, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_9<- MRPC(simu.data_M9,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_9<- mmhc(data.frame(simu.data_M9))

## Estimated graph by existing pc with alpha=0.05
pc.fit_9 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 9
tarmat=matrix(0,nrow=ncol(simu.data_M9),ncol = ncol(simu.data_M9))
colnames(tarmat)=colnames(simu.data_M9)
rownames(tarmat)=colnames(simu.data_M9)
tarmat[1,2]=1
tarmat[1,4]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_9<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_9,main ="(A) Turth Model 9")
plot(MRPC.fit_9$graph, main ="(A) MRPC")
graphviz.plot(mmhc_9, main ="(C) mmhc")
plot(pc.fit_9,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M10), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M10, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_10<- MRPC(simu.data_M10,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_10<- mmhc(data.frame(simu.data_M10))

## Estimated graph by existing pc with alpha=0.05
pc.fit_10 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 10
tarmat=matrix(0,nrow=ncol(simu.data_M10),ncol = ncol(simu.data_M10))
colnames(tarmat)=colnames(simu.data_M10)
rownames(tarmat)=colnames(simu.data_M10)
tarmat[1,2]=1
tarmat[2,3]=1

Truth_10<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_10,main ="(A) Turth Model 10")
plot(MRPC.fit_10$graph, main ="(B) MRPC")
graphviz.plot(mmhc_10, main ="(C) mmhc")
plot(pc.fit_10,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M11), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M11, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_11<- MRPC(simu.data_M11,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_11<- mmhc(data.frame(simu.data_M11))

## Estimated graph by existing pc with alpha=0.05
pc.fit_11 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 11
tarmat=matrix(0,nrow=ncol(simu.data_M11),ncol = ncol(simu.data_M11))
colnames(tarmat)=colnames(simu.data_M11)
rownames(tarmat)=colnames(simu.data_M11)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[2,4]=1
Truth_11<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_11,main ="(A) Turth Model 11")
plot(MRPC.fit_11$graph, main ="(B) MRPC")
graphviz.plot(mmhc_11, main ="(C) mmhc")
plot(pc.fit_11,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M12), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M12, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_12<- MRPC(simu.data_M12,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_12<- mmhc(data.frame(simu.data_M12))

## Estimated graph by existing pc with alpha=0.05
pc.fit_12 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 12
tarmat=matrix(0,nrow=ncol(simu.data_M12),ncol = ncol(simu.data_M12))
colnames(tarmat)=colnames(simu.data_M12)
rownames(tarmat)=colnames(simu.data_M12)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[4,2]=1
Truth_12<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_12,main ="(A) Turth Model 12")
plot(MRPC.fit_12$graph, main ="(A) MRPC")
graphviz.plot(mmhc_12, main ="(C) mmhc")
plot(pc.fit_12,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M13), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M13, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_13<- MRPC(simu.data_M13,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_13<- mmhc(data.frame(simu.data_M13))

## Estimated graph by existing pc with alpha=0.05
pc.fit_13 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 13
tarmat=matrix(0,nrow=ncol(simu.data_M13),ncol = ncol(simu.data_M13))
colnames(tarmat)=colnames(simu.data_M13)
rownames(tarmat)=colnames(simu.data_M13)
tarmat[1,2]=1
tarmat[1,4]=1
tarmat[2,3]=1
Truth_13<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_13,main ="(A) Turth Model 13")
plot(MRPC.fit_13$graph, main ="(A) MRPC")
graphviz.plot(mmhc_13, main ="(C) mmhc")
plot(pc.fit_13,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M14), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M14, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_14<- MRPC(simu.data_M14,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_14<- mmhc(data.frame(simu.data_M14))

## Estimated graph by existing pc with alpha=0.05
pc.fit_14 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 14
tarmat=matrix(0,nrow=ncol(simu.data_M14),ncol = ncol(simu.data_M14))
colnames(tarmat)=colnames(simu.data_M14)
rownames(tarmat)=colnames(simu.data_M14)
tarmat[1,2]=1
tarmat[1,4]=1
tarmat[2,3]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_14<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_14,main ="(A) Turth Model 14")
plot(MRPC.fit_14$graph, main ="(B) MRPC")
graphviz.plot(mmhc_14, main ="(C) mmhc")
plot(pc.fit_14,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M15), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M15, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_15<- MRPC(simu.data_M15,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)
                             
## Estimated graph by mmhc
mmhc_15<- mmhc(data.frame(simu.data_M15))

## Estimated graph by existing pc with alpha=0.05
pc.fit_15 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 15
tarmat=matrix(0,nrow=ncol(simu.data_M15),ncol = ncol(simu.data_M15))
colnames(tarmat)=colnames(simu.data_M15)
rownames(tarmat)=colnames(simu.data_M15)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[3,4]=1

Truth_15<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_15,main ="(A) Turth Model 15")
plot(MRPC.fit_15$graph, main ="(B) MRPC")
graphviz.plot(mmhc_15, main ="(C) mmhc")
plot(pc.fit_15,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M16), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M16, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_16<- MRPC(simu.data_M16,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_16<- mmhc(data.frame(simu.data_M16))

## Estimated graph by existing pc with alpha=0.05
pc.fit_16 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 16
tarmat=matrix(0,nrow=ncol(simu.data_M16),ncol = ncol(simu.data_M16))
colnames(tarmat)=colnames(simu.data_M16)
rownames(tarmat)=colnames(simu.data_M16)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[4,3]=1

Truth_16<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_16,main ="(A) Turth Model 16")
plot(MRPC.fit_16$graph, main ="(B) MRPC")
graphviz.plot(mmhc_16, main ="(C) mmhc")
plot(pc.fit_16,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M17), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M17, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_17<- MRPC(simu.data_M17,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_17<- mmhc(data.frame(simu.data_M17))

## Estimated graph by existing pc with alpha=0.05
pc.fit_17 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 17
tarmat=matrix(0,nrow=ncol(simu.data_M17),ncol = ncol(simu.data_M17))
colnames(tarmat)=colnames(simu.data_M17)
rownames(tarmat)=colnames(simu.data_M17)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[2,4]=1
tarmat[3,4]=1
tarmat[4,3]=1

Truth_17<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_17,main ="(A) Turth Model 17")
plot(MRPC.fit_17$graph, main ="(B) MRPC")
graphviz.plot(mmhc_17, main ="(C) mmhc")
plot(pc.fit_17,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M18), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M18, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_18<- MRPC(simu.data_M18,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_18<- mmhc(data.frame(simu.data_M18))

## Estimated graph by existing pc with alpha=0.05
pc.fit_18 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 18
tarmat=matrix(0,nrow=ncol(simu.data_M18),ncol = ncol(simu.data_M18))
colnames(tarmat)=colnames(simu.data_M18)
rownames(tarmat)=colnames(simu.data_M18)
tarmat[1,2]=1
tarmat[3,2]=1

Truth_18<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_18,main ="(A) Turth Model 18")
plot(MRPC.fit_18$graph, main ="(B) MRPC")
graphviz.plot(mmhc_18, main ="(C) mmhc")
plot(pc.fit_18,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M19), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M19, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_19<- MRPC(simu.data_M19,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_19<- mmhc(data.frame(simu.data_M19))

## Estimated graph by existing pc with alpha=0.05
pc.fit_19 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 19
tarmat=matrix(0,nrow=ncol(simu.data_M19),ncol = ncol(simu.data_M19))
colnames(tarmat)=colnames(simu.data_M19)
rownames(tarmat)=colnames(simu.data_M19)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[2,4]=1

Truth_19<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_19,main ="(A) Turth Model 19")
plot(MRPC.fit_19$graph, main ="(B) MRPC")
graphviz.plot(mmhc_19, main ="(C) mmhc")
plot(pc.fit_19,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M20), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M20, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_20<- MRPC(simu.data_M20,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_20<- mmhc(data.frame(simu.data_M20))

## Estimated graph by existing pc with alpha=0.05
pc.fit_20 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 20
tarmat=matrix(0,nrow=ncol(simu.data_M20),ncol = ncol(simu.data_M20))
colnames(tarmat)=colnames(simu.data_M20)
rownames(tarmat)=colnames(simu.data_M20)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[4,2]=1

Truth_20<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_20,main ="(A) Turth Model 20")
plot(MRPC.fit_20$graph, main ="(B) MRPC")
graphviz.plot(mmhc_20, main ="(C) mmhc")
plot(pc.fit_20,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M21), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M21, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_21<- MRPC(simu.data_M21,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_21<- mmhc(data.frame(simu.data_M21))

## Estimated graph by existing pc with alpha=0.05
pc.fit_21 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 21
tarmat=matrix(0,nrow=ncol(simu.data_M21),ncol = ncol(simu.data_M21))
colnames(tarmat)=colnames(simu.data_M21)
rownames(tarmat)=colnames(simu.data_M21)
tarmat[1,2]=1
tarmat[1,4]=1
tarmat[3,2]=1

Truth_21<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_21,main ="(A) Turth Model 21")
plot(MRPC.fit_21$graph, main ="(B) MRPC")
graphviz.plot(mmhc_21, main ="(C) mmhc")
plot(pc.fit_21,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M22), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M22, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_22<- MRPC(simu.data_M22,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_22<- mmhc(data.frame(simu.data_M22))

## Estimated graph by existing pc with alpha=0.05
pc.fit_22 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 22
tarmat=matrix(0,nrow=ncol(simu.data_M22),ncol = ncol(simu.data_M22))
colnames(tarmat)=colnames(simu.data_M22)
rownames(tarmat)=colnames(simu.data_M22)
tarmat[1,2]=1
tarmat[1,4]=1
tarmat[3,2]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_22<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_22,main ="(A) Turth Model 22")
plot(MRPC.fit_22$graph, main ="(B) MRPC")
graphviz.plot(mmhc_22, main ="(C) mmhc")
plot(pc.fit_22,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M23), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M23, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_23<- MRPC(simu.data_M23,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_23<- mmhc(data.frame(simu.data_M23))

## Estimated graph by existing pc with alpha=0.05
pc.fit_23 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 23
tarmat=matrix(0,nrow=ncol(simu.data_M23),ncol = ncol(simu.data_M23))
colnames(tarmat)=colnames(simu.data_M23)
rownames(tarmat)=colnames(simu.data_M23)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[3,4]=1

Truth_23<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_23,main ="(A) Turth Model 23")
plot(MRPC.fit_23$graph, main ="(B) MRPC")
graphviz.plot(mmhc_23, main ="(C) mmhc")
plot(pc.fit_23,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M24), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M24, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_24<- MRPC(simu.data_M24,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_24<- mmhc(data.frame(simu.data_M24))

## Estimated graph by existing pc with alpha=0.05
pc.fit_24 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 24
tarmat=matrix(0,nrow=ncol(simu.data_M24),ncol = ncol(simu.data_M24))
colnames(tarmat)=colnames(simu.data_M24)
rownames(tarmat)=colnames(simu.data_M24)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[3,4]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_24<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_24,main ="(A) Turth Model 24")
plot(MRPC.fit_24$graph, main ="(B) MRPC")
graphviz.plot(mmhc_24, main ="(C) mmhc")
plot(pc.fit_24,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M25), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M25, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_25<- MRPC(simu.data_M25,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_25<- mmhc(data.frame(simu.data_M25))

## Estimated graph by existing pc with alpha=0.05
pc.fit_25 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 25
tarmat=matrix(0,nrow=ncol(simu.data_M25),ncol = ncol(simu.data_M25))
colnames(tarmat)=colnames(simu.data_M25)
rownames(tarmat)=colnames(simu.data_M25)
tarmat[1,2]=1
tarmat[1,3]=1

Truth_25<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_25,main ="(A) Turth Model 25")
plot(MRPC.fit_25$graph, main ="(B) MRPC")
graphviz.plot(mmhc_25, main ="(C) mmhc")
plot(pc.fit_25,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M26), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M26, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_26<- MRPC(simu.data_M26,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_26<- mmhc(data.frame(simu.data_M26))

## Estimated graph by existing pc with alpha=0.05
pc.fit_26 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 26
tarmat=matrix(0,nrow=ncol(simu.data_M26),ncol = ncol(simu.data_M26))
colnames(tarmat)=colnames(simu.data_M26)
rownames(tarmat)=colnames(simu.data_M26)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,4]=1

Truth_26<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_26,main ="(A) Turth Model 26")
plot(MRPC.fit_26$graph, main ="(B) MRPC")
graphviz.plot(mmhc_26, main ="(C) mmhc")
plot(pc.fit_26,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M27), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M27, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_27<- MRPC(simu.data_M27,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_27<- mmhc(data.frame(simu.data_M27))

## Estimated graph by existing pc with alpha=0.05
pc.fit_27 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 27
tarmat=matrix(0,nrow=ncol(simu.data_M27),ncol = ncol(simu.data_M27))
colnames(tarmat)=colnames(simu.data_M27)
rownames(tarmat)=colnames(simu.data_M27)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[4,2]=1

Truth_27<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_27,main ="(A) Turth Model 27")
plot(MRPC.fit_27$graph, main ="(B) MRPC")
graphviz.plot(mmhc_27, main ="(C) mmhc")
plot(pc.fit_27,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M28), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M28, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_28<- MRPC(simu.data_M28,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_28<- mmhc(data.frame(simu.data_M28))

## Estimated graph by existing pc with alpha=0.05
pc.fit_28 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 28
tarmat=matrix(0,nrow=ncol(simu.data_M28),ncol = ncol(simu.data_M28))
colnames(tarmat)=colnames(simu.data_M28)
rownames(tarmat)=colnames(simu.data_M28)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1

Truth_28<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_28,main ="(A) Turth Model 28")
plot(MRPC.fit_28$graph, main ="(B) MRPC")
graphviz.plot(mmhc_28, main ="(C) mmhc")
plot(pc.fit_28,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M29), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M29, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_29<- MRPC(simu.data_M29,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_29<- mmhc(data.frame(simu.data_M29))

## Estimated graph by existing pc with alpha=0.05
pc.fit_29 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 29
tarmat=matrix(0,nrow=ncol(simu.data_M29),ncol = ncol(simu.data_M29))
colnames(tarmat)=colnames(simu.data_M29)
rownames(tarmat)=colnames(simu.data_M29)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_29<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_29,main ="(A) Turth Model 29")
plot(MRPC.fit_29$graph, main ="(B) MRPC")
graphviz.plot(mmhc_29, main ="(C) mmhc")
plot(pc.fit_29,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M30), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M30, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_30<- MRPC(simu.data_M30,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_30<- mmhc(data.frame(simu.data_M30))

## Estimated graph by existing pc with alpha=0.05
pc.fit_30 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 30
tarmat=matrix(0,nrow=ncol(simu.data_M30),ncol = ncol(simu.data_M30))
colnames(tarmat)=colnames(simu.data_M30)
rownames(tarmat)=colnames(simu.data_M30)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[3,4]=1

Truth_30<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_30,main ="(A) Turth Model 30")
plot(MRPC.fit_30$graph, main ="(B) MRPC")
graphviz.plot(mmhc_30, main ="(C) mmhc")
plot(pc.fit_30,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M31), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M31, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_31<- MRPC(simu.data_M31,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_31<- mmhc(data.frame(simu.data_M31))

## Estimated graph by existing pc with alpha=0.05
pc.fit_31 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 31
tarmat=matrix(0,nrow=ncol(simu.data_M31),ncol = ncol(simu.data_M31))
colnames(tarmat)=colnames(simu.data_M31)
rownames(tarmat)=colnames(simu.data_M31)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[4,3]=1

Truth_31<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_31,main ="(A) Turth Model 31")
plot(MRPC.fit_31$graph, main ="(B) MRPC")
graphviz.plot(mmhc_31, main ="(C) mmhc")
plot(pc.fit_31,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M32), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M32, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_32<- MRPC(simu.data_M32,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_32<- mmhc(data.frame(simu.data_M32))

## Estimated graph by existing pc with alpha=0.05
pc.fit_32 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 32
tarmat=matrix(0,nrow=ncol(simu.data_M32),ncol = ncol(simu.data_M32))
colnames(tarmat)=colnames(simu.data_M32)
rownames(tarmat)=colnames(simu.data_M32)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[3,4]=1
tarmat[4,3]=1

Truth_32<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_32,main ="(A) Turth Model 31")
plot(MRPC.fit_32$graph, main ="(B) MRPC")
graphviz.plot(mmhc_32, main ="(C) mmhc")
plot(pc.fit_32,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M33), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M33, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_33<- MRPC(simu.data_M33,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_33<- mmhc(data.frame(simu.data_M33))

## Estimated graph by existing pc with alpha=0.05
pc.fit_33 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 33
tarmat=matrix(0,nrow=ncol(simu.data_M33),ncol = ncol(simu.data_M33))
colnames(tarmat)=colnames(simu.data_M33)
rownames(tarmat)=colnames(simu.data_M33)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1

Truth_33<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_33,main ="(A) Turth Model 33")
plot(MRPC.fit_33$graph, main ="(B) MRPC")
graphviz.plot(mmhc_33, main ="(C) mmhc")
plot(pc.fit_33,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M34), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M34, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_34<- MRPC(simu.data_M34,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_34<- mmhc(data.frame(simu.data_M34))

## Estimated graph by existing pc with alpha=0.05
pc.fit_34 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 34
tarmat=matrix(0,nrow=ncol(simu.data_M34),ncol = ncol(simu.data_M34))
colnames(tarmat)=colnames(simu.data_M34)
rownames(tarmat)=colnames(simu.data_M34)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[2,4]=1
Truth_34<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_34,main ="(A) Turth Model 34")
plot(MRPC.fit_34$graph, main ="(B) MRPC")
graphviz.plot(mmhc_34, main ="(C) mmhc")
plot(pc.fit_34,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M35), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M35, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_35<- MRPC(simu.data_M35,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_35<- mmhc(data.frame(simu.data_M35))

## Estimated graph by existing pc with alpha=0.05
pc.fit_35 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 35
tarmat=matrix(0,nrow=ncol(simu.data_M35),ncol = ncol(simu.data_M35))
colnames(tarmat)=colnames(simu.data_M35)
rownames(tarmat)=colnames(simu.data_M35)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[4,2]=1
Truth_35<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_35,main ="(A) Turth Model 35")
plot(MRPC.fit_35$graph, main ="(B) MRPC")
graphviz.plot(mmhc_35, main ="(C) mmhc")
plot(pc.fit_35,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M36), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M36, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_36<- MRPC(simu.data_M36,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_36<- mmhc(data.frame(simu.data_M36))

## Estimated graph by existing pc with alpha=0.05
pc.fit_36 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 36
tarmat=matrix(0,nrow=ncol(simu.data_M36),ncol = ncol(simu.data_M36))
colnames(tarmat)=colnames(simu.data_M36)
rownames(tarmat)=colnames(simu.data_M36)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[2,3]=1
tarmat[3,2]=1

Truth_36<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_36,main ="(A) Turth Model 36")
plot(MRPC.fit_36$graph, main ="(B) MRPC")
graphviz.plot(mmhc_36, main ="(C) mmhc")
plot(pc.fit_36,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M37), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M37, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_37<- MRPC(simu.data_M37,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_37<- mmhc(data.frame(simu.data_M37))

## Estimated graph by existing pc with alpha=0.05
pc.fit_37 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 37
tarmat=matrix(0,nrow=ncol(simu.data_M37),ncol = ncol(simu.data_M37))
colnames(tarmat)=colnames(simu.data_M37)
rownames(tarmat)=colnames(simu.data_M37)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[2,4]=1
tarmat[4,2]=1

Truth_37<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_37,main ="(A) Turth Model 37")
plot(MRPC.fit_37$graph, main ="(B) MRPC")
graphviz.plot(mmhc_37, main ="(C) mmhc")
plot(pc.fit_37,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M38), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M38, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_38<- MRPC(simu.data_M38,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_38<- mmhc(data.frame(simu.data_M38))

## Estimated graph by existing pc with alpha=0.05
pc.fit_38 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 38
tarmat=matrix(0,nrow=ncol(simu.data_M38),ncol = ncol(simu.data_M38))
colnames(tarmat)=colnames(simu.data_M38)
rownames(tarmat)=colnames(simu.data_M38)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[3,4]=1

Truth_38<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_38,main ="(A) Turth Model 38")
plot(MRPC.fit_38$graph, main ="(B) MRPC")
graphviz.plot(mmhc_38, main ="(C) mmhc")
plot(pc.fit_38,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M39), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M39, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_39<- MRPC(simu.data_M39,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_39<- mmhc(data.frame(simu.data_M39))

## Estimated graph by existing pc with alpha=0.05
pc.fit_39 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 39
tarmat=matrix(0,nrow=ncol(simu.data_M39),ncol = ncol(simu.data_M39))
colnames(tarmat)=colnames(simu.data_M39)
rownames(tarmat)=colnames(simu.data_M39)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[4,3]=1

Truth_39<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_39,main ="(A) Turth Model 39")
plot(MRPC.fit_39$graph, main ="(B) MRPC")
graphviz.plot(mmhc_39, main ="(C) mmhc")
plot(pc.fit_39,  main ="(D) pc")

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

#Classical correlation 
suffStat_C= list(C = cor(simu.data_M40), n = n)
#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.data_M40, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_40<- MRPC(simu.data_M40,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_40<- mmhc(data.frame(simu.data_M40))

## Estimated graph by existing pc with alpha=0.05
pc.fit_40 <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth for model 40
tarmat=matrix(0,nrow=ncol(simu.data_M40),ncol = ncol(simu.data_M40))
colnames(tarmat)=colnames(simu.data_M40)
rownames(tarmat)=colnames(simu.data_M40)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[1,4]=1
tarmat[2,3]=1
tarmat[3,2]=1
tarmat[3,4]=1
tarmat[4,3]=1

Truth_40<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_40,main ="(A) Turth Model 40")
plot(MRPC.fit_40$graph, main ="(B) MRPC")
graphviz.plot(mmhc_40, main ="(C) mmhc")
plot(pc.fit_40,  main ="(D) pc")


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

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_0_Outlier_beta0.005 <- MRPC(simu.outdata_M0,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest and beta=0
MRPC.fit_0_Outlier_beta0 <- MRPC(simu.outdata_M0,suffStat_C,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_0_Outlier<- mmhc(data.frame(simu.outdata_M0))

## Estimated graph by existing pc with alpha=0.05
pc.fit_0_Outlier <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(simu.data_M0),ncol = ncol(simu.data_M0))
colnames(tarmat)=colnames(simu.data_M0)
rownames(tarmat)=colnames(simu.data_M0)
tarmat[1,2]=1

Truth_0<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(3,2))
plot(Truth_0,main ="(A) Turth Model 0 with Outliers")
plot(MRPC.fit_0_Outlier_beta0.005$graph, main ="(B) MRPC (beta=0.005)")
plot(MRPC.fit_0_Outlier_beta0$graph, main ="(C) MRPC(beta=0)")
graphviz.plot(mmhc_0_Outlier, main ="(D) mmhc")
plot(pc.fit_0_Outlier,  main ="(E) pc")

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

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_1_Outlier_beta0.005 <- MRPC(simu.outdata_M1,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest and beta=0
MRPC.fit_1_Outlier_beta0 <- MRPC(simu.outdata_M1,suffStat_C,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_1_Outlier<- mmhc(data.frame(simu.outdata_M1))

## Estimated graph by existing pc with alpha=0.05
pc.fit_1_Outlier <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(simu.data_M1),ncol = ncol(simu.data_M1))
colnames(tarmat)=colnames(simu.data_M1)
rownames(tarmat)=colnames(simu.data_M1)
tarmat[1,2]=1
tarmat[2,3]=1
Truth_1<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(3,2))
plot(Truth_1,main ="(A) Turth Model 1 with Outliers")
plot(MRPC.fit_1_Outlier_beta0.005$graph, main ="(B) MRPC (beta=0.005)")
plot(MRPC.fit_1_Outlier_beta0$graph, main ="(C) MRPC(beta=0)")
graphviz.plot(mmhc_1_Outlier, main ="(D) mmhc")
plot(pc.fit_1_Outlier,  main ="(E) pc")

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

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_2_Outlier_beta0.005 <- MRPC(simu.outdata_M2,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest and beta=0
MRPC.fit_2_Outlier_beta0 <- MRPC(simu.outdata_M2,suffStat_C,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_2_Outlier<- mmhc(data.frame(simu.outdata_M2))

## Estimated graph by existing pc with alpha=0.05
pc.fit_2_Outlier <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(simu.data_M2),ncol = ncol(simu.data_M2))
colnames(tarmat)=colnames(simu.data_M2)
rownames(tarmat)=colnames(simu.data_M2)
tarmat[1,2]=1
tarmat[3,2]=1
Truth_2<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(3,2))
plot(Truth_2,main ="(A) Turth Model 2 with Outliers")
plot(MRPC.fit_2_Outlier_beta0.005$graph, main ="(B) MRPC (beta=0.005)")
plot(MRPC.fit_2_Outlier_beta0$graph, main ="(C) MRPC(beta=0)")
graphviz.plot(mmhc_2_Outlier, main ="(D) mmhc")
plot(pc.fit_2_Outlier,  main ="(E) pc")

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

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_3_Outlier_beta0.005 <- MRPC(simu.outdata_M3,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest and beta=0
MRPC.fit_3_Outlier_beta0 <- MRPC(simu.outdata_M3,suffStat_C,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_3_Outlier<- mmhc(data.frame(simu.outdata_M3))

## Estimated graph by existing pc with alpha=0.05
pc.fit_3_Outlier <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(simu.data_M3),ncol = ncol(simu.data_M3))
colnames(tarmat)=colnames(simu.data_M3)
rownames(tarmat)=colnames(simu.data_M3)
tarmat[1,2]=1
tarmat[1,3]=1
Truth_3<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(3,2))
plot(Truth_3,main ="(A) Turth Model 3 with Outliers")
plot(MRPC.fit_3_Outlier_beta0.005$graph, main ="(B) MRPC (beta=0.005)")
plot(MRPC.fit_3_Outlier_beta0$graph, main ="(C) MRPC(beta=0)")
graphviz.plot(mmhc_3_Outlier, main ="(D) mmhc")
plot(pc.fit_3_Outlier,  main ="(E) pc")

## Model 4
Outdata=matrix(c(runif(n = 300, min = 80, max = 100)),ncol = 3)
simu.data_M4[c(901:1000),2]=Outdata[,1] #Adding outliers
simu.outdata_M4=simu.data_M4 #New data with outliers

n<-nrow (simu.outdata_M4)    #Number of row
V<-colnames(simu.outdata_M4) #Column names

#Classical correlation (Beta=0)
Rcor_C=RobustCor(simu.outdata_M4, 0)  
suffStat_C= list(C = Rcor_C$RR, n = n)

#Robust correlation (Beta=0.005)
Rcor_R=RobustCor(simu.outdata_M4, 0.005) 
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_4_Outlier_beta0.005 <- MRPC(simu.outdata_M4,suffStat_R,GV=1,FDR=0.05,
                             indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by MRPC using gaussCItest and beta=0
MRPC.fit_4_Outlier_beta0 <- MRPC(simu.outdata_M4,suffStat_C,GV=1,FDR=0.05,
                         indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_4_Outlier<- mmhc(data.frame(simu.outdata_M4))

## Estimated graph by existing pc with alpha=0.05
pc.fit_4_Outlier <- pc(suffStat_C,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(simu.data_M4),ncol = ncol(simu.data_M4))
colnames(tarmat)=colnames(simu.data_M4)
rownames(tarmat)=colnames(simu.data_M4)
tarmat[1,2]=1
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,2]=1

Truth_4<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(3,2))
plot(Truth_4,main ="(A) Turth Model 4 with Outliers")
plot(MRPC.fit_4_Outlier_beta0.005$graph, main ="(B) MRPC (beta=0.005)")
plot(MRPC.fit_4_Outlier_beta0$graph, main ="(C) MRPC(beta=0)")
graphviz.plot(mmhc_4_Outlier, main ="(D) mmhc")
plot(pc.fit_4_Outlier,  main ="(E) pc")


#Continuous phenotypes without genetic variant
#N=1000 with moderate signals

#parameters
N=1000
signal=0.5
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1

T1=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T4=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T2=Case_1P(N=N,P1=T1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T3=Case_1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T5=Case_1P(N=N,P1=T2,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T6=Case_2P(N=N,P1=T1,P2=T5,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
T7=Case_1P(N=N,P1=T6,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T8=Case_2P(N=N,P1=T1,P2=T5,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
#Data
Data_test1=cbind(T1,T2,T3,T4,T5,T6,T7,T8)

n<-nrow (Data_test1)    #Number of row
V<-colnames(Data_test1) #Column names

#Classical correlation
suffStat_C1= list(C = cor(Data_test1), n = n)

#Robust correlation (Beta=0.005)
Rcor_R1=RobustCor(Data_test1, 0.005) 
suffStat_R1= list(C = Rcor_R1$RR, n = n)

## Estimated graph by MRPC using gaussCItest and beta=0.005
MRPC.fit_test1<- MRPC(Data_test1,suffStat_R1,GV=0,FDR=0.05,
                      indepTest ='gaussCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_test1=mmhc(data.frame(Data_test1))

## Estimated graph by existing pc with alpha=0.05
pc.fit_test1<- pc(suffStat_C1,indepTest = gaussCItest,alpha=0.05,labels =V,verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(Data_test1),ncol = ncol(Data_test1))
colnames(tarmat)=colnames(Data_test1)
rownames(tarmat)=colnames(Data_test1)
tarmat[1,2]=1
tarmat[1,6]=1
tarmat[1,8]=1
tarmat[2,3]=1
tarmat[2,5]=1
tarmat[5,6]=1
tarmat[5,8]=1
tarmat[6,7]=1

Truth_test1<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(1,4))
plot(Truth_test1,main=" (A) Truth")
plot(MRPC.fit_test1$graph,main=" (B) MRPC")
graphviz.plot(mmhc_test1,main=" (C) mmhc")
plot(pc.fit_test1,main=" (D) pc")


#Discrete phenotypes with genetic variant example 1,
#N=1000 with strong signals

#parameters
N=1000
p=0.45
signal=1
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1
#Data
V1 <- c(sample(c(0, 1, 2),
               size = N,
               replace = TRUE,
               prob = c((1 - p)^2,
                        2*p*(1 - p),
                        p^2)))
T1=Case_1P(N=N,P1=V1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T2=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T3=Case_2P(N=N,P1=T1,P2=T2,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
T4=Case_1P(N=N,P1=T3,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T5=Case_1P(N=N,P1=T4,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)

#Data
D1=cbind(V1,T1,T2,T3,T4,T5)

#category
D1[which(D1<0)]=0
D1[which(D1>0 & D1<1)]=1
D1[which(D1>1)]=2

#
n<-nrow (D1)    #Number of row
V<-colnames(D1) #Column names
#suffStat for discrete data
suffStat1= list(dm = D1, adaptDF = FALSE)

## Estimated graph by existing MRPC
MRPC_fit_D1=MRPC(D1,suffStat1,GV=1,FDR=0.05,indepTest ='disCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_D1=mmhc(data.frame(D1))

## Estimated graph by existing pc with alpha=0.05
pc.fit_D1=pc(suffStat1,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.05, labels = V, verbose = TRUE)

#plot(M1)
#Truth
tarmat=matrix(0,nrow=ncol(D1),ncol = ncol(D1))
colnames(tarmat)=colnames(D1)
rownames(tarmat)=colnames(D1)
tarmat[1,2]=1
tarmat[2,4]=1
tarmat[3,4]=1
tarmat[4,5]=1
tarmat[5,6]=1
Truth_D1<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_D1, main = " (A) Truth") 
plot(MRPC_fit_D1$graph,main = " (B) MRPC")
graphviz.plot(mmhc_D1,main=" (C) mmhc")
plot(pc.fit_D1,main = " (D) pc")

##Discrete phenotypes with genetic variant example 2
#Parameters
N=1000
p=0.45
signal=1
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1
#Data
V1 <- c(sample(c(0, 1, 2),
               size = N,
               replace = TRUE,
               prob = c((1 - p)^2,
                        2*p*(1 - p),
                        p^2)))
T1=Case_1P(N=N,P1=V1,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T3=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T2=Case_2P(N=N,P1=T1,P2=T3,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
T4=Case_1P(N=N,P1=T3,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T5=Case_1P(N=N,P1=T4,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)

D2=cbind(V1,T1,T2,T3,T4,T5)

#category
D2[which(D2<0)]=0
D2[which(D2>0 & D2<1)]=1
D2[which(D2>1)]=2

#
n<-nrow (D2)    #Number of row
V<-colnames(D2) #Column names

suffStat2= list(dm = D2, adaptDF = FALSE)

## Estimated graph by existing MRPC
MRPC_fit_D2=MRPC(D2,suffStat2,GV=1,FDR=0.05,indepTest ='disCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_D2=mmhc(data.frame(D2))

## Estimated graph by existing pc with alpha=0.05
pc.fit_D2=pc(suffStat2,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.05, labels = V, verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(D2),ncol = ncol(D2))
colnames(tarmat)=colnames(D2)
rownames(tarmat)=colnames(D2)
tarmat[1,2]=1
tarmat[2,3]=1
tarmat[4,3]=1
tarmat[4,5]=1
tarmat[5,6]=1
Truth_D2<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_D2, main = " (A) Truth") 
plot(MRPC_fit_D2$graph,main = " (B) MRPC")
graphviz.plot(mmhc_D2,main=" (C) mmhc")
plot(pc.fit_D2,main = " (D) pc")



##Discrete phenotypes without genetic variant example 1
#Parameters
N=1000
#p=0.45
signal=1
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1
#Data
T1=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T2=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T3=Case_2P(N=N,P1=T1,P2=T2,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
T4=Case_1P(N=N,P1=T3,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T5=Case_1P(N=N,P1=T4,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)

#Data
D3=cbind(T1,T2,T3,T4,T5)
#category
D3[which(D3<0)]=0
D3[which(D3>0 & D3[,]<1)]=1
D3[which(D3>1)]=2

n<-nrow (D3)    #Number of row
V<-colnames(D3) #Column names
suffStat3= list(dm = D3, adaptDF = FALSE)

## Estimated graph by MRPC 
MRPC_fit_D3=MRPC(D3,suffStat3,GV=0,FDR=0.05,indepTest ='disCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc
mmhc_D3=mmhc(data.frame(D3))

## Estimated graph by pc
pc.fit_D3=pc(suffStat3,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.05, labels = V, verbose = TRUE)


#Truth
tarmat=matrix(0,nrow=ncol(D3),ncol = ncol(D3))
colnames(tarmat)=colnames(D3)
rownames(tarmat)=colnames(D3)
tarmat[1,3]=1
tarmat[2,3]=1
tarmat[3,4]=1
tarmat[4,5]=1

Truth_D3<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_D3, main = " (A) Truth") 
plot(MRPC_fit_D3$graph,main = " (B) MRPC")
graphviz.plot(mmhc_D3,main=" (C) mmhc")
plot(pc.fit_D3,main = " (D) pc")


##Discrete phenotypes without genetic variant example 2
#Parameters
N=1000
signal=1
b0.1 = 0
b1.1 =signal
b1.2 =signal
b1.3 =signal
sd.1 = 1
#Data
T1=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T3=Case_NP(N=N,b0.1=b0.1,sd.1=sd.1)
T2=Case_2P(N=N,P1=T1,P2=T3,b0.1=b0.1,b1.1=b1.1,b1.2=b1.2,sd.1=sd.1)
T4=Case_1P(N=N,P1=T3,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)
T5=Case_1P(N=N,P1=T4,b0.1=b0.1,b1.1=b1.1,sd.1=sd.1)

D4=cbind(T1,T2,T3,T4,T5)

#category
D4[which(D4<0)]=0
D4[which(D4>0 & D4<1)]=1
D4[which(D4>1)]=2

#
n<-nrow (D4)    #Number of row
V<-colnames(D4) #Column names
suffStat4= list(dm = D4, adaptDF = FALSE)

## Estimated graph by MRPC 
MRPC_fit_D4=MRPC(D4,suffStat4,GV=0,FDR=0.05,indepTest ='disCItest',labels=V,verbose = TRUE)

## Estimated graph by mmhc 
mmhc_D4=mmhc(data.frame(D4))

## Estimated graph by pc with alpha = 0.05
pc.fit_D4=pc(suffStat4,
           ## independence test: G^2 statistic
           indepTest = disCItest, alpha = 0.05, labels = V, verbose = TRUE)

#Truth
tarmat=matrix(0,nrow=ncol(D4),ncol = ncol(D4))
colnames(tarmat)=colnames(D4)
rownames(tarmat)=colnames(D4)
tarmat[1,2]=1
tarmat[3,2]=1
tarmat[3,4]=1
tarmat[4,5]=1

Truth_D4<-as(tarmat, "graphNEL")

## Plot of the results
## Show estimated graph
par(mfrow=c(2,2))
plot(Truth_D4, main = " (A) Truth") 
plot(MRPC_fit_D4$graph,main = " (B) MRPC")
graphviz.plot(mmhc_D4,main=" (C) mmhc")
plot(pc.fit_D4,main = " (D) pc")

}
