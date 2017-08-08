\name{ModiSkeleton}
\alias{ModiSkeleton}

\title{
Draw Undirected Graph using the Modified Skeleton Function
}
\description{
This Modiskeleton () function is used for estimate the undirected graph (skeleton) using the modified skeleton function. The main difference is that the ModiSkeleton () function is used alpha (significance level) as a function of the prior outcomes (reject or accept hypothesis) but the existing \link{skeleton} function used alpha as a fixed pre-assigned significance level. See details in \bold{Choice of Significance Level} below.
}
\usage{
ModiSkeleton(data,suffStat,FDR, indepTest, labels, p,
method = c("stable","original", "stable.fast"),
m.max = Inf,fixedGaps = NULL, fixedEdges = NULL,
NAdelete = TRUE, verbose = FALSE)
}

\arguments{
Since, this is the modification of the existing \link{skeleton} function, therefore, most of the terms and sentences are taken from \link{skeleton} (Spirtes et al., 2000). Also, most of the arguments are identical to \link{MRPC} and \link{pc}.

\item{data}{
Data matrix, where row is the samples and column is the genes. The column order from  Single-nucleotide polymorphism (SNPs), indels,copy number variation (CNV) or Expression quantitative trait loci (eQTL) to genes. If we consider one SNPs/CNV/indels/CNV/eQTL, then the first column of the input matrix is SNPs/CNV/indels/CNV/eQTL and the remaining columns are the gene expression data and so on.
}
\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistic consists of the correlation matrix of the data and the sample size. Note that, \link{MRPC} used robust correlation and existing \link{pc} used classical pearson correlation.
}
\item{FDR}{
Need to specify pre-assigned level.If FDR=0.05, that ensures FDR and mFDR remains below 0.05.
}
\item{indepTest}{A \link{function} for testing conditional independence. It is used to tests the conditional independence of x and y given S, called as indepTest(x,y,S,suffStat). Where, x and y are variables, and S is a (possibly empty) vector of variables. suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. The different indepTest is used for different data types, for example, Gaussian data = gaussCItest, Discrete data = disCItest and Binary data = binCItest. See help(gaussCItest)

The \link{ci.test} (Marco Scutari, 2010) is also used for testing conditional independence and return value of indepTest is the p-value. If none is specified, the default test statistic is the mutual information for categorical variables, the Jonckheere-Terpstra test for ordered factors and the linear correlation for continuous variables.See help(ci.test)

Remember that need to specify the which indepTest would like for Independence testing. For example, if you would like to use gaussCItest you would type indepTest='gaussCItest' into the function otherwise indepTest='citest'. Note that, we used \link{gaussCItest} to compare our \link{MRPC} with the existing \link{pc}, because of \link{ci.test} is not robust. See details in example.
}
\item{labels}{
character vector of variable (or node) names.
}
\item{p}{
(optional) number of variables (or nodes). May be specified if labels are not, in which case labels is set to 1:p.
}
  \item{method}{
(optional) Character string specifying method; the default, "stable" provides an order-independent skeleton.
}
  \item{m.max}{
(optional) Maximal size of the conditioning sets that are considered in the conditional independence tests.
}
  \item{fixedGaps}{
(optional) A logical matrix of dimension p*p. If entry [x,y] or [y,x] (or both) are TRUE, the edge x---y is removed before starting the algorithm. Therefore, this edge is guaranteed to be absent in the resulting graph.
}
  \item{fixedEdges}{
(optional) A logical matrix of dimension p*p. If entry [x,y] or [y,x] (or both) are TRUE, the edge x---y is never considered for removal. Therefore, this edge is guaranteed to be present in the resulting graph.
}
  \item{NAdelete}{
(optional) If indepTest returns NA and this option is TRUE, the corresponding edge is deleted. If this option is FALSE, the edge is not deleted.
}
  \item{verbose}{
(optional) If TRUE, detailed output is provided. Default is FALSE for no output details
}
}
\details{
We incorporated sequential hypothesis testing to draw the undirected graph (skeleton) by \link{ModiSkeleton} function  (similar as the pc algorithm by \link{skeleton} function but, difference is that the \link{ModiSkeleton} function is used alpha (significance level) as a function of the prior outcomes (reject or accept hypothesis) but the \link{skeleton} function used alpha as a fixed pre-assigned significance level and \link{ModiSkeleton} take a robust correlation approach (Badsha et al., 2013) reduces the impact of outliers but the \link{skeleton} used classical correlation which is highly influnced by outliers. See details in in \link{ModiSkeleton} and \link{skeleton}). The other's term is the same as existing \link{skeleton}.

The skeleton function starts with a complete undirected graph. Then a series of conditional independence tests is done and edges are deleted in the following way.

\bold{First}, all pairs of nodes are tested for marginal independence. If two nodes x and y are judged to be marginally independent at level alpha (we adjusted the significance level at every test, see below for the details \bold{Choice of Significance Level}), the edge between them is deleted and the empty set is saved as separation sets S[x, y] and S[y, x]. After all pairs have been tested for marginal independence and some edges might have been removed.

\bold{Second}, all pairs of nodes (x, y) are tested for conditional independence given any single node or all possible subsets of the remaining nodes. If there is any node z such that x and y are conditionally independent given z, the edge between x and y is removed and node z is saved as separation sets (sepset) S[x, y] and S[y, x]. The algorithm continues in this way by increasing the size of the conditioning set step by step. The algorithm stops if all adjacency sets in the current graph are smaller than the size of the conditioning set. The result is the skeleton in which every edge is still undirected.

Now,each triple of vertices (x, z, y) such that the pairs (x, z) and (y, z) are each adjacent in the skeleton but (x, y) are not (such a triple is called an 'unshielded triple'), is oriented based on the information saved in the conditioning sepset S[x, y] and S[y, x]. More precisely, an unshielded triple x---z---y is oriented as x--->z<---y (called v-structure) if z is not in S[x, y] = S[y, x] i.e., x and y are conditionally dependent in given z.

\bold{Choice of Significance Level}:
One of the most important part of our method is that choice of significance level. To draw the undirected skeleton, the significance level (number in (0,1)) treated as a tuning parameter for the marginal and conditional independence tests. At each step of these test, we must decide whether to reject or accept hypothesis, without having access to the number of hypotheses (potentially infinite) or the future p-values. Hypotheses with p-values below a significance level, typically 0.05, are considered to be statistically significant. While this ritual controls type I errors for single testing problems, in case of testing multiple hypotheses we need to adjust the significance levels to control other metrics such as family-wise error rate (FWER) or false discovery rate (FDR) or marginal FDR (mFDR). FWER is the probability of rejecting at least one true null hypothesis, which is very conservative for multiple testing. On the other hand, FDR is the expected proportion of false positives among rejected null hypotheses (Benjamini and Hochberg, 1995) and mFDR is the expected number of false discoveries to the expected number of discoveries. FDR controls a property of one realized set of tests, while mFDR is the ratio of two expectations over many realizations.

To draw the undirected skeleton, existing \link{skeleton} algorithm (Spirtes et al., 2000) used alpha as a fixed pre-assigned significance level, which does not control the FDR and mFDR. Kalisch and Buhlmann, 2007 also used alpha as a fixed pre-assigned value based on Structural Hamming Distance (SHD) suggested by Tsamardinos et al., (2006) for estimating high-dimensional DAG with the \link{pc} algorithm. But, problem is that to calculate SHD, we have to need the true graph, which is unknown in our study.

Another problem of using fixed pre-assigned significance level (alpha) is that if the p-values close to the value of alpha, then it's difficult to accept or reject null hypotheses. For example, if p-value=0.02312 and if we used fixed alpha=0.05, then it's consider to be statistically significant. But, if we used alpha=0.01, then it's consider to be statistically insignificant (see examples in \link{MRPC} for model 0 and model 3). So, now, question is that how can we choose the appropriate value of alpha.
To avoid these problem, we need to adjust the alpha value at each step of the test. Benjamini and Hochberg, 1995 also proposed a sequential testing procedure to control FDR in multiple testing problems assuming that all the p-values are given a priori but it does not address the scenario described above.

Therefore, the \link{ModiSkeleton} function used LOND algorithm for alpha (level of significance) as a function of the prior outcomes (reject or accept) that control FDR and mFDR as online manner (Javanmard and Montanari, 2015), which ensures that FDR remains below a pre-assigned leve alpha.

}

\value{
An object containing an estimate of the skeleton of the underlying DAG as follow:
\itemize{
\item \code{obj}: object of class "graph" of the directed and undirected graph that was estimated.
\item \code{test}: current test number
\item \code{alpha}: level of significance that updated in every step of sequential FDR test.
\item \code{R}: rejection number.
}
}
\references{
1.Badsha, M. B. Mollah, M. N. Jahan, N. and Kurata, H. (2013). Robust complementary hierarchical clustering for gene expression data analysis by beta-divergence. J Biosci Bioeng 116(3): 397-407.

2. Benjamini, Y. and Y Hochberg(1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing, J. R. Statist. Soc. B, B, 57, 289-300.

3. Javanmard and Montanari (March 5, 2015) On Online Control of False Discovery Rate. arXiv:150206197 [statME].

4. Kalisch, M. and Buhlmann, P. (2007) Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm, Journal of Machine Learning Research, 8, 613-636.

5.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

6.Tsamardinos, I. Brown, L.E. and Aliferis, C.F. (2006). The Max-Min Hill-Climbing Bayesian Network Structure Learning Algorithm. JMLR 65, 31-78.

7. Spirtes, P. Glymour, C and Scheines, R (2000). Causation, Prediction, and Search, 2nd edition. The MIT Press.
}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
 \link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based PC (MRPC) algorithm; \link{EdgeOrientation} for orientation rules for edges in MRPC; \link{SimulatedData} simulated data generating function.
}
\examples{
## Load packages
library(pcalg)  #library for existing pc
library(bnlearn)#library for ci.test
## Load predefined simulated data
## Data pre-processing

#The 1st column of the input matrix will be the
#genotype of the expression quantitative trait loci
#(eQTL)/Copy number variation (CNVs) and the remaining
#columns are the gene expression data.
# We used pre-assigned level alpha=0.05 that ensures
#FDR and mFDR remains below 0.05.

## Model 0
## Estimate Skeleton
n<-nrow (simu.data_M0)    #Number of row
V<-colnames(simu.data_M0) #Column names

Rcor_R=RobustCor(simu.data_M0, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

Skel.fit_0 <- ModiSkeleton(simu.data_M0,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Plot of the results
## Show estimated skeleton
plot(Skel.fit_0$graph,  main ="Estimated Skeleton for Model 0")

## Model 1
## Estimate Skeleton
n<-nrow (simu.data_M1)    #Number of row
V<-colnames(simu.data_M1) #Column names

Rcor_R=RobustCor(simu.data_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

Skel.fit_1 <- ModiSkeleton(simu.data_M1,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Plot of the results
## Show estimated skeleton
plot(Skel.fit_1$graph,  main ="Estimated Skeleton for Model 1")

## Model 2
## Estimate Skeleton
n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M2) #Column names

Rcor_R=RobustCor(simu.data_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

Skel.fit_2 <- ModiSkeleton(simu.data_M2,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Plot of the results
## Show estimated skeleton
plot(Skel.fit_2$graph,  main ="Estimated Skeleton for Model 2")

## Model 3
## Estimate Skeleton
n<-nrow(simu.data_M3)     #Number of row
V<-colnames(simu.data_M3) #Column names

Rcor_R=RobustCor(simu.data_M3, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

Skel.fit_3 <- ModiSkeleton(simu.data_M3,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Plot of the results
## Show estimated skeleton
plot(Skel.fit_3$graph,  main ="Estimated Skeleton for Model 3")

## Model 4
## Estimate Skeleton
n<-nrow (simu.data_M4)    #Number of row
V<-colnames(simu.data_M4) #Column names

Rcor_R=RobustCor(simu.data_M4, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

Skel.fit_4<- ModiSkeleton(simu.data_M4,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Plot of the results
## Show estimated skeleton
plot(Skel.fit_4$graph,  main ="Estimated Skeleton for Model 4")

}

