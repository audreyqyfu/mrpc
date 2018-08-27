\name{EdgeOrientation}
\alias{EdgeOrientation}

\title{
Edges Orientation Rules for the MRPC Algorithm
}
\description{
This function perform the last step of \link{MRPC} algorithm, where we determine the edge direction for undirected graph. The functions first determine the edges between genetic variants with the others based on MR, then orient v-structures, then orient the remaining edges, weather MR is applicable or not. MR is a new way for edge direction determination based four different cases. See below for the details.
}
\usage{
EdgeOrientation(gInput,GV=GV,suffStat, FDR, indepTest = indepTest, verbose = FALSE)

}

\arguments{
\item{gInput}{Object containing skeleton, marginal and conditional independence information.
}
\item{GV}{Need to specify for how many genetic variants (SNPs/indels/CNV/eQTL) in your input data. For example, if your dataset have one genetic variants (first column), then GV=1, if 2 (1st and 2nd Column), then GV=2..and so on.
}
\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistic consists of the correlation matrix of the data and the sample size.
}

\item{FDR}{
Need to specify pre-assigned level.If FDR=0.05, that ensures FDR and mFDR remains below 0.05.
}
\item{indepTest}{A \link{function} for testing conditional independence. It is used to tests the conditional independence of x and y given S, called as indepTest(x,y,S,suffStat). Where, x and y are variables, and S is a (possibly empty) vector of variables. suffStat is a list, see the argument above. The return value of indepTest is the p-value of the test for conditional independence. The different indepTest is used for different data types, for example, Gaussian data = gaussCItest, Discrete data = disCItest and Binary data = binCItest. See help(gaussCItest)}
\item{verbose}{(optional) 1: detailed output is provided; 0: No output is provided}
}
\details{
The orientation of the edges direction based on Mendelian randomization using the four different cases. Here, we consider x is a genetic variants, y and z are the gene expression data. The 1st column of the input matrix will be the genetic variants and the remaining column are the gene expression data.

Four different cases are as follows:

\bold{Case-1}: Relation between x (genetic variants) and other's nodes,then genetic variants will regulates the other's nodes (genes) and direction will be genetic variants-->other's nodes. Note that, if dataset have more than one genetic variants and two genetic variants have edges, then direction will be genetic variants<-->genetic variants, which indicating that which indicating that we have evidence that the two genetic variants are not independent, but we do not have enough information to determine which genetic variants is the regulator and which is the target.

\bold{Case-2}: If, y and z are adjacent and, x and z conditionally independent given y, then gene y will regulates the expression of gene z and the edge direction will be y-->z.


\bold{Case-3}: If, y and z are adjacent and, x and z conditionally dependent given y, then gene z will regulates the expression of gene y and the edge direction will be z-->y.

\bold{Case-4}: If, y and z are adjacent with x and y conditionally dependent given z with x and z conditionally dependent given y, then the edge direction will be y<-->z.

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

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{SimulatedData} for simulated data generating function.

}

\examples{
## Load packages
library(pcalg)  #library for existing pc
## Load predefined data
## Data pre-processing

#The 1st column of the input matrix will be the
#genetic variants and the remaining columns are the gene expression #data.

## Model 1
Truth<-MRPCtruth$M1   #Truth for model 1
data<-simu.data_M1    #data load for model 1
n<-nrow (data)    #Number of row
V<-colnames(data) #Column names

Rcor_R<-RobustCor(data, 0.005) #Robust correlation (Beta=0.005)
suffStat_R<- list(C = Rcor_R$RR, n = n)

## Estimate skeleton
Skel.fit<- ModiSkeleton(data,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation<- EdgeOrientation(Skel.fit,GV=1,suffStat_R,FDR=0.05,indepTest
='gaussCItest',verbose=1)

## Plot of the results
#if (require(Rgraphviz)){
par(mfrow=c(1,2))
plot(Truth,main ="(A) Truth")
plot(Edge_orientation,main ="(B) MRPC ")
#}
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

