\name{EdgesOrientation}
\alias{EdgesOrientation}

\title{
Edges Orientation Rules for the MRPC Algorithm
}
\description{
This function perform the last step of \link{MRPC} algorithm, where we determine the edge direction for undirected graph. The functions first determine the edges between genetic variants with the others based on MR, then orient v-structures, then orient the remaining edges, weather MR is applicable or not. MR is a new way for edge direction determination based four different cases. See below for the details.
}
\usage{
EdgesOrientation(gInput,NQ, suffStat, FDR, verbose = FALSE)

}

\arguments{
\item{gInput}{Object containing skeleton, marginal and conditional independence information.
}

\item{NQ}{Need to specify for how many genetic variants (SNPs/indels/CNV/eQTL) in your input data. For example, if your dataset have one genetic variants (first column), then NQ=1, if 2 (1st and 2nd Column), then NQ=2..and so on.
}
\item{suffStat}{A \link{list} of sufficient statistics, containing all necessary elements for the conditional independence tests in the function indepTest for gaussCItest. The sufficient statistic consists of the correlation matrix of the data and the sample size.
}

\item{FDR}{
Need to specify pre-assigned level.If FDR=0.05, that ensures FDR and mFDR remains below 0.05.
}
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
A class of objects that containing the information for directed graph.
}


\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
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

## Model 0
n<-nrow (simu.data_M0)    #Number of row
V<-colnames(simu.data_M0) #Column names

Rcor_R=RobustCor(simu.data_M0, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)


## Estimate skeleton

Skel.fit_0 <- ModiSkeleton(simu.data_M0,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation_0 <- EdgesOrientation(Skel.fit_0,NQ=1, suffStat_R,FDR=0.05,verbose=1)

## Plot of the results

plot(Edge_orientation_0$obj,main ="Model 0")

## Model 1
n<-nrow (simu.data_M1)    #Number of row
V<-colnames(simu.data_M1) #Column names

Rcor_R=RobustCor(simu.data_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)
## Estimate skeleton

Skel.fit_1 <- ModiSkeleton(simu.data_M1,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation_1 <- EdgesOrientation(Skel.fit_1,NQ=1,suffStat_R,FDR=0.05, verbose=1)

## Plot of the results

plot(Edge_orientation_1$obj,main ="Model 1")

## Model 2
n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M2) #Column names

Rcor_R=RobustCor(simu.data_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate skeleton

Skel.fit_2 <- ModiSkeleton(simu.data_M2,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation_2 <- EdgesOrientation(Skel.fit_2,NQ=1,suffStat_R,FDR=0.05, verbose=1)

## Plot of the results

plot(Edge_orientation_2$obj,main ="Model 2")

## Model 3
n<-nrow (simu.data_M3)    #Number of row
V<-colnames(simu.data_M3) #Column names

Rcor_R=RobustCor(simu.data_M3, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimate skeleton

Skel.fit_3 <- ModiSkeleton(simu.data_M3,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation_3 <- EdgesOrientation(Skel.fit_3,NQ=1,suffStat_R,FDR=0.05, verbose=1)

## Plot of the results

plot(Edge_orientation_3$obj,main ="Model 3")

## Model 4
n<-nrow (simu.data_M4)     #Number of row
V<-colnames(simu.data_M4)  #Column names

Rcor_R=RobustCor(simu.data_M4, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)
## Estimate skeleton

Skel.fit_4 <- ModiSkeleton(simu.data_M4,suffStat_R,FDR=0.05,indepTest=
'gaussCItest',labels=V,verbose=TRUE)

## Orientation edges
Edge_orientation_4 <- EdgesOrientation(Skel.fit_4,NQ=1,suffStat_R,FDR=0.05, verbose=1)

## Plot of the results

plot(Edge_orientation_4$obj,main ="Model 4")

}

