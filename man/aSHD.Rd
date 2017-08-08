\name{aSHD}
\alias{aSHD}

\title{
adjusted Structural Hamming Distance (aSHD) 
}
\description{
The SHD, as is implemented in the R package pcalg (Kalisch et al., 2012) and bnlearn(Scutari, 2010), counts how many differences exist between two directed graphs. This distance is 1 if an edge exists in one graph but missing in the other, or if the direction of an edge is different in the two graphs. The larger this distance, the more different the two graphs are. We adjusted the SHD to reduce the penalty on the wrong direction of an edge to 0.5. For example, between two graphs V-->T1<--T2 and V-->T1-->T2, the SHD is 1 and our aSHD is 0.5.
}
\usage{
aSHD(g1, g2, GV)
}

\arguments{
  \item{g1}{
First graph object
}
  \item{g2}{
Second graph object
}
\item{GV}{Need to specify for how many genetic variants (SNPs/indels/CNV/eQTL) in your input data for the consistence of MRPC with others inference, because only the MRPC implemented the gMR.For example, if your dataset have one genetic variants (first column), then GV=1, if 2 (1st and 2nd Column), then GV=2..and so on.
}
}

\references{
1.Kalisch, M., Machler, M., Colombo, D., Maathuis, M.H. & Buhlmann, P. Causal Inference Using Graphical Models with the R Package pcalg. J. Stat.Softw. 47, 26 (2012).

2.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
## True model
n<-nrow (simu.data_M1)     #Number of row
V<-colnames(simu.data_M1)  #Column names

Rcor_R=RobustCor(simu.data_M1, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_1_GT <- MRPC(simu.data_M1,suffStat_R,GV=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

#Inferred graph 

n<-nrow (simu.data_M2)    #Number of row
V<-colnames(simu.data_M1) #Column names

Rcor_R=RobustCor(simu.data_M2, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by MRPC using gaussCItest
MRPC.fit_2_GT <- MRPC(simu.data_M2,suffStat_R,GV=1,FDR=0.05,
indepTest ='gaussCItest', labels= V, verbose = TRUE)

#
Distance=aSHD(MRPC.fit_1_GT$graph,MRPC.fit_2_GT$graph,GV=1)

}

