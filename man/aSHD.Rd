\name{aSHD}
\alias{aSHD}

\title{
Adjusted Structural Hamming Distance (aSHD) 
}
\description{
The SHD, as is implemented in the R package pcalg (Kalisch et al., 2012) and bnlearn(Scutari, 2010), counts how many differences exist between two directed graphs. This distance is 1 if an edge exists in one graph but missing in the other, or if the direction of an edge is different in the two graphs. The larger this distance, the more different the two graphs are. We adjusted the SHD to reduce the penalty on the wrong direction of an edge to 0.5. For example, between two graphs V-->T1<--T2 and V-->T1-->T2, the SHD is 1 and our aSHD is 0.5.
}
\usage{
aSHD(g1, g2, GV,edge.presence, edge.direction)
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
\item{edge.presence}{
specify the weights for edge presence
}
  \item{edge.direction}{
specify the weights for edge direction
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
#Truth graph (V1-->T1-->T2-->T3)
tarmat_s1=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s1)=c("V1","T1","T2","T3")
rownames(tarmat_s1)=colnames(tarmat_s1)
#Create an adjacency matrix for truth
tarmat_s1[1,2]=1
tarmat_s1[2,3]=1
tarmat_s1[3,4]=1
#Graph object for truth
Truth<-as(tarmat_s1, "graphNEL")


#Inferred graph (V1-->T1<--T2-->T3)
tarmat_s2=matrix(0,nrow=4,ncol = 4)
colnames(tarmat_s2)=c("V1","T1","T2","T3")
rownames(tarmat_s2)=colnames(tarmat_s2)
#Create an adjacency matrix for truth
tarmat_s2[1,2]=1
tarmat_s2[3,2]=1
tarmat_s2[3,4]=1
#Graph objects for inferred
Inferred<-as(tarmat_s2, "graphNEL")
#
Distance=aSHD(Truth,Inferred,GV=1,edge.presence=1.0, edge.direction=0.5)

}

