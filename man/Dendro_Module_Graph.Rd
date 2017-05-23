\name{Dendro_Module_Graph}
\alias{Dendro_Module_Graph}
\title{
Plot Nodes Dendrogram with Modules and Graph
}
\description{
Clustering dendrogram of nodes, with dissimilarity based on topological overlap, together with assigned module colors and graphical visualization.
}
\usage{
Dendro_Module_Graph(Adj_directed,minModuleSize,NQ)
}

\arguments{
\item{Adj_directed}{
Adjacency matrix from directed graph
}
  \item{minModuleSize}{
Need to specify minimum module size.
}
\item{NQ}{
Need to specify for how many SNPs/indels/CNV/eQTL in your input data matrix. For example, If your dataset have one SNPs/indels/CNV/eQTL (first column), then NQ=1, if 2 SNPs/indels/CNV/eQTL (1st and 2nd Column), then NQ=2..and so on.
}

}

\value{
An object of \link{class} that containing the graph object as follows:
\itemize{
\item \code{obj}: object of class "graph" of the directed and undirected graph that was estimated.
\item \code{GroupMods}: list of modules with length of the nodes.
\item \code{GroupModsColors}: list of color of the modules with length of the nodes.
}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based PC (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{EdgesOrientation} for orientation rules determination for edges in MRPC; \link{SimulatedData} for simulated data generating function.
}
\examples{
## Model 0
n<-nrow (simu.data_M0)    #Number of row
V<-colnames(simu.data_M0) #Column names

Rcor_R=RobustCor(simu.data_M0, 0.005) #Robust correlation (Beta=0.005)
suffStat_R= list(C = Rcor_R$RR, n = n)

## Estimated graph by existing pc

## Estimated graph by MRPC using gaussCItest
MRPC.fit_0_GT <- MRPC(simu.data_M0,suffStat_R,NQ=1,FDR=0.05,
indepTest ='gaussCItest',labels=V,verbose = TRUE)

#Adjacency matrix from directed graph
Adj_directed=as(MRPC.fit_0_GT$obj@graph,"matrix")
#Plot of the graph
Dendro_Module_Graph(Adj_directed,minModuleSize = 1,NQ=1)
}
