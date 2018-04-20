\name{DendroModuleGraph }
\alias{DendroModuleGraph }
\title{
Visualization of Nodes by Dendrogram with Modules and Inferred Graph
}
\description{
Clustering dendrogram of nodes, with dissimilarity based on topological overlap, together with assigned module colors and graphical visualization.
}
\usage{
DendroModuleGraph (Adj_directed,minModuleSize,GV)
}

\arguments{
\item{Adj_directed}{
Adjacency matrix from directed graph
}
  \item{minModuleSize}{
Need to specify minimum module size.
}
\item{GV}{
Need to specify for how many SNPs/indels/CNV/eQTL in your input data matrix. For example, If your dataset have one SNPs/indels/CNV/eQTL (first column), then GV=1, if 2 SNPs/indels/CNV/eQTL (1st and 2nd Column), then GV=2..and so on.
}

}

\value{
An object of list of \link{class} that containing the graph object as follows:
\itemize{
\item \code{obj}: object of class "graph" of the directed and undirected graph that was estimated.
\item \code{GroupMods}: list of modules with length of the nodes.
\item \code{Adjmatrixdirected}: Adjacency matrix from directed graph that can help to extract (directed) edges along with their nodes.

}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based PC (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function; \link{EdgeOrientation} for orientation rules determination for edges in MRPC; \link{SimulatedData} for simulated data generating function.
}
\examples{

#Adjacency matrix from directed example graph
Adj_directed=as(ExampleMRPC$complex$cont$withGV$graph,"matrix")
#Plot of the graph
DendroModuleGraph(Adj_directed,minModuleSize = 5,GV=14)
}
