\name{DendroModuleGraph }
\alias{DendroModuleGraph }
\title{
Visualization of Nodes by Dendrogram with Modules and Inferred Graph
}
\description{
Clustering dendrogram of nodes with dissimilarity based on topological overlap, together with assigned module colors and graphical visualization.
}
\usage{
DendroModuleGraph (Adj_directed,minModuleSize,GV)
}

\arguments{
\item{Adj_directed}{
Adjacency matrix from directed graph
}
  \item{minModuleSize}{
Minimum module size.
}
\item{GV}{
Number of  SNPs/indels/CNV/eQTL in the input data matrix. For example, If the data has one SNPs/indels/CNV/eQTL (first column), then GV = 1, if 2 SNPs/indels/CNV/eQTL (1st and 2nd Column), then GV = 2, and so on.
}

}

\value{
A list containing the graph objects as follows:
\itemize{
\item \code{obj}: An object of class "graph" of the estimated graph.
\item \code{GroupMods}: A list of modules containing the number of nodes in each module.
\item \code{Adjmatrixdirected}: The adjacency matrix of the directed graph used to extract directed edges along with their nodes.

}
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based PC (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using the modified skeleton function; \link{EdgeOrientation} for orientation rules used to determine the edges in MRPC; \link{SimulatedData} for simulated data generating function.
}
\examples{

# Adjacency matrix from directed example graph
Adj_directed <- as(ExampleMRPC$complex$cont$withGV$graph,
                   "matrix")

# Plot of the graph
DendroModuleGraph(Adj_directed,
                  minModuleSize = 5,
                  GV = 14)
}
