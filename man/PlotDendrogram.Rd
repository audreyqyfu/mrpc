\name{PlotDendrogram}
\alias{PlotDendrogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clustering of nodes by dendrogram with modules
}
\description{
Clustering and visualization dendrogram of nodes with dissimilarity based on topological overlap, together with assigned modules colors.
}
\usage{
PlotDendrogram(Adj_directed,minModuleSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Adj_directed}{
Adjacency matrix from directed graph
}
\item{minModuleSize}{
Minimum module size.
}
}

\value{
A list containing the graph objects as follows:
\itemize{
\item \code{PlotDendrogramObj}: An object of class "graph" of the estimated graph.
\item \code{dynamicColors}: A list of colors with corresponding nodes.
\item \code{GroupMods}: Dynamic tree cut to identify modules whose expression profiles are very similar.
\item \code{GroupModsColors}: A table for number of nodes with corresponding colors.
\item \code{Adj_symmetric_matrix}: A symmetric matrix from ddjacency matrix of directed graph.
}

}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based PC (MRPC) algorithm.
}
\examples{
# Adjacency matrix from directed example graph
Adj_directed <- as(data_examples$complex$cont$withGV$graph,
                   "matrix")

# Plot of dendrogram with modules colors of nodes
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                  minModuleSize = 5)
}

