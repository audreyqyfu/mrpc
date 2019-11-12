\name{PlotGraphWithModules}
\alias{PlotGraphWithModules}

\title{
Plot a graph with nodes in modules indicated by colors 
}
\description{
Visualization of a graph with nodes in modules inferred from the clustering dendrogram by \link{PlotDendrogram}.
}
\usage{
PlotGraphWithModules(Adj_directed, PlotDendrogramObj, 
                     GV = GV, node.size = 8, arrow.size = 5,
                     label.size = 3,alpha = 1,...)

}

\arguments{
  \item{Adj_directed}{
Adjacency matrix of a graph.
}
\item{PlotDendrogramObj}{
The graphical objects from \link{PlotDendrogram}.
}
\item{GV}{
The number of genetic variants (SNPs/indels/CNVs/eQTL) in the input data matrix. For example, if the data has one SNPs/indels/CNV/eQTL in the first column, then GV = 1, if 2 SNPs/indels/CNVs/eQTL in the 1st and 2nd Column, then GV = 2, and so on. If no GV then GV = 0.
}
\item{node.size}{
The size of the nodes in the graph. Defaults to 8.
}
\item{arrow.size}{
The size of the arrows for directed network edges, in points. Defaults to 5.
}
\item{label.size}{
The size of the node labels in points, as a numeric value, a vector of numeric values, or as a vertex attribute containing numeric values. Defaults to 3.

}
\item{alpha}{
The level of transparency of the edges and nodes. Defaults to 1 (no transparency).
}
\item{...}{
Other arguments passed to \link{ggnet2}.
}

}

\value{
\itemize{
\item \code{PlotGraphWithModulesObj}: An object of class "graph" of the graph.
}

}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{PlotDendrogram}
}
\examples{
# Adjacency matrix from a graph in the example
Adj_directed <- as(data_examples$complex$cont$withGV$graph,
                   "matrix")

# A clustering dendrogram with nodes grouped in colored modules
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                                    minModuleSize = 5)
                  
# A graph object with nodes in modules
PlotGraphWithModulesObj <- PlotGraphWithModules(Adj_directed,
                                                PlotDendrogramObj,
                                                GV=14,
                                                node.size=8,
                                                arrow.size = 5,
                                                label.size = 3,
                                                alpha = 1) 

# Plot the graph with nodes in different colors
plot(PlotGraphWithModulesObj) 
}

