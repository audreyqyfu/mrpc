\name{PlotGraphWithModules}
\alias{PlotGraphWithModules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of inferred graph with modules colors 
}
\description{
Visualization of inferred graph with modules colors from clustering dendrogram of nodes \link{PlotDendrogram}.
}
\usage{
PlotGraphWithModules(Adj_directed,PlotDendrogramObj,GV=GV,node.size=8,arrow.size = 5,
label.size = 3,alpha = 1,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Adj_directed}{
Adjacency matrix from directed graph
}
\item{PlotDendrogramObj}{
Plot objects from \link{PlotDendrogram}
}
\item{GV}{
The number of genetic variants (SNPs/indels/CNV/eQTL) in the input data matrix. For example, if the data has one SNPs/indels/CNV/eQTL, first column, then GV = 1, if 2 SNPs/indels/CNV/eQTL, 1st and 2nd Column, then GV = 2, if no GV then GV = 0, and so on.
}
\item{node.size}{
The size of the nodes. Defaults to 8.
}
\item{arrow.size}{
The size of the arrows for directed network edges, in points. Defaults to 5.
}
\item{label.size}{
The size of the node labels, in points, as a numeric value, a vector of numeric values, or as a vertex attribute containing numeric values. Defaults to 3.

}
\item{alpha}{
The level of transparency of the edges and nodes. Defaults to 1 (no transparency).
}
\item{...}{
Other arguments passed to the function.
}

}

\value{
\itemize{
\item \code{PlotGraphWithModulesObj}: An object of class "graph" of the estimated graph.
}

}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{PlotDendrogram} visualization of nodes by dendrogram with modules.
}
\examples{
# Adjacency matrix from directed example graph
Adj_directed <- as(data_examples$complex$cont$withGV$graph,
                   "matrix")

# Plot of dendrogram with modules colors of nodes
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                  minModuleSize = 5)
                  
#Visualization of inferred graph with modules colors
PlotGraphWithModulesObj<-PlotGraphWithModules(Adj_directed,PlotDendrogramObj,GV=14,node.size=8
,arrow.size = 5,label.size = 3,alpha = 1) 
#plot 
plot(PlotGraphWithModulesObj) 
}

