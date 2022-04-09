\name{PlotDendrogram}
\alias{PlotDendrogram}

\title{
Plot a dendrogram and display node groups in colored modules
}
\description{
Generate a dendrogram of nodes with dissimilarity based on topological overlap, and group nodes into modules indicated by colors.
}
\usage{
PlotDendrogram(Adj_directed, minModuleSize, groupLabels = " ", 
               dendroLabels = FALSE, hclustHang = 0.03, 
               dendroAddGuide = FALSE, dendroGuideHang = 0.05, 
               dendroMain = "Dendrogram with modules of nodes in colors", ...)
}

\arguments{
  \item{Adj_directed}{
Adjacency matrix from directed graph
}
\item{minModuleSize}{
Minimum module size.
}
\item{groupLabels}{
Argument for \link{plotDendroAndColors}. Labels for the colorings given in colors. The labels will be printed to the left of the color rows in the plot. 
}
\item{dendroLabels}{
Argument for \link{plotDendroAndColors}. Dendrogram labels. 
}
\item{hclustHang}{
Argument \code{hang} for \link{plot.hclust}. The fraction of the plot height by which labels should hang below the rest of the plot. 
}
\item{dendroAddGuide}{
Argument \code{addGuide} for \link{plotDendroAndColors}. Logical: should vertical "guide lines" be added to the dendrogram plot? The lines make it easier to identify color codes with individual samples. 
}
\item{dendroGuideHang}{
Argument \code{guideHang} for \link{plotDendroAndColors}. The fraction of the dendrogram height to leave between the top end of the guide line and the dendrogram merge height.  
}
\item{dendroMain}{
Argument \code{main} for \link{plot.hclust}. Title of the plot.  
}
\item{...}{
Additional plotting arguments for \link{plotDendroAndColors} and \link{plot.hclust}.  
}
}

\value{
A list containing the graph objects as follows:
\itemize{
\item \code{PlotDendrogramObj}: An object of class "graph" of the estimated graph.
\item \code{dynamicColors}: A list of colors with corresponding nodes.
\item \code{GroupMods}: Dynamic tree cut to identify modules whose phenotype profiles are very similar.
\item \code{GroupModsColors}: A table for number of nodes with corresponding colors.
\item \code{Adj_symmetric_matrix}: A symmetric matrix from ddjacency matrix of directed graph.
}

}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\seealso{
\link{MRPC}.
}
\examples{
\dontrun{
# Adjacency matrix from directed example graph
Adj_directed <- as(data_examples$complex$cont$withGV$graph,
                   "matrix")

# Plot of dendrogram with modules colors of nodes
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                                    minModuleSize = 5)
                                    
}
}

