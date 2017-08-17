\name{ExampleMRPC}
\alias{ExampleMRPC}

\title{
Gaussian Example Data for Simple and Complex Model
}
\description{
Example of the numeric data matrix and components of a graph
}
\usage{
data(ExampleMRPC)
}

\details{
Contain a numeric data matrix and a graph. The data with 1000 samples in row and five variables in column for simple model. First two column is the genetic variants and remaning columns are the genes expression.On the other hand, for complex model 1000 samples in row and 22 variables in column. First 14 column is the genetic variants and remaning columns are the genes expression.   
}
\value{
A list that containing the numeric data matrix and components of a graph.
\itemize{
\item \code{data}: Data matrix.
\item \code{graph}: Components of a graph.
}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
library(pcalg) #For graph plot
data ("ExampleMRPC ")
#If you need the data and a graph for simple model invloving 5 nodes
Data=ExampleMRPC$simple$data
plot(ExampleMRPC$simple$graph)
#If you need the data and a graph for simple model invloving 22 nodes
Data=ExampleMRPC$complex$data
plot(ExampleMRPC$complex$graph)
}
