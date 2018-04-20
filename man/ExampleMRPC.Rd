\name{ExampleMRPC}
\alias{ExampleMRPC}

\title{
Example Data for Simple and Complex Model
}
\description{
Example of the numeric data (Continuous and Discrete) matrix and components of a graph
}
\usage{
data(ExampleMRPC)
}

\details{
Contain a numeric data matrix and a graph for both continuous and discrete. 

\bold{For continuous data with genetic information}: 1000 samples in row and 6 variables in column. First two column is the genetic variants and remaning columns are the genes expression.

\bold{Continuous data without genetic information}: 1000 samples in row and 8 variables in column.

\bold{Discrete data with genetic information}: 1000 samples in row and 6 variables in column. First column is the genetic variants and remaning columns are the genes expression.

\bold{Discrete data without genetic information}: 1000 samples in row and 5 variables in column.

\bold{Continuous data with genetic information for complex model}: 1000 samples in row and 22 variables in column. First 14 column is the genetic variants and remaning columns are the genes expression.   
}
\value{
A list that containing the numeric data matrix and components of a graph.
\itemize{
\item \code{simple}: Simple model.
\item \code{complex}: Complex model.
\item \code{cont}: Continuous.
\item \code{disc}: Discrete.
\item \code{withGV}: With genetic information.
\item \code{withoutGV}: Without genetic information.
\item \code{data}: Data matrix.
\item \code{graph}: Components of a graph.
}
}

\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
library(pcalg) #For graph plot
#Continuous data with genetic information
# load the data
data ("ExampleMRPC")  
# extract the sample size
n <- nrow (ExampleMRPC$simple$cont$withGV$data)    
# extract the node/column names
V <- colnames(ExampleMRPC$simple$cont$withGV$data)    
# calculate robust correlation matrix
Rcor_R <- RobustCor(ExampleMRPC$simple$cont$withGV$data,
Beta=0.005) 

suffStat_R <- list(C = Rcor_R$RR, n = n)
data.mrpc.cont.withGV <- 
MRPC(data=ExampleMRPC$simple$cont$withGV$data,
suffStat=suffStat_R,GV=2, FDR=0.05, indepTest =
'gaussCItest',labels=V, verbose = TRUE)


par (mfrow=c(1,2))
# plot the truth graph
plot (ExampleMRPC$simple$cont$withGV$graph, main="truth") 
# plot the inferred causal graph
plot (data.mrpc.cont.withGV$graph, main="inferred")    


#Continuous data without genetic information
 # load the data
data ("ExampleMRPC")    
 # extract the sample size
n <- nrow (ExampleMRPC$simple$cont$withoutGV$data) 
# extract the node/column names
V <- colnames(ExampleMRPC$simple$cont$withoutGV$data)   # calculate robust correlation matrix
Rcor_R <- RobustCor(ExampleMRPC$simple$cont$withoutGV$data,
Beta=0.005) 
suffStat_R <- list(C = Rcor_R$RR, n = n)
data.mrpc.cont.withoutGV <- 
MRPC(data=ExampleMRPC$simple$cont$withoutGV$data,
suffStat=suffStat_R, GV=0,FDR=0.05, 
indepTest='gaussCItest', labels=V, verbose = TRUE)


par (mfrow=c(1,2))
# plot the truth graph
plot (ExampleMRPC$simple$cont$withoutGV$graph, main="truth") 
# plot the inferred causal graph
plot (data.mrpc.cont.withoutGV$graph, main="inferred")    
#Discrete data with genetic information
# load the data
data ("ExampleMRPC")     
# extract the sample size
n <- nrow (ExampleMRPC$simple$disc$withGV$data)  
 # extract the node/column names
V <- colnames(ExampleMRPC$simple$disc$withGV$data)
# calculate robust correlation matrix
Rcor_R <- RobustCor(ExampleMRPC$simple$disc$withGV$data, Beta=0.005) 
suffStat_R <- list(C = Rcor_R$RR, n = n)
data.mrpc.disc.withGV <-
MRPC(data=ExampleMRPC$simple$disc$withGV$data, 
suffStat=suffStat_R,GV=1,FDR=0.05, indepTest = 
'gaussCItest',labels=V, verbose = TRUE)

par (mfrow=c(1,2))
# plot the truth graph
plot (ExampleMRPC$simple$disc$withGV$graph, main="truth") 
# plot the inferred causal graph
plot (data.mrpc.disc.withGV$graph, main="inferred") 

#Discrete data without genetic information
 # load the data
data ("ExampleMRPC")    
 # extract the sample size
n <- nrow (ExampleMRPC$simple$disc$withoutGV$data)  
# extract the node/column names
V <- colnames(ExampleMRPC$simple$disc$withoutGV$data)   
# calculate robust correlation matrix
Rcor_R <- RobustCor(ExampleMRPC$simple$disc$withoutGV$data, Beta=0.005) 
suffStat_R <- list(C = Rcor_R$RR, n = n)
data.mrpc.disc.withoutGV <- 
MRPC(data=ExampleMRPC$simple$disc$withoutGV$data, 
suffStat=suffStat_R, GV=1,FDR=0.05, 
indepTest='gaussCItest', labels=V, verbose = TRUE)

par (mfrow=c(1,2))
# plot the truth graph
plot (ExampleMRPC$simple$disc$withoutGV$graph, main="truth") 
# plot the inferred causal graph
plot (data.mrpc.disc.withoutGV$graph, main="inferred") 

#Continuous data with genetic information for complex model
# load the data
data ("ExampleMRPC")   
# graph without clustering
plot(ExampleMRPC$complex$cont$withGV$graph) 
Adj_directed <- 
as(ExampleMRPC$complex$cont$withGV$graph,"matrix")

# dendrogram and graph with clustering
DendroModuleGraph(Adj_directed, minModuleSize = 5, GV=14) 
}
