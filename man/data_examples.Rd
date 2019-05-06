\name{data_examples}
\alias{data_examples}

\title{
Example data under simple and complex models
}
\description{
Example data under the simple and complex graphs.  Data may be continuous or discrete.
}
\usage{
data(data_examples)
}

\details{
For each model, the graph and a simulated data matrix are available for both continuous and discrete data. 

\bold{For continuous data with genetic information}: 1000 samples in row and 6 variables in column. First two columns are the genetic variants and remaning columns are gene expression.

\bold{Continuous data without genetic information}: 1000 samples in row and 8 variables in column.

\bold{Discrete data with genetic information}: 1000 samples in row and 6 variables in column. First column is the genetic variant and remaning columns are the gene expression.

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
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{

# Continuous data with genetic information
# load the data
data("data_examples")  

# extract the sample size
n <- nrow(data_examples$simple$cont$withGV$data)  

# extract the node/column names
V <- colnames(data_examples$simple$cont$withGV$data)  

# calculate robust correlation matrix
Rcor_R <- RobustCor(data_examples$simple$cont$withGV$data,
                    Beta = 0.005) 

suffStat_R <- list(C = Rcor_R$RR, n = n)

data.mrpc.cont.withGV <- MRPC(data = data_examples$simple$cont$withGV$data,
                              suffStat = suffStat_R, GV = 2,
                              FDR = 0.05, indepTest = 'gaussCItest',
                              labels = V, verbose = TRUE)

par(mfrow = c(1, 2))
# plot the true graph
plot(data_examples$simple$cont$withGV$graph,
     main = "truth") 
# plot the inferred graph
plot(data.mrpc.cont.withGV,
     main = "inferred")    


# Continuous data without genetic information
# load the data
data("data_examples")    

# extract the sample size
n <- nrow(data_examples$simple$cont$withoutGV$data) 

# extract the node/column names
V <- colnames(data_examples$simple$cont$withoutGV$data)

# calculate robust correlation matrix
Rcor_R <- RobustCor(data_examples$simple$cont$withoutGV$data,
                    Beta = 0.005) 
                    
suffStat_R <- list(C = Rcor_R$RR, n = n)

data.mrpc.cont.withoutGV <- MRPC(data = data_examples$simple$cont$withoutGV$data,
                                 suffStat = suffStat_R, GV = 0,
                                 FDR = 0.05, indepTest = 'gaussCItest',
                                 labels = V, verbose = TRUE)

par(mfrow = c(1, 2))
# plot the true graph
plot(data_examples$simple$cont$withoutGV$graph,
     main = "truth") 
# plot the inferred graph
plot(data.mrpc.cont.withoutGV,
     main = "inferred") 

# Discrete data with genetic information
# load the data
data("data_examples")     

# extract the sample size
n <- nrow(data_examples$simple$disc$withGV$data)  

# extract the node/column names
V <- colnames(data_examples$simple$disc$withGV$data)

# calculate robust correlation matrix
Rcor_R <- RobustCor(data_examples$simple$disc$withGV$data,
                    Beta = 0.005)
                    
suffStat_R <- list(C = Rcor_R$RR, n = n)

data.mrpc.disc.withGV <- MRPC(data = data_examples$simple$disc$withGV$data,
                              suffStat = suffStat_R, GV = 1,
                              FDR = 0.05, indepTest = 'gaussCItest',
                              labels = V, verbose = TRUE)

par (mfrow = c(1, 2))
# plot the true graph
plot(data_examples$simple$disc$withGV$graph,
     main = "truth") 
# plot the inferred causal graph
plot(data.mrpc.disc.withGV,
     main = "inferred") 

# Discrete data without genetic information
# load the data
data("data_examples")    

# extract the sample size
n <- nrow (data_examples$simple$disc$withoutGV$data)  

# extract the node/column names
V <- colnames(data_examples$simple$disc$withoutGV$data)   

# calculate robust correlation matrix
Rcor_R <- RobustCor(data_examples$simple$disc$withoutGV$data,
                    Beta = 0.005) 
                    
suffStat_R <- list(C = Rcor_R$RR, n = n)

data.mrpc.disc.withoutGV <- MRPC(data = data_examples$simple$disc$withoutGV$data,
                                 suffStat = suffStat_R, GV = 1,
                                 FDR = 0.05, indepTest = 'gaussCItest',
                                 labels = V, verbose = TRUE)

par(mfrow = c(1, 2))
# plot the true graph
plot(data_examples$simple$disc$withoutGV$graph,
     main = "truth") 
# plot the inferred graph
plot(data.mrpc.disc.withoutGV,
     main = "inferred") 

# Continuous data with genetic information for complex model
# load the data
data("data_examples")   

# graph without clustering
plot(data_examples$complex$cont$withGV$graph) 

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
