\name{data_examples}
\alias{data_examples}

\title{
Example data under simple and complex models
}
\description{
Example data under the simple and complex graphs. Data may be continuous or discrete.
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
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\examples{

# Continuous data with genetic varitant (GV)
# load the data
data("data_examples")  
data <- data_examples$simple$cont$withGV$data
# Extract the sample size
n <- nrow(data)  

# Extract the node/column names
V <- colnames(data)  

# Calculate Pearson correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Infer the graph by MRPC
data.mrpc.cont.withGV <- MRPC(data = data,
                              suffStat = suffStat_C, 
                              GV = 2,
                              FDR = 0.05, 
                              indepTest = 'gaussCItest',
                              labels = V,
                              FDRcontrol = 'LOND',
                              verbose = FALSE)
# Plot the results
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
data <- data_examples$simple$cont$withoutGV$data
  
# Extract the sample size
n <- nrow(data) 

# Extract the node/column names
V <- colnames(data)

# Calculate Pearson correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Infer the graph by MRPC
data.mrpc.cont.withoutGV <- MRPC(data = data,
                                 suffStat = suffStat_C, 
                                 GV = 0,
                                 FDR = 0.05,
                                 indepTest = 'gaussCItest',
                                 labels = V, 
                                 FDRcontrol = 'LOND',
                                 verbose = FALSE)
# Plot the results
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
data <- data_examples$simple$disc$withGV$data
  
# Extract the sample size
n <- nrow(data)  

# Extract the node/column names
V <- colnames(data)

suffStat_C <- list (dm = data, adaptDF = FALSE, n.min = 1000)

# Infer the graph by MRPC
data.mrpc.disc.withGV <- MRPC(data = data,
                              suffStat = suffStat_C, 
                              GV = 1,
                              FDR = 0.05, 
                              indepTest = 'disCItest',
                              labels = V, 
                              FDRcontrol = 'LOND',
                              verbose = FALSE)
# Plot the results
par (mfrow = c(1, 2))
# plot the true graph
plot(data_examples$simple$disc$withGV$graph,
     main = "truth") 
# Plot the inferred causal graph
plot(data.mrpc.disc.withGV,
     main = "inferred") 

# Discrete data without genetic information
# load the data
data("data_examples")    
data <- data_examples$simple$disc$withoutGV$data
  
# Extract the sample size
n <- nrow (data)  

# Extract the node/column names
V <- colnames(data)   

suffStat_C <- list (dm = data, adaptDF = FALSE, n.min = 1000)

# Infer the graph by MRPC
data.mrpc.disc.withoutGV <- MRPC(data = data,
                                 suffStat = suffStat_C, 
                                 GV = 1,
                                 FDR = 0.05, 
                                 indepTest = 'disCItest',
                                 labels = V, 
                                 FDRcontrol = 'LOND',
                                 verbose = FALSE)
# Plot the results
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

# Graph without clustering
plot(data_examples$complex$cont$withGV$graph) 

# Adjacency matrix from directed example graph
Adj_directed <- as(data_examples$complex$cont$withGV$graph,
                   "matrix")

# Plot of dendrogram with modules colors of nodes
PlotDendrogramObj <- PlotDendrogram(Adj_directed,
                                    minModuleSize = 5)
                  
# Visualization of inferred graph with modules colors
PlotGraphWithModulesObj <- PlotGraphWithModules(Adj_directed,
                                                PlotDendrogramObj,
                                                GV = 14,
                                                node.size = 8,
                                                arrow.size = 5,
                                                label.size = 3,
                                                alpha = 1) 
# plot 
plot(PlotGraphWithModulesObj)

# Run MRPC on the complex data set with ADDIS as the FDR control method.
data <- data_examples$complex$cont$withGV$data
n <- nrow (data)        # Number of rows
V <- colnames(data)     # Column names

# Calculate Pearson correlation
suffStat_C <- list(C = cor(data),
                   n = n)

# Infer the graph by MRPC
MRPC.addis <- MRPC(data,
                   suffStat = suffStat_C,
                   GV = 14,
                   FDR = 0.05,
                   indepTest = 'gaussCItest',
                   labels = V,
                   FDRcontrol = 'ADDIS',
                   tau = 0.5,
                   lambda = 0.25,
                   verbose = FALSE)
                   
# Plot the true and inferred graphs.
par(mfrow = c(1, 2))
plot(data_examples$complex$cont$withGV$graph,
     main = 'True graph')
plot(MRPC.addis,
     main = 'Inferred graph')
}
