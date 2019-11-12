\name{SimulateData}
\alias{SimulateData}

\title{
Simulate data under certain graphs
}
\description{
This function simulates data using linear models for several graphs: the five basic topologies and three topologies that are common in biology, namely the multi-parent graph, the star graph and the layered graph. See references for details.
}
\usage{
SimulateData(N, p, model, b0.1, b1.1, b1.2, b1.3, sd.1)
}

\arguments{

\item{N}{The number of observations.}
\item{p}{Population frequency of the reference allele. Real number between 0 to 1, which is the number of a particular allele is present.}
\item{model}{The model for which data will be simulated. For example, if you want to generate data for model 0 you would type 'model0' into the function.}
\item{b0.1}{Intercept of b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
\item{b1.1}{Slope of P1 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
\item{b1.2}{Slope of P2 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
\item{b1.3}{Slope of P3 for b0.1 + b1.1*P1 + b1.2*P2 + b1.3*P3, where P1, P2, and P3 are the parents of the corresponding node.}
\item{sd.1}{Standard deviation for corresponding data generated nodes.}

}

\details{
The first column of the input matrix is the genotype of the expression quantitative trait loci (eQTL)/Copy number variation (CNVs) and the remaining columns are the node expression data.
}

\value{
Matrix
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).
  
2. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899. 
}

\seealso{
\link{MRPC}; \link{SimulateDataNP}, which simulates data for a node with no parent; \link{SimulateData1P} for a node with one parent; \link{SimulateData2P} for a node with two parents.

}
\examples{

# Data pre-processing

# If you use only one genotype of the expression quantitative trait loci
# (eQTL)/Copy number variation (CNV), the 1st column of
# the input matrix will be #eQTL/CNV and the remaining
# columns are the gene expression data.


## Model 0
simu_data_M0 <- SimulateData(N = 10^3, 
                             p = 0.45,
                             'model0', 
                             b0.1 = 0,
                             b1.1 = 1, 
                             b1.2 = 1,
                             b1.3 = 1, 
                             sd.1 = 1)
                             
## Model 1
simu_data_M1 <- SimulateData(N = 10^3, 
                             p = 0.45,
                             'model1', 
                             b0.1 = 0,
                             b1.1 = 1, 
                             b1.2 = 1,
                             b1.3 = 1, 
                             sd.1 = 1)
                             
## Model 2
simu_data_M2 <- SimulateData(N = 10^3, 
                             p = 0.45,
                             'model2', 
                             b0.1 = 0,
                             b1.1 = 1, 
                             b1.2 = 1,
                             b1.3 = 1, 
                             sd.1 = 1)
                             
## Model 3
simu_data_M3 <- SimulateData(N = 10^3, 
                             p = 0.45,
                             'model3', 
                             b0.1 = 0,
                             b1.1 = 1, 
                             b1.2 = 1,
                             b1.3 = 1, 
                             sd.1 = 1)
                             
## Model 4
simu_data_M4 <- SimulateData(N = 10^3, 
                             p = 0.45,
                             'model4', 
                             b0.1 = 0,
                             b1.1 = 1, 
                             b1.2 = 1,
                             b1.3 = 1, 
                             sd.1 = 1)
                             
## Multiple Parent Model
simu_data_multiparent <- SimulateData(N = 10^3, 
                                      p = 0.45,
                                      'multiparent', 
                                      b0.1 = 0,
                                      b1.1 = 1, 
                                      b1.2 = 1,
                                      b1.3 = 1, 
                                      sd.1 = 1)

## Star Model
simu_data_starshaped <- SimulateData(N = 10^3, 
                                     p = 0.45,
                                     'starshaped', 
                                     b0.1 = 0,
                                     b1.1 = 1, 
                                     b1.2 = 1,
                                     b1.3 = 1, 
                                     sd.1 = 1)

## Layered Model
simu_data_layered <- SimulateData(N = 10^3, 
                                  p = 0.45,
                                  'layered', 
                                  b0.1 = 0,
                                  b1.1 = 1, 
                                  b1.2 = 1,
                                  b1.3 = 1, 
                                  sd.1 = 1)
}
