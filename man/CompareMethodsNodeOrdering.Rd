\name{CompareMethodsNodeOrdering}
\alias{CompareMethodsNodeOrdering}

\title{Comparision of inference accuracy using the same data but with different node orderings.}

\description{
Investigate the impact of the same data but with different node orderings across three methods: \link{MRPC} (Badsha and Fu, 2019; Badsha et al., 2018), \link{mmhc} implemented in \code{bnlearn} (Scutari, 2010) and \link{pc} implemented in \code{pcalg} (Kalisch et al., 2012). 
}
\usage{
CompareMethodsNodeOrdering(N, model, signal, n_data, n_nodeordering)
}

\arguments{
  \item{N}{The number of observation.}
  \item{model}{Either 'truth1' or 'truth2' to specify the model you want to generate data from.}
  \item{signal}{The signal strength which is the coefficient of the parent nodes in the linear model.}
  \item{n_data}{The number of data sets you want to generate for testing the results of the three methods.}
  \item{n_nodeordering}{The number of times you want to reorder the nodes for testing the results of the three methods.}
}
\details{
The output is a matrix where each row indicates the data sets and the columns are the results of the different node orderings by five methods MRPC, pc, pc.stable, mmpc, and mmhc respectively. 
}
\value{
Matrix
}
\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Badsha MB, Martin EA and Fu AQ (2018). MRPC: An R package for accurate inference of causal graphs. arXiv:1806.01899.

3. Kalisch M, Machler M, Colombo D, Maathuis MH and Buhlmann P (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47, 26.

4. Scutari M (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.


}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
# We want to generate 2 different data sets for truth1 (V1-->T1-->T2-->T3) 
# with signal = 1, N = 100 and 6 different node orderings. Therefore, we will 
# get 2 by 30 output matrix, where first and second row are the data 1 and 2 respectively.
# Columns 1:6, 7:12, 13:18, 19:24 and 25:30 are for MRPC, pc, pc.stable, mmpc and
# mmhc respectively with node order 1, 2, ..., 6.

library(MRPC)    #  MRPC
library(pcalg)   #  pc
library(bnlearn) #  pc.stable,mmpc and mmhc

# Run
Output <- CompareMethodsNodeOrdering(N = 100, 
                                    'truth1',
                                     signal = 1, 
                                     n_data = 2,
                                     n_nodeordering = 6)
}

