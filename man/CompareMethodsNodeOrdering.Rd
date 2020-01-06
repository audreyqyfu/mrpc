\name{CompareMethodsNodeOrdering}
\alias{CompareMethodsNodeOrdering}

\title{Comparison of inference accuracy using the same data but with different node orderings.}

\description{
Investigate the performance of five methods on the same data but with different node orderings: \link{MRPC} (Badsha and Fu, 2019; Badsha et al., 2018), \link{pc}, implemented in \code{pcalg} (Kalisch et al., 2012), and \link{pc.stable}, \link{mmpc}, and \link{mmhc}, the last three all implemented in \code{bnlearn} (Scutari, 2010). See details in Badsha et al., 2018.
}
\usage{
CompareMethodsNodeOrdering(N, signal, model, n_data, n_nodeordering)
}

\arguments{
  \item{N}{The number of observations.}
  \item{signal}{The signal strength which is the coefficient of the parent nodes in the linear model.}
    \item{model}{Either 'truth1' or 'truth2' to specify the model to generate data from.}
  \item{n_data}{The number of independent data sets to generate.}
  \item{n_nodeordering}{The number of times to reorder the nodes.}
}
\details{
We generated different data sets from each of the two graphs (V1-->T1-->T2-->T3 and V1-->T1<--T2-->T3), where V1 is the genetic varitant node and T1, T2 and T3 are the phenotype nodes. For each data set, we reordered the columns of the data matrix with different permutations of the T nodes, and thus generated new permuted data sets. We then applied all the methods to each of the data sets (restricting edge direction wherever necessary and possible), obtained different inferred graphs, and counted the number of unique graphs among the inferred graphs for each data set. 

The output is the number of unique graphs inferred by each method across all permutations.The columns are indicates which methods (MRPC, pc, pc.stable, mmpc and mmhc), and each rows are indicates the independent data sets under different node orderings. See details in Badsha et al., 2018.

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
# We will generate 10 different data sets from truth1 (V1-->T1-->T2-->T3) 
# with signal = 1.0, N = 100 and 6 different T nodes orderings. Therefore, we
# will get 10 by 5 output matrix, where the rows are 10 independent data 
# sets and columns are for five methods MRPC, pc, pc.stable, mmpc 
# and mmhc respectively.

library(MRPC)    #  MRPC
library(pcalg)   #  pc
library(bnlearn) #  pc.stable, mmpc, and mmhc

# Run
Output <- CompareMethodsNodeOrdering(N = 100,
                                     signal = 1.0, 
                                     model = 'truth1',
                                     n_data = 10,
                                     n_nodeordering = 6)
}

