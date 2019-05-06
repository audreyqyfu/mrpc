\name{CompareMethodsNodeOrdering}
\alias{CompareMethodsNodeOrdering}

\title{Comparision of inference accuracy using the same data but with different node orderings.}

\description{
Investigate the impact of the same data but with different node orderings across three methods: \link{MRPC} (Badsha et al., 2018), \link{mmhc} implemented in \code{bnlearn} (Scutari, 2010) and \link{pc} implemented in \code{pcalg} (Kalisch et al., 2012). 
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
The output is a matrix where each row indicates the data sets and the columns are the results of the different node orderings by three methods MRPC, mmhc, and pc respectively. 
}
\value{
Matrix
}
\references{
1.Md. Bahadur Badsha, Audrey Qiuyan Fu: Learning causal biological networks with generalized Mendelian randomization.bioRxiv (2017), 171348. 

2.Marco Scutari (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.

3.Kalisch, M., Machler, M., Colombo, D., Maathuis, M.H. & Buhlmann, P. Causal Inference Using Graphical Models with the R Package pcalg. J. Stat.Softw. 47, 26 (2012).
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
# We want to generate 2 different data sets for truth1 with
# signal = 1, N = 100 and 6 different node orderings. Therefore,
# you will get 2 by 18 output matrix, where first and second row 
# are the data 1 and 2 respectively. Column 1:6 is the result of 
# MRPC with node order 1, 2, ..., 6, then column 7:12 is for mmhc and 
# column 13:18 is for pc.

library(bnlearn) # for mmhc
library(pcalg)   # for pc

# Run
Output <- CompareMethodsNodeOrdering(N = 100, 'truth1',
                         signal = 1, n_data = 2,
                         n_nodeordering = 6)

}

