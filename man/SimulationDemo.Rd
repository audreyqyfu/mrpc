\name{SimulationDemo}
\alias{SimulationDemo}

\title{Same data but different node ordering
}
\description{
Investigate the impact of same data but different node ordering across three methods: MRPC (Batsha et al., 2018), the mmhc method (implemented in bnlearn; Scutari, 2010) and the pc method (implemented in pcalg; Kalisch et al., 2012). 
}
\usage{
SimulationDemo(N, model,signal,n_data,n_nodeordering)
}

\arguments{
  \item{N}{Number of observation}
  \item{model}{To specify the truth you want that shown in table 1 and 2. For example, if you want to truth 1 you need type 'truth1' into the function.}
  \item{signal}{signal strengths (which is the coefficient of parent nodes in the linear model}
  \item{n_data}{how many data you want to generate for testing the results}
  \item{n_nodeordering}{how many times you want to order the node for testing the results}
}
\details{
The output results as a matrix each row indicates the data sets and col's are the different node ordering by three methods MRPC, mmhc and pc respectively. 
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
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\examples{
#We want to generate 2 different data sets for truth1 with
#signal=1,N=1000 and 6 different node ordering. Therefore,
#you will get 2 by 18 output matrix, where first and second row 
#are the data 1 and 2 respectively. Column 1:6 is the result of 
#MRPC with node #order 1,2...6,then column 7:12 is for mmhc and 
#column 13:18 is for pc, #respectively

library(bnlearn) # for mmhc
library(pcalg)   #for pc
#Run
Output=SimulationDemo(N=1000,'truth1',signal=1,n_data=2, n_nodeordering=6)

}

