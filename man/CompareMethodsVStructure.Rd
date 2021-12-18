\name{CompareMethodsVStructure}
\alias{CompareMethodsVStructure}

\title{
Comparison of inference accuracy of different methods on data with and without a v-structure
}
\description{
This function compares inference accuracy on graphs with and without a v-structure in terms of recall and precision by six methods \link{MRPC}, \link{pc}, \link{pc.stable}, \link{mmpc}, \link{mmhc}, and \link{hc}, across multiple data sets. See details in  Badsha and Fu, 2019 and Badsha et al., 2021.
}
\usage{
CompareMethodsVStructure(N, signal, model, includeGV, ita)
}

\arguments{

\item{N}{Number of observations.}
\item{signal}{
The coefficient of parent nodes in the linear model. For example, strong = 1.0, moderate = 0.5, and weak = 0.2.
}
  \item{model}{
The graph from which the data is generated.  Specifically, two graphs are considered here: 'model 1' (V1->T1->T2), which does not contain a v-structure, and 'model 2' (V1->T1<-T2), which is a v-structure.
}
\item{includeGV}{
If TRUE, include edges involving genetic variants (GVs) when comparing the true and inferred graphs. If FALSE, exclude such edges.  
}
\item{ita}{Number of independent data sets to simulate.}

}
\details{
The output is a matrix, where the rows are the six methods: MRPC, pc, pc.stable, mmpc, mmhc, and hc, and the columns are the mean of recall, sd of recall, mean of precision, and sd of precision, respectively. Mean and sd are calculated across all the simulated data sets. For methods from the bnlearn package (pc.stable, mmpc, mmhc, and hc), we apply the blacklist argument to exclude edges pointing at the genetic variant, and therefore evaluate recall and precision including the edges involving these edges (i.e., includeGV = TRUE).
}
\value{
Matrix
}
\references{
1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10(460).

2. Badsha MB, Martin EA and Fu AQ (2021). An R Package for Inference of Causal Graphs. Frontiers in Genetics, 10(651812).
  
3. Kalisch M, Machler M, Colombo D, Maathuis MH and Buhlmann P (2012). Causal Inference Using Graphical Models with the R Package pcalg. Journal of Statistical Software, 47, 26.

4. Scutari M (2010). Learning Bayesian Networks with the bnlearn R Package. Journal of Statistical Software, 35(3), 1-22.
}

\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}

\seealso{
\link{RecallPrecision}: Performance evaluation in terms of recall and precision.
}
\examples{
  
 \donttest{
# For demonstration purposes, only 10 data sets 
# with a sample size of 100 and signal is 1, are simulated here.

# Comparison of inference accuracy on model 1 without a v-structure 
Result1 <- CompareMethodsVStructure(N = 100, 
                                    signal = 1.0,
                                    'model1', 
                                    includeGV = TRUE,
                                    ita = 10)

# Comparison of inference accuracy with a v-structure
Result2 <- CompareMethodsVStructure(N = 100, 
                                    signal = 1.0,
                                    'model2',
                                    includeGV = TRUE,
                                    ita = 10)

}
  }
