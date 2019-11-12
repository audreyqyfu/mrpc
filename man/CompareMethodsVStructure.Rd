\name{CompareMethodsVStructure}
\alias{CompareMethodsVStructure}

\title{
Comparison of inference accuracy with and without a v-structure
}
\description{
This is the example comparison of inference accuracy with and without a v-structure by recall and precision of five methods \link{MRPC}, \link{pc}, \link{pc.stable}, \link{mmpc}, and \link{mmhc}, across 1000 data sets. Signal values indicate the signal strength: the larger the value the stronger the signal. The sample size in each simulated data set is 1000. See Badsha and Fu, 2019 for the comparison of inference accuracy with sample sizes of 50, 200, and 500 .  
}
\usage{
CompareMethodsVStructure(N, signal, model, includeGV, ita)
}

\arguments{

\item{N}{Number of observation}
\item{signal}{
The coefficient of parent nodes in the linear model. For example, strong = 1.0, moderate = 0.5, and weak = 0.2 signal strengths.
}
  \item{model}{
The model data will be simulated from. For example, if you want comparison of inference accuracy without a v-structure you would type 'model1' into the function. If you want comparison of inference accuracy with a v-structure you would type 'model2' into the function.
}
\item{includeGV}{
If TRUE, include edges involving genetic variants (GV) when comparing the true and inferred graphs. If FALSE, exclude edges involving genetic variants (GV) when comparing the true and inferred graphs.
}
\item{ita}{Number of different data sets}

}
\details{
The output is a matrix, where the rows are the five methods MRPC, pc, pc.stable, mmpc, and mmhc, and the columns are mean Recall, sd Recall, mean Precision, and sd Precision, respectively.
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

\seealso{
\link{RecallPrecision}: Performance Evaluation by Recall and Precision in MRPC.
}
\examples{
# Comparison of inference accuracy without a v-structure 
# across 10 data sets where the sample size in each simulated 
# data set is set to 100 and signal = 1.0. 
# We excluded edges involving genetic variants.

Result1 <- CompareMethodsVStructure(N = 100, 
                                    signal = 1.0,
                                    'model1', 
                                    includeGV = FALSE,
                                    ita = 10)

# Comparison of inference accuracy with a v-structure
# across 10 data sets where the sample size in each simulated 
# data set is set to 100 and signal = 1.0. 
# We excluded edges involving genetic variants.

Result2 <- CompareMethodsVStructure(N = 100, 
                                    signal = 1.0,
                                    'model2',
                                    includeGV = FALSE, 
                                    ita = 10)

}
