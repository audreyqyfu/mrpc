\name{InferenceAccuracyExample}
\alias{InferenceAccuracyExample}

\title{
Comparison of inference accuracy with and without a v-structure
}
\description{
This is the example comparison of inference accuracy with and without a v-structure by recall and precision of three methods (MPRC, mmhc and pc) across 1000 data sets. Signal values of indicate the signal strength: the larger the stronger the signal.The sample size in each simulated data set is set to 1000. The comparison of inference accuracy for the sample size 50, 200 and 500, see Badsha and Fu, 2018.  
}
\usage{
InferenceAccuracyExample(N=1000,signal,model,ita=1000)
}

\arguments{
\item{N}{Number of observation}
\item{signal}{
Which is the coefficient of parent nodes in the linear model. For example, three different (strong=1.0, moderate=0.5 and weak=0.2) signal strengths.
}
  \item{model}{
To specify the model you want. For example, if you want comparison of inference accuracy without a v-structure you would type 'model1' into the function. If you want comparison of inference accuracy with a v-structure you would type 'model2' into the function.
}
\item{ita}{Number of different data sets}
}
\details{
The output as a matrix, where rows are the three methods (MRPC,mmhc and pc), and colums are mean Recall, sd Recall, mean Precision and sd Precision, respectively.
}
\value{
Matrix
}
\references{
1. Badsha MB, Fu AQ (2018). Learning causal biological networks with the principle of Mendelian randomization." bioRxiv. doi:10.1101/171348.
}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}

\seealso{
\link{Recall_Precision}: Performance Evaluation by Recall and Precision in MRPC.
}
\examples{
#Comparison of inference accuracy without a v-structure 
#across 1000 data sets and sample size in each simulated 
#data set is set to 1000. 
Result1=InferenceAccuracyExample(N=1000,signal=1.0,'model1',ita=1000)
#Comparison of inference accuracy with a v-structure
#across 1000 data sets and sample size in each simulated 
#data set is set to 1000. 
Result2=InferenceAccuracyExample(N=1000,signal=1.0,'model2',ita=1000)

}
