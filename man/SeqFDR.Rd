\name{SeqFDR}
\alias{SeqFDR}

\title{
Sequential FDR
}
\description{
Sequential FDR method that controls the FDR and mFDR in an online manner.
}
\usage{
SeqFDR(m, FDR, a=2, R)
}

\arguments{

  \item{m}{
The number of current the test.
}
  \item{FDR}{
FDR level.
}
  \item{a}{
A constant.
}
  \item{R}{
All of the decisions from the tests that have already been performed.
}
}
\details{
We used the LOND (significance Levels based On Number of Discoveries) algorithm that controls FDR in an online manner (Javanmard and Montanari, 2015). The significance level (i.e., the type I error rate) is dynamic and based on the total number of discoveries made so far. 

}
\value{
The value of alpha.
}
\references{
1. Javanmard A and Montanari A (2015). On Online Control of False Discovery Rate. arXiv:150206197.
}
\author{
Md Bahadur Badsha (mbbadshar@gmail.com)
}


\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function.
}
