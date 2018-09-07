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
We used the LOND (significance Levels based On Number of Discoveries) algorithm that controls FDR and mFDR in an online manner (Javanmard and Montanari, 2015). Where the significance level, alpha, is based on the total number of discoveries made so far. Which is similar to the algorithm called alpha-investing rules introduced by (Foster and Staine, 2007) to control only mFDR in an online manner.

}
\value{
The value of alpha.
}
\references{
[1] Javanmard and Montanari (March 5, 2015): On Online Control of False Discovery Rate. arXiv:150206197 [statME].

[2] Foster, D. P. and Stine, R. A. (2007): Alpha-investing: A procedure for sequential control of expected false discoveries. http://gosset.wharton.upenn.edu/research/edc.pdf,.

}
\author{
Md Bahadur Badsha <mdbadsha@uidaho.edu>, Evan A Martin <mart9986@vandals.uidaho.edu>, and Audrey Qiuyan Fu <audreyf@uidaho.edu>
}


\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function.
}
