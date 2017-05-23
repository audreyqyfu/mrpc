\name{SeqFDR}
\alias{SeqFDR}

\title{
Sequential FDR
}
\description{
Sequential FDR methd that control FDR and mFDR as online manner 
}
\usage{
SeqFDR(m, FDR, a=2, R)
}

\arguments{

  \item{m}{
Number of current test
}
  \item{FDR}{
FDR level
}
  \item{a}{
Constant value =2
}

  \item{R}{
Rejection number
}
}
\details{
We used LOND (significance Levels based On Number of Discoveries) algorithm that control FDR and mFDR on online manner (Javanmard and Montanari, 2015). Where, the significance levels (alpha) based on total number of discoveries made so far. Which is the similar algorithm called alpha-investing rules introduces by Foster and Staine, 2007 to control only mFDR in online manner.

}
\value{
Integer alpha value
}
\references{
[1] Javanmard and Montanari (March 5, 2015): On Online Control of False Discovery Rate. arXiv:150206197 [statME].

[2] Foster, D. P. and Stine, R. A. (2007): Alpha-investing: A procedure for sequential control of expected false discoveries. http://gosset.wharton.upenn.edu/research/edc.pdf,.

}
\author{
Md. Bahadur Badsha (mdbadsha@uidaho.edu)
}


\seealso{
\link{MRPC} for estimating a DAG using the Mendelian Randomization (MR) based (MRPC) algorithm; \link{ModiSkeleton} for estimating a skeleton using modified skeleton function.
}
