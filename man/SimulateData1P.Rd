\name{SimulateData1P}
\alias{SimulateData1P}

\title{
Simulate data for a node with one parent
}
\description{
Simulate data for a node with one parent
}
\usage{
SimulateData1P(N, P1, b0.1, b1.1, sd.1)
}
\arguments{
  \item{N}{
Number of observations
}
\item{P1}{Data vector of the parent node P1.}

\item{b0.1}{Intercept of b0.1 + b1.1*P1, where P1 is the parent of the corresponding node.}

\item{b1.1}{Slope of P1 for b0.1 + b1.1*P1, where P1 is the parent of the corresponding node.}

\item{sd.1}{Standard deviation for corresponding data generated nodes.}
}

\value{
Vector
}

\author{
Md Bahadur Badsha (mdbadsha@uidaho.edu)
}
\seealso{
\link{SimulateData} for simulated data generating function.
}
\examples{
Data1P <- SimulateData1P(N = 10^3, 
                         P1 = 1, 
                         b0.1 = 0, 
                         b1.1 = 1, 
                         sd.1 = 1)
}
